import subprocess
from pathlib import Path
import shlex
import time
import os 

class Executor():

    def __init__(
        self,
        ssh_host: str, 
        ssh_user: str, 
        ssh_key_path: str, 
        remote_work_dir: str, 
        local_work_dir: Path = Path("."),
        benchmark_dir: str = "projects/benchmark", 
        script: str = "run_juwels_booster_likwid.sh",
        sync_excludes: tuple = (".git/", "projects/benchmark/runs/", "projects/hermes/"),
    ):
        self.ssh_host = ssh_host
        self.ssh_user = ssh_user
        self.ssh_key_path = ssh_key_path
        self.remote_work_dir = remote_work_dir.rstrip("/")
        self.local_work_dir = Path(local_work_dir)
        self.benchmark_dir = benchmark_dir.strip("/")
        self.script = script
        self.sync_excludes = sync_excludes

    def ensure_connection(self) -> None: 
        control_dir = os.path.expanduser("~/.ssh/controlmasters")
        os.makedirs(control_dir, mode=0o700, exist_ok=True)
        self.control_path = f"{control_dir}/%r@%h:%p"

        check = subprocess.run(
            ["ssh", "-O", "check", "-o", f"ControlPath={self.control_path}",
            f"{self.ssh_user}@{self.ssh_host}"],
            capture_output=True, text=True,
        )
        if check.returncode == 0:
            print(f"Reusing existing SSH connection to {self.ssh_host}.")
            return

        print(f"No active SSH connection to {self.ssh_host}. Authentication needed.")

        establish = subprocess.run([
            "ssh", "-M", "-S", self.control_path, "-o", "ControlPersist=6h",
            f"{self.ssh_user}@{self.ssh_host}", "true",
        ])

        if establish.returncode != 0:
            raise RuntimeError("Failed to establish SSH control connection.")
        print(f"SSH connection to {self.ssh_host} established.")


    def ssh_base(self) -> list:
        return ["ssh", "-o", f"ControlPath={self.control_path}", "-i", self.ssh_key_path, f"{self.ssh_user}@{self.ssh_host}"]
 
    def run(self, cmd: list, timeout: int = 600) -> subprocess.CompletedProcess:
        print(f"\n$ {' '.join(shlex.quote(c) for c in cmd)}")
        return subprocess.run(cmd, text=True, capture_output=True, timeout=timeout)
 
    def rsync(self) -> None: 
        cmd = ["rsync", "-avz", "--delete", "-e", f"ssh -i {self.ssh_key_path} -o ControlPath={self.control_path}"]
        for pattern in self.sync_excludes:
            cmd += ["--exclude", pattern]
        cmd += [f"{self.local_work_dir}/", f"{self.ssh_user}@{self.ssh_host}:{self.remote_work_dir}/"]
        res = self.run(cmd, timeout=1800)
        if res.returncode != 0:
            raise RuntimeError(f"rsync to remote failed:\n{res.stderr}")

    def submit_job(self, experiment_id: str, env: dict) -> str:
        export_vars = {"RUN_ID": experiment_id, **env}
        export_str = "ALL," + ",".join(f"{k}={v}" for k, v in export_vars.items())
        remote_scripts_dir = f"{self.remote_work_dir}/{self.benchmark_dir}/scripts"
        remote_cmd = (
            f"cd {shlex.quote(remote_scripts_dir)} && "
            f"sbatch --export={shlex.quote(export_str)} {shlex.quote(self.script)}"
        )
        res = self.run(self.ssh_base() + [remote_cmd])
        if res.returncode != 0:
            raise RuntimeError(f"sbatch submission failed:\n{res.stderr}")
        for token in res.stdout.split():  
            if token.isdigit():
                return token
        raise RuntimeError(f"Could not parse job id from sbatch output:\n{res.stdout}")

    def wait_for_job(self, job_id: str, poll_interval: int = 30, timeout: int = 7200) -> str:
        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline:
            squeue = self.run(self.ssh_base() + [f"squeue -j {job_id} -h -o %T"])
            state = squeue.stdout.strip()
            if not state:
                break
            print(f"Job {job_id}: {state}")
            time.sleep(poll_interval)
        else:
            raise TimeoutError(f"Job {job_id} did not finish within {timeout}s")

        for attempt in range(5):
            sacct = self.run(
                self.ssh_base() + [f"sacct -j {job_id} --format=State --noheader --parsable2"]
            )
            first_line = sacct.stdout.strip().splitlines()[0] if sacct.stdout.strip() else ""
            if first_line.strip():
                return first_line.strip()
            print(f"  sacct returned nothing yet for job {job_id}, retrying ({attempt + 1}/5)...")
            time.sleep(5)

        return "UNKNOWN"

    def get_results(self, experiment_id: str, res_dir: Path) -> Path:
        remote_run_dir =  f"{self.remote_work_dir}/{self.benchmark_dir}/runs/{experiment_id}/"
        res_dir.mkdir(parents=True, exist_ok=True)
        cmd = [
            "rsync", "-avz", "-e", f"ssh -i {self.ssh_key_path} -o ControlPath={self.control_path}",
            f"{self.ssh_user}@{self.ssh_host}:{remote_run_dir}",
            f"{res_dir}/",
        ]
        res = self.run(cmd, timeout=600)
        if res.returncode != 0:
            raise RuntimeError(f"rsync fetch failed:\n{res.stderr}")
        return res_dir

    def run_benchmark(self, experiment_id: str, env: dict, res_dir: Path, 
                      poll_interval: int = 30, timeout: int = 7200) -> tuple:
        self.rsync()
        job_id = self.submit_job(experiment_id, env)
        print(f"Submitted job {job_id} (experiment_id={experiment_id})")
        state = self.wait_for_job(job_id, poll_interval=poll_interval, timeout=timeout)
        res_dir = self.get_results(experiment_id, res_dir)
        print(state)
        return state == "COMPLETED", res_dir
 
