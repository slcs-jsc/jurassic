import argparse
import datetime
import os
import subprocess
import json
import yaml
from dataclasses import asdict, dataclass
from pathlib import Path
from statistics import mean

from projects.hermes.scripts.checkpoint_manager import CheckpointManager
from projects.hermes.scripts.execution import Executor
from projects.hermes.scripts.parsing import parse_run_dir, summarize
from projects.hermes.scripts.plot_results import plot_agent_progress, plot_metrics

from run_agent import AIAgent

with open(Path(os.environ.get("HERMES_OPT_CONFIG", Path(__file__).parent / "config.yaml"))) as f: 
    cfg = yaml.safe_load(f)

timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
RESULTS_FILE = (Path(cfg["experiment_dir"]) / f"run_{timestamp}"/ cfg["json_out"]).resolve()
RESULTS_DIR = (Path(cfg["experiment_dir"]) / f"run_{timestamp}"/ cfg["res_dir"]).resolve()

@dataclass
class Experiment:
    iteration: int
    run_id: str
    accepted: bool
    score: float | None
    reason: str 
    diff: str
    response: str = ""

# Helper Functions
def run_command( command: list[str], *, timeout: int = 3600) -> subprocess.CompletedProcess[str]: 
    print(f"\n$ {' '.join(command)}") 
    return subprocess.run(command, text=True, capture_output=True, timeout=timeout)

def check_validation(res_dir: Path, dry_run=False) -> tuple[bool, str]:
    if dry_run:
        return True, "Correctness verified (dry run)."

    status_file = res_dir / "validation_status.txt"
    if not status_file.exists():
        return True, "No validation_status.txt found. Correctness not verified."

    content = status_file.read_text().strip()
    _, _, value = content.partition("=")
    if value == "skipped":
        return True, "Validation explicitly skipped for this run (SKIP_VALIDATION=1). Correctness not verified."
    try: 
        return_code =  int(value)
    except ValueError: 
        return False, f"Could not parse validation_status.txt: {content!r}"
    
    log_path = res_dir / "validation.log"
    log_tail = log_path.read_text()[-2000:] if log_path.exists() else "(no validation.log found)"
    return return_code == 0, f"validation exit_code={return_code}\n{log_tail}"

def score(configs: list, thread_count: int = cfg["thread_count_for_scoring"]) -> float | None:
    values = [
        c["timers"]["TIMER_GROUP_ANALYSIS"]
        for c in configs
        if c["omp_threads"] == thread_count and "TIMER_GROUP_ANALYSIS" in c.get("timers", {})
    ]
    return mean(values) if values else None

def log_experiment(experiment: Experiment) -> None: 
    RESULTS_FILE.parent.mkdir(parents=True, exist_ok=True)
    with open(RESULTS_FILE, "a") as f:
        f.write(json.dumps(asdict(experiment)) + "\n")

def build_prompt(benchmark_summary: str, diff: str, last_res: str) -> str: 
    prompt = [f"Last attempt: {last_res}\n"]
    if diff.strip():
        prompt.append(f"Diff of the last change:\n```diff\n{diff}\n```\n")
    prompt.append(f"Current benchmark results:\n{benchmark_summary}\n")
    prompt.append("Propose and apply the next optimization.")
    return "\n".join(prompt)

def optimize(args: argparse.Namespace, dry_run=False):
    checkpoints = CheckpointManager(
        repo_dir=Path(cfg["repo_dir"]), 
        worktree_dir=Path(cfg["worktree_dir"]), 
    )
    checkpoints.setup()

    executor = Executor(        
        ssh_host=args.ssh_host,
        ssh_user=args.ssh_user,
        ssh_key_path=args.ssh_key_path,
        remote_work_dir=args.remote_work_dir,
        local_work_dir=Path(cfg["worktree_dir"]),
        benchmark_dir=args.benchmark_dir,
        script=args.run_script,
        sync_excludes=tuple(args.sync_excludes) if args.sync_excludes else
            (".git/", "projects/benchmark/runs/", "projects/hermes/"),
    )
    executor.ensure_connection() # keep ssh connection alive via ControlMaster

    worktree_abs = Path(cfg["worktree_dir"]).resolve()
    os.environ["HERMES_WRITE_SAFE_ROOT"] = str(worktree_abs)
    os.environ["TERMINAL_CWD"] = str(worktree_abs)
    os.chdir(worktree_abs)

    cwd_notice = (
        f"\n\nYour working directory for ALL file reads and writes is exactly:\n"
        f"{worktree_abs}\n"
        f"Use this absolute path directly, prefixed to every file path you read or write "
        f"(e.g. {worktree_abs}/src/mathops.c). Do not rely on os.getcwd(), search_files, or "
        f"execute_code to discover this -- those tools may report a different, incorrect "
        f"location in this environment. Treat the path above as ground truth."
    )

    agent = AIAgent(
        model=cfg["model_name"],
        quiet_mode=True,        # set False for debugging
        ephemeral_system_prompt=cfg["instructions"] + cwd_notice,
        disabled_toolsets=["terminal", "computer_use"], 
        skip_memory=True, 
        max_iterations=cfg["max_agent_iterations"]
    )
    conversation_history = []

    # Peform baselin run
    success, res_dir = executor.run_benchmark(
        f"run_{timestamp}_iter_000_baseline", cfg["benchmark_env"], RESULTS_DIR / "iter_000_baseline", 
        poll_interval=args.poll_interval, timeout=args.job_timeout,
    )
    if not success: 
        raise RuntimeError(f"Could not execute baseline. Abort optimization loop.") 

    baseline_configs = parse_run_dir(res_dir)
    best_score = score(baseline_configs)
    best_configs = baseline_configs
    best_summary = summarize(baseline_configs)

    print(f"Baseline score (TIMER_GROUP_ANALYSIS @ {cfg['thread_count_for_scoring']} threads): {best_score}")

    last_res = "None yet. This is the first attempt"
    last_diff = ""

    for i in range(1, cfg["num_iterations"] + 1):
        experiment_id = f"run_{timestamp}_iteration_{i:03d}"
        prompt = build_prompt(best_summary, last_diff, last_res)
        result = agent.run_conversation(
            user_message=prompt, 
            conversation_history=conversation_history, 
            task_id=experiment_id
        )

        print("RAW RESULT:", json.dumps(result, indent=2, default=str))

        conversation_history = result["messages"]
        agent_response = result["final_response"]
    
        diff = checkpoints.diff_last_good()
        if not diff.strip():
            last_res = "No changes where made. Try a different directive"
            last_diff = ""
            log_experiment(Experiment(i, experiment_id, False, None, last_res, last_diff, agent_response))
            continue

        success, res_dir = executor.run_benchmark(
            experiment_id, cfg["benchmark_env"], RESULTS_DIR / experiment_id, 
            poll_interval=args.poll_interval, timeout=args.job_timeout,
        )

        skipped_log = res_dir / "skipped_profiling.txt"
        if skipped_log.exists():
            # Either validation failed or none of the requested LIKWID_GROUPS is available
            checkpoints.discard()
            last_res = f"Profiling was skipped on the cluster ({skipped_log.read_text().strip()}). Change reverted."
            last_diff = diff
            log_experiment(Experiment(i, experiment_id, False, None, last_res, last_diff, agent_response))
            continue

        if not success: 
            checkpoints.discard()
            last_res =  f"Benchmark job failed to complete, change reverted."
            last_diff = diff
            log_experiment(Experiment(i, experiment_id, False, None, last_res, last_diff, agent_response))
            continue
        
        correct, test_out = check_validation(res_dir, dry_run)
        if not correct: 
            checkpoints.discard()
            last_res = f"Correctness check failed, change reverted:\n{test_out[-2000:]}"
            last_diff = diff
            log_experiment(Experiment(i, experiment_id, False, None, last_res, last_diff, agent_response))
            continue

        new_configs = parse_run_dir(res_dir)
        new_score = score(new_configs)

        # If candidate > baseline Then baseline = candidate 
        # else rollback changes 
        if new_score is not None and best_score is not None and new_score < best_score:
            checkpoints.commit(i, {"score": new_score})
            last_res = f"Accepted. Score improved from {best_score:.3f}s to {new_score:.3f}s"
            best_score = new_score
            best_configs = new_configs
            best_summary = summarize(new_configs)
            log_experiment(Experiment(i, experiment_id, True, new_score, last_res, diff, agent_response))
        else:
            checkpoints.rollback()
            last_res = f"Rejected. Score {new_score} did not improve on previous score {best_score}. Change reverted"
            log_experiment(Experiment(i, experiment_id, False, new_score, last_res, diff, agent_response))

        last_diff = diff
        plot_metrics(new_configs, res_dir)

    plot_agent_progress(RESULTS_FILE, res_dir)
    print(f"\nFinished. Best score: {best_score}")


def parse_args() -> argparse.Namespace:
    executor_defaults = cfg.get("executor") or {}
    parser = argparse.ArgumentParser(description="LLM-optimization loop for JURASSIC")
    parser.add_argument("--ssh-host", default=executor_defaults.get("ssh_host", "juwels.fz-juelich.de"),
                        help="Cluster login node hostname")
    parser.add_argument("--ssh-user", default=executor_defaults.get("ssh_user") or os.environ.get("HERMES_SSH_USER"), 
                        help="Cluster username or set HERMES_SSH_USER")
    parser.add_argument("--ssh-key-path", default=executor_defaults.get("ssh_key_path") or os.environ.get("HERMES_SSH_KEY"), 
                        help="Path to the SSH private key (or set HERMES_SSH_KEY)")
    parser.add_argument("--remote-work-dir", default=executor_defaults.get("remote_work_dir", cfg["repo_dir"]), 
                        help="Absolute path to the jurassic repo root on the cluster")
    parser.add_argument("--benchmark-dir", default=executor_defaults.get("benchmark_dir", "projects/benchmark"), 
                        help="Path to the benchmark subdir, relative to the remote repo root")
    parser.add_argument("--run-script", default=executor_defaults.get("run_script", "run_juwels_booster_likwid.sh"), 
                        help="sbatch script name, run from <benchmark-subdir>/scripts")
    parser.add_argument("--sync-exclude", dest="sync_excludes", action="append",
                        default=executor_defaults.get("sync_excludes"),
                        help="rsync --exclude pattern when pushing the worktree to the cluster")
    parser.add_argument("--poll-interval", type=int, default=executor_defaults.get("poll_interval", 30), 
                        help="Seconds between SLURM squeue polls")
    parser.add_argument("--job-timeout", type=int, default=executor_defaults.get("job_timeout", 7200),
                        help="Seconds to wait for a benchmark job to finish before giving up")

    args = parser.parse_args()

    if args.remote_work_dir is None:
        parser.error("--remote_work_dir (or executor.remote_work_dir in config.yaml) is required")
    if args.ssh_user is None:
        parser.error("--ssh-user (or executor.ssh_user in config.yaml, or HERMES_SSH_USER) is required")
    if args.ssh_key_path is None:
        parser.error("--ssh-key-path (or executor.ssh_key_path in config.yaml, or HERMES_SSH_KEY) is required")

    return args

if __name__ == "__main__":
    optimize(parse_args())
