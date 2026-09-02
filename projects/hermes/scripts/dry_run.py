import os
import glob
import subprocess
import sys
import tempfile
from pathlib import Path
import yaml
 
REAL_CONFIG = Path(os.environ.get("HERMES_OPT_CONFIG", Path(__file__).parent / "config.yaml"))
SCRATCH_REPO = Path(tempfile.mkdtemp(prefix="optloop_dry_run_"))

DRY_RUN_INSTRUCTIONS = """\
This is a test, not a real optimization task. There is no real project
name -- do not search for, reference, or assume any specific named codebase. Only
look at files already present in your current working directory.
 
You will find a single C file with one loop marked TODO. Add an OpenMP
`#pragma omp parallel for` directive to that loop and apply the change immediately
by calling the file-write tool in this turn -- do not just describe the change.

Always use the full absolute path when reading or writing files (e.g.
/path/to/your/working/dir/src/mathops.c), not a bare relative path like
./mathops.c or src/mathops.c. Determine the absolute path by listing your working
directory first if you're unsure of it.
"""

def cleanup_stale_scratch_dirs() -> None:
    for stale in glob.glob(str(Path(tempfile.gettempdir()) / "optloop_dry_run_*")):
        subprocess.run(["rm", "-rf", stale])
        print(f"Removed stale scratch dir: {stale}")

def make_scratch_repo() -> None:
    # Throwaway git repo with a small C dummy file
    src = SCRATCH_REPO / "src"
    src.mkdir(parents=True)
    (src / "mathops.c").write_text(
        "#include <stdio.h>\n\n"
        "void compute(double *out, int n) {\n"
        "    for (int i = 0; i < n; i++) {\n"
        "        out[i] = i * 0.5;  // TODO: parallelize this loop\n"
        "    }\n"
        "}\n"
    )
    subprocess.run(["git", "init", "-q"], cwd=SCRATCH_REPO, check=True)
    subprocess.run(["git", "add", "-A"], cwd=SCRATCH_REPO, check=True)
    subprocess.run(
        ["git", "-c", "user.email=dryrun@test", "-c", "user.name=dryrun",
         "commit", "-q", "-m", "initial"],
        cwd=SCRATCH_REPO, check=True,
    )
    print(f"Scratch repo: {SCRATCH_REPO}")
 
def make_dry_config() -> Path:
    # Uses real config.yaml, with repo_dir/worktree_dir/experiment_dir
    with open(REAL_CONFIG) as f:
        cfg = yaml.safe_load(f)
    cfg["repo_dir"] = str(SCRATCH_REPO)
    cfg["worktree_dir"] = str(SCRATCH_REPO / "_worktree")
    cfg["experiment_dir"] = str(SCRATCH_REPO / "_experiment")
    cfg["num_iterations"] = 2
    cfg["max_agent_iterations"] = 8
    cfg["instructions"] = DRY_RUN_INSTRUCTIONS
    dry_config_path = SCRATCH_REPO / "dry_config.yaml"
    with open(dry_config_path, "w") as f:
        yaml.safe_dump(cfg, f)
    return dry_config_path
 
def make_fake_results(res_dir: Path) -> None:
    # Minimal, real-format LIKWID output, so the actual parser can be tested
    res_dir.mkdir(parents=True, exist_ok=True)
    (res_dir / "likwid_available_groups.txt").write_text(
        "Group name\tDescription\n" + "-" * 40 + "\nFLOPS_DP\tDouble Precision MFLOP/s\n"
    )
    csv = (
        "STRUCT,Info,3\n"
        "CPU name:,Fake CPU,\n"
        "CPU type:,Fake,\n"
        "CPU clock:,2.5 GHz,\n"
        "TABLE,Group 1 Metric,FLOPS_DP,1\n"
        "Metric,HWThread 0,\n"
        "Runtime (RDTSC) [s],1.0,\n"
    )
    (res_dir / "log.omp48.FLOPS_DP.csv").write_text(csv)
    (res_dir / "log.omp48.FLOPS_DP.txt").write_text(
        "TIMER_GROUP_INPUT = 0.1 s\nTIMER_GROUP_ANALYSIS = 0.9 s\nTIMER_TOTAL = 1.0 s\n"
    )
 
def fake_run_benchmark(self, experiment_id, env, res_dir, poll_interval=30, timeout=7200):
    print(f"  [dry-run] pretending to submit+run {experiment_id} (no SSH, no sbatch)")
    make_fake_results(res_dir)
    return True, res_dir
 
def main():
    cleanup_stale_scratch_dirs()
    make_scratch_repo()
    dry_config_path = make_dry_config()
    os.environ["HERMES_OPT_CONFIG"] = str(dry_config_path)  
 
    sys.path.insert(0, str(Path(__file__).parent))
    import projects.hermes.scripts.execution as execution
    execution.Executor.run_benchmark = fake_run_benchmark  # no real SSH/sbatch
 
    import projects.hermes.scripts.optimize as optimize
    sys.argv = ["optimize.py", "--ssh-user", "dry-run-unused", "--ssh-key-path", "dry-run-unused"]
    args = optimize.parse_args()
 
    try:
        optimize.optimize(args, dry_run=True)
    finally:
        print(f"\nDone. Scratch repo + logs: {SCRATCH_REPO}")
        print(f"Can be deleted via: rm -rf {SCRATCH_REPO}")
 
 
if __name__ == "__main__":
    main()
