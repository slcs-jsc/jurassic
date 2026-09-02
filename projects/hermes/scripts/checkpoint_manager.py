import subprocess
from pathlib import Path

class CheckpointManager:

    def __init__(self, repo_dir: Path, worktree_dir: Path, branch: str = "agent-optimization"):
        self.repo_dir = Path(repo_dir)
        self.worktree_dir = Path(worktree_dir)
        self.branch = branch

    def cmd(self, *args, cwd: Path | None = None) -> subprocess.CompletedProcess:
        cwd = cwd or self.worktree_dir
        return subprocess.run(
            ["git", "-C", str(cwd), *args], text=True, capture_output=True, check=True
        )

    def setup(self) -> None:
        branch_exists = self.cmd("branch", "--list", self.branch, cwd=self.repo_dir).stdout

        if not branch_exists.strip():
            self.cmd("branch", self.branch, cwd=self.repo_dir)
        if not self.worktree_dir.exists():
            self.cmd("worktree", "add", str(self.worktree_dir), self.branch, cwd=self.repo_dir)
            head = self.cmd("rev-parse", "HEAD").stdout.strip()
            self.cmd("tag", "-f", "last_good", head)
            print(f"Created new branch: {self.branch}")

    def commit(self, iteration: int, profile:dict) -> str:

        self.cmd("add", "-A")
        msg = f"iteration-{iteration:03d}: " + ", ".join(f"{k}={v}" for k, v in profile.items())
        self.cmd("commit", "-m", msg, "--allow-empty")
        head = self.cmd("rev-parse", "HEAD").stdout.strip()
        self.cmd("tag", "-f", "last_good", head)
        return head

    def discard(self) -> None:
        self.cmd("reset", "--hard", "HEAD")
        self.cmd("clean", "-fd")

    def rollback(self) -> None:
        self.cmd("reset", "--hard", "last_good")
        self.cmd("clean", "-fd")

    def diff_last_good(self) -> str:
        return self.cmd("diff", "last_good").stdout
    