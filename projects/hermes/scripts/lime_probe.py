#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
import time
import threading
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from lime_test import explain_baseline_profile, make_single_shot_agent_call  

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--profile", required=True, type=Path)
    p.add_argument("--target-function", action="append", required=True, dest="target_functions")
    p.add_argument("--instructions-file", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--model", required=True)
    p.add_argument("--max-agent-iterations", type=int, default=2)
    p.add_argument("--num-samples", type=int, default=200)
    p.add_argument("--num-workers", type=int, default=1)
    p.add_argument("--min-delay-s", type=float, default=1.0)
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args()

def build_agent(instructions_text: str, model: str, max_agent_iterations: int):
    from run_agent import AIAgent

    return AIAgent(
        model=model,
        quiet_mode=True,
        ephemeral_system_prompt=instructions_text,
        disabled_toolsets=["terminal", "computer_use", "execute_code"],
        skip_memory=True,
        max_iterations=max_agent_iterations,
    )

def build_prompt(instructions_text: str, profile_text: str) -> str:
    return (
        f"{instructions_text}\n\n"
        "Here is the current performance profile:\n\n"
        f"{profile_text}\n\n"
        "Which function should be optimized first, and why? "
        "If you would change the code, name the exact function(s) you would edit."
    )

def throttle(fn, min_delay_s: float):
    lock = threading.Lock()
    last_call = {"t": 0.0}
 
    def wrapped(*args, **kwargs):
        with lock:
            wait = last_call["t"] + min_delay_s - time.monotonic()
            if wait > 0:
                time.sleep(wait)
            last_call["t"] = time.monotonic()
        return fn(*args, **kwargs)
 
    return wrapped

def main() -> None:
    args = parse_args()

    if not args.profile.is_file():
        sys.exit(f"--profile not found: {args.profile}")
    if not args.instructions_file.is_file():
        sys.exit(f"--instructions-file not found: {args.instructions_file}")

    profile_text = args.profile.read_text().rstrip("\n")
    instructions_text = args.instructions_file.read_text()

    if args.dry_run:
        rows = profile_text.split("\n")
        print(f"Profile: {args.profile} ({len(rows)} rows)")
        print("-" * 60)
        print(profile_text)
        print("-" * 60)
        example_perturbed = "\n".join(rows[1:])
        print(build_prompt(instructions_text, example_perturbed))
        return

    build_agent_fn = lambda: build_agent(instructions_text, args.model, args.max_agent_iterations)
    agent_call_fn = make_single_shot_agent_call(build_agent_fn=build_agent_fn, base_instructions=instructions_text)
    agent_call_fn = throttle(agent_call_fn, args.min_delay_s)

    def build_prompt_fn(perturbed_profile_text: str) -> str:
        return build_prompt(instructions_text, perturbed_profile_text)

    explain_baseline_profile(
        profile_text=profile_text,
        build_prompt_fn=build_prompt_fn,
        agent_call_fn=agent_call_fn,
        target_functions=args.target_functions,
        out_dir=args.out_dir,
        num_samples=args.num_samples,
        num_workers=args.num_workers,
    )


if __name__ == "__main__":
    main()