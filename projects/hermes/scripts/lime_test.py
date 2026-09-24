from __future__ import annotations

import json
import re
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Callable, Iterable, Sequence

import numpy as np


def default_focus_score(agent_response_text: str, target_functions: Sequence[str]) -> float:
    if not target_functions:
        return 0.0

    text = agent_response_text.lower()
    per_fn_scores = []
    for fn in target_functions:
        fn_l = fn.lower()
        count = len(re.findall(re.escape(fn_l), text))
        per_fn_scores.append(min(1.0, 0.7 * min(count, 1) + 0.2 * min(max(count - 1, 0), 1)))

    return float(max(per_fn_scores)) if per_fn_scores else 0.0


def diff_aware_focus_score(agent_response_text: str, target_functions: Sequence[str]) -> float:
    if not target_functions:
        return 0.0

    diff_lines = [
        line for line in agent_response_text.splitlines()
        if line.startswith(("+", "-")) and not line.startswith(("+++", "---"))
    ]
    if not diff_lines:
        return 0.0

    diff_text = "\n".join(diff_lines).lower()
    hits = sum(1 for fn in target_functions if fn.lower() in diff_text)
    return min(1.0, hits / len(target_functions))


def make_single_shot_agent_call(
    agent,
    base_instructions: str,
) -> Callable[[str], str]:
    call_counter = {"n": 0}

    def call(prompt_text: str) -> str:
        call_counter["n"] += 1
        task_id = f"lime_probe_{call_counter['n']:04d}"

        result = agent.run_conversation(
            user_message=prompt_text,
            conversation_history=[],
            task_id=task_id,
        )
        response = result["final_response"]
        if not response:
            print(f"[lime_explain] WARNING: empty/failed response for {task_id}")
            return ""
        return response

    return call


def make_lime_classifier(
    build_prompt_fn: Callable[[str], str],
    agent_call_fn: Callable[[str], str],
    score_fn: Callable[[str], float],
    cache: dict,
    num_workers: int = 8,
    log_every: int = 10,
):
    calls_made = {"n": 0}

    def classify_one(text: str) -> float:
        if text in cache:
            return cache[text]
        prompt = build_prompt_fn(text)
        response = agent_call_fn(prompt)
        score = score_fn(response)
        cache[text] = score
        return score

    def classifier_fn(perturbed_texts: Iterable[str]) -> np.ndarray:
        texts = list(perturbed_texts)
        to_run = [t for t in texts if t not in cache]

        if to_run:
            with ThreadPoolExecutor(max_workers=num_workers) as pool:
                for text, score in zip(to_run, pool.map(classify_one, to_run)):
                    cache[text] = score

        calls_made["n"] += len(to_run)
        if log_every and calls_made["n"] % log_every < len(to_run):
            print(f"[lime_explain] agent calls so far: {calls_made['n']} "
                  f"(unique cached prompts: {len(cache)})")

        return np.array([[1.0 - cache[t], cache[t]] for t in texts])

    return classifier_fn


def explain_baseline_profile(
    profile_text: str,
    build_prompt_fn: Callable[[str], str],
    agent_call_fn: Callable[[str], str],
    target_functions: Sequence[str],
    out_dir: Path | str,
    score_fn: Callable[[str], float] | None = None,
    num_samples: int = 200,
    num_features: int | None = None,
    num_workers: int = 8,
) -> list[tuple[str, float]]:
    from lime.lime_text import LimeTextExplainer

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if score_fn is None:
        score_fn = lambda resp: default_focus_score(resp, target_functions)

    n_rows = len([r for r in profile_text.split("\n") if r.strip()])
    if num_features is None:
        num_features = min(15, n_rows)

    if n_rows > 25:
        print(f"[lime_explain] WARNING: profile has {n_rows} rows, "
              f"num_samples={num_samples} may be too low for a stable fit.")

    cache: dict[str, float] = {}
    classifier_fn = make_lime_classifier(
        build_prompt_fn=build_prompt_fn,
        agent_call_fn=agent_call_fn,
        score_fn=score_fn,
        cache=cache,
        num_workers=num_workers,
    )

    explainer = LimeTextExplainer(
        class_names=["ignored", "focused"],
        split_expression=lambda x: x.split("\n"),
        bow=False,
    )

    print(f"[lime_explain] explaining baseline profile ({n_rows} rows), "
          f"num_samples={num_samples}, target_functions={list(target_functions)}")

    exp = explainer.explain_instance(
        text_instance=profile_text,
        classifier_fn=classifier_fn,
        num_features=num_features,
        num_samples=num_samples,
    )

    weights = exp.as_list()

    (out_dir / "row_importance.json").write_text(
        json.dumps(
            {
                "target_functions": list(target_functions),
                "num_samples": num_samples,
                "num_unique_agent_calls": len(cache),
                "weights": weights,
            },
            indent=2,
        )
    )
    exp.save_to_file(str(out_dir / "lime_baseline_profile.html"))
    (out_dir / "baseline_profile.txt").write_text(profile_text)

    print(f"[lime_explain] done. {len(cache)} unique agent calls made "
          f"(num_samples={num_samples} requested).")
    print("[lime_explain] row importance (most influential first):")
    for row, weight in weights:
        print(f"  {weight:+.4f}  {row!r}")
    print(f"[lime_explain] wrote {out_dir / 'lime_baseline_profile.html'} "
          f"and {out_dir / 'row_importance.json'}")

    return weights