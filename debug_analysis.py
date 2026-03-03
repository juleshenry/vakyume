#!/usr/bin/env python3
"""Debug script: pretty-prints the analysis.json for a project."""

import argparse
import json
import os
import sys


def main():
    parser = argparse.ArgumentParser(description="Print analysis.json summary")
    parser.add_argument(
        "project",
        nargs="?",
        default="projects/VacuumTheory",
        help="Project directory (default: projects/VacuumTheory)",
    )
    parser.add_argument(
        "--full", action="store_true", help="Print full JSON, not just summary"
    )
    args = parser.parse_args()

    path = os.path.join(args.project, "reports", "analysis.json")
    if not os.path.exists(path):
        print(f"Not found: {path}")
        sys.exit(1)

    with open(path) as f:
        data = json.load(f)

    solved = data.get("solved", [])
    inconsistent = data.get("inconsistent", {})
    failed = data.get("failed", [])

    print(f"=== {path} ===")
    print(f"Solved:       {len(solved)}")
    print(f"Inconsistent: {len(inconsistent)}")
    print(f"Failed:       {len(failed)}")
    print()

    if args.full:
        print(json.dumps(data, indent=2))
        return

    # --- Inconsistent detail ---
    if inconsistent:
        print("--- INCONSISTENT FAMILIES ---")
        for fam, info in sorted(inconsistent.items()):
            scores = info.get("scores", {})
            broken = info.get("broken", [])
            trusted = info.get("trusted", [])
            mismatches = info.get("mismatches", {})

            score_str = "  ".join(f"{v}={s}" for v, s in sorted(scores.items()))
            print(f"\n  {fam}")
            print(f"    scores:  {score_str}")
            print(f"    broken:  {broken}")
            print(f"    trusted: {trusted}")

            if mismatches:
                for var, trials in mismatches.items():
                    if not trials:
                        continue
                    print(f"    mismatches for '{var}':")
                    for i, trial in enumerate(trials[:3]):
                        if "inputs" in trial:
                            print(f"      trial {i}: inputs={trial['inputs']}")
                            print(f"               output={trial.get('output')}")
                            for m in trial.get("mismatches", []):
                                target = m.get("target", "?")
                                expected = m.get("expected", "?")
                                got = m.get("got", "?")
                                err = m.get("error")
                                if err:
                                    print(f"               vs {target}: ERROR {err}")
                                else:
                                    print(
                                        f"               vs {target}: expected={expected} got={got}"
                                    )
                        elif "error" in trial:
                            print(f"      trial {i}: ERROR {trial['error']}")

    # --- Failed detail ---
    if failed:
        print("\n--- FAILED FAMILIES ---")
        for entry in failed:
            name = entry.get("file", "?")
            err = entry.get("error", "?")
            print(f"  {name}: {err[:200]}")

    # --- Solved list (compact) ---
    if solved:
        print(f"\n--- SOLVED ({len(solved)}) ---")
        for s in sorted(solved):
            print(f"  {s}")


if __name__ == "__main__":
    main()
