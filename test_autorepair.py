import os
import sys
from vakyume.master import run_pipeline


def main():
    project_dir = "dummy_project"
    # We want to test --overwrite as requested
    # Note: run_pipeline signature is run_pipeline(project_dir=".", max_rounds=DEFAULT_MAX_ROUNDS, auto_repair=True, overwrite=False)
    print(f"Running pipeline on {project_dir} with overwrite=True...")
    analysis = run_pipeline(
        project_dir=project_dir, max_rounds=2, auto_repair=True, overwrite=True
    )
    print("\nPipeline finished.")
    print(f"Solved: {analysis.get('solved', [])}")
    print(f"Inconsistent: {list(analysis.get('inconsistent', {}).keys())}")
    print(f"Failed: {[f['file'] for f in analysis.get('failed', [])]}")


if __name__ == "__main__":
    main()
