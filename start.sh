#!/usr/bin/env bash
# Fresh clone -> everything: installs uv if needed, the environment, runs the
# test suite, verifies every project, and regenerates each project's docs.
#
#   ./start.sh            # everything (CPU heavy: whole-project verification)
#   QUICK=1 ./start.sh    # setup + fast tests only
#   JOBS=4 ./start.sh     # limit verification worker processes
set -euo pipefail
cd "$(dirname "$0")"

JOBS="${JOBS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 2)}"

if ! command -v uv >/dev/null 2>&1; then
  echo "==> installing uv"
  curl -LsSf https://astral.sh/uv/install.sh | sh
  export PATH="$HOME/.local/bin:$PATH"
fi

echo "==> installing dependencies (Python >= 3.11)"
uv sync

echo "==> fast tests"
uv run pytest -q

if [[ "${QUICK:-0}" == "1" ]]; then
  echo "==> QUICK=1: skipping project verification"
  exit 0
fi

for project in projects/*/; do
  project="${project%/}"
  echo "==> verifying $project ($JOBS workers)"
  uv run vakyume report "$project" --jobs "$JOBS"
  echo "==> docs for $project"
  uv run python -c "import sys; from vakyume.gen_docs import generate_equation_report as g; print(g(sys.argv[1]))" "$project"
done

echo
echo "Done. Per-project results: projects/*/docs/STATUS.md"
echo "Try: uv run vakyume solve projects/VacuumTheory 2-1 rho=1.2 D=0.1 v=3 mu=1.8e-5"
