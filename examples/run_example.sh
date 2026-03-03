#!/bin/bash

# Demonstrates the Vakyume pipeline on a small example project.
#
# Prerequisites:
#   pip install vakyume   (or run from the repo root)
#   Ollama running locally (for the repair stage)

set -euo pipefail

EXAMPLE_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$EXAMPLE_DIR/sample_project"

# Create the project skeleton expected by vakyume
mkdir -p "$PROJECT_DIR/notes"

# Write a tiny notes file
cat <<'EOF' > "$PROJECT_DIR/notes/15_heat_exchangers.py"
# 15-1 heat transfer rate
"""
Q := heat transfer rate
U := overall heat transfer coefficient
A := heat transfer area
delta_T := temperature difference
"""
Q = U * A * delta_T
EOF

echo "Running Vakyume pipeline on sample project..."
python3 -m vakyume run "$PROJECT_DIR" --max-rounds 3

echo ""
echo "Done.  Check $PROJECT_DIR/shards/ for generated solver shards."
