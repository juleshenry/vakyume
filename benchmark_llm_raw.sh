#!/usr/bin/env bash
#
# benchmark_llm_raw.sh
#
# Benchmark: how many equation families can each model solve
# with NO scaffolds, maxrounds 5? (LLM-only repair)
#
# Models tested:
#   1. phi3:latest
#   2. llama3.2:latest
#
# Results saved to benchmark_results/

set -euo pipefail

PROJECT="projects/VacuumTheory"
MAX_ROUNDS=5
RESULTS_DIR="benchmark_results"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

mkdir -p "$RESULTS_DIR"

# ── Helper: run one benchmark ──────────────────────────────────────────────
run_benchmark() {
    local model="$1"
    local label="$2"

    echo ""
    echo "============================================================"
    echo "  BENCHMARK: $label  (model=$model, max_rounds=$MAX_ROUNDS, no scaffolds)"
    echo "============================================================"
    echo ""

    # Wipe shards so we start fresh (--overwrite does this)
    python3 vakyume.py \
        --llm-model "$model" \
        run "$PROJECT" \
        --max-rounds "$MAX_ROUNDS" \
        --overwrite \
        2>&1 | tee "$RESULTS_DIR/${label}_${TIMESTAMP}.log"

    # Extract the final counts from the pipeline output
    local log_file="$RESULTS_DIR/${label}_${TIMESTAMP}.log"
    local solved=$(grep -E "^Solved:" "$log_file" | tail -1 | awk '{print $2}')
    local inconsistent=$(grep -E "^Inconsistent:" "$log_file" | tail -1 | awk '{print $2}')
    local failed=$(grep -E "^Failed:" "$log_file" | tail -1 | awk '{print $2}')

    # Also grab the analysis.json for exact numbers
    if [ -f "$PROJECT/reports/analysis.json" ]; then
        cp "$PROJECT/reports/analysis.json" "$RESULTS_DIR/${label}_analysis_${TIMESTAMP}.json"
    fi

    echo ""
    echo "------------------------------------------------------------"
    echo "  RESULT: $label"
    echo "    Solved:       ${solved:-?}"
    echo "    Inconsistent: ${inconsistent:-?}"
    echo "    Failed:       ${failed:-?}"
    echo "------------------------------------------------------------"
    echo ""

    # Append to summary
    echo "$label | Solved: ${solved:-?} | Inconsistent: ${inconsistent:-?} | Failed: ${failed:-?}" \
        >> "$RESULTS_DIR/summary_${TIMESTAMP}.txt"
}

# ── Summary header ─────────────────────────────────────────────────────────
echo "Vakyume LLM-raw benchmark (no scaffolds, max_rounds=$MAX_ROUNDS)" \
    > "$RESULTS_DIR/summary_${TIMESTAMP}.txt"
echo "Started: $(date)" >> "$RESULTS_DIR/summary_${TIMESTAMP}.txt"
echo "---" >> "$RESULTS_DIR/summary_${TIMESTAMP}.txt"

# ── Run benchmarks ─────────────────────────────────────────────────────────
run_benchmark "phi3:latest"     "phi3"
run_benchmark "llama3.2:latest" "llama3.2"

# ── Print summary ──────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo "  FINAL SUMMARY"
echo "============================================================"
cat "$RESULTS_DIR/summary_${TIMESTAMP}.txt"
echo ""
echo "Full logs:     $RESULTS_DIR/*_${TIMESTAMP}.log"
echo "Analysis JSON: $RESULTS_DIR/*_analysis_${TIMESTAMP}.json"
echo "Summary:       $RESULTS_DIR/summary_${TIMESTAMP}.txt"
