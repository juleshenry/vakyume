#!/usr/bin/env bash
#
# benchmark_llm_raw.sh
#
# Benchmark: how many equation families can each model solve
# with NO scaffolds, maxrounds 5? (LLM-only repair)
#
# Models tested:
#   1. phi3:latest       (best known: 148/148 pipeline, 144/148 pytest)
#   2. llama3.2:latest
#
# Two scores per model:
#   - Pipeline score:  solved/inconsistent/failed (golden-tuple verification
#                      with 90% noise tolerance)
#   - Pytest score:    strict cross-validation via Verify.verify() (random
#                      inputs, no noise tolerance). Captures the honest
#                      pass rate. Run with --tb=no -q for concise output.
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

    # Clean reports from previous run so this model starts fresh
    rm -rf "$PROJECT/reports"

    # Wipe shards so we start fresh (--overwrite does this)
    # || true: allow the pipeline to exit non-zero (e.g. internal pytest
    # failures or unsolved shards) without killing the benchmark.
    python3 vakyume.py \
        --llm-model "$model" \
        run "$PROJECT" \
        --max-rounds "$MAX_ROUNDS" \
        --overwrite \
        2>&1 | tee "$RESULTS_DIR/${label}_${TIMESTAMP}.log" || true

    # Extract the final counts from the pipeline output.
    # vakyume.py prints authoritative counts at the end as:
    #   Solved: N
    #   Inconsistent: N
    #   Failed: N
    # We anchor on "^Solved:" etc. to avoid false matches from SymPy/pytest
    # lines that also contain "solved" or "failed".
    local log_file="$RESULTS_DIR/${label}_${TIMESTAMP}.log"
    local solved=$(grep -E '^Solved: ' "$log_file" | tail -1 | awk '{print $2}')
    local inconsistent=$(grep -E '^Inconsistent: ' "$log_file" | tail -1 | awk '{print $2}')
    local failed=$(grep -E '^Failed: ' "$log_file" | tail -1 | awk '{print $2}')

    # Also grab the analysis.json for exact numbers
    if [ -f "$PROJECT/reports/analysis.json" ]; then
        cp "$PROJECT/reports/analysis.json" "$RESULTS_DIR/${label}_analysis_${TIMESTAMP}.json"
    fi

    # Run pytest for strict cross-validation score
    echo ""
    echo "  Running pytest strict cross-validation..."
    local pytest_log="$RESULTS_DIR/${label}_pytest_${TIMESTAMP}.log"
    python3 -m pytest tests/test_shards.py -q --tb=no 2>&1 | tee "$pytest_log" || true

    # Extract pytest pass/fail counts
    # Format: "144 passed, 4 failed" or "148 passed"
    local pytest_passed=$(grep -oE '[0-9]+ passed' "$pytest_log" | awk '{print $1}')
    local pytest_failed=$(grep -oE '[0-9]+ failed' "$pytest_log" | awk '{print $1}')
    pytest_passed=${pytest_passed:-0}
    pytest_failed=${pytest_failed:-0}
    local pytest_total=$((pytest_passed + pytest_failed))

    echo ""
    echo "------------------------------------------------------------"
    echo "  RESULT: $label"
    echo "    Pipeline:"
    echo "      Solved:       ${solved:-?}"
    echo "      Inconsistent: ${inconsistent:-?}"
    echo "      Failed:       ${failed:-?}"
    echo "    Pytest (strict):"
    echo "      Passed:       ${pytest_passed}/${pytest_total}"
    echo "------------------------------------------------------------"
    echo ""

    # Append to summary
    {
        echo "$label"
        echo "  Pipeline:  Solved: ${solved:-?} | Inconsistent: ${inconsistent:-?} | Failed: ${failed:-?}"
        echo "  Pytest:    ${pytest_passed}/${pytest_total} passed (strict cross-validation)"
        echo ""
    } >> "$RESULTS_DIR/summary_${TIMESTAMP}.txt"
}

# ── Summary header ─────────────────────────────────────────────────────────
{
    echo "Vakyume LLM-raw benchmark (no scaffolds, max_rounds=$MAX_ROUNDS)"
    echo "Started: $(date)"
    echo "---"
} > "$RESULTS_DIR/summary_${TIMESTAMP}.txt"

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
echo "Full logs:      $RESULTS_DIR/*_${TIMESTAMP}.log"
echo "Pytest logs:    $RESULTS_DIR/*_pytest_${TIMESTAMP}.log"
echo "Analysis JSON:  $RESULTS_DIR/*_analysis_${TIMESTAMP}.json"
echo "Summary:        $RESULTS_DIR/summary_${TIMESTAMP}.txt"
