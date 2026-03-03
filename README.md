# Vakyume
# A Python Library for Vacuum Design

Inspired by the 1986 edition of *Process Vacuum System Design and Operation* by James L. Ryans and Daniel L. Roper, this project offers functions related to engineering vacuums.

We utilize a "One-Odd-Out" (OOO) verification methodology: for any equation family (e.g., $PV=nRT$), if solvers exist for every variable ($P, V, n, R, T$), they must be mathematically consistent. If one solver produces a result that breaks the equality when plugged into the others, it is flagged for LLM-assisted repair.

---

## Table of Contents

- [Quick Start](#quick-start)
- [PDF to C++ Pipeline](#pdf-to-c-pipeline)
- [CLI Reference](#cli-reference)
- [PDF Scraper Tips](#pdf-scraper-tips)
- [Project Structure](#project-structure)
- [Setting Up a New Project](#setting-up-a-new-project)
- [Architecture](#architecture)
  - [Module Separation of Concerns](#module-separation-of-concerns)
  - [Module Dependency Graph](#module-dependency-graph)
- [Algorithm Deep Dive](#algorithm-deep-dive)
  - [Stage 1: Parse Equation Notes](#stage-1-parse-equation-notes)
  - [Stage 2: Generate Solvers (SymPy)](#stage-2-generate-solvers-sympy)
  - [Stage 3: Scaffold Pattern Matching](#stage-3-scaffold-pattern-matching)
  - [Stage 4: LLM Delegation](#stage-4-llm-delegation)
  - [Stage 5: One-Odd-Out Verification](#stage-5-one-odd-out-verification)
  - [Stage 6: Repair Loop](#stage-6-repair-loop)
  - [Stage 7: Reconstruct Library](#stage-7-reconstruct-library)
  - [Stage 8: C++ Transpile](#stage-8-c-transpile)
- [Scaffold Algorithm Detail](#scaffold-algorithm-detail)
  - [Pattern Handler Registry](#pattern-handler-registry)
  - [Branch 1: Transcendental (Target in Exponent)](#branch-1-transcendental-target-in-exponent)
  - [Branch 2: Ratio Power -- Target in Numerator](#branch-2-ratio-power----target-in-numerator)
  - [Branch 3: Ratio Power -- Target in Denominator](#branch-3-ratio-power----target-in-denominator)
  - [Branch 4: Ratio Power -- Re-match Numerator](#branch-4-ratio-power----re-match-numerator)
  - [Branch 5: Bare Power](#branch-5-bare-power)
  - [Branch 6: Fractional Exponent with Multiple Difference Terms](#branch-6-fractional-exponent-with-multiple-difference-terms)
  - [Branch 7: Fractional Exponent with Single Occurrence](#branch-7-fractional-exponent-with-single-occurrence)
  - [Fallbacks F1 and F2](#fallbacks-f1-and-f2)
- [LLM Delegation: Why Algebraic Step-by-Step?](#llm-delegation-why-algebraic-step-by-step)
  - [The Problem with Unconstrained LLM Solving](#the-problem-with-unconstrained-llm-solving)
  - [Scaffold as Algebraic Guardrail](#scaffold-as-algebraic-guardrail)
  - [Prompt Architecture](#prompt-architecture)
  - [Bypass vs. Fill-in vs. Full Delegation](#bypass-vs-fill-in-vs-full-delegation)
- [The `@kwasak` Decorator](#the-kwasak-decorator)
- [Resume & Testing](#resume--testing)
- [Examples](#examples)
- [Accomplished & Discoveries](#accomplished--discoveries)

---

## Quick Start

### 1. Install Dependencies
```bash
python3 -m pip install sympy scipy timeout-decorator numpy httpx ollama pymupdf
```

### 2. Run the Orchestrator
The main entry point is `vakyume.py`. It uses a project-centric directory structure.

#### Example: Building Models (Physics Kinematics)
To scrape a PDF, verify equations, and generate C++ code in one command:
```bash
python3 vakyume.py build projects/BuildingModels \
    --pdf projects/BuildingModels/BuildingModelsToDescribeOurWorld.pdf \
    --max-rounds 5 \
    --verbose
```

This runs the full pipeline:
1. **Scrape**: Extracts formulas from the PDF into `projects/BuildingModels/notes/`.
2. **Verify**: Shards the formulas and runs One-Odd-Out algebraic consistency checks.
3. **Repair**: Uses LLM-assisted repair for any inconsistent or broken solvers.
4. **Reconstruct**: Merges verified shards into `projects/BuildingModels/certified.py`.
5. **C++**: Transpiles the Python library into `projects/BuildingModels/vakyume.cpp`.

For more granular control, see the [PDF to C++ Pipeline](#pdf-to-c-pipeline) section.

## PDF to C++ Pipeline

End-to-end: take a textbook PDF, extract equations, generate verified solvers,
assemble a Python library, and compile to a C++ binary.

### 1. Create a Project and Extract Equations
```bash
mkdir -p projects/MyTextbook/notes

# Scrape all chapters from the PDF into notes files
python3 vakyume.py scrape textbook.pdf -o projects/MyTextbook/notes

# Or use the interactive wizard for guided extraction
python -m vakyume.pdf_scraper textbook.pdf --wizard
```

### 2. Generate Solvers, Verify, and Repair
```bash
python3 vakyume.py run projects/MyTextbook
```

This parses the notes, generates SymPy solvers for every variable, runs
One-Odd-Out cross-validation, and auto-repairs broken solvers via LLM.

### 3. Reconstruct the Python Library
```bash
python3 vakyume.py reconstruct projects/MyTextbook
```

Produces `certified.py` (flat single-file) and `py/` (modular package).

### 4. Transpile to C++ and Compile
```bash
python3 vakyume.py make-cpp projects/MyTextbook
```

Produces `cpp/` (headers + main) and `bin/` (compiled binary).

### All-in-One
```bash
# Resume an existing project through to C++
python3 vakyume.py build projects/VacuumTheory
```

## PDF Scraper Tips

The PDF scraper (`vakyume/pdf_scraper.py`) uses PyMuPDF for text extraction
and a local LLM for equation identification.  Here are tips for getting the
best results.

### Model Selection

- **phi3:latest** (3.8B) is the default.  Fast and usually sufficient for
  well-formatted textbooks.  Occasionally produces minor formatting artifacts.
- **llama3** (8B) tends to produce cleaner structured output with fewer
  artifacts, at the cost of slower extraction.  Recommended for difficult PDFs.
- Use `--wizard` or `-m MODEL` to switch models.

### What Gets Extracted (and What Doesn't)

The scraper targets **explicitly numbered equations** -- those marked with
labels like (3.1), (9.2), etc.  It deliberately ignores:
- Checkpoint questions and worked examples
- Unnumbered inline equations
- Derivative/integral notation (not algebraic -- vakyume cannot solve these)
- Vector equations (scalar/magnitude forms are extracted instead)

### Common LLM Artifacts and How They're Handled

| Artifact | Scraper Fix |
|----------|-------------|
| `^` instead of `**` | Regex replacement (`^2` → `**2`) |
| Unicode superscripts (`²`, `³`) | Replaced with `**2`, `**3` |
| Function notation (`x(t) =`) | Stripped to `x =` |
| `* *` spacing from PDF text | Collapsed to `**` |
| Trailing equation numbers `(3.1)` | Stripped from equation text |
| Derivative/integral notation | Rejected by validity filter |

### Improving Extraction Quality

1. **Run chapter-by-chapter first** (`-c 3`) to spot-check output before
   running the full book.
2. **Inspect the output notes files** -- the scraper writes one `.py` per
   chapter under the output directory, plus a `_extraction_summary.json`.
3. **Edit notes by hand** if the LLM gets an equation wrong.  The notes
   format is simple enough to fix manually.
4. **Re-run individual chapters** without affecting the rest: use `-c N` to
   overwrite just that chapter's output file.
5. **Chapters with no equations** (e.g., "Scientific Method", "Vectors",
   "Calculus") are auto-skipped.  See `SKIP_CHAPTERS` in `pdf_scraper.py`.

### Notes Format Reference

Each notes file contains one or more equation blocks:
```python
# 3-3 Position with constant acceleration
"""
x := position
x_0 := initial position
v_0 := initial velocity
a := acceleration
t := time
"""
x = x_0 + v_0 * t + 0.5 * a * t**2
```

Rules:
- Header: `# <chapter>-<eq_number> <title>`
- Docstring: triple-quoted, one `var := description` per line
- Equation: `lhs = rhs` using Python math syntax (`**`, `sqrt()`, `log()`)
- Variables should be lowercase with underscores for subscripts (`v_0`, `F_g`)

## Project Structure
Each project (e.g., `projects/VacuumTheory/`) contains:
- `notes/`: Python-like equation definitions (Input).
- `shards/`: Individual generated solvers (Intermediate).
- `reports/`: Verification and analysis logs.
- `repair_prompts/`: LLM interaction logs.
- `certified.py`: Flat single-file certified Python library.
- `py/`: Modular Python package (one file per chapter class + `__init__.py`).
- `cpp/`: Modular C++ source (`main.cpp` + per-chapter `.hpp` headers).
- `bin/vacuum_theory_binary`: Compiled C++ test binary.

## Setting Up a New Project

To start a new project from a textbook PDF:

1. **Create the Project Directory**:
   ```bash
   mkdir -p projects/MyNewProject/notes
   ```
2. **Extract Equations from the PDF**:
   ```bash
   # Quick: scrape all chapters at once
   python3 vakyume.py scrape path/to/textbook.pdf -o projects/MyNewProject/notes

   # Guided: interactive wizard with model and chapter selection
   python -m vakyume.pdf_scraper path/to/textbook.pdf --wizard
   ```
   See [PDF Scraper Tips](#pdf-scraper-tips) for guidance on getting clean output.

3. **Run the Full Pipeline**:
   ```bash
   python3 vakyume.py build projects/MyNewProject
   ```
   This will generate solvers, verify, reconstruct, and compile to C++.
   
4. **Manual Notes (Optional)**:
   You can also manually add `.py` files to the `notes/` directory. Use the following format:
   ```python
   # 1-1 Ideal Gas Law
   """
   p := pressure
   V := volume
   n := moles
   R := gas constant
   T := temperature
   """
   p * V = n * R * T
   ```

5. **Verify and Generate**:
   ```bash
   python3 vakyume.py build projects/MyNewProject
   ```

---

## Architecture

### Module Separation of Concerns

Each module in the `vakyume/` package has a single, well-defined responsibility.
No module reaches outside its domain -- the pipeline orchestrator is the only
component that wires them together.

| Module | Responsibility | Depends On |
|--------|----------------|------------|
| `config.py` | Constants, exceptions, import headers | (nothing) |
| `kwasak.py` | `@kwasak` missing-variable dispatch decorator | (nothing) |
| `parser.py` | Tokenize equations, build normal form, drive SymPy, write shard files | `config` |
| `scaffold.py` | Regex pattern matching on equation forms, generate algebraic derivation scaffolds | (nothing -- pure string transforms) |
| `verifier.py` | One-Odd-Out cross-validation engine (harmony checks, random trials, invariant detection) | `kwasak` |
| `llm.py` | Ollama/Phi-3 interface, prompt templates, code extraction from LLM responses | (nothing -- pure I/O) |
| `repair.py` | Orchestrate scaffold generation + LLM calls, post-process and write repaired shards | `config`, `llm`, `parser`, `scaffold` |
| `utils.py` | AST-based source code extraction helpers | (nothing) |
| `pipeline.py` | Pipeline orchestrator: wires parse -> verify -> repair -> certify | all of the above |
| `reconstruct.py` | Scan shard directories, extract functions via AST, assemble single importable module | `utils` |
| `cpp_gen.py` | Deterministic Python AST to C++ transpiler, LLM fallback, compile + test | `llm` |
| `pdf_scraper.py` | PyMuPDF text extraction, LLM-based equation identification, notes output, interactive wizard | `llm` |

### Module Dependency Graph

```
                        vakyume.py  (CLI entry point)
                             |
                             v
                      ┌─────────────┐
                      │ pipeline.py │  <-- orchestrator, the ONLY module that
                      │             │      imports from all others
                      └──────┬──────┘
                             |
            ┌────────────────┼────────────────────┐
            |                |                    |
            v                v                    v
     ┌────────────┐   ┌────────────┐       ┌──────────────┐
     │ parser.py  │   │ verifier.py│       │  repair.py   │
     │            │   │            │       │              │
     │ Tokenize   │   │ OOO cross- │       │ Scaffold +   │
     │ Normal form│   │ validation │       │ LLM prompt   │
     │ SymPy solve│   │ Harmony    │       │ Post-process │
     │ Write shard│   │ Invariant  │       │ Write shard  │
     └─────┬──────┘   └─────┬──────┘       └───┬──────┬───┘
           |                |                   |      |
           v                v                   v      v
      ┌─────────┐     ┌──────────┐      ┌──────────┐ ┌───────┐
      │config.py│     │kwasak.py │      │scaffold.py│ │llm.py │
      │         │     │          │      │           │ │       │
      │Constants│     │ @kwasak  │      │ Pattern   │ │Ollama │
      │Imports  │     │ dispatch │      │ matching  │ │Phi-3  │
      │Exceptions│    │          │      │ Algebraic │ │Prompts│
      └─────────┘     └──────────┘      │ scaffolds │ │Extract│
                                        └───────────┘ └───────┘

     ┌──────────────┐   ┌────────────┐
     │reconstruct.py│   │ cpp_gen.py │
     │              │   │            │
     │ Shard -> lib │   │ Python AST │
     │ AST extract  │   │ -> C++ src │
     │ @kwasak wire │   │ LLM fallbk │
     └──────────────┘   └────────────┘
```

**Key separation principles:**
- `scaffold.py` is a **pure function** module -- it takes an equation string and target variable, returns a scaffold string. No I/O, no LLM calls, no file access.
- `llm.py` is a **pure I/O** module -- it sends prompts to Ollama and returns raw text. No algebra, no file writes.
- `repair.py` is the **decision maker** that bridges scaffold and LLM -- it decides whether to bypass the LLM entirely (complete scaffold), prompt the LLM with a partial scaffold (fill-in-the-blank), or delegate fully to the LLM (no scaffold).
- `verifier.py` knows nothing about repair or scaffolds -- it only reports scores and mismatches.
- `pipeline.py` is the **only module that imports from all others**. It runs the verify-repair loop and manages state.

---

## Algorithm Deep Dive

Below is the end-to-end processing pipeline. Each stage maps directly to one
or more modules, maintaining the separation of concerns described above.

### Stage 1: Parse Equation Notes

**Module:** `parser.py` -- `shard_from_chapters()`

```
   notes/*.py                            Example input:
   ┌─────────────────────────────┐       # 1-7 ideal gas law
   │  For each .py file:         │       """
   │                             │       p := pressure
   │  1. Parse header comment    │       V := volume
   │     "# X-Y Title"           │       n := moles
   │     -> eqn_number = "1_7"   │       R := gas constant
   │                             │       T := temperature
   │  2. Parse docstring for     │       """
   │     variable definitions    │       p * V = n * R * T
   │     "var := description"    │
   │                             │
   │  3. Extract equation body   │
   │     "lhs = rhs"             │
   │                             │
   │  4. Tokenize: split on      │
   │     operators *()/-+        │
   │     -> tokens = [p,V,n,R,T] │
   │                             │
   │  5. Build normal form:      │
   │     (rhs) - (lhs) = 0       │
   │     "n*R*T - p*V"           │
   │     Replace ln() -> log()   │
   └─────────────────────────────┘
```

The tokenizer (`Solver.tokenize`) dilates operator characters (adding spaces
around `*()/-+`) then splits, filtering out math keywords like `ln` and `log`.
The normal form is always `rhs - lhs`, so that solving for any variable means
finding the root of a single expression equal to zero.

### Stage 2: Generate Solvers (SymPy)

**Module:** `parser.py` -- `Solver.get_solns_vanilla_nf()`

For each variable token V in the equation, attempt symbolic isolation:

```
   For each variable V in tokens:
   ┌───────────────────────────────────────┐
   │  sympy.solve(sympify(normal_form),    │
   │              Symbol(V))               │
   │                                       │
   │  Timeout: 1 second (MAX_COMP_TIME)    │
   │                                       │
   │  ┌─── Success ───┐  ┌─── Failure ──┐ │
   │  │                │  │              │ │
   │  │ Write shard:   │  │ Write stub:  │ │
   │  │ result = []    │  │ raise        │ │
   │  │ V = <soln_1>   │  │ Unsolved     │ │
   │  │ result.append  │  │ Exception    │ │
   │  │ V = <soln_2>   │  │ ("Pending    │ │
   │  │ result.append  │  │  LLM/Manual  │ │
   │  │ return result  │  │  Repair")    │ │
   │  └────────────────┘  └──────────────┘ │
   └───────────────────────────────────────┘
```

**Why the 1-second timeout?** SymPy's `solve()` can hang indefinitely on
transcendental or high-degree polynomial equations. The timeout ensures the
pipeline keeps moving -- unsolvable shards are deferred to the scaffold +
LLM repair path.

**Multi-valued solutions:** When SymPy returns multiple roots (e.g., a
quadratic has two solutions), all are stored in a `result` list. The
verifier's `are_similar()` checks if ANY element matches.

Each shard is written as an independent `.py` file:
```
shards/<FamilyName>/eqn_X_Y__variable.py
```

A **main shard** (`eqn_X_Y.py`) is also generated containing a class that
imports all variant solvers and exposes a `@kwasak`-decorated dispatch method.

### Stage 3: Scaffold Pattern Matching

**Module:** `scaffold.py` -- `generate_derivation_scaffold()`

This is the core algebraic intelligence of Vakyume. When SymPy cannot isolate
a variable (or produces a placeholder), the scaffold module inspects the
equation's **syntactic form** using regex pattern matching and generates a
partially or fully completed algebraic derivation.

The scaffold is NOT a SymPy alternative -- it is a **regex-driven algebraic
template engine** that recognizes common equation patterns from engineering
textbooks and produces step-by-step derivations with concrete variable names
and coefficients.

See [Scaffold Algorithm Detail](#scaffold-algorithm-detail) for the full
breakdown of each pattern handler.

### Stage 4: LLM Delegation

**Module:** `repair.py` -- `attempt_repair_shard()`, `llm.py`

LLM delegation is **never unconstrained**. The system always attempts to
provide algebraic structure before involving the LLM. The decision tree:

```
   Broken shard detected (UnsolvedException or failed verification)
                    |
                    v
   ┌────────────────────────────────┐
   │  Generate scaffold             │
   │  scaffold.py ->                │
   │  generate_derivation_scaffold()│
   └───────────────┬────────────────┘
                   |
          ┌────────┴────────┐
          |                 |
   Scaffold complete?    Scaffold partial/empty?
   (has assignment +     (ends with "V =" or
    optionally return)    returns None)
          |                 |
          v                 v
   ┌──────────────┐  ┌──────────────────────┐
   │ BYPASS LLM   │  │ Has partial scaffold?│
   │              │  └──────┬───────────────┘
   │ Write scaffold│        |
   │ directly as   │   ┌────┴────┐
   │ repaired shard│   |         |
   │              │  yes        no
   │ (Zero LLM    │   |         |
   │  cost)       │   v         v
   └──────────────┘ ┌─────────┐ ┌──────────────┐
                    │ FILL-IN │ │ FULL LLM     │
                    │         │ │ DELEGATION   │
                    │ Prompt: │ │              │
                    │ "Complete│ │ Prompt:      │
                    │  this   │ │ "Solve for V"│
                    │  func"  │ │ + equation   │
                    │ +scaffold│ │ + example    │
                    │ +example │ │ + header     │
                    │ +errors │ │ + errors     │
                    └─────────┘ └──────────────┘
```

See [LLM Delegation: Why Algebraic Step-by-Step?](#llm-delegation-why-algebraic-step-by-step)
for the detailed rationale and prompt architecture.

### Stage 5: One-Odd-Out Verification

**Module:** `verifier.py` -- `Verify` class

The OOO methodology cross-validates all solvers in an equation family:

```
   For an equation with N variables, and N solver variants:

   For each source variable V (3 trials):
   ┌────────────────────────────────────────────────────┐
   │  1. Generate random inputs for all OTHER N-1 vars  │
   │     (uniform random in [1.1, 10.0])                │
   │                                                    │
   │  2. Compute V using V's solver                     │
   │     source_value = solver_V(other_inputs)          │
   │                                                    │
   │  3. HARMONY CHECK: plug ALL N values into the      │
   │     original equation  lhs - rhs  and verify       │
   │     |residual| < 1e-4                              │
   │     (Uses a dynamically-generated subshard module)  │
   │                                                    │
   │  4. CROSS-CHECK: for every OTHER variable W:       │
   │     - Remove W from the full value set              │
   │     - Call solver_W with the remaining N-1 values   │
   │     - Check: are_similar(known_W, computed_W)?      │
   │       (tolerance 1e-6 on real and imaginary parts)  │
   │     - Multi-valued: ANY solution matching counts    │
   │                                                    │
   │  5. Score = number of successful cross-checks      │
   │     Perfect score = N (all variants agree)          │
   │     Broken = score < N                             │
   └────────────────────────────────────────────────────┘
```

**Invariant detection:** Before running trials, the verifier probes each
solver 3 times with different random inputs. If a solver always returns the
exact same value regardless of inputs, that variable has **cancelled out** of
the equation (e.g., appears as `p_i / p_i` in a ratio). These invariant
variants are excluded from cross-checking and given automatic perfect scores
so they don't poison the family's overall result.

**Harmony subshard:** The verifier dynamically generates a small Python module
(`harmony_check_<md5>.py`) containing a `check_harmony()` function that
computes `lhs - rhs` for the original equation. This is used as ground truth
independent of any solver.

### Stage 6: Repair Loop

**Module:** `pipeline.py` -- `run_pipeline()`

The pipeline runs a verify-repair loop with convergence:

```
   ┌──────────────────────────────────────────────────┐
   │  For round = 1 to max_rounds:                    │
   │                                                  │
   │    1. Verify ALL families (or only broken ones   │
   │       in --repair-only mode)                     │
   │                                                  │
   │    2. Analyze: categorize into                   │
   │       solved / inconsistent / failed             │
   │                                                  │
   │    3. If nothing broken -> STOP (all certified)  │
   │                                                  │
   │    4. Pick ONE family to repair                  │
   │       (alphabetically first inconsistent)        │
   │                                                  │
   │    5. For each broken variant in that family:    │
   │       attempt_repair_shard()                     │
   │       (up to 3 attempts per family per round)    │
   │                                                  │
   │    6. Re-verify the repaired family immediately  │
   │                                                  │
   │    7. If family passes -> mark as solved,        │
   │       continue to next round                     │
   │       If still broken -> retry with updated      │
   │       mismatch diagnostics                       │
   └──────────────────────────────────────────────────┘
```

**Focus strategy:** Only one family is repaired per round. This prevents
wasting LLM calls on families that might be fixed by a single repair, and
ensures the pipeline converges efficiently.

### Stage 7: Reconstruct Library

**Module:** `reconstruct.py` -- `reconstruct_from_shards()`

After all families are certified, assemble shards into importable Python
modules:

1. Scan all family directories in `shards/`
2. Group methods by class name
3. Extract function source using AST parsing (handles both standalone functions
   and class methods)
4. Emit `certified.py` (flat single-file) with proper class structure,
   `@kwasak` decorators, and import headers
5. Emit `py/` package: one module per chapter class (e.g., `basic.py`,
   `rotary.py`) plus an `__init__.py` that re-exports all classes
6. Clean up `# [.pyeqn]` tags to plain comments

### Stage 8: C++ Transpile

**Module:** `cpp_gen.py`

Optional stage that converts the Python library to C++:

1. **Deterministic AST transpiler** walks the Python AST and maps:
   - `log/sqrt/exp/pow` to `std::log/std::sqrt/std::exp/std::pow`
   - `**` operator to `std::pow(base, exp)`
   - `e/pi` to `M_E/M_PI`
   - `I` (imaginary) to `std::complex<double>(0.0, 1.0)`
   - `result.append()` to `result.push_back()`
   - `LambertW` to a shared Halley-iteration implementation (`lambertw.hpp`)
   - All variables typed as `double`; complex mode if `I` is detected

2. **Modular output**: each chapter class gets its own `.hpp` header under
   `cpp/`, with a `main.cpp` that `#include`s all headers and runs the test
   suite. Shared utilities (e.g., `lambertw.hpp`) are extracted to avoid
   redefinition errors.

3. **LLM fallback** for constructs the transpiler cannot handle (shells out
   to `ollama run phi3`)

4. **Compile and test**: `g++ -std=c++17 -I cpp/`, output binary goes to
   `bin/vacuum_theory_binary`

---

## Scaffold Algorithm Detail

The scaffold module (`scaffold.py`) is the algebraic brain of Vakyume. It
uses a **priority-ordered handler registry** where each handler inspects the
equation form using regex and returns either a complete scaffold string or
`None` to fall through to the next handler.

### Pattern Handler Registry

```python
_SCAFFOLD_HANDLERS = [
    _scaffold_exponent_target,                      # Branch 1: transcendental
    _scaffold_ratio_power_target_in_numerator,      # Branch 2: (target/C)**n
    _scaffold_ratio_power_target_in_denominator,    # Branch 3: (A/target)**n
    _scaffold_ratio_power_solve_numerator,          # Branch 4: (target/B)**n rematch
    _scaffold_bare_power,                           # Branch 5: target**n
    _scaffold_frac_exp_multi_diff,                  # Branch 6: frac exp, 2+ diffs
    _scaffold_frac_exp_single_occurrence,           # Branch 7: frac exp, 1 occurrence
]
```

After all handlers have been tried, two fallbacks apply:
- **F1** (single occurrence, no pattern match): stub `target =` for LLM
- **F2** (multiple occurrences, no pattern match): complete `brentq` solver

### Decision Flowchart

```
   generate_derivation_scaffold(pyeqn, target_var, header)
                         |
                         v
                  Parse "lhs = rhs"
                  Count occurrences of target_var in equation
                         |
                    ┌────┴────┐
                    |         |
                 count=0   count>=1
                    |         |
                    v         v
                 return ""   Try handlers in priority order
                             |
   ┌─────────────────────────┼─────────────────────────────┐
   |                         |                             |
   v                         v                             v
 Branch 1               Branch 2-5                   Branch 6-7
 Target in              Power/ratio                  Fractional
 exponent?              patterns                     exponent with
 **(...target...)        (target/C)**n               (target-X) terms
   |                     target**n                       |
   |                     (A/target)**n                   |
   v                         |                           v
 TRANSCENDENTAL              v                    Can linearize?
 -> brentq              Algebraic                    /       \
 numerical              inversion                  yes       no
 solver                 step-by-step              /           \
 (COMPLETE)             (COMPLETE)          Linear solve    brentq
                                           (COMPLETE)    (COMPLETE)
                                                |
                                           All handlers
                                           returned None?
                                                |
                                           ┌────┴────┐
                                           |         |
                                        count=1   count>=2
                                           |         |
                                           v         v
                                     Fallback F1  Fallback F2
                                     "target ="   brentq solver
                                     (PARTIAL)    (COMPLETE)
```

### Branch 1: Transcendental (Target in Exponent)

**Detection:** Regex `\*\*\s*\([^)]*target[^)]*\)` -- the target variable
appears inside an exponent expression like `X ** (a * target + b)`.

**Strategy:** These equations cannot be solved algebraically (they are
transcendental). The scaffold generates a **complete `brentq` numerical
solver** that:
1. Defines a residual function `_res(target_val)` = rhs(target_val) - lhs
2. Scans `x = 1.01, 1.02, ... 1000.0` looking for a sign change
3. Once a sign change `[lo, hi]` is found, calls `scipy.optimize.brentq`
4. Returns the root

**Example:** For equation `y = A * B ** (C * T)`, solving for `T`:
```python
def eqn_X_Y__T(self, y, A, B, C, **kwargs):
    # [.pyeqn] y = A * B ** (C * T)
    # T appears in the exponent -- use numerical solver
    from scipy.optimize import brentq
    def _res(T_val):
        return A * B ** (C * T_val) - y
    lo, hi = None, None
    prev = _res(1.01)
    for i in range(1, 100000):
        x = 1.01 + i * 0.01
        try:
            cur = _res(x)
        except Exception:
            continue
        if prev * cur < 0:
            lo, hi = x - 0.01, x
            break
        prev = cur
    if lo is None:
        raise UnsolvedException("No sign change found for T")
    T = brentq(_res, lo, hi)
    return [T]
```

This scaffold is **complete** -- the LLM is bypassed entirely.

### Branch 2: Ratio Power -- Target in Numerator

**Detection:** Regex `(target / C) ** n` where `C` is a constant/variable
and `n` is a numeric exponent.

**Strategy:** Algebraic inversion by clearing the exponent.

**Sub-case 2a:** With `- 1` modifier (e.g., `outer * ((target/C)**n - 1) = lhs`):
```
Step 1: (target / C)**n - 1 = lhs / outer
Step 2: (target / C)**n     = lhs / outer + 1
Step 3: target / C          = (lhs / outer + 1) ** (1/n)
Step 4: target              = C * (lhs / outer + 1) ** (1/n)
```

**Sub-case 2b:** Simple form (e.g., `outer * (target/C)**n = lhs`):
```
Step 1: (target / C)**n = lhs / outer
Step 2: target / C      = (lhs / outer) ** (1/n)
Step 3: target          = C * (lhs / outer) ** (1/n)
```

Both sub-cases produce **complete** scaffolds with the final assignment.

### Branch 3: Ratio Power -- Target in Denominator

**Detection:** Regex `(A / target) ** n`.

**Strategy:** Same exponent-clearing algebra as Branch 2, but inverted:
```
Step 1: (A / target)**n = ratio_expr
Step 2: A / target      = ratio_expr ** (1/n)
Step 3: target          = A / (ratio_expr ** (1/n))
```

Also handles the `- 1` modifier variant. Produces **complete** scaffolds.

### Branch 4: Ratio Power -- Re-match Numerator

**Detection:** Same regex as Branch 2 (`(target / B) ** n`), but this
handler fires when Branch 2 did not match (different structural context).

This covers cases where the equation structure prevents Branch 2's
`- 1` detection from triggering but the core ratio-power form is present.
Produces **complete** scaffolds.

### Branch 5: Bare Power

**Detection:** Regex `target ** n` (no surrounding ratio), with single
occurrence only (`total == 1`).

**Strategy:** Isolate the power term and take the inverse root.

**Sub-case 5a:** Power inside parentheses with additive terms
(e.g., `outer * (additive + coeff * target**n * cofactors) = lhs`):
```
Step 1: lhs / outer  = additive + coeff * target**n * cofactors
Step 2: (lhs/outer - additive) = coeff * target**n * cofactors
Step 3: target**n    = (lhs/outer - additive) / (coeff * cofactors)
Step 4: target       = ((lhs/outer - additive) / (coeff * cofactors)) ** (1/n)
```

**Sub-case 5b:** Simple `outer * target**n = lhs`:
```
Step 1: target**n = lhs / outer
Step 2: target    = (lhs / outer) ** (1/n)
```

Both produce **complete** scaffolds.

### Branch 6: Fractional Exponent with Multiple Difference Terms

**Detection:** Fractional exponent (`**0.X`) AND target appears in 2+
`(target - X)` difference terms within the exponent's base expression.

This is the pattern seen in thermodynamic correction equations like:
```
S_p = S_Th * ((P - p_s) * (460 + T_i) / ((P - p_c) * (460 + T_e))) ** 0.6
```

**Strategy:** Clear the fractional exponent by raising both sides to the
inverse power, then the equation becomes **linear** in the target variable
because each occurrence is in a `(target - constant)` factor.

```
Step 1: (lhs / outer) ** (1/0.6) = inner_expr     [clear exponent]
Step 2: R * den_other * (target - dd) = num_other * (target - nd)
                                                    [cross-multiply]
Step 3: target * (R*den_other - num_other) = R*den_other*dd - num_other*nd
                                                    [collect target]
Step 4: target = (R*den_other*dd - num_other*nd) / (R*den_other - num_other)
```

**Bare-target fallback:** If the target also appears as a bare variable (not
inside a difference term) in the denominator, linearization fails and the
handler falls back to a `brentq` numerical solver.

Produces **complete** scaffolds in the linear case.

### Branch 7: Fractional Exponent with Single Occurrence

**Detection:** Fractional exponent (`**0.X`) AND target appears exactly once,
inside a `(target - X)`, `(X - target)`, or `(C + target)` term.

**Strategy:** Clear the exponent (raise to inverse power to get `R`), then
isolate the difference/additive term:

For `(target - X)` in numerator:
```
R = (lhs / outer) ** (1/exp)     [clear exponent]
(target - X) = R * den / num_remaining
target = X + R * den / num_remaining
```

For `(X - target)` in numerator:
```
target = X - R * den / num_remaining
```

For `(target - X)` or `(X - target)` in denominator:
```
(target - X) = num / (R * den_remaining)
target = X + num / (R * den_remaining)
```

For `(C + target)` additive form:
```
(C + target) = factor
target = factor - C
```

All sub-patterns produce **complete** scaffolds.

### Fallbacks F1 and F2

When no pattern handler matches:

| Fallback | Condition | Output | LLM needed? |
|----------|-----------|--------|-------------|
| **F1** | Target appears exactly once, no pattern match | `target =` (incomplete stub) | Yes -- LLM must fill in the rearrangement |
| **F2** | Target appears 2+ times, no pattern match | Complete `brentq` numerical solver | No -- bypasses LLM |

**F2 rationale:** When a variable appears multiple times in an equation and
no algebraic pattern is recognized, symbolic isolation is unlikely to succeed.
A numerical solver is the pragmatic choice and is deterministically generated
without LLM involvement.

### Summary: Scaffold Completeness by Branch

| Branch | Form | Scaffold Type | LLM Needed? |
|--------|------|---------------|-------------|
| 1 | `X ** (...target...)` | Complete (brentq) | No |
| 2 | `(target / C) ** n` | Complete (algebraic) | No |
| 3 | `(A / target) ** n` | Complete (algebraic) | No |
| 4 | `(target / B) ** n` rematch | Complete (algebraic) | No |
| 5 | `target ** n` | Complete (algebraic) | No |
| 6 | `**0.X` + 2+ `(target-X)` | Complete (algebraic or brentq) | No |
| 7 | `**0.X` + 1 `(target-X)` | Complete (algebraic) | No |
| F1 | 1 occurrence, no pattern | Partial (`target =`) | Yes |
| F2 | 2+ occurrences, no pattern | Complete (brentq) | No |

**Key insight:** The majority of scaffolds are complete and bypass the LLM
entirely. The LLM is only needed for F1 -- single-occurrence variables in
equations whose form is not recognized by any pattern handler.

---

## LLM Delegation: Why Algebraic Step-by-Step?

### The Problem with Unconstrained LLM Solving

Naively asking an LLM "solve this equation for X" produces unreliable results:

1. **Sign errors:** LLMs frequently flip signs during multi-step algebra,
   especially when distributing negatives through parentheses.
2. **Missing solutions:** Polynomial equations have multiple roots; LLMs
   tend to return only one.
3. **Domain errors:** LLMs may produce expressions with `log(-x)` or
   `sqrt(-y)` without using `cmath` for complex-safe math.
4. **Hallucinated simplifications:** LLMs may "simplify" `(a-b)/(c-d)` to
   `a/c - b/d`, which is algebraically wrong.
5. **SymPy dependency:** Despite being told not to, LLMs often import
   `sympy.solve()` in their output, defeating the purpose of generating a
   standalone solver.

### Scaffold as Algebraic Guardrail

The scaffold approach **constrains the LLM's creative freedom** by
pre-computing as much of the algebra as possible:

```
WITHOUT scaffold (unreliable):
   "Solve S_p = S_Th * ((P-p_s)*(460+T_i)/((P-p_c)*(460+T_e)))**0.6 for P"
   -> LLM must: parse equation, identify pattern, derive algebra, write code
   -> Error rate: HIGH (5+ algebraic steps, each can go wrong)

WITH scaffold (reliable):
   "Complete this function:
    def eqn_8_15__P(self, S_p, S_Th, p_s, T_i, p_c, T_e, **kwargs):
        # Step 1: (S_p/S_Th) ** (1/0.6) = ((P-p_s)*(460+T_i)/((P-p_c)*(460+T_e)))
        R = (S_p / (S_Th)) ** (1.666667)
        # Step 2: R * (460+T_e) * (P - p_c) = (460+T_i) * (P - p_s)
        # Step 3: P * (R*(460+T_e) - (460+T_i)) = R*(460+T_e)*p_c - (460+T_i)*p_s
        # Step 4: P = (R*(460+T_e)*p_c - (460+T_i)*p_s) / (R*(460+T_e) - (460+T_i))
        P = (R * (460+T_e) * p_c - (460+T_i) * p_s) / ..."
   -> LLM only needs to: copy the expression from Step 4
   -> Error rate: LOW (expression is pre-derived)
```

In practice, when the scaffold is complete, the LLM is **bypassed entirely**
-- the scaffold itself becomes the repaired shard. The LLM is only invoked
when the scaffold is partial (Fallback F1).

### Prompt Architecture

When the LLM must be involved, prompts are carefully structured:

**System prompt constraints (always enforced):**
- "Output ONLY code"
- "No markdown"
- "Do NOT use sympy"
- "Use cmath for complex-safe math"
- "Return as `return [value]`"
- "If the variable is in the exponent (transcendental), use `scipy.optimize.brentq`
  with a sign-change scan"

**User prompt structure (three tiers):**

| Tier | When | Prompt Contents |
|------|------|-----------------|
| **Scaffold fill-in** | Partial scaffold exists | Working sibling shard (example) + scaffold with step-by-step algebra + mismatch diagnostics |
| **Full delegation** | No scaffold | Equation string + working sibling shard + function header + mismatch diagnostics |
| **Repair with context** | Previous attempt failed | Same as above + specific input/output/expected triples from failed verification trials |

**Mismatch diagnostics** are critical for iterative repair. When a previous
LLM attempt produced incorrect code, the next prompt includes up to 2 concrete
examples:
```
Previous attempt failed with these errors:
  Inputs: {'S_Th': 3.14, 'p_s': 2.71, ...}
  Got: [42.0]
  Expected P=7.38, got [42.0]
```

This gives the LLM concrete numerical feedback on where its previous answer
was wrong, dramatically improving repair success rates.

### Bypass vs. Fill-in vs. Full Delegation

The three modes of LLM interaction, in order of preference:

```
                     ┌────────────────────────────────────────────────┐
                     │           SCAFFOLD COMPLETENESS               │
                     │                                                │
                     │  COMPLETE            PARTIAL         NONE      │
                     │  (has assignment     (ends with     (handler   │
                     │   + return)          "target =")    returned   │
                     │                                     None)      │
                     │      |                  |              |       │
                     │      v                  v              v       │
                     │  ┌─────────┐     ┌───────────┐  ┌──────────┐  │
                     │  │ BYPASS  │     │ FILL-IN   │  │ FULL     │  │
                     │  │         │     │           │  │ DELEGATE │  │
                     │  │ Write   │     │ LLM sees  │  │          │  │
                     │  │ scaffold│     │ steps 1-3 │  │ LLM sees │  │
                     │  │ directly│     │ and fills │  │ equation │  │
                     │  │         │     │ step 4    │  │ + example│  │
                     │  │ LLM     │     │           │  │ + header │  │
                     │  │ cost: 0 │     │ LLM       │  │          │  │
                     │  │         │     │ cost: LOW │  │ LLM      │  │
                     │  │ Error   │     │           │  │ cost:HIGH│  │
                     │  │ rate: 0 │     │ Error     │  │          │  │
                     │  │         │     │ rate: LOW │  │ Error    │  │
                     │  └─────────┘     └───────────┘  │ rate:MED │  │
                     │                                  └──────────┘  │
                     └────────────────────────────────────────────────┘
```

**Post-processing (all LLM outputs):**
1. Strip markdown fences
2. Find target `def` by name (fallback: first `def`)
3. Normalize function signature (`self, params, **kwargs`)
4. Ensure `return [value]` format
5. Validate syntax with `ast.parse()`
6. Preserve `# [.pyeqn]` comment
7. Prepend standard import header

If `ast.parse()` fails, the repair attempt is logged and the shard remains
broken for the next round.

---

## The `@kwasak` Decorator

The `@kwasak` decorator enables kwarg-solver dispatch.  Given a function
`eqn_X_Y` with variant solvers `eqn_X_Y__V1`, `eqn_X_Y__V2`, etc., calling
with one keyword argument omitted automatically dispatches to the correct
variant solver:

```python
@kwasak
def ideal_gas(P=None, V=None, n=None, R=None, T=None):
    """P * V = n * R * T"""
    pass

# Omit P -> dispatches to ideal_gas__P(V, n, R, T)
ideal_gas(V=1.0, n=1.0, R=8.314, T=300)  # returns P
```

---

## Resume & Testing
The orchestrator naturally supports resuming. If a shard file already exists in the project's `shards/` directory, it is tested but not re-scraped unless you use the `--overwrite` flag.

## Examples

### Vacuum Theory
Equations extracted from *Process Vacuum System Design and Operation* (1986).
```bash
python3 vakyume.py build projects/VacuumTheory
```

### Building Models
Physics kinematics equations extracted from *Building Models To Describe Our World*.
```bash
python3 vakyume.py build projects/BuildingModels
```

---

## Accomplished & Discoveries

- **One-Odd-Out (OOO) Methodology**: We've proven that verifying $f(x, y) \to z$ and $f(x, z) \to y$ consistency is highly effective at catching algebraic errors in LLM-generated or SymPy-isolated functions.
- **Scaffold-First Architecture**: The majority of "LLM-repaired" shards are actually solved by the scaffold pattern matcher with zero LLM calls. The LLM is a fallback, not the primary solver.
- **Project-Centric Structure**: Successfully transitioned to a formal Python package (`vakyume/`) with a project-based directory structure (e.g., `projects/VacuumTheory/`).
- **Resumability**: The pipeline automatically skips existing shards, allowing for long-running verification or repair tasks to be resumed.
- **C++ Library Generation**: A deterministic AST transpiler converts Python solvers to C++ (`log`, `sqrt`, `exp`, `pow`, arithmetic). Falls back to local LLMs (Phi-3 via Ollama) for constructs the transpiler cannot handle. Enforces `double` precision throughout. Output is modular: per-chapter `.hpp` headers under `cpp/`, compiled binary under `bin/`.
- **Modular Output Structure**: Both Python and C++ outputs are organized into modular packages. Python: `certified.py` (flat) + `py/` (per-class modules with `__init__.py`). C++: `cpp/` (per-class `.hpp` headers + `main.cpp`) + `bin/vacuum_theory_binary`.
- **Parallelization**: Parallelizing Ollama calls (using 2-4 workers) significantly speeds up the verification of large equation sets (~500+ functions).
