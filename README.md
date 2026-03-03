# Vakyume
# A Python Library for Vacuum Design

Inspired by the 1986 edition of *Process Vacuum System Design and Operation* by James L. Ryans and Daniel L. Roper, this project offers functions related to engineering vacuums.

We utilize a "One-Odd-Out" (OOO) verification methodology: for any equation family (e.g., $PV=nRT$), if solvers exist for every variable ($P, V, n, R, T$), they must be mathematically consistent. If one solver produces a result that breaks the equality when plugged into the others, it is flagged for LLM-assisted repair.

## Quick Start

### 1. Install Dependencies
```bash
python3 -m pip install sympy scipy timeout-decorator numpy httpx ollama
```

### 2. Run the Orchestrator
The main entry point is `vakyume.py`. It uses a project-centric directory structure:
```bash
python3 vakyume.py projects/VacuumTheory
```

## CLI Reference

Vakyume exposes four subcommands.  Run `python3 vakyume.py --help` for an overview.

### `run` -- Shard, Verify & Repair

Run the core pipeline: generate solver shards from equation notes, verify them
with One-Odd-Out, and auto-repair broken solvers via LLM.

```
python3 vakyume.py run <project> [options]
```

| Flag | Description |
|------|-------------|
| `--pdf PATH` | Extract equations from a textbook PDF (OCR stage) |
| `--cpp` | Generate C++ library after certification |
| `--max-rounds N` | Maximum repair rounds (default 10) |
| `--overwrite` | Wipe and regenerate all shards |
| `--verbose` | Show detailed verification output |
| `--repair-only` | Skip shard generation; only repair broken families |

```bash
# Basic run
python3 vakyume.py run projects/VacuumTheory

# With C++ generation and verbose output
python3 vakyume.py run projects/VacuumTheory --cpp --verbose

# Resume and only repair broken equations
python3 vakyume.py run projects/VacuumTheory --repair-only --max-rounds 5
```

### `reconstruct` -- Rebuild Importable Python Library

Reassemble individual solver shards into a single, importable Python module
(`vakyume_lib.py`).

```
python3 vakyume.py reconstruct <project> [options]
```

| Flag | Description |
|------|-------------|
| `-o, --output PATH` | Custom output file (default `<project>/vakyume_lib.py`) |
| `--stdout` | Print to stdout instead of writing a file |

```bash
python3 vakyume.py reconstruct projects/VacuumTheory
python3 vakyume.py reconstruct projects/VacuumTheory -o mylib.py
python3 vakyume.py reconstruct projects/VacuumTheory --stdout | head
```

### `make-cpp` -- Generate & Compile C++ Library

Convert the reconstructed Python library (or any certified file) to C++,
compile it with `g++`, and run the built-in test suite.

The converter uses a **deterministic AST transpiler** for standard arithmetic
and math functions (`log`, `sqrt`, `exp`, `pow`), falling back to a local LLM
(Phi-3 via Ollama) only for constructs the transpiler cannot handle.

```
python3 vakyume.py make-cpp <project> [options]
```

| Flag | Description |
|------|-------------|
| `-i, --input FILE` | Python source to convert (default: `vakyume_lib.py`, then `vakyume_certified.py`) |

```bash
# Convert the reconstructed library (default)
python3 vakyume.py make-cpp projects/VacuumTheory

# Convert a specific file
python3 vakyume.py make-cpp projects/VacuumTheory -i vakyume_certified.py
```

Output:
- `<project>/vakyume.cpp` -- generated C++ source
- `<project>/vakyume_test` -- compiled test binary

### `build` -- Full Pipeline (All-in-One)

Run the complete pipeline in a single command: scrape equations, verify &
repair, reconstruct the Python library, and generate the C++ library.

```
python3 vakyume.py build <project> [options]
```

| Flag | Description |
|------|-------------|
| `--pdf PATH` | Extract equations from a textbook PDF |
| `--max-rounds N` | Maximum repair rounds (default 10) |
| `--overwrite` | Wipe and regenerate all shards |
| `--verbose` | Show detailed verification output |

```bash
# Full build from scratch
python3 vakyume.py build projects/MyProject --overwrite

# Build from a PDF
python3 vakyume.py build projects/MyProject --pdf textbook.pdf

# Resume an existing project through to C++
python3 vakyume.py build projects/VacuumTheory
```

Pipeline stages:
1. **Scraping** -- generate solver shards from `notes/*.py`
2. **Verification** -- One-Odd-Out cross-checks all solvers
3. **Repair** -- LLM-assisted fix of broken solvers (up to `max-rounds`)
4. **Reconstruction** -- assemble shards into `vakyume_lib.py`
5. **C++ generation** -- transpile to C++, compile, and test

## Project Structure
Each project (e.g., `projects/VacuumTheory/`) contains:
- `notes/`: Python-like equation definitions (Input).
- `shards/`: Individual generated solvers (Intermediate).
- `reports/`: Verification and analysis logs.
- `repair_prompts/`: LLM interaction logs.
- `vakyume_certified.py`: The finalized, verified Python library.
- `vakyume_lib.py`: Reconstructed single-file Python library.
- `vakyume.cpp`: The generated C++ library.
- `vakyume_test`: Compiled C++ test binary.

## Setting Up a New Project

To start a new project from a textbook PDF:

1. **Create the Project Directory**:
   ```bash
   mkdir -p projects/MyNewProject/notes
   ```
2. **Run Vakyume with the PDF**:
   ```bash
   python3 vakyume.py build projects/MyNewProject --pdf path/to/textbook.pdf
   ```
   This will extract formulas, generate solvers, verify, reconstruct, and compile to C++.
   
   *Note: For best results, ensure `pdftotext` is installed on your system.*

3. **Manual Notes (Optional)**:
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

4. **Verify and Generate**:
   ```bash
   python3 vakyume.py build projects/MyNewProject
   ```

## Resume & Testing
The orchestrator naturally supports resuming. If a shard file already exists in the project's `shards/` directory, it is tested but not re-scraped unless you use the `--overwrite` flag.

## Vacuum Theory Example
The default `chapters/` directory contains equations extracted from the 1986 edition of *Process Vacuum System Design and Operation*. These files (`01_vacuum_theory.py`, etc.) serve as a primary example of how the pipeline processes a complex engineering textbook.

## Accomplished & Discoveries

- **One-Odd-Out (OOO) Methodology**: We've proven that verifying $f(x, y) \to z$ and $f(x, z) \to y$ consistency is highly effective at catching algebraic errors in LLM-generated or SymPy-isolated functions.
- **Project-Centric Structure**: Successfully transitioned to a formal Python package (`vakyume/`) with a project-based directory structure (e.g., `projects/VacuumTheory/`).
- **Resumability**: The pipeline automatically skips existing shards, allowing for long-running verification or repair tasks to be resumed.
- **C++ Library Generation**: A deterministic AST transpiler converts Python solvers to C++ (`log`, `sqrt`, `exp`, `pow`, arithmetic). Falls back to local LLMs (Phi-3 via Ollama) for constructs the transpiler cannot handle. Enforces `double` precision throughout.
- **Parallelization**: Parallelizing Ollama calls (using 2-4 workers) significantly speeds up the verification of large equation sets (~500+ functions).

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


## Architecture

### Core Components

| Module | Purpose |
|--------|---------|
| `vakyume.py` | CLI entry point (`run`, `reconstruct`, `make-cpp`, `build`) |
| `vakyume/master.py` | Pipeline orchestrator: scraping, verification, repair |
| `vakyume/verifier.py` | One-Odd-Out cross-validation engine |
| `vakyume/llm.py` | LLM interface (Phi-3 / Ollama) |
| `vakyume/kwasak.py` | `@kwasak` decorator for kwarg-solver dispatch |
| `vakyume/reconstruct.py` | Shard reassembly into a single importable library |
| `vakyume/cpp_gen.py` | AST-based Python-to-C++ transpiler |
| `vakyume/config.py` | Shared constants and import headers |

### Algorithmic Flowchart

Below is the end-to-end processing pipeline that Vakyume executes when you
run `python3 vakyume.py build <project>`.

```
                         notes/*.py
                    (Python-like equations)
                             |
                             v
               ┌─────────────────────────┐
               │  1. PARSE EQUATION FILE │
               │                         │
               │  For each notes/*.py:   │
               │  - Parse header comment │
               │    "# X-Y Title"        │
               │  - Parse docstring for  │
               │    variable definitions │
               │    "var := description"  │
               │  - Extract equation     │
               │    body  "lhs = rhs"    │
               └────────────┬────────────┘
                            |
                            v
               ┌─────────────────────────┐
               │  2. GENERATE SOLVERS    │
               │                         │
               │  For each variable V    │
               │  in the equation:       │
               │                         │
               │  ┌───────────────────┐  │
               │  │ 2a. SymPy Solve   │  │
               │  │  sympy.solve(     │  │
               │  │    Eq(lhs, rhs),  │  │
               │  │    Symbol(V))     │  │
               │  └────────┬──────────┘  │
               │           |             │
               │       success?          │
               │        /     \          │
               │      yes      no        │
               │       |        |        │
               │       |        v        │
               │       │  ┌───────────┐  │
               │       │  │ 2b. Build │  │
               │       │  │ Scaffold  │  │
               │       │  │ (pattern  │  │
               │       │  │ matching) │  │
               │       │  └─────┬─────┘  │
               │       │        |        │
               │       │        v        │
               │       │  ┌───────────┐  │
               │       │  │ 2c. LLM   │  │
               │       │  │ Fill-in   │  │
               │       │  │ (Phi-3)   │  │
               │       │  └─────┬─────┘  │
               │       │        |        │
               │       v        v        │
               │     Write shard file:   │
               │     shards/<Family>/    │
               │       eqn_X_Y__V.py    │
               └────────────┬────────────┘
                            |
                            v
               ┌─────────────────────────┐
               │  3. VERIFY (OOO)        │
               │                         │
               │  For each equation      │
               │  family with N vars:    │
               │                         │
               │  For each variable V:   │
               │    Generate random      │
               │    inputs for other     │
               │    N-1 variables        │
               │    ┌─────────────────┐  │
               │    │ Compute V via   │  │
               │    │ its solver      │  │
               │    └────────┬────────┘  │
               │             |           │
               │    ┌─────────────────┐  │
               │    │ Plug V + inputs │  │
               │    │ into ALL other  │  │
               │    │ solvers & check │  │
               │    │ round-trip      │  │
               │    └────────┬────────┘  │
               │             |           │
               │         pass/fail?      │
               │          /      \       │
               │        pass    FAIL     │
               │         |       |       │
               │         |       v       │
               │         | Flag variant  │
               │         | as broken     │
               └─────────┬──┬────────────┘
                         |  |
                    pass |  | broken
                         |  v
                         |  ┌─────────────────────┐
                         |  │  4. REPAIR (LLM)    │
                         |  │                     │
                         |  │  Show LLM:          │
                         |  │  - The equation     │
                         |  │  - A working solver │
                         |  │    (as example)     │
                         |  │  - The error info   │
                         |  │                     │
                         |  │  LLM writes fixed   │
                         |  │  solver code        │
                         |  │                     │
                         |  │  Repeat verify +    │
                         |  │  repair up to       │
                         |  │  --max-rounds N     │
                         |  └──────────┬──────────┘
                         |             |
                         v             v
               ┌─────────────────────────┐
               │  5. RECONSTRUCT         │
               │                         │
               │  Assemble all certified │
               │  shards into one file:  │
               │  vakyume_lib.py         │
               │                         │
               │  Each equation family   │
               │  becomes a @kwasak      │
               │  decorated method with  │
               │  variant sub-functions: │
               │                         │
               │  @kwasak                │
               │  def eqn_X_Y(...):      │
               │      ...                │
               │  def eqn_X_Y__V1(...):  │
               │      return <solution>  │
               │  def eqn_X_Y__V2(...):  │
               │      return <solution>  │
               └────────────┬────────────┘
                            |
                            v
               ┌─────────────────────────┐
               │  6. C++ TRANSPILE       │
               │  (optional --cpp)       │
               │                         │
               │  AST-based converter:   │
               │  - Parse Python AST     │
               │  - Map math.log -> log  │
               │  - Map ** -> pow()      │
               │  - Enforce double types │
               │  - Fallback to LLM for  │
               │    complex constructs   │
               │                         │
               │  Output: vakyume.cpp    │
               │  Compile: g++ -> test   │
               └─────────────────────────┘
```

### Scaffold Pattern Matching (Step 2b Detail)

When SymPy cannot isolate a variable, the scaffold builder
(`_generate_derivation_scaffold`) tries pattern handlers in order:

| Priority | Pattern | Strategy |
|----------|---------|----------|
| 1 | Target in exponent: `X ** (... target ...)` | Numerical `brentq` solver |
| 2 | `(target / C) ** n` | Invert: `C * (lhs/outer) ** (1/n)` |
| 3 | `(A / target) ** n` | Invert: `A / (ratio) ** (1/n)` |
| 4 | `(target / B) ** n` (re-match) | Invert for numerator |
| 5 | `target ** n` (bare power) | Isolate and apply `(1/n)` root |
| 6 | Fractional exponent, 2+ `(target - X)` terms | Clear exponent, cross-multiply, solve linear |
| 7 | Fractional exponent, single `(target - X)` | Clear exponent, isolate difference |
| F1 | Single occurrence (no pattern match) | Stub for LLM to fill in |
| F2 | Multiple occurrences (no pattern match) | Numerical `brentq` solver |

Each handler either returns a complete scaffold string or `None` to fall
through to the next handler.

### The `@kwasak` Decorator

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
