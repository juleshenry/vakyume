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
- **Sympy Solvers**: Initial algebraic isolation of variables.
- **Numerical Fallbacks**: Uses `scipy.optimize.newton` or `fsolve` for transcendental equations where algebraic isolation is impossible.
- **One-Odd-Out Verification**: Random "fuzzed" inputs are used to verify that moving between any two variables in an equation family returns to the source of truth.

---

# Kwarg Solver Decorator
(Existing documentation...)

As long as one keyword argument is not given, its value is calculated

E.g. 
```
@kwarg_solver
def einstein(...):
	"""e = m * c ** 2"""
	return ...

einstein(e = 1000) # returns m, (1000 / 8.98755179 e16), ~1.11265 e -14
einstein(m = 1000) # returns e, 1000 * 8.98755179 e16, ~8.98755179 e19
```

Transcription Phase
- [x] Transcription of Formulae
- [x] Develop universal format for these notes
- [x] Confirm adherence to strict format
- [x] Filling in physics constants

Implementation Phase
- [x] Use sympy to arbitrary solve all equations for one variable 🐍📐🎊
- [x] Implement solver-decorator

Integration Phase
- [x] Separate the python library that converts equations to dynamic classes
- [x] Use c++-imports to speed up python library
- [x] Make class decorator to ensure invariants on `kwarg_solver`
- [x] Metacode Motherlode, the class that does it all

Namely, kwarg_solver is a decorator that requires for x kwargs default to zero, there are x auxiliary functions of the form `vanilla__a, vanilla__b, etc.., vanilla__z`. Implement this as a class decorator. 

Finally, metacode motherlode : 

UniversalSolver({"eqn_name" : "name of method", "eqn": "normal form of equation"})
-> 
UniversalSolver({"eqn_name" : "Einstein", "eqn": "0 = m  * 8.98755179e16 - e"})
->
`
class Einstein:
    @kwarg_solver
    def einstein(s, e: float = None, m: float = None, **kwargs):
        return  # decorator skips return
    def einstein__m(s, e: float):
        return e / 8.98755179e16
    def einstein__e(s, m: float) -> float:
        return m * 8.98755179e16
`

You get returned a generated class code that solves the equation for parameters
Einstein().einstein(e = 1000)# Instantly returns ~1.11265 e -14

# Structure:

- 1. Original notes, Python equations
- 2. Python class creation
- 3. Testing of all existing methods, filling in unresolved equation-stubs or deferring (need help here!)
- 4. Conversion to C++
- 5. Testing in C++ all calls
- 6. Speed tests over compute space

# Process of creation:
- 1. Sympy attempt
- 2. LLM failover
- 3. Can say, with LLM, flow is: MAKE, EVAL(CODE), TEST_CODE 
- 4. a loop could be created to run in perpetuity until failovers are all working.
- 5. could even make another function to reiterate over a particular failing test for efficiency

# Final Goal:
python3 vakyume.py build your_project
[x] loading formulae
[x] solving via sympy
[x] now leveraging LLM's to solve remaining equations...
[x] spitting out Python libraries
[x] spitting out C++ library derived from the Python...
