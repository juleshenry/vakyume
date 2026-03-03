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

## Project Structure
Each project (e.g., `projects/VacuumTheory/`) contains:
- `notes/`: Python-like equation definitions (Input).
- `shards/`: Individual generated solvers (Intermediate).
- `reports/`: Verification and analysis logs.
- `repair_prompts/`: LLM interaction logs.
- `vakyume_certified.py`: The finalized, verified Python library.
- `vakyume.cpp`: The generated C++ library (if `--cpp` is used).

## Setting Up a New Project

To start a new project from a textbook PDF:

1. **Create the Project Directory**:
   ```bash
   mkdir -p projects/MyNewProject/notes
   ```
2. **Run Vakyume with the PDF**:
   ```bash
   python3 vakyume.py projects/MyNewProject --pdf path/to/textbook.pdf
   ```
   This will attempt to extract formulas from the PDF into `projects/MyNewProject/notes/extracted_notes.py`. 
   
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
   python3 vakyume.py projects/MyNewProject --cpp
   ```

## Resume & Testing
The orchestrator naturally supports resuming. If a shard file already exists in the project's `shards/` directory, it is tested but not re-scraped unless you use the `--overwrite` flag.

## Vacuum Theory Example
The default `chapters/` directory contains equations extracted from the 1986 edition of *Process Vacuum System Design and Operation*. These files (`01_vacuum_theory.py`, etc.) serve as a primary example of how the pipeline processes a complex engineering textbook.

## Accomplished & Discoveries

- **One-Odd-Out (OOO) Methodology**: We've proven that verifying $f(x, y) \to z$ and $f(x, z) \to y$ consistency is highly effective at catching algebraic errors in LLM-generated or SymPy-isolated functions.
- **Project-Centric Structure**: Successfully transitioned to a formal Python package (`vakyume/`) with a project-based directory structure (e.g., `projects/VacuumTheory/`).
- **Resumability**: The pipeline automatically skips existing shards, allowing for long-running verification or repair tasks to be resumed.
- **C++ Library Generation**: High-performance C++ code is generated using local LLMs (Phi-3 via Ollama), strictly enforcing `double` precision and `std::complex` for imaginary numbers.
- **Parallelization**: Parallelizing Ollama calls (using 2-4 workers) significantly speeds up the verification of large equation sets (~500+ functions).

## Examples

### Vacuum Theory
Equations extracted from *Process Vacuum System Design and Operation* (1986).
```bash
python3 vakyume.py projects/VacuumTheory --cpp
```

### Building Models
Physics kinematics equations extracted from *Building Models To Describe Our World*.
```bash
python3 vakyume.py projects/BuildingModels --cpp
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
python3 vakyume.py your_text_book.pdf 
[x] loading formulae
[x] solving via sympy
[x]  now leveraging LLM's to solve remaining equations...
[ ] spitting out Python libraries
[ ]  spitting out C++ library derived from the Python...
