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
The main entry point is `vakyume_master.py`. It handles the entire lifecycle:
```bash
python3 vakyume_master.py
```
This script will:
1. **Verify Shards**: Check all equations in `shards/` for mathematical consistency using the `tru.py` engine.
2. **Analyze Results**: Identify "solved" vs "inconsistent" equations, saved to `reports/analysis.json`.
3. **Generate Repair Prompts**: Create targeted prompts in `repair_prompts/` for any inconsistent equations.
4. **Assemble Certified Library**: Create `vakyume_certified.py` containing only the 100% verified equations.

## Workflow for LLM-Assisted Repair

If you have inconsistent shards (check the output of the master run), follow this process to fix them:

1. **Locate Prompt**: Find the corresponding file in `repair_prompts/repair_<shard_name>.txt`.
2. **Run Repair**: Feed the content of the prompt to an LLM (e.g., via Ollama or ChatGPT).
   ```bash
   # Example using Ollama
   PROMPT=$(cat repair_prompts/repair_SelectingPump_eqn_8_06.py.txt)
   ollama run phi3 "$PROMPT"
   ```
3. **Apply Fix**: Paste the corrected methods back into the shard file in `shards/`.
4. **Re-Verify**: Run `python3 vakyume_master.py` again. The repaired equation will now be promoted to the certified library if it passes consistency checks.

## Key Files
- `vakyume_master.py`: Root orchestrator for verification and certification.
- `tru.py`: The mathematical n-way consistency verification engine (The "Fuzzer").
- `shards/`: Individual Python files for each equation family.
- `vakyume_certified.py`: The production-ready, 100% verified subset of the library.
- `reports/analysis.json`: Detailed report on the consistency of every equation.

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
