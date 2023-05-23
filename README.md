# Vakyume
# A Python Library for Vacuum Design

Inspired by the 1986 edition of *Process Vacuum System Design and Operation* by James L. Ryans and Daniel L. Roper, this project offers functions related to engineering vacuums. 

We revisit the history of low pressure devices, recreating the example calculations programatically. Furthermore, a supplementary decorator is developed to calculate arbitrary missing values for a given formula. 

In a nutshell, we demonstrate how to take a textbook's equations and turn it into a complete equation solver for all permutations. 

Later, we convert the Python code to C++ for optimization.


# Kwarg Solver Decorator

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
- [ ] Inflation adjusting money calculations in chapter 8
- [ ] Codify Test Examples from Book

Implementation Phase
- [x] Use sympy to arbitrary solve all equations for one variable 🐍📐🎊
- [x] Implement solver-decorator
- [ ] Confirm worked examples in book via test cases
- [ ] Recreate some graphs in book
- [ ] Create arbitrary solver demo that scrolls through various calculations, for fun!

Integration Phase
- [ ] Separate the python library that converts equations to dynamic classes
- [ ] Use c++-imports to speed up python library
- [ ] Make class decorator to ensure invariants on `kwarg_solver`
- [ ] Metacode Motherlode, the class that does it all

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

Structure:

- 1. Original notes, Python equations
- 2. Python class creation
- 3. Testing of all existing methods, filling in unresolved equation-stubs or deferring
- 4. Conversion to C++
- 5. Testing in C++ all calls
- 6. Speed tests over compute space

