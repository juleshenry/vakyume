# Equation Permuter
The purpose of this library is to convert formulae to a python library that, if possible, solves for every variable, given a single equation.

The equations must be saved in a `.pyeqn` file. A .pyeqn file is markdown for formulae that uses Pythonic syntax to describe an equation.

You can validate your `.pyeqn` file with the pyeqn_validator.py script

# Installation

To use this pyproject.toml file
Make sure you have Poetry installed.
Run poetry install to install the specified dependencies.

# Usage

Equation Permutation: 

```python3 go.py --infile "tests/linear_eqn.pyeqn" --outfile "permuted_equations.py"```

Pyeqn Validator:
```python3 pyeqn_validator.py --infile "tests/linear_eqn.pyeqn"```

