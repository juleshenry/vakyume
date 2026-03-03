class IntractableSolution(Exception):
    """Exception raised when only numerical methods are allowed."""

    pass


class OllamaOffline(Exception):
    """Exception raised when the Ollama is offline."""

    pass


class UnsolvedException(Exception):
    """Exception raised when equation solving is unavailable."""

    pass


# Code generation constants
TAB = " " * 4
MAX_COMP_TIME_SECONDS = 1
COOLDOWN_SECONDS = 2
FUNKTORZ = "*()/-+"

# ── Shared import headers for generated code ────────────────────────────────
# SHARD_IMPORT_HEADER:  used for individual solver shard files (no kwasak).
# LIBRARY_IMPORT_HEADER: used for assembled/reconstructed libraries (with kwasak).

SHARD_IMPORT_HEADER = (
    "from cmath import log, sqrt, exp\n"
    "from math import e, pi\n"
    "from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp\n"
    "from scipy.optimize import newton\n"
    "import numpy as np\n"
    "from vakyume.config import UnsolvedException\n\n"
)

LIBRARY_IMPORT_HEADER = (
    "from cmath import log, sqrt, exp\n"
    "from math import e, pi\n"
    "from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp\n"
    "from scipy.optimize import newton\n"
    "from vakyume.kwasak import kwasak\n"
    "from vakyume.config import UnsolvedException\n"
    "import numpy as np\n"
)
