"""Vakyume: verified solvers for every variable of textbook equations.

    >>> from vakyume import Library
    >>> lib = Library.load("projects/VacuumTheory")
    >>> lib["2-1"].solve(rho=1.2, D=0.1, v=3.0, mu=1.8e-5)   # -> Re
    >>> lib.FluidFlowVacuumLines.eqn_2_1(Re=2.0e4, rho=1.2, D=0.1, v=3.0)  # -> mu

The original scrape/shard/verify/repair pipeline is still importable
(``from vakyume import run_pipeline``) with the ``legacy`` extra installed.
"""

from importlib import import_module

from .equations import Chapter, Equation, Library, Variable, load_toml, parse_notes
from .solve import AmbiguousSolution, NoSolution, Solution

__all__ = [
    "AmbiguousSolution",
    "Chapter",
    "Equation",
    "Library",
    "NoSolution",
    "Solution",
    "Variable",
    "load_toml",
    "parse_notes",
]

# Legacy names, imported on first use so the new core does not pull in
# timeout-decorator, ollama, etc.
_LEGACY = {
    "run_pipeline": ".pipeline",
    "Solver": ".parser",
    "Verify": ".verifier",
    "reconstruct_from_shards": ".reconstruct",
    "reconstruct_cli": ".reconstruct",
    "UnsolvedException": ".config",
    "TAB": ".config",
    "MAX_COMP_TIME_SECONDS": ".config",
    "COOLDOWN_SECONDS": ".config",
    "OPERATORS": ".config",
    "FUNKTORZ": ".config",
    "llm_config": ".config",
}


def __getattr__(name: str):
    if name == "run_cpp_gen":
        return import_module(".cpp_gen", __name__).main
    if name in _LEGACY:
        return getattr(import_module(_LEGACY[name], __name__), name)
    raise AttributeError(f"module 'vakyume' has no attribute {name!r}")
