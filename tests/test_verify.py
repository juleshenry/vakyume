import pytest

from vakyume.equations import Equation, Library
from vakyume.verify import verify_equation, verify_library

from .conftest import PROJECTS


def test_correct_equation_passes():
    r = verify_equation(Equation("t", "C = C_1 * D ** 4 / (mu * L) * P + C_2 * D ** 3 / L"), n=5)
    assert r.status == "pass"
    assert r.samples == 5
    assert all(t.samples == 5 for t in r.targets)


def test_nothing_passes_by_default():
    assert verify_equation(Equation("t", "F = *x"), n=3).status == "invalid"
    assert verify_equation(Equation("t", "Reff = 1 / (1/R1 + ...)"), n=3).status == "unsupported"
    # log of a negative number everywhere: no real sample point exists
    r = verify_equation(Equation("t", "y = log(-x ** 2 - 1)"), n=3)
    assert r.status == "unsampled"


def test_broken_solver_is_caught():
    e = Equation("t", "z = x + y")
    solver = e.solver("z")
    _ = solver.pieces  # prepare
    solver._piece_fns[0] = lambda x, y: [x - y]  # deliberate sign error
    r = verify_equation(e, n=5)
    z = next(t for t in r.targets if t.target == "z")
    # wrong candidates are rejected by the residual check, so the solver
    # returns nothing and the round trip fails -- a complete solver that
    # misses the hidden value is a failure, not a pass.
    assert z.status == "fail"


@pytest.mark.slow
@pytest.mark.parametrize("project", ["VacuumTheory", "BuildingModels"])
def test_whole_project(project):
    """Heavy: verifies every equation. Run with `pytest -m slow`."""
    results = verify_library(Library.load(PROJECTS / project), n=8)
    failing = [(r.id, r.status, r.note) for r in results if r.status == "fail"]
    assert not failing
