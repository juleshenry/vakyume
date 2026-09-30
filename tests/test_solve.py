import math

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from vakyume import AmbiguousSolution, NoSolution
from vakyume.equations import Equation, Library

from .conftest import PROJECTS

VT = Library.load(PROJECTS / "VacuumTheory")


def eq(text: str) -> Equation:
    return Equation("t", text)


def test_kwasak_style_dispatch():
    e = eq("Re = rho * D * v / mu")
    assert e.solve(rho=2.0, D=3.0, v=4.0, mu=6.0) == pytest.approx(4.0)
    assert e(Re=4.0, rho=2.0, D=3.0, v=4.0) == pytest.approx(6.0)
    with pytest.raises(TypeError):
        e.solve(rho=2.0, D=3.0)


def test_all_real_roots_and_caller_selects():
    e = eq("q = D ** 4 * k")
    sol = e.solve(q=16.0, k=1.0, select="all")
    assert sol.complete and sol.roots == pytest.approx((-2.0, 2.0))
    with pytest.raises(AmbiguousSolution):
        e.solve(q=16.0, k=1.0)
    assert e.solve(q=16.0, k=1.0, where=lambda r: r > 0) == pytest.approx(2.0)
    assert e.solve(q=16.0, k=1.0, select="min") == pytest.approx(-2.0)
    with pytest.raises(NoSolution):
        e.solve(q=-16.0, k=1.0)
    assert len(e.solve(q=-16.0, k=1.0, select="all", complex=True).roots) == 4


def test_math_domain_only():
    # log needs a positive argument: x = -e is not a root even though it
    # would satisfy log(|x|) = 1.
    e = eq("y = log(x)")
    assert e.solve(y=1.0, select="all").roots == pytest.approx((math.e,))


def test_odd_integer_power_keeps_negative_real_root():
    e = eq("y = x ** 3")
    assert e.solve(y=-8.0) == pytest.approx(-2.0)


@pytest.mark.parametrize(
    "eq_id,target,method",
    [
        ("2-34", "D", "polynomial"),  # quartic: was a 10 KB closed form
        ("8-7", "P_1", "closed"),  # float exponent 0.286: was "Pending Repair"
        ("10-19", "P", "closed"),  # rational power of a rational function
        ("10-20", "P", "closed"),
        ("2-5", "D", "polynomial"),  # complex branches filtered out
    ],
)
def test_regression_equations_are_complete(eq_id, target, method):
    s = VT[eq_id].solver(target)
    assert s.complete
    assert s.method == method


def test_transcendental_falls_back_to_numeric():
    e = VT["7-14b"]
    s = e.solver("del_T_1")
    assert not s.complete and "numeric" in s.method
    knowns = dict(Q_condensor_heat_duty=1000.0, U=2.0, del_T_2=10.0)
    lmtd = (30.0 - 10.0) / math.log(3.0)
    knowns["A"] = 1000.0 / (2.0 * lmtd)
    roots = e.solve_for("del_T_1", knowns, select="all").roots
    assert any(r == pytest.approx(30.0) for r in roots)


@settings(max_examples=40, deadline=None)
@given(
    rho=st.floats(1e-2, 1e3),
    D=st.floats(1e-2, 1e3),
    v=st.floats(1e-2, 1e3),
    mu=st.floats(1e-2, 1e3),
)
def test_roundtrip_property(rho, D, v, mu):
    e = VT["2-1"]
    Re = e.solve(rho=rho, D=D, v=v, mu=mu)
    assert e.solve(Re=Re, rho=rho, v=v, mu=mu) == pytest.approx(D, rel=1e-9)


@settings(max_examples=25, deadline=None)
@given(
    C_1=st.floats(0.1, 10),
    C_2=st.floats(0.1, 10),
    D=st.floats(0.1, 10),
    mu=st.floats(0.1, 10),
    L=st.floats(0.1, 10),
    P_p=st.floats(0.1, 10),
)
def test_quartic_roundtrip_property(C_1, C_2, D, mu, L, P_p):
    e = VT["2-34"]
    C = e.solve(C_1=C_1, C_2=C_2, D=D, mu=mu, L=L, P_p=P_p)
    roots = e.solve(C=C, C_1=C_1, C_2=C_2, mu=mu, L=L, P_p=P_p, select="all").roots
    assert any(r == pytest.approx(D, rel=1e-7) for r in roots)
