import sympy as sp

from vakyume.equations import Library, dump_toml, load_toml, parse_notes, parse_relation


def test_parses_physics_names_as_symbols():
    rel = parse_relation("E = I * S + N * beta + gamma")
    assert rel.kind == "algebraic"
    assert {s.name for s in (rel.lhs - rel.rhs).free_symbols} == {"E", "I", "S", "N", "beta", "gamma"}


def test_float_exponents_become_rationals_and_pi():
    rel = parse_relation("y = 3.141592653589793 * (P_2 / P_1) ** 0.286")
    assert rel.rhs.has(sp.pi)
    assert sp.Rational(143, 500) in rel.rhs.atoms(sp.Rational)
    assert not rel.rhs.atoms(sp.Float)


def test_ln_and_approx_equal():
    rel = parse_relation("W ~= ln(P) * V ** (3/5)")
    assert rel.kind == "algebraic"
    assert rel.rhs.has(sp.log)


def test_unsupported_kinds():
    assert parse_relation("Reff = 1 / (1 / R1 + 1 / R2 + ...)").kind == "aggregate"
    assert parse_relation("p_v / p_s <= P_0_v / P_D").kind == "inequality"
    assert parse_relation("L = r ^ p").kind == "vector"
    assert parse_relation("PS = -V * dP / dT + Q_0").kind == "differential"
    assert parse_relation("m * d2x / dt2 = k * x").kind == "differential"
    assert parse_relation("W = U(B) - U(A)").kind == "functional"


def test_delta_names_are_not_derivatives():
    assert parse_relation("x = delta_P / delta_h").kind == "algebraic"


def test_extraction_losses_are_invalid():
    assert parse_relation("F = *x").kind == "invalid"
    rel = parse_relation("x = 0 * (x + v*t)/(t*v/c**2)")
    assert rel.kind == "invalid" and "cancel" in rel.error


NOTES = '''# Chapter 7 : Things
# 7-1 First
"""
P := pressure, torr
1 <= P <= 10 #torr
"""
W ~= 0.026 * P ** 0.34
# 7-2 Chained 10-3 inside an equation must not retag
f_m = 24 * Q_r / h = 0.0266 * Q_r
q = x * 10-3
"""
stray unbalanced docstring
# 7-3 After the stray quotes
y = 2 * x
"""
x := after-the-fact note
"""

def eqn_7_3(x):
    y = 2 * x
    return y
'''


def test_parse_notes_is_robust(tmp_path):
    p = tmp_path / "07_things.py"
    p.write_text(NOTES)
    ch = parse_notes(p)
    assert ch.name == "Things" and ch.title == "Things" and ch.number == 7
    ids = [e.id for e in ch.equations]
    assert ids == ["7-1", "7-2", "7-2.2", "7-2.3", "7-3"]
    first = ch.equations[0]
    assert (first.variable("P").min, first.variable("P").max) == (1.0, 10.0)
    assert first.variable("P").desc == "pressure, torr"
    assert ch.equations[-1].variable("x").desc == "after-the-fact note"


def test_toml_roundtrip(tmp_path):
    p = tmp_path / "07_things.py"
    p.write_text(NOTES)
    ch = parse_notes(p)
    t = tmp_path / "07_things.toml"
    t.write_text(dump_toml(ch))
    back = load_toml(t)
    assert [(e.id, e.text, e.title) for e in back.equations] == [
        (e.id, e.text, e.title) for e in ch.equations
    ]
    assert back.equations[0].variable("P").max == 10.0


def test_library_attribute_access(tmp_path):
    (tmp_path / "notes").mkdir()
    (tmp_path / "notes" / "02_fluid_flow.py").write_text("# 2-1 Reynolds\nRe = rho * D * v / mu\n")
    lib = Library.load(tmp_path)
    assert lib.FluidFlow.eqn_2_1 is lib["2-1"]
