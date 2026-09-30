"""Residual-based verification. The equation itself is the only oracle.

For each algebraic equation:

1. Build consistent sample points: draw all variables but a random *pivot*,
   solve for the pivot, and keep the point only if it satisfies the
   equation (checked independently, see below). Which solver produced the
   point does not matter -- only that the equation holds there.
2. For every target variable and every sample point: hide the target, solve
   for it, and check
   * **soundness** -- every returned root satisfies the equation, re-checked
     with SymPy at 30 digits (independent of the lambdified solver code);
   * **round trip** -- the hidden value is among the returned roots. A miss
     is a failure if the solver claims completeness, and a ``partial``
     otherwise.

Nothing passes by default: an equation with no sample points is
``unsampled``, non-algebraic ones are ``unsupported``/``invalid``.
"""

from __future__ import annotations

import json
import traceback
import zlib
from collections.abc import Callable
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np
import sympy as sp

from .equations import ALGEBRAIC, Equation, Library

CHECK_TOL = 1e-8  # scaled residual at 30 digits
ROUNDTRIP_RTOL = 1e-6

# status ordering, worst first
STATUSES = ("fail", "invalid", "unsampled", "partial", "unsupported", "pass")


@dataclass
class TargetResult:
    target: str
    method: str = ""
    complete: bool = False
    samples: int = 0
    unsound: int = 0  # roots that do not satisfy the equation
    missed: int = 0  # hidden value not among the roots
    errors: list[str] = field(default_factory=list)

    @property
    def status(self) -> str:
        if self.errors or self.unsound or (self.missed and self.complete):
            return "fail"
        if self.samples == 0:
            return "unsampled"
        if self.missed:
            return "partial"
        return "pass"


@dataclass
class EquationResult:
    id: str
    chapter: str
    text: str
    kind: str
    status: str
    note: str = ""
    samples: int = 0
    targets: list[TargetResult] = field(default_factory=list)

    def to_dict(self) -> dict:
        d = asdict(self)
        for t, tr in zip(d["targets"], self.targets):
            t["status"] = tr.status
        return d


# ── independent residual check ───────────────────────────────────────────────


def exact_residual(eq: Equation, point: dict[str, float]) -> float:
    """Scaled residual evaluated by SymPy at 30 significant digits."""
    rel = eq.relation
    assert rel.lhs is not None and rel.rhs is not None
    subs = {sp.Symbol(k): sp.Float(repr(float(v)), 30) for k, v in point.items()}
    terms = list(sp.Add.make_args(rel.lhs)) + list(sp.Add.make_args(rel.rhs))
    try:
        a = complex(rel.lhs.evalf(30, subs=subs))
        b = complex(rel.rhs.evalf(30, subs=subs))
        scale = sum(abs(complex(t.evalf(30, subs=subs))) for t in terms)
    except (TypeError, ValueError, ZeroDivisionError, OverflowError):
        return float("inf")
    err = abs(a - b)
    if not np.isfinite(err) or not np.isfinite(scale):
        return float("inf")
    if err == 0:
        return 0.0
    return err / scale if scale > 0 else float("inf")


# ── sampling ─────────────────────────────────────────────────────────────────

Draw = Callable[[str, "tuple[float | None, float | None] | None"], float]


def rng_draw(rng: np.random.Generator) -> Draw:
    """Default sampler: log-uniform magnitudes in [1e-2, 1e3], or within the
    variable's hint. Positive values are a sampling choice (realistic inputs),
    not a restriction on what the solvers may return."""

    def draw(_name: str, hint: tuple[float | None, float | None] | None) -> float:
        if hint and hint[0] is not None and hint[1] is not None and hint[1] > hint[0]:
            lo, hi = hint
            if lo > 0:
                return float(np.exp(rng.uniform(np.log(lo), np.log(hi))))
            return float(rng.uniform(lo, hi))
        return float(10 ** rng.uniform(-2, 3))

    return draw


def sample_points(
    eq: Equation,
    n: int,
    draw: Draw,
    rng: np.random.Generator,
    max_attempts: int | None = None,
) -> list[dict[str, float]]:
    points: list[dict[str, float]] = []
    names = list(eq.symbols)
    for _ in range(max_attempts or n * 20):
        if len(points) >= n:
            break
        pivot = names[int(rng.integers(len(names)))]
        point = {v: draw(v, eq.variable(v).hint) for v in names if v != pivot}
        try:
            sol = eq.solve_for(pivot, point, select="all")
        except Exception:
            continue
        if not sol.roots:
            continue
        point[pivot] = float(sol.roots[int(rng.integers(len(sol.roots)))])
        if exact_residual(eq, point) <= CHECK_TOL:
            points.append(point)
    return points


# ── verification ─────────────────────────────────────────────────────────────


def _seed(eq: Equation, seed: int) -> int:
    return seed ^ zlib.crc32(f"{eq.chapter}/{eq.id}".encode())


def verify_target(eq: Equation, target: str, points: list[dict[str, float]]) -> TargetResult:
    res = TargetResult(target)
    try:
        solver = eq.solver(target)
        res.method, res.complete = solver.method, solver.complete
    except Exception as exc:
        res.errors.append(f"prepare: {exc!r}")
        return res
    for point in points:
        knowns = {k: v for k, v in point.items() if k != target}
        hidden = point[target]
        try:
            sol = eq.solve_for(target, knowns, select="all")
        except Exception as exc:
            res.errors.append(f"{knowns}: {exc!r}")
            continue
        res.samples += 1
        for r in sol.roots:
            if exact_residual(eq, {**knowns, target: r}) > CHECK_TOL:
                res.unsound += 1
        if not _recovered(eq, knowns, target, hidden, sol.roots):
            res.missed += 1
    return res


def _recovered(eq: Equation, knowns: dict[str, float], target: str, hidden: float, roots) -> bool:
    """Is ``hidden`` among ``roots``? Relative closeness first; for
    ill-conditioned points (terms nearly cancelling) the closest root counts
    as the same one if the equation also holds at their midpoint."""
    if not roots:
        return False
    closest = min(roots, key=lambda r: abs(r - hidden))
    if abs(closest - hidden) <= ROUNDTRIP_RTOL * max(abs(hidden), 1e-300):
        return True
    mid = (closest + hidden) / 2
    return exact_residual(eq, {**knowns, target: mid}) <= CHECK_TOL * 10


def verify_equation(eq: Equation, n: int = 12, seed: int = 0, draw: Draw | None = None) -> EquationResult:
    base = EquationResult(eq.id, eq.chapter, eq.text, eq.kind, status="unsupported")
    if eq.kind != ALGEBRAIC:
        base.status = "invalid" if eq.kind == "invalid" else "unsupported"
        base.note = eq.relation.error
        return base
    rng = np.random.default_rng(_seed(eq, seed))
    try:
        points = sample_points(eq, n, draw or rng_draw(rng), rng)
    except Exception:
        base.status, base.note = "fail", traceback.format_exc(limit=2)
        return base
    base.samples = len(points)
    if not points:
        base.status, base.note = "unsampled", "no consistent real sample point found"
        return base
    base.targets = [verify_target(eq, t, points) for t in eq.symbols]
    statuses = {t.status for t in base.targets}
    base.status = next(s for s in STATUSES if s in statuses)
    return base


def verify_library(
    lib: Library,
    n: int = 12,
    seed: int = 0,
    only: Callable[[Equation], bool] | None = None,
    progress: Callable[[EquationResult], None] | None = None,
) -> list[EquationResult]:
    results = []
    for eq in lib.equations:
        if only and not only(eq):
            continue
        r = verify_equation(eq, n=n, seed=seed)
        results.append(r)
        if progress:
            progress(r)
    return results


# ── reporting ────────────────────────────────────────────────────────────────


def summarize(results: list[EquationResult]) -> dict[str, int]:
    counts = {s: 0 for s in STATUSES}
    for r in results:
        counts[r.status] += 1
    methods: dict[str, int] = {}
    for r in results:
        for t in r.targets:
            key = t.method.split("+")[0] or "?"
            methods[key] = methods.get(key, 0) + 1
    return {**counts, **{f"solvers_{k}": v for k, v in sorted(methods.items())}}


def write_reports(results: list[EquationResult], project_dir: Path, title: str) -> tuple[Path, Path]:
    docs = project_dir / "docs"
    docs.mkdir(parents=True, exist_ok=True)
    json_path = docs / "status.json"
    json_path.write_text(json.dumps([r.to_dict() for r in results], indent=1) + "\n")

    s = summarize(results)
    lines = [
        f"# Equation status: {title}",
        "",
        "Generated by `vakyume report`. Every root is substituted back into the",
        "original equation (30-digit SymPy check); no solver is trusted as an oracle.",
        "",
        "| status | equations |",
        "|---|---|",
        *[f"| {k} | {s[k]} |" for k in STATUSES],
        "",
        "Solvers by method: "
        + ", ".join(f"{k.removeprefix('solvers_')} {v}" for k, v in s.items() if k.startswith("solvers_")),
        "",
        "| id | status | equation | per-variable (method: pass/partial/fail) | note |",
        "|---|---|---|---|---|",
    ]
    for r in results:
        per = ", ".join(f"{t.target}:{t.method}{'' if t.complete else '~'}/{t.status}" for t in r.targets)
        text = r.text.replace("|", "\\|")
        note = r.note.replace("|", "\\|").replace("\n", " ")[:80]
        lines.append(f"| {r.id} | {r.status} | `{text}` | {per} | {note} |")
    lines += ["", "`~` marks a solver that cannot prove it found every root (numeric scan)."]
    md_path = docs / "STATUS.md"
    md_path.write_text("\n".join(lines) + "\n")
    return md_path, json_path
