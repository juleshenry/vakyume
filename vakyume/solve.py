"""Tiered solver: deterministic SymPy first, numeric fallback, residual as judge.

For each (equation, target) pair a :class:`TargetSolver` is prepared once:

1. **Peel.** Invert the outer layers around the target (``+``, ``*``,
   non-integer rational powers, ``exp``, ``log``, ``a**x``) so that e.g.
   ``S = k * B(P)**(3/5)`` becomes ``B(P) = (S/k)**(5/3)``.
2. **Closed form / polynomial.** If the target is then isolated we have a
   closed form. If the residual's numerator is a polynomial in the target we
   keep its coefficients: degree <= 2 gets exact roots, higher degrees are
   solved numerically with ``numpy.roots``. Both are *complete*: every real
   root is found.
3. **SymPy solve.** Otherwise ``sympy.solve`` runs in a worker process with a
   timeout. Its answers are candidates only; the result is marked incomplete
   and a numeric scan is merged in at call time.

At call time every candidate root is substituted back into the original
equation and kept only if the scaled residual is below tolerance. That
check -- not agreement with other solvers -- decides what is returned.

Only the math domain is applied (principal branches, real values). Physical
constraints are the caller's job: see ``select``/``where``/``bounds``.
"""

from __future__ import annotations

import hashlib
import os
import pickle
import warnings
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import sympy as sp
from scipy.optimize import brentq, minimize_scalar

ANALYZER_VERSION = "1"
SOLVE_TIMEOUT_S = float(os.environ.get("VAKYUME_SOLVE_TIMEOUT", "10"))
MAX_CLOSED_OPS = 400  # reject giant closed forms (e.g. the 10 KB quartic)
REL_TOL = 1e-9  # scaled residual tolerance for accepting a root
_LAMBDIFY_MODULES = ["numpy", "scipy"]


class NoSolution(ValueError):
    pass


class AmbiguousSolution(ValueError):
    def __init__(self, msg: str, roots: tuple[Any, ...]):
        super().__init__(msg)
        self.roots = roots


@dataclass(frozen=True)
class Solution:
    target: str
    roots: tuple[Any, ...]
    complete: bool  # True: every real root (in the math domain) is listed
    method: str  # closed | polynomial | sympy | numeric (joined with '+')


# ── preparation ──────────────────────────────────────────────────────────────


@dataclass
class _Piece:
    """One branch after peeling: solve ``inner(x) = outer(knowns)``."""

    kind: str  # closed | polynomial | sympy | numeric
    exprs: list[sp.Expr]  # closed/sympy: root exprs; polynomial: coeffs (high->low)
    complete: bool


_INVERTIBLE_FUNCS = (sp.exp, sp.log)


def _peel(lhs: sp.Expr, rhs: sp.Expr, x: sp.Symbol) -> list[tuple[sp.Expr, sp.Expr]]:
    """Isolate ``x`` as far as possible. Returns [(outer, inner)] pairs with
    ``inner`` containing x and ``outer`` free of it."""
    if lhs.has(x) and rhs.has(x):
        return [(sp.Integer(0), rhs - lhs)]
    outer, inner = (rhs, lhs) if lhs.has(x) else (lhs, rhs)
    pairs = [(outer, inner)]
    for _ in range(50):
        changed = False
        nxt = []
        for L, R in pairs:
            if R == x:
                nxt.append((L, R))
                continue
            if R.is_Add or R.is_Mul:
                with_x = [a for a in R.args if a.has(x)]
                if len(with_x) == 1:
                    rest = R.func(*[a for a in R.args if not a.has(x)])
                    L = L - rest if R.is_Add else L / rest
                    nxt.append((L, with_x[0]))
                    changed = True
                    continue
            elif R.is_Pow:
                b, e = R.args
                if b.has(x) and not e.has(x) and e.is_Rational and not e.is_Integer:
                    # principal branch: b**(p/q) with q > 1 needs b >= 0,
                    # so the inverse is unique.
                    nxt.append((L ** (1 / e), b))
                    changed = True
                    continue
                if b.has(x) and not e.has(x) and e.is_Integer and e == -1:
                    nxt.append((1 / L, b))
                    changed = True
                    continue
                if e.has(x) and not b.has(x):
                    nxt.append((sp.log(L) / sp.log(b), e))
                    changed = True
                    continue
            elif isinstance(R, sp.exp):
                nxt.append((sp.log(L), R.args[0]))
                changed = True
                continue
            elif isinstance(R, sp.log) and len(R.args) == 1:
                nxt.append((sp.exp(L), R.args[0]))
                changed = True
                continue
            nxt.append((L, R))
        pairs = nxt
        if not changed:
            break
    return pairs


def _polynomial_coeffs(expr: sp.Expr, x: sp.Symbol) -> list[sp.Expr] | None:
    num, _den = sp.fraction(sp.together(expr))
    num = sp.expand(num)
    if not num.has(x) or not num.is_polynomial(x):
        return None
    try:
        return sp.Poly(num, x).all_coeffs()
    except sp.PolynomialError:
        return None


def _sympy_solve_worker(expr: sp.Expr, x: sp.Symbol) -> list[sp.Expr]:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sols = sp.solve(expr, x, dict=False, check=True)
    return [s for s in sols if not isinstance(s, (dict, tuple))]


def _cache_dir() -> Path:
    base = os.environ.get("VAKYUME_CACHE_DIR") or os.path.join(
        os.environ.get("XDG_CACHE_HOME", os.path.expanduser("~/.cache")), "vakyume"
    )
    return Path(base) / f"solve-v{ANALYZER_VERSION}"


def _sympy_solve_with_timeout(expr: sp.Expr, x: sp.Symbol) -> list[sp.Expr] | None:
    """``sympy.solve`` in a killable worker. None on timeout or error.

    Results are cached on disk keyed by the expression, so the expensive part
    runs once per equation, not once per process.
    """
    key = hashlib.sha256(f"{sp.srepr(expr)}|{x.name}".encode()).hexdigest()
    path = _cache_dir() / f"{key}.pkl"
    try:
        return pickle.loads(path.read_bytes())
    except (OSError, pickle.PickleError, EOFError, AttributeError):
        pass

    import multiprocessing as mp

    from pebble import ProcessPool

    # fork where available: the worker only runs pure-Python SymPy, and spawn
    # would re-import the caller's __main__ (breaking unguarded scripts).
    ctx = mp.get_context("fork" if "fork" in mp.get_all_start_methods() else "spawn")
    result: list[sp.Expr] | None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)  # "fork with threads" notice
        pool = ProcessPool(max_workers=1, context=ctx)
    with pool:
        fut = pool.schedule(_sympy_solve_worker, args=(expr, x), timeout=SOLVE_TIMEOUT_S)
        try:
            result = fut.result()
        except Exception:  # timeout, NotImplementedError, anything SymPy throws
            result = None

    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(pickle.dumps(result))
    except OSError:
        pass
    return result


def _single_monotone_occurrence(expr: sp.Expr, x: sp.Symbol) -> bool:
    """x occurs once, only under operations whose real inverse is unique
    (or fully enumerated by SymPy)."""
    if expr.count(x) != 1:
        return False
    node = expr
    while node != x:
        if not (node.is_Add or node.is_Mul or node.is_Pow or isinstance(node, _INVERTIBLE_FUNCS)):
            return False
        node = next(a for a in node.args if a.has(x))
    return True


def prepare_pieces(lhs: sp.Expr, rhs: sp.Expr, x: sp.Symbol) -> list[_Piece]:
    pieces = []
    for outer, inner in _peel(lhs, rhs, x):
        if inner == x:
            pieces.append(_Piece("closed", [outer], complete=True))
            continue
        coeffs = _polynomial_coeffs(inner - outer, x)
        if coeffs is not None:
            degree = len(coeffs) - 1
            if degree <= 2:
                roots = sp.roots(sp.Poly(coeffs, x), multiple=True)
                if len(roots) == degree:
                    pieces.append(_Piece("closed", roots, complete=True))
                    continue
            pieces.append(_Piece("polynomial", coeffs, complete=True))
            continue
        sols = _sympy_solve_with_timeout(inner - outer, x)
        sols = [s for s in (sols or []) if sp.count_ops(s) <= MAX_CLOSED_OPS]
        complete = bool(sols) and _single_monotone_occurrence(inner, x)
        pieces.append(_Piece("sympy" if sols else "numeric", sols, complete=complete))
    return pieces


# ── runtime ──────────────────────────────────────────────────────────────────


def _terms(e: sp.Expr) -> list[sp.Expr]:
    return list(sp.Add.make_args(e))


class TargetSolver:
    def __init__(
        self,
        lhs: sp.Expr,
        rhs: sp.Expr,
        target: str,
        variables: Iterable[str],
        hint: tuple[float | None, float | None] | None = None,
    ):
        self.target = target
        self.hint = hint
        self._syms = {s.name: s for s in (lhs - rhs).free_symbols}
        self.x = self._syms[target]
        self.knowns = [v for v in variables if v != target]
        args = [self._syms[n] for n in self.knowns] + [self.x]
        self._lhs, self._rhs = lhs, rhs
        self._args = args
        # residual pieces for the acceptance check
        self._f_lhs = sp.lambdify(args, lhs, modules=_LAMBDIFY_MODULES)
        self._f_rhs = sp.lambdify(args, rhs, modules=_LAMBDIFY_MODULES)
        self._f_terms = sp.lambdify(args, _terms(lhs) + _terms(rhs), modules=_LAMBDIFY_MODULES)
        self._pieces: list[_Piece] | None = None
        self._piece_fns: list[Callable[..., Any]] = []

    # lazily, because sympy.solve may be expensive
    @property
    def pieces(self) -> list[_Piece]:
        if self._pieces is None:
            self._pieces = prepare_pieces(self._lhs, self._rhs, self.x)
            known_syms = [self._syms[n] for n in self.knowns]
            self._piece_fns = [
                sp.lambdify(known_syms, p.exprs, modules=_LAMBDIFY_MODULES) if p.exprs else None
                for p in self._pieces
            ]
        return self._pieces

    @property
    def complete(self) -> bool:
        return all(p.complete for p in self.pieces)

    @property
    def method(self) -> str:
        kinds = [p.kind for p in self.pieces]
        if not self.complete and "numeric" not in kinds:
            kinds.append("numeric")
        return "+".join(dict.fromkeys(kinds))

    # residual ----------------------------------------------------------------

    def scaled_residual(self, vals: list[Any], root: Any) -> float:
        with np.errstate(all="ignore"), warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                a = self._f_lhs(*vals, root)
                b = self._f_rhs(*vals, root)
                scale = sum(abs(t) for t in self._f_terms(*vals, root))
            except (ZeroDivisionError, OverflowError, ValueError, TypeError):
                return float("inf")
        err = abs(complex(a) - complex(b))
        if not np.isfinite(err) or not np.isfinite(scale):
            return float("inf")
        if err == 0:
            return 0.0
        return err / scale if scale > 0 else float("inf")

    def _accept(self, vals: list[Any], root: Any) -> bool:
        return self.scaled_residual(vals, root) <= REL_TOL

    def _polish(self, vals: list[float], r: float) -> float:
        """A few secant/Newton steps on the residual; keeps the best."""

        def f(t: float) -> float:
            with np.errstate(all="ignore"):
                return float(self._f_rhs(*vals, t) - self._f_lhs(*vals, t))

        best, best_res = r, self.scaled_residual(vals, r)
        x0 = r
        for _ in range(4):
            h = 1e-7 * max(abs(x0), 1e-12)
            try:
                fx, d = f(x0), (f(x0 + h) - f(x0 - h)) / (2 * h)
            except (ZeroDivisionError, OverflowError, ValueError):
                break
            if not np.isfinite(fx) or not np.isfinite(d) or d == 0:
                break
            x0 = x0 - fx / d
            if abs(x0 - r) > 1e-6 * max(abs(r), 1e-300):
                break  # refine only; never walk to a different root
            res = self.scaled_residual(vals, x0)
            if res < best_res:
                best, best_res = x0, res
        return best

    # candidates --------------------------------------------------------------

    def _candidates(self, vals: list[float]) -> list[complex]:
        cvals = [np.complex128(v) for v in vals]
        out: list[complex] = []
        for piece, fn in zip(self.pieces, self._piece_fns):
            if fn is None:
                continue
            with np.errstate(all="ignore"), warnings.catch_warnings():
                warnings.simplefilter("ignore")
                try:
                    values = [complex(v) for v in fn(*cvals)]
                except (ZeroDivisionError, OverflowError, ValueError, TypeError):
                    continue
            if piece.kind == "polynomial":
                coeffs = np.array(values, dtype=complex)
                if not np.all(np.isfinite(coeffs)) or not np.any(coeffs):
                    continue
                values = list(np.roots(coeffs))
            out += [v for v in values if np.isfinite(v)]
        return out

    def _f_real(self, vals: list[float]) -> Callable[[Any], Any]:
        def f(t: Any) -> Any:
            return self._f_rhs(*vals, t) - self._f_lhs(*vals, t)

        return f

    def _scan(self, vals: list[float]) -> list[float]:
        """Sign-change + local-minimum scan over a symmetric log grid.

        Without bounds this cannot promise completeness, which is why the
        result of an incomplete method is flagged ``complete=False``.
        """
        pos = np.logspace(-12, 12, 1201)
        grids = [-pos[::-1], [0.0], pos]
        if self.hint and self.hint[0] is not None and self.hint[1] is not None:
            lo, hi = self.hint
            grids.append(np.linspace(lo, hi, 401))
        grid = np.unique(np.concatenate(grids))
        f = self._f_real(vals)
        with np.errstate(all="ignore"), warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                y = np.asarray(f(grid), dtype=float) * np.ones_like(grid)
            except (ZeroDivisionError, OverflowError, ValueError, TypeError):
                y = np.array([_safe_float(f, t) for t in grid])
        ok = np.isfinite(y)
        roots: list[float] = list(grid[ok & (y == 0)])
        s = np.sign(y)
        idx = np.nonzero(ok[:-1] & ok[1:] & (s[:-1] * s[1:] < 0))[0]
        for i in idx:
            try:
                roots.append(
                    brentq(f, grid[i], grid[i + 1], xtol=1e-300, rtol=4 * np.finfo(float).eps, maxiter=200)
                )
            except (ValueError, RuntimeError, ZeroDivisionError, OverflowError):
                pass
        # touching roots (even multiplicity): local minima of |f|
        a = np.where(ok, np.abs(y), np.inf)
        mins = np.nonzero((a[1:-1] < a[:-2]) & (a[1:-1] < a[2:]))[0] + 1
        for i in mins:
            lo, hi = grid[i - 1], grid[i + 1]
            try:
                r = minimize_scalar(
                    lambda t: abs(_safe_float(f, t)),
                    bounds=(lo, hi),
                    method="bounded",
                    options={"xatol": 1e-14 * max(abs(grid[i]), 1e-300)},
                )
                roots.append(float(r.x))
            except (ValueError, RuntimeError):
                pass
        return roots

    # public ------------------------------------------------------------------

    def __call__(self, knowns: dict[str, float], complex: bool = False) -> Solution:
        vals = [float(knowns[n]) for n in self.knowns]
        found: list[Any] = []
        for c in self._candidates(vals):
            if not complex:
                if abs(c.imag) > 1e-6 * max(1.0, abs(c.real)):
                    continue
                c = self._polish(vals, float(c.real))
            if self._accept(vals, c):
                found.append(c)
        if not self.complete:
            for r in self._scan(vals):
                if self._accept(vals, r):
                    found.append(r)
        return Solution(self.target, _dedupe(found), self.complete, self.method)


def _safe_float(f: Callable[[Any], Any], t: float) -> float:
    try:
        with np.errstate(all="ignore"):
            v = complex(f(t))
        return v.real if v.imag == 0 else float("nan")
    except (ZeroDivisionError, OverflowError, ValueError, TypeError):
        return float("nan")


def _dedupe(roots: list[Any], rel: float = 1e-9) -> tuple[Any, ...]:
    out: list[Any] = []
    for r in sorted(roots, key=lambda z: (np.real(z), np.imag(z))):
        if out and abs(r - out[-1]) <= rel * max(abs(r), abs(out[-1]), 1e-300):
            continue
        out.append(r)
    return tuple(float(r) if isinstance(r, (float, np.floating)) else r for r in out)


def select_root(
    sol: Solution,
    select: str = "unique",
    where: Callable[[Any], bool] | None = None,
    bounds: tuple[float | None, float | None] | None = None,
) -> Any:
    roots = sol.roots
    if bounds is not None:
        lo, hi = bounds
        roots = tuple(r for r in roots if (lo is None or r >= lo) and (hi is None or r <= hi))
    if where is not None:
        roots = tuple(r for r in roots if where(r))
    if select == "all":
        return Solution(sol.target, roots, sol.complete, sol.method)
    if not roots:
        raise NoSolution(f"no real solution for {sol.target} (method: {sol.method})")
    if select == "min":
        return min(roots)
    if select == "max":
        return max(roots)
    if select == "unique":
        if len(roots) > 1:
            raise AmbiguousSolution(
                f"{len(roots)} roots for {sol.target}: {roots}; pass select='all'|'min'|'max', "
                "where=..., or bounds=... to choose",
                roots,
            )
        return roots[0]
    raise ValueError(f"unknown select={select!r}")
