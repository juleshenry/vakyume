"""Equation records: the single source of truth for every solver.

An equation is stored once, as text (``lhs = rhs``) plus metadata. Solvers
for each variable are derived from it on demand (see :mod:`vakyume.solve`),
so there are no per-variable generated files to drift out of sync.

Two on-disk formats are supported:

* ``equations/*.toml`` -- the structured format (preferred).
* ``notes/*.py``       -- the legacy notes DSL, parsed here so it can be
  converted once with ``vakyume convert-notes``.
"""

from __future__ import annotations

import json
import math
import re
import tomllib
from collections.abc import Callable
from dataclasses import dataclass, field
from functools import cached_property
from pathlib import Path
from tokenize import TokenError
from typing import TYPE_CHECKING, Any

import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations

if TYPE_CHECKING:
    from .solve import TargetSolver

# Kinds of relation. Only ``algebraic`` equations get solvers; the rest are
# kept as reference material so nothing is silently dropped.
ALGEBRAIC = "algebraic"
KINDS = (
    ALGEBRAIC,
    "inequality",  # 1 <= P <= 10, p_v / p_s <= P_0_v / P_D
    "aggregate",  # sums with ellipses: 1/R1 + 1/R2 + ...
    "differential",  # dP / dt, d2x / dt2
    "functional",  # U(B) - U(A)
    "vector",  # L = r ^ p
    "invalid",  # unparseable or degenerate (e.g. operands lost in extraction)
)

FUNCTIONS: dict[str, Any] = {
    "ln": sp.log,
    "log": sp.log,
    "log10": lambda x: sp.log(x, 10),
    "exp": sp.exp,
    "sqrt": sp.sqrt,
    "abs": sp.Abs,
    "Abs": sp.Abs,
    "sin": sp.sin,
    "cos": sp.cos,
    "tan": sp.tan,
    "asin": sp.asin,
    "acos": sp.acos,
    "atan": sp.atan,
    "sinh": sp.sinh,
    "cosh": sp.cosh,
    "tanh": sp.tanh,
}
CONSTANTS: dict[str, Any] = {"pi": sp.pi}

_IDENT_RE = re.compile(r"[A-Za-z_]\w*")
_INEQUALITY_RE = re.compile(r"<=|>=|<|>")
# dP / dt, d2x / dt2, dL/dt, d r / dt -- but not delta_P / delta_h
_DERIVATIVE_RE = re.compile(r"\bd\s?\d?[A-Za-z](?:_\w+)?\s*/\s*d\s?[A-Za-z]\d?\b")


class ParseError(ValueError):
    pass


@dataclass
class Variable:
    name: str
    desc: str = ""
    unit: str = ""
    # Optional physical range. Used to bias sampling and to seed the numeric
    # root search; never used to discard roots from a solve.
    min: float | None = None
    max: float | None = None

    @property
    def hint(self) -> tuple[float | None, float | None] | None:
        if self.min is None and self.max is None:
            return None
        return (self.min, self.max)


@dataclass
class Relation:
    """Result of parsing an equation's text."""

    kind: str
    lhs: sp.Expr | None = None
    rhs: sp.Expr | None = None
    error: str = ""


@dataclass
class Equation:
    id: str
    text: str
    title: str = ""
    chapter: str = ""
    variables: dict[str, Variable] = field(default_factory=dict)
    kind_override: str | None = None

    # ── parsing ──────────────────────────────────────────────────────────

    @cached_property
    def relation(self) -> Relation:
        return parse_relation(self.text, self.kind_override)

    @property
    def kind(self) -> str:
        return self.relation.kind

    @property
    def is_solvable(self) -> bool:
        return self.kind == ALGEBRAIC

    @cached_property
    def symbols(self) -> tuple[str, ...]:
        """Variable names in order of first appearance in the text."""
        rel = self.relation
        if rel.lhs is None or rel.rhs is None:
            return ()
        free = {s.name for s in (rel.lhs - rel.rhs).free_symbols}
        seen: list[str] = []
        for name in _IDENT_RE.findall(self.text):
            if name in free and name not in seen:
                seen.append(name)
        return tuple(seen)

    @property
    def method_name(self) -> str:
        """Legacy-style attribute name, e.g. ``eqn_2_1``."""
        return "eqn_" + re.sub(r"\W", "_", self.id)

    def variable(self, name: str) -> Variable:
        return self.variables.get(name) or Variable(name)

    # ── solving ──────────────────────────────────────────────────────────

    def solver(self, target: str) -> TargetSolver:
        from .solve import TargetSolver

        cache: dict[str, TargetSolver] = self.__dict__.setdefault("_solvers", {})
        if target not in cache:
            if not self.is_solvable:
                raise ParseError(
                    f"Equation {self.id} is {self.kind!r}, not solvable: {self.relation.error or self.text}"
                )
            if target not in self.symbols:
                raise KeyError(f"{target!r} is not a variable of {self.id}: {self.symbols}")
            rel = self.relation
            assert rel.lhs is not None and rel.rhs is not None
            cache[target] = TargetSolver(
                rel.lhs,
                rel.rhs,
                target,
                self.symbols,
                hint=self.variable(target).hint,
            )
        return cache[target]

    def solve_for(
        self,
        target: str,
        knowns: dict[str, float] | None = None,
        /,
        *,
        select: str = "unique",
        where: Callable[[float], bool] | None = None,
        bounds: tuple[float | None, float | None] | None = None,
        complex: bool = False,
        **kw: float,
    ) -> Any:
        """Solve for ``target`` given every other variable.

        ``select``: ``"unique"`` (return the root, raise if there are 0 or
        several), ``"all"`` (return the :class:`Solution`), ``"min"`` or
        ``"max"``. ``where`` and ``bounds`` filter roots before selecting;
        they are how a caller applies physical context.
        """
        from .solve import select_root

        values = {**(knowns or {}), **kw}
        missing = [s for s in self.symbols if s != target and s not in values]
        if missing:
            raise TypeError(f"{self.id}: missing values for {missing}")
        extra = set(values) - set(self.symbols)
        if extra:
            raise TypeError(f"{self.id}: unknown variables {sorted(extra)}")
        sol = self.solver(target)(values, complex=complex)
        return select_root(sol, select=select, where=where, bounds=bounds)

    def solve(
        self,
        knowns: dict[str, float] | None = None,
        /,
        *,
        select: str = "unique",
        where: Callable[[float], bool] | None = None,
        bounds: tuple[float | None, float | None] | None = None,
        complex: bool = False,
        **kw: float,
    ) -> Any:
        """kwasak-style: pass every variable but one; solve for that one."""
        values = {**(knowns or {}), **kw}
        missing = [s for s in self.symbols if values.get(s) is None]
        if len(missing) != 1:
            raise TypeError(f"{self.id}: pass all but exactly one of {self.symbols}; missing {missing}")
        target = missing[0]
        values.pop(target, None)
        return self.solve_for(target, values, select=select, where=where, bounds=bounds, complex=complex)

    __call__ = solve


# ── relation parsing ─────────────────────────────────────────────────────────


def _rationalize(expr: sp.Expr) -> sp.Expr:
    """Float exponents -> Rationals, and the literal 3.14159... -> pi.

    ``(P_2/P_1)**0.286`` is trivially invertible as ``**(143/500)`` but
    SymPy often gives up on the float form.
    """
    floats = {f: sp.pi for f in expr.atoms(sp.Float) if abs(float(f) - math.pi) < 1e-12}
    if floats:
        expr = expr.xreplace(floats)
    return expr.replace(
        lambda e: e.is_Pow and e.exp.is_Float,
        lambda e: sp.Pow(e.base, sp.nsimplify(e.exp, rational=True)),
    )


def _parse_side(text: str) -> sp.Expr:
    local: dict[str, Any] = {}
    for name in _IDENT_RE.findall(text):
        if name in FUNCTIONS:
            local[name] = FUNCTIONS[name]
        elif name in CONSTANTS:
            local[name] = CONSTANTS[name]
        else:
            # Everything else is a variable -- including E, I, S, N, beta ...
            # which SymPy would otherwise treat as built-ins.
            local[name] = sp.Symbol(name)
    return parse_expr(
        text,
        local_dict=local,
        global_dict={"Integer": sp.Integer, "Float": sp.Float, "Rational": sp.Rational, "Symbol": sp.Symbol},
        transformations=standard_transformations,
    )


def parse_relation(text: str, kind_override: str | None = None) -> Relation:
    body = text.split("#", 1)[0].strip()
    body = body.replace("math.", "").replace("~=", "=")

    if kind_override and kind_override != ALGEBRAIC:
        return Relation(kind_override, error="kind set explicitly")
    if "..." in body:
        return Relation("aggregate", error="open-ended sum ('...')")
    if _INEQUALITY_RE.search(body):
        return Relation("inequality", error="inequality / validity range")
    if "^" in body:
        return Relation("vector", error="'^' (cross product) is not scalar algebra")
    if _DERIVATIVE_RE.search(body):
        return Relation("differential", error="derivative notation (d?/d?)")

    parts = re.split(r"(?<![=!<>])=(?!=)", body)
    if len(parts) != 2:
        return Relation("invalid", error=f"expected exactly one '=' in {body!r}")
    try:
        lhs, rhs = (_rationalize(_parse_side(p.strip())) for p in parts)
    except TypeError as exc:
        if "not callable" in str(exc):
            return Relation("functional", error=f"unknown function call: {exc}")
        return Relation("invalid", error=f"parse error: {exc}")
    except (SyntaxError, TokenError, ValueError, AttributeError) as exc:
        return Relation("invalid", error=f"parse error: {exc}")

    written = {n for n in _IDENT_RE.findall(body) if n not in FUNCTIONS and n not in CONSTANTS}
    present = {s.name for s in (lhs - rhs).free_symbols}
    if not present:
        return Relation("invalid", error="no variables")
    lost = sorted(written - present)
    if lost:
        # e.g. "x = 0 * (x + v*t)" after gamma was lost in extraction
        return Relation("invalid", lhs, rhs, error=f"variables cancel out: {lost}")
    return Relation(ALGEBRAIC, lhs, rhs)


# ── chapters & libraries ─────────────────────────────────────────────────────


@dataclass
class Chapter:
    name: str  # class-style name, e.g. FluidFlowVacuumLines
    title: str = ""
    number: int | None = None
    source: str = ""
    equations: list[Equation] = field(default_factory=list)

    def __getattr__(self, attr: str) -> Equation:
        # lib.FluidFlowVacuumLines.eqn_2_1(Re=..., rho=..., ...)
        if attr.startswith("eqn_"):
            for eq in self.equations:
                if eq.method_name == attr:
                    return eq
        raise AttributeError(attr)


class Library:
    def __init__(self, chapters: list[Chapter], root: Path | None = None):
        self.chapters = chapters
        self.root = root
        self._by_id = {eq.id: eq for ch in chapters for eq in ch.equations}

    @classmethod
    def load(cls, project_dir: str | Path) -> Library:
        root = Path(project_dir)
        eq_dir, notes_dir = root / "equations", root / "notes"
        if eq_dir.is_dir() and any(eq_dir.glob("*.toml")):
            chapters = [load_toml(p) for p in sorted(eq_dir.glob("*.toml"))]
        elif notes_dir.is_dir():
            chapters = [parse_notes(p) for p in sorted(notes_dir.glob("*.py")) if not p.name.startswith("_")]
        else:
            raise FileNotFoundError(f"No equations/ or notes/ directory in {root}")
        return cls(chapters, root)

    @property
    def equations(self) -> list[Equation]:
        return [eq for ch in self.chapters for eq in ch.equations]

    def __getitem__(self, eq_id: str) -> Equation:
        return self._by_id[eq_id]

    def __getattr__(self, name: str) -> Chapter:
        for ch in self.__dict__.get("chapters", []):
            if ch.name == name:
                return ch
        raise AttributeError(name)


# ── legacy notes DSL ─────────────────────────────────────────────────────────

_TAG_RE = re.compile(r"#\s*(?:Equation\s+)?(\d{1,2}-(?:\d{1,2}[a-z]?|\?|[A-Za-z]\w*))\b[\s,:.-]*(.*)")
_CHAPTER_RE = re.compile(r"#\s*Chapter\s*(\d+)\s*[:.-]?\s*(.*)", re.I)
_VARNOTE_RE = re.compile(r"^([A-Za-z_]\w*)\s*:=\s*(.*)$")
_NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"
_RANGE_RE = re.compile(rf"^({_NUM})\s*<=?\s*([A-Za-z_]\w*)\s*<=?\s*({_NUM})")


def _class_name(stem: str) -> str:
    parts = stem.split("_")
    if parts and parts[0].isdigit():
        parts = parts[1:]
    return "".join(p.capitalize() for p in parts)


def parse_notes(path: str | Path) -> Chapter:
    """Parse a legacy ``notes/NN_chapter.py`` file.

    Lines are classified one at a time rather than by tracking docstring
    state, because real notes files have unbalanced triple quotes and
    equations inside docstrings. Fixes over the old regex parser: tags are
    only read from comment lines (a ``10-3`` inside an equation no longer
    renumbers it), ``~=`` is an equality, chained ``a = b = c`` is split,
    indented code (pasted ``def`` bodies) is skipped, and range lines like
    ``1 <= P <= 10`` become hints for ``P``.
    """
    path = Path(path)
    stem = path.stem
    number = int(stem.split("_")[0]) if stem.split("_")[0].isdigit() else None
    chapter = Chapter(name=_class_name(stem), number=number, source=str(path))

    cur_id, cur_title = "", ""
    notes: dict[str, Variable] = {}
    block: list[Equation] = []
    anon = 0
    count_for_id: dict[str, int] = {}

    def flush() -> None:
        # Variable notes may come before or after the equation in a block.
        for eq in block:
            used = set(_IDENT_RE.findall(eq.text))
            eq.variables = {k: Variable(**vars(v)) for k, v in notes.items() if k in used}
        block.clear()

    def add(text: str) -> None:
        nonlocal anon
        eq_id = cur_id
        if not eq_id:
            anon += 1
            eq_id = f"{number or 0}-u{anon}"
        n = count_for_id.get(eq_id, 0) + 1
        count_for_id[eq_id] = n
        if n > 1:
            eq_id = f"{eq_id}.{n}"
        eq = Equation(id=eq_id, text=text, title=cur_title, chapter=chapter.name)
        chapter.equations.append(eq)
        block.append(eq)

    for raw in path.read_text().splitlines():
        if raw[:1].isspace() or raw.startswith(("def ", "import ", "from ", "return")):
            continue
        s = raw.replace('"""', "").strip()
        if not s:
            continue
        if s.startswith("#"):
            if m := _CHAPTER_RE.match(s):
                chapter.title = m.group(2).strip()
            elif m := _TAG_RE.match(s):
                flush()
                notes = {}
                cur_id, cur_title = m.group(1), m.group(2).strip().strip(",")
            continue
        if m := _VARNOTE_RE.match(s):
            name, desc = m.group(1), m.group(2).strip()
            notes.setdefault(name, Variable(name)).desc = desc
            continue
        if m := _RANGE_RE.match(s):
            lo, name, hi = float(m.group(1)), m.group(2), float(m.group(3))
            v = notes.setdefault(name, Variable(name))
            v.min, v.max = lo, hi
            continue
        body = s.split("#", 1)[0]
        if "=" not in body or ":=" in body:
            continue
        # chained a = b = c  ->  a = b, b = c
        pieces = re.split(r"(?<![=!<>~])=(?!=)", body)
        if len(pieces) > 2 and not _INEQUALITY_RE.search(body):
            for a, b in zip(pieces, pieces[1:]):
                add(f"{a.strip()} = {b.strip()}")
        else:
            add(s)
    flush()

    # Make ids unique across the chapter (e.g. two "8-?" tags).
    seen: dict[str, int] = {}
    for eq in chapter.equations:
        if eq.id in seen:
            seen[eq.id] += 1
            eq.id = f"{eq.id}#{seen[eq.id]}"
        else:
            seen[eq.id] = 1
    return chapter


# ── TOML format ──────────────────────────────────────────────────────────────


def _toml_str(s: str) -> str:
    return json.dumps(s, ensure_ascii=False)


def dump_toml(chapter: Chapter) -> str:
    out = []
    if chapter.source:
        out.append(f"# Converted from {Path(chapter.source).name} by `vakyume convert-notes`.")
    out.append(f"chapter = {_toml_str(chapter.name)}")
    if chapter.title:
        out.append(f"title = {_toml_str(chapter.title)}")
    if chapter.number is not None:
        out.append(f"number = {chapter.number}")
    for eq in chapter.equations:
        out += ["", "[[equation]]", f"id = {_toml_str(eq.id)}"]
        if eq.title:
            out.append(f"title = {_toml_str(eq.title)}")
        out.append(f"eq = {_toml_str(eq.text)}")
        if eq.kind_override:
            out.append(f"kind = {_toml_str(eq.kind_override)}")
        if eq.variables:
            out.append("[equation.variables]")
            for name, v in eq.variables.items():
                fields = []
                if v.desc:
                    fields.append(f"desc = {_toml_str(v.desc)}")
                if v.unit:
                    fields.append(f"unit = {_toml_str(v.unit)}")
                if v.min is not None:
                    fields.append(f"min = {v.min!r}")
                if v.max is not None:
                    fields.append(f"max = {v.max!r}")
                out.append(f"{name} = {{ {', '.join(fields)} }}")
    return "\n".join(out) + "\n"


def load_toml(path: str | Path) -> Chapter:
    path = Path(path)
    data = tomllib.loads(path.read_text())
    chapter = Chapter(
        name=data.get("chapter") or _class_name(path.stem),
        title=data.get("title", ""),
        number=data.get("number"),
        source=str(path),
    )
    for rec in data.get("equation", []):
        variables = {name: Variable(name, **fields) for name, fields in (rec.get("variables") or {}).items()}
        chapter.equations.append(
            Equation(
                id=str(rec["id"]),
                text=rec["eq"],
                title=rec.get("title", ""),
                chapter=chapter.name,
                variables=variables,
                kind_override=rec.get("kind"),
            )
        )
    return chapter
