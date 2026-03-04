"""Post-scrape validation for Vakyume equation notes.

Validates that scraped note files contain well-formed equations that the
parser/shard-generator can actually handle.  Catches the common failure
modes produced by LLM-based PDF extraction:

  * Missing operands        (``F = * x``)
  * Truncated expressions   (``E = U +``)
  * Unrecognised functions  (``sqrt``, ``U(B)``)
  * Derivative / calculus   (``dL / dt``, ``d2x / dt2``)
  * Cross-product / XOR     (``L = r ^ p``)
  * Concatenated variables  (``mgH`` instead of ``m * g * H``)
  * Invalid equation tags   (``# 8-?``, ``# 5-unknown``)
  * Too few solvable vars   (single-variable identities)

Each equation is classified as **valid** or **invalid** with a human-
readable reason.  The pipeline uses this to skip invalid equations during
shard generation and to print a summary so the user knows what to fix.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field

from .config import OPERATORS

# ---------------------------------------------------------------------------
#  Regex for the equation-number tag the parser expects
# ---------------------------------------------------------------------------
_EQN_TAG_RE = re.compile(r"\d{1,2}-\d{1,2}\w*")

# Operators that may legally appear at the start of the RHS / end of a side
_BINARY_OPS = set("*+/-")

# Functions the parser actually knows how to handle
_KNOWN_FUNCS = {"ln", "log"}

# Derivative-like tokens that signal calculus notation
_DERIV_TOKENS = re.compile(
    r"\bd\s*\d*[a-zA-Z]\s*/\s*d\s*t\s*\d*\b"  # d2x/dt2, dL/dt, dr/dt …
)

# Cross-product / XOR operator
_XOR_RE = re.compile(r"\s\^\s")

# Function-call syntax the parser can't handle (but not ln/log)
_FUNC_CALL_RE = re.compile(r"\b(?!ln\b|log\b)([A-Za-z_]\w*)\s*\(")

# Recognised math functions that aren't variables but also aren't solvable
# in the current parser (sqrt is imported via cmath but the *parser*
# tokeniser doesn't exclude it, so it becomes a variable).
_UNSUPPORTED_FUNCS = {"sqrt", "sin", "cos", "tan", "exp", "abs"}


# ---------------------------------------------------------------------------
#  Result types
# ---------------------------------------------------------------------------
@dataclass
class EquationResult:
    """Validation result for a single equation line."""

    file: str
    line_no: int
    eqn_tag: str
    raw_line: str
    valid: bool
    reasons: list[str] = field(default_factory=list)


@dataclass
class FileResult:
    """Validation results for an entire notes file."""

    path: str
    equations: list[EquationResult] = field(default_factory=list)

    @property
    def valid_count(self) -> int:
        return sum(1 for eq in self.equations if eq.valid)

    @property
    def invalid_count(self) -> int:
        return sum(1 for eq in self.equations if not eq.valid)


@dataclass
class ValidationReport:
    """Aggregate validation report across all notes files."""

    files: list[FileResult] = field(default_factory=list)

    @property
    def total_equations(self) -> int:
        return sum(len(f.equations) for f in self.files)

    @property
    def total_valid(self) -> int:
        return sum(f.valid_count for f in self.files)

    @property
    def total_invalid(self) -> int:
        return sum(f.invalid_count for f in self.files)

    def valid_tags_for_file(self, filepath: str) -> set[str]:
        """Return the set of valid equation tags in *filepath*."""
        for fr in self.files:
            if fr.path == filepath:
                return {eq.eqn_tag for eq in fr.equations if eq.valid}
        return set()


# ---------------------------------------------------------------------------
#  Per-equation checks
# ---------------------------------------------------------------------------


def _check_equation(
    raw_line: str,
    eqn_tag: str,
    file: str,
    line_no: int,
) -> EquationResult:
    """Run all checks on a single equation line and return a result."""
    reasons: list[str] = []
    stripped = raw_line.strip()

    # Remove inline comments
    code = stripped.split("#")[0].strip()

    # ── 1. Tag validity ──────────────────────────────────────────────────
    if not eqn_tag:
        reasons.append("missing equation tag")
    elif not _EQN_TAG_RE.fullmatch(eqn_tag.replace("_", "-")):
        reasons.append(f"invalid equation tag '{eqn_tag}' (must match N-N pattern)")

    # ── 2. Exactly one '=' ───────────────────────────────────────────────
    parts = code.split("=")
    if len(parts) != 2:
        reasons.append(f"expected exactly one '=' but found {len(parts) - 1}")
        return EquationResult(file, line_no, eqn_tag, stripped, False, reasons)

    lhs, rhs = parts[0].strip(), parts[1].strip()

    # ── 3. Empty sides ───────────────────────────────────────────────────
    if not lhs:
        reasons.append("empty left-hand side")
    if not rhs:
        reasons.append("empty right-hand side")
    if reasons:
        return EquationResult(file, line_no, eqn_tag, stripped, False, reasons)

    # ── 4. Truncated expression (ends/starts with binary op) ─────────────
    if rhs[-1] in _BINARY_OPS:
        reasons.append(f"truncated: RHS ends with '{rhs[-1]}'")
    if lhs[-1] in _BINARY_OPS:
        reasons.append(f"truncated: LHS ends with '{lhs[-1]}'")
    # Leading operator on RHS (excluding unary minus)
    if rhs[0] in (_BINARY_OPS - {"-"}):
        reasons.append(f"missing operand: RHS starts with '{rhs[0]}'")

    # ── 5. Missing operand around '*' (e.g. "= * x", "- * x") ───────────
    # Collapse whitespace and check for '*' preceded by an operator or '='
    # which indicates the left operand is missing.  We exclude '**' (power).
    compressed = re.sub(r"\s+", "", code)
    if re.search(r"(?<=[=+\-/(^])\*(?!\*)", compressed) or re.match(
        r"\*(?!\*)", compressed
    ):
        reasons.append("missing operand before '*'")

    # ── 6. Derivative / calculus notation ────────────────────────────────
    if _DERIV_TOKENS.search(code):
        reasons.append("derivative/calculus notation (not algebraic)")

    # ── 7. Cross-product / XOR ───────────────────────────────────────────
    if _XOR_RE.search(code):
        reasons.append("'^' operator (cross product / XOR, not algebraic)")

    # ── 8. Unsupported function calls ────────────────────────────────────
    for m in _FUNC_CALL_RE.finditer(code):
        fname = m.group(1)
        if fname in _UNSUPPORTED_FUNCS:
            reasons.append(f"unsupported function '{fname}()' (use ** 0.5 for sqrt)")
        elif not fname.startswith("eqn_"):
            # Arbitrary function call like U(B) — parser can't handle
            reasons.append(f"function-call syntax '{fname}()' not supported")

    # ── 9. Concatenated variables ────────────────────────────────────────
    # A token that's 3+ alpha chars with mixed case and no underscores
    # is likely concatenated vars (e.g. "mgH" = m*g*H).
    # Only flag if it's not a known identifier pattern.
    for token in re.findall(r"\b([a-zA-Z_]\w*)\b", code):
        if (
            len(token) >= 3
            and token not in _KNOWN_FUNCS
            and token not in _UNSUPPORTED_FUNCS
            and not token.startswith("_")
            and re.search(r"[a-z]", token)
            and re.search(r"[A-Z]", token)
            # Skip normal camelCase identifiers that could be variable names
            # Only flag when lowercase letters are mixed among uppercase in a
            # physics-equation-like way (e.g. "mgH" but not "v_0x")
            and not "_" in token
            and re.fullmatch(r"[a-z]+[A-Z][a-zA-Z]*", token)
        ):
            reasons.append(
                f"possible concatenated variables '{token}' (should use * between vars)"
            )

    # ── 10. Too few tokens (trivial identity like a = a_A) ───────────────
    # Extract solvable tokens the same way the parser does
    tokens = _extract_tokens(code)
    if len(tokens) < 2:
        reasons.append(
            f"only {len(tokens)} variable(s) — not enough for shard generation"
        )

    valid = len(reasons) == 0
    return EquationResult(file, line_no, eqn_tag, stripped, valid, reasons)


def _extract_tokens(eqn: str) -> list[str]:
    """Mimic ``Solver.get_tokens`` without instantiating the class."""
    excluded = {"ln", "log"}
    for m in excluded:
        eqn = f"{m}".join(o.strip() for o in eqn.split(m))

    eqn = eqn.split("#")[0]
    dilated = ""
    for i, ch in enumerate(eqn):
        if (
            0 < i < len(eqn) - 1
            and eqn[i + 1] != "*"
            and eqn[i - 1] != "*"
            and ch in OPERATORS
        ):
            dilated += f" {ch} "
        else:
            dilated += ch
    raw_tokens = dilated.split()

    tokens = set()
    for t in raw_tokens:
        clean = t.strip().replace("(", "").replace(")", "").split("**")[0].strip()
        if clean.isidentifier() and clean not in excluded:
            tokens.add(clean)
    return sorted(tokens)


# ---------------------------------------------------------------------------
#  File-level validation
# ---------------------------------------------------------------------------


def validate_notes_file(filepath: str) -> FileResult:
    """Validate every equation in a single notes file."""
    result = FileResult(path=filepath)
    basename = os.path.basename(filepath)

    with open(filepath, "r") as f:
        lines = f.readlines()

    eqn_tag = ""
    for line_no_0, line in enumerate(lines, start=1):
        stripped = line.strip()

        # Track the current equation tag
        if tag_match := re.findall(r"\d{1,2}-\d{1,2}\w*", line):
            eqn_tag = tag_match[0].replace("-", "_")

        # Detect equation lines (same logic as parser.py)
        if (
            " = " in line
            and not stripped.startswith("#")
            and not stripped.startswith('"""')
        ):
            eq_result = _check_equation(stripped, eqn_tag, basename, line_no_0)
            result.equations.append(eq_result)

    return result


# ---------------------------------------------------------------------------
#  Directory-level validation
# ---------------------------------------------------------------------------


def validate_notes_dir(notes_dir: str) -> ValidationReport:
    """Validate all notes files in a directory.  Returns a full report."""
    report = ValidationReport()

    if not os.path.isdir(notes_dir):
        return report

    for fname in sorted(os.listdir(notes_dir)):
        if fname.endswith(".py") and not fname.startswith("__"):
            fpath = os.path.join(notes_dir, fname)
            report.files.append(validate_notes_file(fpath))

    return report


# ---------------------------------------------------------------------------
#  Pretty-printed summary (for terminal output)
# ---------------------------------------------------------------------------


def print_validation_report(report: ValidationReport) -> None:
    """Print a human-readable validation summary to stdout."""
    if report.total_equations == 0:
        print("  No equations found in notes files.")
        return

    print(
        f"\n  Notes validation: {report.total_valid} valid, "
        f"{report.total_invalid} invalid "
        f"(out of {report.total_equations} equations)"
    )

    if report.total_invalid == 0:
        print("  All equations are well-formed — ready for shard generation.")
        return

    print()
    for fr in report.files:
        if fr.invalid_count == 0:
            continue
        print(f"  {os.path.basename(fr.path)}:")
        for eq in fr.equations:
            if eq.valid:
                continue
            tag_str = eq.eqn_tag if eq.eqn_tag else "???"
            print(f"    [{tag_str}] L{eq.line_no}: {eq.raw_line}")
            for reason in eq.reasons:
                print(f"           -> {reason}")
        print()

    print(
        f"  {report.total_invalid} equation(s) will be SKIPPED during shard generation."
    )
    print("  Fix the notes files above and re-run to include them.\n")
