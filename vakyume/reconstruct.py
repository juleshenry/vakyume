"""
Shard Reconstruction CLI

Rebuilds chapter-libraries (shards/) back into a coherent, importable Python
module so that downstream code can do:

    from vakyume import LiquidRing
    pump = LiquidRing()
    pump.eqn_10_1(D_r=5.0, w=1800)   # -> solves for sig_R
"""

import ast
import os
import re
import sys
import textwrap

from .config import TAB, LIBRARY_IMPORT_HEADER


# ── Helpers ──────────────────────────────────────────────────────────────────

IMPORT_HEADER = LIBRARY_IMPORT_HEADER


def _family_class_name(family_dir_name: str) -> str:
    """LiquidRing_eqn_10_1 -> LiquidRing"""
    return family_dir_name.split("_eqn_")[0]


def _extract_standalone_function(source: str, func_name: str) -> str | None:
    """Return the source text of a top-level `def func_name(...)` from *source*."""
    try:
        tree = ast.parse(source)
    except SyntaxError:
        return None

    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name == func_name:
            lines = source.splitlines(keepends=True)
            # walk backwards from node line to pick up decorators
            start = node.lineno - 1
            while start > 0 and lines[start - 1].strip().startswith("@"):
                start -= 1
            end = node.end_lineno
            return "".join(lines[start:end])
    return None


def _extract_class_methods(source: str, class_name: str) -> list[tuple[str, str]]:
    """Return [(method_name, method_source), ...] for every *method* in *class_name*.

    Attribute assignments (e.g. ``eqn_10_1__D_r = eqn_10_1__D_r``) are
    intentionally skipped — in shard files they wire up imports that are
    irrelevant after reconstruction.
    """
    try:
        tree = ast.parse(source)
    except SyntaxError:
        return []

    results = []
    for node in tree.body:
        if isinstance(node, ast.ClassDef) and node.name == class_name:
            lines = source.splitlines(keepends=True)
            for item in node.body:
                if isinstance(item, ast.FunctionDef):
                    start = item.lineno - 1
                    # walk backwards for decorators (within class body)
                    while start > node.lineno and lines[start - 1].strip().startswith(
                        "@"
                    ):
                        start -= 1
                    end = item.end_lineno
                    method_src = "".join(lines[start:end])
                    results.append((item.name, method_src))
                # Assign nodes (e.g. eqn_X_Y__var = eqn_X_Y__var) are shard-
                # import wiring and are skipped; the actual function bodies
                # are picked up from the standalone solver shard files.
    return results


# ── Core reconstruction ─────────────────────────────────────────────────────


def reconstruct_from_shards(shards_dir: str) -> str:
    """
    Scan *shards_dir* and reconstruct a single Python module containing all
    chapter classes with their equation methods.

    Returns the generated source code as a string.
    """
    if not os.path.isdir(shards_dir):
        raise FileNotFoundError(f"Shards directory not found: {shards_dir}")

    # class_name -> { method_name -> method_source }
    class_methods: dict[str, dict[str, str]] = {}
    # track ordering of classes by first-seen chapter number
    class_order: dict[str, str] = {}

    family_dirs = sorted(
        d
        for d in os.listdir(shards_dir)
        if os.path.isdir(os.path.join(shards_dir, d)) and "_eqn_" in d
    )

    for family_dir_name in family_dirs:
        class_name = _family_class_name(family_dir_name)
        if class_name not in class_methods:
            class_methods[class_name] = {}
            class_order[class_name] = family_dir_name

        family_path = os.path.join(shards_dir, family_dir_name)
        seen = class_methods[class_name]

        # Process all .py files in the family directory (sorted so main shard
        # comes before the individual solver files alphabetically, but we
        # handle both patterns regardless of order).
        py_files = sorted(
            f
            for f in os.listdir(family_path)
            if f.endswith(".py") and not f.startswith("__")
        )

        for py_file in py_files:
            file_path = os.path.join(family_path, py_file)
            with open(file_path, "r") as fh:
                content = fh.read()

            try:
                tree = ast.parse(content)
            except SyntaxError:
                continue

            for node in tree.body:
                if isinstance(node, ast.FunctionDef):
                    # Top-level function (individual solver shard)
                    if node.name not in seen:
                        src = _extract_standalone_function(content, node.name)
                        if src:
                            seen[node.name] = src

                elif isinstance(node, ast.ClassDef) and node.name == class_name:
                    # Main shard file — extract methods from class body
                    for method_name, method_src in _extract_class_methods(
                        content, class_name
                    ):
                        if method_name not in seen:
                            seen[method_name] = method_src

    # ── Assemble output ──────────────────────────────────────────────────

    out_lines = [IMPORT_HEADER, ""]

    # Sort classes by their first-seen family name (gives chapter ordering)
    sorted_classes = sorted(class_methods.keys(), key=lambda c: class_order[c])

    for class_name in sorted_classes:
        methods = class_methods[class_name]
        if not methods:
            continue

        out_lines.append(f"class {class_name}:")

        # Separate kwasak dispatchers, attribute assignments, and solver methods.
        # We want a clean ordering:
        #   1. @kwasak dispatcher + attribute assignments for that equation
        #   2. solver methods (__var)
        # grouped by equation number.

        # Group methods by equation family
        eqn_groups: dict[str, dict[str, str]] = {}
        for method_name, method_src in methods.items():
            # Extract equation identifier: eqn_X_Y from eqn_X_Y or eqn_X_Y__var
            m = re.match(r"(eqn_\d+_\d+\w*?)(?:__.*)?$", method_name)
            if m:
                eqn_key = m.group(1)
            else:
                eqn_key = "__other__"
            if eqn_key not in eqn_groups:
                eqn_groups[eqn_key] = {}
            eqn_groups[eqn_key][method_name] = method_src

        for eqn_key in sorted(eqn_groups.keys()):
            group = eqn_groups[eqn_key]

            # Emit the @kwasak dispatcher first (the bare eqn_X_Y method)
            if eqn_key in group:
                src = group[eqn_key]
                out_lines.append(_indent_method(src))

            # Then solver methods (eqn_X_Y__var), sorted
            for mname in sorted(group.keys()):
                if mname == eqn_key:
                    continue
                src = group[mname]
                out_lines.append(_indent_method(src))

        out_lines.append("")  # blank line between classes

    return _cleanup_source("\n".join(out_lines))


def _cleanup_source(source: str) -> str:
    """Post-process the assembled source to strip reconstruction artifacts.

    • ``# [.pyeqn] sig_R = ...`` → ``# sig_R = ...``  (keep equation, drop tag)
    """
    source = re.sub(r"# \[\.pyeqn\] ", "# ", source)
    return source


def _indent_method(source: str) -> str:
    """
    Ensure *source* is indented as a class method (one level of indentation).

    The source may already be indented (from a class body) or unindented (from
    a standalone function).  We normalize to exactly one TAB level.
    """
    lines = source.splitlines()
    if not lines:
        return ""

    # Detect current indentation of the first non-blank line
    first_line = lines[0]
    current_indent = len(first_line) - len(first_line.lstrip())

    result = []
    for line in lines:
        if line.strip():
            # Remove existing indent, add one TAB
            stripped = (
                line[current_indent:] if len(line) >= current_indent else line.lstrip()
            )
            result.append(TAB + stripped)
        else:
            result.append("")
    return "\n".join(result)


# ── CLI entry point ─────────────────────────────────────────────────────────


def reconstruct_cli(project_dir: str, output: str | None = None, stdout: bool = False):
    """
    Run reconstruction for a project.  Writes the output module to
    ``<project_dir>/vakyume_lib.py`` by default.
    """
    shards_dir = os.path.join(project_dir, "shards")
    source = reconstruct_from_shards(shards_dir)

    if stdout:
        print(source)
        return

    out_path = output or os.path.join(project_dir, "vakyume_lib.py")
    with open(out_path, "w") as f:
        f.write(source)

    # Count classes and methods for summary
    try:
        tree = ast.parse(source)
        n_classes = sum(1 for n in tree.body if isinstance(n, ast.ClassDef))
        n_methods = sum(
            sum(1 for item in n.body if isinstance(item, ast.FunctionDef))
            for n in tree.body
            if isinstance(n, ast.ClassDef)
        )
    except SyntaxError:
        n_classes = n_methods = "?"

    print(f"Reconstructed library written to {out_path}")
    print(f"  {n_classes} classes, {n_methods} methods")
