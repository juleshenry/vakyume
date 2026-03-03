"""Source-code extraction utilities.

Helper functions for extracting function and method source code from Python
files, used by the pipeline to assemble certified libraries.
"""

import re


def get_standalone_method_source(file_path: str, method_name: str) -> str:
    """Extract the source of a top-level function from *file_path*."""
    with open(file_path, "r") as f:
        lines = f.readlines()

    target_idx = -1
    for idx, line in enumerate(lines):
        if line.startswith(f"def {method_name}("):
            target_idx = idx
            break

    if target_idx == -1:
        return ""

    start_idx = target_idx
    while start_idx > 0 and lines[start_idx - 1].strip().startswith("@"):
        start_idx -= 1

    end_idx = target_idx + 1
    while end_idx < len(lines):
        if lines[end_idx].strip():
            if not lines[end_idx].startswith(" "):
                stripped = lines[end_idx].lstrip()
                if stripped.startswith(("def ", "class ", "@")):
                    break
        end_idx += 1
    return "".join(lines[start_idx:end_idx])


def get_method_source_from_class(
    file_path: str,
    class_name: str,
    method_name: str,
) -> str:
    """Extract a method's source from a class definition in *file_path*."""
    with open(file_path, "r") as f:
        lines = f.readlines()

    class_idx = -1
    for idx, line in enumerate(lines):
        if line.startswith(f"class {class_name}:") or line.startswith(
            f"class {class_name}("
        ):
            class_idx = idx
            break

    if class_idx == -1:
        return ""

    target_idx = -1
    for idx in range(class_idx + 1, len(lines)):
        line = lines[idx]
        if re.search(rf"def\s+{method_name}\s*\(", line):
            target_idx = idx
            break
        if line.strip() and not line.startswith(" "):
            break

    if target_idx == -1:
        return ""

    start_idx = target_idx
    while start_idx > class_idx + 1 and lines[start_idx - 1].strip().startswith("@"):
        start_idx -= 1

    line_with_indent = lines[target_idx]
    indent_len = len(line_with_indent) - len(line_with_indent.lstrip())

    end_idx = target_idx + 1
    while end_idx < len(lines):
        line = lines[end_idx]
        if line.strip():
            curr_indent = len(line) - len(line.lstrip())
            if curr_indent <= indent_len:
                stripped = line.lstrip()
                if stripped.startswith(("def ", "class ", "@")):
                    break
                if curr_indent == 0:
                    break
        end_idx += 1

    return "".join(lines[start_idx:end_idx])
