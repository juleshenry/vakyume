"""
Unit test for method-block extraction from Python source text.

This is a pure-Python utility test -- no external dependencies.
"""

import re


def extract_method_blocks(code_text):
    """Extract ``def eqn_*`` blocks from *code_text*, including decorators.

    Returns a dict mapping method name -> source text (with original indentation).
    """
    lines = code_text.splitlines()
    blocks = {}
    i = 0
    while i < len(lines):
        line = lines[i]
        stripped = line.lstrip()
        if stripped.startswith("def "):
            match = re.search(r"def\s+(eqn_\w+)\s*\(", stripped)
            if match:
                name = match.group(1)
                start = i
                # Walk backwards to include preceding decorators
                while start > 0 and (
                    lines[start - 1].lstrip().startswith("@")
                    or not lines[start - 1].strip()
                ):
                    if lines[start - 1].lstrip().startswith("@"):
                        start -= 1
                    else:
                        j = start - 1
                        found_decorator = False
                        while j >= 0:
                            if lines[j].lstrip().startswith("@"):
                                found_decorator = True
                                break
                            if lines[j].strip():
                                break
                            j -= 1
                        if found_decorator:
                            start = j
                        else:
                            break

                indent = len(line) - len(line.lstrip())
                end = i + 1
                while end < len(lines):
                    if lines[end].strip():
                        curr_indent = len(lines[end]) - len(lines[end].lstrip())
                        if curr_indent <= indent and not lines[end].lstrip().startswith(
                            ")"
                        ):
                            if lines[end].lstrip().startswith("def ") or lines[
                                end
                            ].lstrip().startswith("@"):
                                break
                    end += 1

                block_lines = lines[start:end]
                blocks[name] = "\n".join(block_lines)
                i = end
                continue
        i += 1
    return blocks


# ── Tests ───────────────────────────────────────────────────────────────────

SAMPLE_CODE = """\
from math import log
class Test:
    @staticmethod
    def eqn_1__x(y):
        return [y + 1]

    @staticmethod
    def eqn_1__y(x):
        return [x - 1]

    @staticmethod
    def eqn_1__broken(z):
        this is a syntax error
"""


def test_extract_finds_all_methods():
    methods = extract_method_blocks(SAMPLE_CODE)
    assert set(methods.keys()) == {"eqn_1__x", "eqn_1__y", "eqn_1__broken"}


def test_extracted_method_is_executable():
    methods = extract_method_blocks(SAMPLE_CODE)
    block = methods["eqn_1__x"]
    # De-indent to top level
    lines = block.splitlines()
    min_indent = min(len(l) - len(l.lstrip()) for l in lines if l.strip())
    deindented = "\n".join(l[min_indent:] for l in lines)
    deindented = deindented.replace("@staticmethod\n", "")

    locs = {}
    exec(deindented, {}, locs)
    func = locs["eqn_1__x"]
    assert func(10) == [11]


def test_broken_method_fails_exec():
    methods = extract_method_blocks(SAMPLE_CODE)
    block = methods["eqn_1__broken"]
    lines = block.splitlines()
    min_indent = min(len(l) - len(l.lstrip()) for l in lines if l.strip())
    deindented = "\n".join(l[min_indent:] for l in lines)
    deindented = deindented.replace("@staticmethod\n", "")

    with __import__("pytest").raises(SyntaxError):
        exec(deindented, {}, {})
