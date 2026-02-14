import sys
import os
import re
import ast


def extract_method_blocks(code_text):
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
                # Include preceding decorators
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


def test_granular():
    code = """
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
    methods = extract_method_blocks(code)
    print(f"Found methods: {list(methods.keys())}")

    imports = "from math import log\n"

    for name, block in methods.items():
        # De-indent
        lines = block.splitlines()
        if not lines:
            continue
        min_indent = min(len(l) - len(l.lstrip()) for l in lines if l.strip())
        deindented = "\n".join(l[min_indent:] for l in lines)

        # Remove @staticmethod if we want to exec as top-level
        deindented = deindented.replace("@staticmethod", "")

        full_code = imports + deindented
        print(f"--- Executing {name} ---")
        try:
            locs = {}
            exec(full_code, {}, locs)
            func = locs.get(name)
            if func:
                print(f"Success! {name}(10) = {func(10)}")
            else:
                print(f"Failed: {name} not found in locals after exec")
        except Exception as e:
            print(f"Failed: {type(e).__name__}: {e}")


if __name__ == "__main__":
    test_granular()
