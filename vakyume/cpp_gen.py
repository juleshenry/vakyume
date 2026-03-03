"""
C++ code generator for Vakyume solver libraries.

Two-tier strategy:
  1. **Deterministic AST transpiler** — walks the Python AST of each solver
     method and emits equivalent C++.  Handles arithmetic, ``log``, ``sqrt``,
     ``exp``, ``pow``, list construction, and the ``result.append(…)`` pattern
     used by every Vakyume shard.
  2. **LLM fallback** (Phi-3 via Ollama) — used only when the AST transpiler
     cannot handle a construct (e.g. SymPy symbolic calls).
"""

import ast
import subprocess
import os
import re
import textwrap
from concurrent.futures import ThreadPoolExecutor


# ── AST helpers ──────────────────────────────────────────────────────────────

# Python names -> C++ equivalents
_MATH_FUNCS = {
    "log": "std::log",
    "sqrt": "std::sqrt",
    "exp": "std::exp",
    "pow": "std::pow",
    "abs": "std::abs",
}

_MATH_CONSTS = {
    "e": "M_E",
    "pi": "M_PI",
}


class _CppTranspileError(Exception):
    """Raised when the deterministic transpiler cannot handle a construct."""


class _ExprTranspiler(ast.NodeVisitor):
    """Convert a Python expression AST node to a C++ expression string."""

    def visit_Constant(self, node):
        if isinstance(node.value, (int, float)):
            v = repr(node.value)
            # Ensure floating-point literal
            if isinstance(node.value, int):
                v += ".0"
            return v
        raise _CppTranspileError(f"unsupported constant: {node.value!r}")

    def visit_Name(self, node):
        if node.id in _MATH_CONSTS:
            return _MATH_CONSTS[node.id]
        return node.id

    def visit_UnaryOp(self, node):
        operand = self.visit(node.operand)
        if isinstance(node.op, ast.USub):
            return f"(-{operand})"
        if isinstance(node.op, ast.UAdd):
            return f"(+{operand})"
        raise _CppTranspileError(f"unsupported unary op: {type(node.op).__name__}")

    def visit_BinOp(self, node):
        left = self.visit(node.left)
        right = self.visit(node.right)
        ops = {
            ast.Add: "+",
            ast.Sub: "-",
            ast.Mult: "*",
            ast.Div: "/",
        }
        op_type = type(node.op)
        if op_type == ast.Pow:
            return f"std::pow({left}, {right})"
        if op_type in ops:
            return f"({left} {ops[op_type]} {right})"
        raise _CppTranspileError(f"unsupported binop: {op_type.__name__}")

    def visit_Call(self, node):
        # Plain function call: log(x), sqrt(x), etc.
        if isinstance(node.func, ast.Name):
            fname = node.func.id
            if fname in _MATH_FUNCS:
                args = ", ".join(self.visit(a) for a in node.args)
                return f"{_MATH_FUNCS[fname]}({args})"
            raise _CppTranspileError(f"unsupported function: {fname}")

        # Attribute call: e.g. result.append(x) — handled at statement level
        raise _CppTranspileError("unsupported call expression")

    def visit_Attribute(self, node):
        # e.g. math.log — we treat the attribute name directly
        if isinstance(node.value, ast.Name):
            combined = f"{node.value.id}.{node.attr}"
            # math.log -> std::log, etc.
            if node.attr in _MATH_FUNCS:
                return _MATH_FUNCS[node.attr]
            if node.attr in _MATH_CONSTS:
                return _MATH_CONSTS[node.attr]
        raise _CppTranspileError(f"unsupported attribute: {ast.dump(node)}")

    def visit_List(self, node):
        # [expr, ...] — used in result = [val]
        elts = ", ".join(self.visit(e) for e in node.elts)
        return f"{{{elts}}}"

    def visit_Subscript(self, node):
        value = self.visit(node.value)
        sl = self.visit(node.slice)
        return f"{value}[{sl}]"

    # fallback
    def generic_visit(self, node):
        raise _CppTranspileError(f"unsupported AST node: {type(node).__name__}")


def _transpile_expr(node):
    """Transpile a single Python expression node to C++."""
    return _ExprTranspiler().visit(node)


def _transpile_method(class_name, method_name, source):
    """
    Attempt a deterministic Python→C++ conversion of a solver method.

    Returns the C++ function string, or raises ``_CppTranspileError`` if
    the method contains constructs the transpiler cannot handle.
    """
    tree = ast.parse(textwrap.dedent(source))
    func_node = None
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef):
            func_node = node
            break
    if func_node is None:
        raise _CppTranspileError("no function found")

    full_name = f"{class_name}_{method_name}"

    # Build parameter list (skip **kwargs, self, cls)
    params = []
    for arg in func_node.args.args:
        name = arg.arg
        if name in ("self", "cls"):
            continue
        params.append(f"double {name}")
    # Skip **kwargs (func_node.args.kwarg)

    sig = f"std::vector<double> {full_name}({', '.join(params)})"

    body_lines = []
    body_lines.append("    std::vector<double> result;")

    tx = _ExprTranspiler()

    for stmt in func_node.body:
        # Skip: docstrings, comments (ast doesn't keep comments, but docstrings are Expr(Constant(str)))
        if (
            isinstance(stmt, ast.Expr)
            and isinstance(stmt.value, ast.Constant)
            and isinstance(stmt.value.value, str)
        ):
            body_lines.append(f"    // {stmt.value.value.strip()}")
            continue

        # result = []
        if (
            isinstance(stmt, ast.Assign)
            and len(stmt.targets) == 1
            and isinstance(stmt.targets[0], ast.Name)
            and stmt.targets[0].id == "result"
            and isinstance(stmt.value, ast.List)
            and len(stmt.value.elts) == 0
        ):
            # Already emitted above
            continue

        # result.append(expr)
        if (
            isinstance(stmt, ast.Expr)
            and isinstance(stmt.value, ast.Call)
            and isinstance(stmt.value.func, ast.Attribute)
            and isinstance(stmt.value.func.value, ast.Name)
            and stmt.value.func.value.id == "result"
            and stmt.value.func.attr == "append"
        ):
            arg_expr = tx.visit(stmt.value.args[0])
            body_lines.append(f"    result.push_back({arg_expr});")
            continue

        # return result
        if isinstance(stmt, ast.Return):
            if stmt.value is None:
                body_lines.append("    return result;")
            elif isinstance(stmt.value, ast.Name) and stmt.value.id == "result":
                body_lines.append("    return result;")
            else:
                val = tx.visit(stmt.value)
                body_lines.append(f"    return {val};")
            continue

        # Simple assignment: var = expr
        if (
            isinstance(stmt, ast.Assign)
            and len(stmt.targets) == 1
            and isinstance(stmt.targets[0], ast.Name)
        ):
            var_name = stmt.targets[0].id
            val = tx.visit(stmt.value)
            body_lines.append(f"    double {var_name} = {val};")
            continue

        raise _CppTranspileError(f"unsupported statement: {type(stmt).__name__}")

    # Only append a trailing return if the body doesn't already end with one
    if not body_lines or not body_lines[-1].strip().startswith("return"):
        body_lines.append("    return result;")

    return sig + " {\n" + "\n".join(body_lines) + "\n}"


# ── LLM fallback ────────────────────────────────────────────────────────────


def convert_with_phi3(class_name, method_name, python_code):
    full_func_name = f"{class_name}_{method_name}"
    prompt = f"""Convert this Python static method from class '{class_name}' to a C++ function.
Function name: {full_func_name}
Return type: std::vector<double> (if all real) or std::vector<std::complex<double>> (if complex I is used).
Math: Use std::log, std::sqrt, std::exp, std::pow from <cmath>.
Complex: Use std::complex<double> and std::complex<double>(0, 1) for I.
Floating point: Use double for all calculations. Use .0 suffix for integers (e.g. 530.0 instead of 530).
Return a vector containing the results.
DO NOT add any logic that is not in the Python code. DO NOT add conversion factors.
ONLY output the C++ code, no preamble, no markdown blocks, no explanation.

Python:
{python_code}
"""
    try:
        result = subprocess.run(
            ["ollama", "run", "phi3", prompt],
            capture_output=True,
            text=True,
            check=True,
        )
        code = result.stdout.strip()
        # Clean up markdown
        code = re.sub(r"```cpp\n?", "", code)
        code = re.sub(r"```\n?", "", code)
        code = code.strip("`").strip()
        # Find the start of the C++ code
        if "std::vector" in code:
            code = code[code.find("std::vector") :]
        return code
    except Exception as e:
        return f"// Error converting {full_func_name}: {e}"


# ── Conversion dispatch ─────────────────────────────────────────────────────


def convert_method(class_name, method_name, source):
    """Try deterministic transpilation first; fall back to LLM."""
    try:
        return _transpile_method(class_name, method_name, source)
    except _CppTranspileError:
        return convert_with_phi3(class_name, method_name, source)


# ── Public helpers ───────────────────────────────────────────────────────────


def get_methods_from_file(file_path):
    with open(file_path, "r") as f:
        tree = ast.parse(f.read())

    classes_methods = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef):
            class_name = node.name
            methods = []
            for item in node.body:
                if isinstance(item, ast.FunctionDef):
                    if "__" in item.name and not item.name.startswith("__"):
                        methods.append(item)
            if methods:
                classes_methods[class_name] = methods
    return classes_methods


def get_method_source(file_path, method_node):
    with open(file_path, "r") as f:
        lines = f.readlines()
    return "".join(lines[method_node.lineno - 1 : method_node.end_lineno])


# ── Main entry point ────────────────────────────────────────────────────────


def main(project_dir=".", input_file=None):
    project_dir = os.path.abspath(project_dir)

    # Resolve the Python source to convert.  Prefer an explicit *input_file*,
    # then fall back to vakyume_lib.py (reconstruct output), and finally to the
    # legacy vakyume_certified.py.
    if input_file:
        certified_file = (
            input_file
            if os.path.isabs(input_file)
            else os.path.join(project_dir, input_file)
        )
    else:
        lib_file = os.path.join(project_dir, "vakyume_lib.py")
        cert_file = os.path.join(project_dir, "vakyume_certified.py")
        certified_file = lib_file if os.path.exists(lib_file) else cert_file

    cpp_output_path = os.path.join(project_dir, "vakyume.cpp")
    test_binary_path = os.path.join(project_dir, "vakyume_test")

    if not os.path.exists(certified_file):
        print(f"Error: {certified_file} not found.")
        return

    print(f"C++ source: {certified_file}")

    classes_methods = get_methods_from_file(certified_file)

    all_tasks = []
    for class_name, methods in classes_methods.items():
        for method in methods:
            source = get_method_source(certified_file, method)
            all_tasks.append(
                (class_name, method.name, source, [arg.arg for arg in method.args.args])
            )

    # Set a limit if you want to test on a subset first, e.g., limit = 10
    limit = None
    tasks_to_run = all_tasks[:limit] if limit else all_tasks

    print(f"Converting {len(tasks_to_run)} functions to C++...")

    # First pass: deterministic transpiler (fast, no LLM).
    # Collect indices that need LLM fallback.
    cpp_functions = [None] * len(tasks_to_run)
    llm_indices = []

    for i, (class_name, method_name, source, _args) in enumerate(tasks_to_run):
        try:
            cpp_functions[i] = _transpile_method(class_name, method_name, source)
        except _CppTranspileError:
            llm_indices.append(i)

    if llm_indices:
        print(
            f"  {len(tasks_to_run) - len(llm_indices)} converted deterministically, "
            f"{len(llm_indices)} require LLM..."
        )
        with ThreadPoolExecutor(max_workers=2) as executor:
            futures = {
                idx: executor.submit(
                    convert_with_phi3,
                    tasks_to_run[idx][0],
                    tasks_to_run[idx][1],
                    tasks_to_run[idx][2],
                )
                for idx in llm_indices
            }
            for idx, fut in futures.items():
                cpp_functions[idx] = fut.result()
    else:
        print(f"  All {len(tasks_to_run)} converted deterministically (no LLM needed).")

    cpp_output = [
        "#include <iostream>",
        "#include <vector>",
        "#include <cmath>",
        "#include <complex>",
        "#include <stdexcept>",
        "",
        "using namespace std::complex_literals;",
        "",
    ]
    cpp_output.extend(f for f in cpp_functions if f)

    test_suite = [
        "",
        "int main() {",
        '    std::cout << "Running test suite..." << std::endl;',
        "    int pass = 0, fail = 0;",
    ]

    for class_name, method_name, source, args in tasks_to_run:
        # Filter out self/cls/kwargs from arg count
        real_args = [a for a in args if a not in ("self", "cls")]
        dummy_args = ", ".join(["1.0"] * len(real_args))
        full_name = f"{class_name}_{method_name}"
        test_suite.append(f'    std::cout << "Testing {full_name}... ";')
        test_suite.append(
            f"    try {{ auto r = {full_name}({dummy_args}); "
            f'std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; }}'
            f' catch(std::exception& e) {{ std::cout << "FAIL: " << e.what() << std::endl; fail++; }}'
            f' catch(...) {{ std::cout << "FAIL" << std::endl; fail++; }}'
        )

    test_suite.append("")
    test_suite.append(
        '    std::cout << "\\n" << pass << " passed, " << fail << " failed." << std::endl;'
    )
    test_suite.append("    return fail > 0 ? 1 : 0;")
    test_suite.append("}")

    cpp_output.extend(test_suite)

    with open(cpp_output_path, "w") as f:
        f.write("\n".join(cpp_output))

    print(f"C++ source created: {cpp_output_path}")
    print("Compiling...")
    compile_res = subprocess.run(
        ["g++", "-std=c++17", cpp_output_path, "-o", test_binary_path, "-lm"],
        capture_output=True,
        text=True,
    )
    if compile_res.returncode != 0:
        print("Compilation FAILED:")
        print(compile_res.stderr)
    else:
        print(f"Compilation SUCCESS. Binary: {test_binary_path}")
        print("Running tests...")
        run_res = subprocess.run([test_binary_path], capture_output=True, text=True)
        print(run_res.stdout)
        if run_res.returncode != 0:
            print(f"Tests exited with code {run_res.returncode}")


if __name__ == "__main__":
    main()
