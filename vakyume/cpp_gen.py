import ast
import subprocess
import os
import re
from concurrent.futures import ThreadPoolExecutor


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


def main(project_dir="."):
    project_dir = os.path.abspath(project_dir)
    certified_file = os.path.join(project_dir, "vakyume_certified.py")
    cpp_output_path = os.path.join(project_dir, "vakyume.cpp")
    test_binary_path = os.path.join(project_dir, "vakyume_test")

    if not os.path.exists(certified_file):
        print(f"Error: {certified_file} not found.")
        return

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

    # Using ThreadPoolExecutor to speed up Ollama calls
    # Note: Too many workers might overwhelm the local Ollama instance
    with ThreadPoolExecutor(max_workers=2) as executor:
        futures = [
            executor.submit(convert_with_phi3, t[0], t[1], t[2]) for t in tasks_to_run
        ]
        cpp_functions = [f.result() for f in futures]

    cpp_output = [
        "#include <iostream>",
        "#include <vector>",
        "#include <cmath>",
        "#include <complex>",
        "",
        "using namespace std::complex_literals;",
        "",
    ]
    cpp_output.extend(cpp_functions)

    test_suite = [
        "",
        "int main() {",
        '    std::cout << "Running test suite..." << std::endl;',
    ]

    for class_name, method_name, source, args in tasks_to_run:
        dummy_args = ", ".join(["1.0"] * len(args))
        full_name = f"{class_name}_{method_name}"
        test_suite.append(f'    std::cout << "Testing {full_name}... ";')
        test_suite.append(
            f'    try {{ {full_name}({dummy_args}); std::cout << "OK" << std::endl; }} catch(...) {{ std::cout << "FAIL" << std::endl; }}'
        )

    test_suite.append('    std::cout << "Test suite finished." << std::endl;')
    test_suite.append("    return 0;")
    test_suite.append("}")

    cpp_output.extend(test_suite)

    with open(cpp_output_path, "w") as f:
        f.write("\n".join(cpp_output))

    print(f"C++ source created: {cpp_output_path}")
    print("Compilation...")
    # Using g++ as it's the standard for C++ (gcc can also work but g++ is safer)
    compile_res = subprocess.run(
        ["g++", cpp_output_path, "-o", test_binary_path], capture_output=True, text=True
    )
    if compile_res.returncode != 0:
        print("Compilation FAILED:")
        print(compile_res.stderr)
    else:
        print(f"Compilation SUCCESS. Binary created at {test_binary_path}")
        print("Running tests...")
        run_res = subprocess.run([test_binary_path], capture_output=True, text=True)
        print(run_res.stdout)


if __name__ == "__main__":
    main()
