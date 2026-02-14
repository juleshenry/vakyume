import os
import re
import importlib.util
import numpy as np
from scipy.optimize import fsolve
import json
import sys

# Add shards directory to sys.path
sys.path.append(os.path.abspath("shards"))
sys.path.append(os.path.abspath("."))


def get_equation(shard_name):
    # Try both with and without .py in the prompt name
    paths = [
        f"repair_prompts/{shard_name}.pyeqn",
        f"repair_prompts/{shard_name}.py.pyeqn",
    ]
    for prompt_path in paths:
        if os.path.exists(prompt_path):
            with open(prompt_path, "r") as f:
                return f.read().strip()
    return None


def extract_variables(equation):
    # Basic variable extractor, avoiding keywords
    vars = re.findall(r"[a-zA-Z_][a-zA-Z0-9_]*", equation)
    reserved = {
        "log",
        "sqrt",
        "exp",
        "pow",
        "e",
        "I",
        "Piecewise",
        "LambertW",
        "Eq",
        "symbols",
        "solve",
    }
    return sorted(list(set(vars) - reserved))


def solve_numerically(equation, target_var, fixed_values):
    # Split equation into LHS and RHS
    if "=" not in equation:
        return None
    lhs_str, rhs_str = equation.split("=")

    def objective(x):
        local_scope = fixed_values.copy()
        local_scope[target_var] = (
            float(x[0]) if isinstance(x, (list, np.ndarray)) else float(x)
        )
        # Use a safe eval with numpy/math functions
        safe_dict = {
            "log": np.log,
            "sqrt": np.sqrt,
            "exp": np.exp,
            "pow": np.power,
            "e": np.e,
            "abs": np.abs,
            "sin": np.sin,
            "cos": np.cos,
            "tan": np.tan,
        }
        safe_dict.update(local_scope)
        try:
            lhs = eval(lhs_str, {"__builtins__": None}, safe_dict)
            rhs = eval(rhs_str, {"__builtins__": None}, safe_dict)
            diff = lhs - rhs
            # Ensure we return a real float for fsolve
            if isinstance(diff, complex):
                return diff.real
            return float(diff)
        except Exception:
            return 1e9

    # Try fsolve with a few starting points to handle nonlinearity
    for start in [1.0, 10.0, 0.1, 100.0, -1.0]:
        try:
            solution, info, ier, msg = fsolve(objective, [start], full_output=True)
            if ier == 1:
                # Double check the solution
                if abs(objective(solution[0])) < 1e-3:
                    return float(solution[0])
        except Exception:
            continue
    return None

    lhs_str, rhs_str = equation.split("=")

    def objective(x):
        local_scope = fixed_values.copy()
        local_scope[target_var] = x
        # Use a safe eval with numpy/math functions
        safe_dict = {
            "log": np.log,
            "sqrt": np.sqrt,
            "exp": np.exp,
            "pow": np.power,
            "e": np.e,
            "abs": np.abs,
            "sin": np.sin,
            "cos": np.cos,
            "tan": np.tan,
        }
        safe_dict.update(local_scope)
        try:
            lhs = eval(lhs_str, {"__builtins__": None}, safe_dict)
            rhs = eval(rhs_str, {"__builtins__": None}, safe_dict)
            return lhs - rhs
        except Exception:
            return 1e9

    # Try fsolve with a few starting points to handle nonlinearity
    for start in [1.0, 10.0, 0.1, 100.0]:
        solution, info, ier, msg = fsolve(objective, start, full_output=True)
        if ier == 1:
            return float(solution[0])
    return None


def verify_shard(shard_path):
    shard_name = os.path.basename(shard_path)
    if shard_name.endswith(".py"):
        shard_name = shard_name[:-3]

    equation = get_equation(shard_name)
    if not equation:
        return {"error": "Equation not found"}

    variables = extract_variables(equation)
    # print(f"DEBUG: variables={variables}")

    # Load the shard
    spec = importlib.util.spec_from_file_location(shard_name, shard_path)
    module = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(module)
    except Exception as e:
        return {"error": f"Failed to load module: {e}"}

    # Determine the class name (usually first class in module)
    class_name = next(
        (
            name
            for name, obj in module.__dict__.items()
            if isinstance(obj, type) and obj.__module__ == shard_name
        ),
        None,
    )
    # print(f"DEBUG: class_name={class_name}")
    if not class_name:
        # Fallback to any class if the module check fails
        class_name = next(
            (
                name
                for name, obj in module.__dict__.items()
                if isinstance(obj, type) and not name.startswith("__")
            ),
            None,
        )

    if not class_name:
        return {"error": "No class found in shard"}

    cls = getattr(module, class_name)
    results = {}

    for target_var in variables:
        method_name = f"eqn_{shard_name.split('_eqn_')[-1]}__{target_var}"
        if not hasattr(cls, method_name):
            continue

        method = getattr(cls, method_name)

        # Test Case
        np.random.seed(42)  # Deterministic tests
        fixed_vars = [v for v in variables if v != target_var]
        test_inputs = {v: np.random.uniform(1.0, 10.0) for v in fixed_vars}

        gold_val = solve_numerically(equation, target_var, test_inputs)
        if gold_val is None:
            results[target_var] = "Gold solve failed"
            continue

        try:
            # Sub-methods return a list
            output_list = method(**test_inputs)
            # Check if gold_val (or something close) is in the list
            passed = any(
                np.isclose(
                    gold_val,
                    float(val.real if hasattr(val, "real") else val),
                    rtol=1e-5,
                )
                for val in output_list
            )

            results[target_var] = {
                "status": "PASS" if passed else "FAIL",
                "gold": gold_val,
                "got": [
                    float(v.real if hasattr(v, "real") else v) for v in output_list
                ],
                "inputs": test_inputs,
            }
        except Exception as e:
            results[target_var] = {"status": "ERROR", "msg": str(e)}

    return results


if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1:
        path = sys.argv[1]
        if os.path.isfile(path):
            print(json.dumps(verify_shard(path), indent=2))
        elif os.path.isdir(path):
            all_results = {}
            for f in os.listdir(path):
                if f.endswith(".py") and "_eqn_" in f:
                    print(f"Verifying {f}...")
                    all_results[f] = verify_shard(os.path.join(path, f))
            with open("reports/gold_audit.json", "w") as out:
                json.dump(all_results, out, indent=2)
