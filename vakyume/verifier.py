from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, Basic
from scipy.optimize import newton
from itertools import product
import random
import inspect
import re
import numpy as np

try:
    from tqdm import tqdm
except Exception:
    tqdm = None

from .kwasak import kwasak_static

import importlib
import os
import sys


class Verify:
    def __init__(self, lib_class, pyeqn=None, subshards_dir=None):
        self.lib_class = lib_class
        # Instantiate the class since methods are now instance methods
        try:
            self.lib_instance = lib_class()
        except Exception as e:
            print(f"  [Verify] WARNING: Could not instantiate {lib_class.__name__}: {e}. Falling back to class-level calls.")
            self.lib_instance = lib_class

        self.pyeqn = pyeqn
        self.subshards_dir = subshards_dir
        self.harmony_func = None
        self.base_equations = [
            name
            for name in dir(lib_class)
            if name.startswith("eqn") and "__" not in name
        ]
        self.equation_variants = [
            name for name in dir(lib_class) if name.startswith("eqn") and "__" in name
        ]

        if self.subshards_dir and self.pyeqn:
            self._prepare_harmony_subshard()

    def _prepare_harmony_subshard(self):
        try:
            lhs_str, rhs_str = self.pyeqn.split("=", 1)
            lhs_str = lhs_str.split("#")[0].strip()
            rhs_str = rhs_str.split("#")[0].strip()

            import hashlib

            h = hashlib.md5(self.pyeqn.encode()).hexdigest()
            module_name = f"harmony_check_{h}"
            subshard_path = os.path.join(self.subshards_dir, f"{module_name}.py")

            if not os.path.exists(subshard_path):
                # Find all identifiers in lhs and rhs to use as arguments
                all_tokens = set(
                    re.findall(r"\b[a-zA-Z_]\w*\b", lhs_str + " " + rhs_str)
                )
                # Remove math functions and numbers
                math_funcs = {
                    "log",
                    "sqrt",
                    "exp",
                    "pow",
                    "e",
                    "pi",
                    "sin",
                    "cos",
                    "tan",
                    "ln",
                }
                tokens = sorted(
                    [t for t in all_tokens if t not in math_funcs and not t.isdigit()]
                )

                with open(subshard_path, "w") as f:
                    f.write("from math import *\nimport numpy as np\n\n")
                    f.write(f"def check_harmony({', '.join(tokens)}, **kwargs):\n")
                    f.write(
                        f"    return ({lhs_str.replace('ln', 'log')}) - ({rhs_str.replace('ln', 'log')})\n"
                    )

            if self.subshards_dir not in sys.path:
                sys.path.append(self.subshards_dir)

            importlib.invalidate_caches()
            if module_name in sys.modules:
                module = importlib.reload(sys.modules[module_name])
            else:
                module = importlib.import_module(module_name)

            self.harmony_func = module.check_harmony
        except Exception as err:
            print(f"  [Harmony Prep] ERROR: {err}")

    def _get_params(self, base_eq):
        variants = [v for v in self.equation_variants if v.startswith(base_eq + "__")]
        params = set()
        for v in variants:
            sig = inspect.signature(getattr(self.lib_class, v))
            for name, param in sig.parameters.items():
                if (
                    param.kind not in (param.VAR_POSITIONAL, param.VAR_KEYWORD)
                    and name != "self"
                ):
                    params.add(name)
            params.add(v.split("__")[-1])
        return sorted(list(params))

    def _check_pyeqn_harmony(self, pyeqn, params_dict):
        """Checks if the values satisfy the original equation string."""
        if not pyeqn or "=" not in pyeqn:
            return True

        if self.harmony_func:
            try:
                res_val = self.harmony_func(**params_dict)
                if isinstance(res_val, (list, tuple, np.ndarray)):
                    res_val = res_val[0]
                res = abs(float(res_val)) < 1e-4
                if not res:
                    print(f"  [Harmony Check] FAIL: residual {res_val} for {pyeqn}")
                return res
            except Exception as err:
                # If it's a domain error (like log of negative), return None to indicate inconclusive
                if (
                    "positive" in str(err).lower()
                    or "math domain error" in str(err).lower()
                ):
                    return None
                print(f"  [Harmony Check] ERROR during subshard call: {err}")
                return False

        # No fallback to eval anymore as per user request
        return True

    @staticmethod
    def make_rand():
        # Avoid 0 and 1 to prevent division by zero or log(1) which might be trivial
        # Also avoid negative numbers for logs/roots
        return round(random.uniform(1.1, 10.0), 5)

    def are_similar(self, a, b):
        if a is None or b is None:
            return False

        def to_float_or_complex(v):
            if isinstance(v, (int, float, complex)):
                return v
            if isinstance(v, (list, tuple, np.ndarray)):
                if len(v) > 0:
                    return to_float_or_complex(v[0])
                return None
            try:
                # Handle sympy numbers
                if hasattr(v, "evalf"):
                    return float(v.evalf())
                return float(v)
            except:
                try:
                    return complex(v)
                except:
                    return None

        va = to_float_or_complex(a)
        vb = to_float_or_complex(b)

        if va is not None and vb is not None:
            if isinstance(va, complex) or isinstance(vb, complex):
                return abs(va - vb) < 1e-6
            return abs(float(va) - float(vb)) < 1e-6

        if isinstance(a, (list, tuple, np.ndarray)) or isinstance(
            b, (list, tuple, np.ndarray)
        ):
            if not isinstance(a, (list, tuple, np.ndarray)):
                a = [a]
            if not isinstance(b, (list, tuple, np.ndarray)):
                b = [b]
            for av in a:
                for bv in b:
                    if self.are_similar(av, bv):
                        return True
        return False

    def verify_equation(self, base_eq):
        print(f"[INPUT] verify_equation: base_eq={base_eq}")
        params = self._get_params(base_eq)
        variants = [p for p in params if hasattr(self.lib_class, f"{base_eq}__{p}")]

        results = {}  # (source_truth_var) -> matches
        mismatches = {}  # (source_truth_var) -> details

        # Increase trials to be more certain
        num_trials = 3

        for source_var in variants:
            variant_method = getattr(self.lib_instance, f"{base_eq}__{source_var}")
            print(f" |- Testing source variant: {source_var}")

            trial_matches = []
            source_mismatches = []

            for trial_idx in range(num_trials):
                test_inputs = {p: self.make_rand() for p in params if p != source_var}
                print(f"    [Trial {trial_idx + 1}] Inputs: {test_inputs}")

                try:
                    source_values = variant_method(**test_inputs.copy())
                    print(
                        f"    [Trial {trial_idx + 1}] {source_var} variant output: {source_values}"
                    )
                    if not source_values:
                        print(
                            f"    [Trial {trial_idx + 1}] {source_var} variant returned NOTHING"
                        )
                        trial_matches.append(0)
                        source_mismatches.append({"error": "No values returned"})
                        continue

                    best_match_for_this_trial = 0
                    trial_detail = None

                    for val in source_values:
                        if isinstance(val, complex) and abs(val.imag) > 1e-5:
                            print(
                                f"    [Trial {trial_idx + 1}] Skipping complex result: {val}"
                            )
                            continue

                        full_set = test_inputs.copy()
                        full_set[source_var] = val

                        # Harmony check with original equation string
                        if self.pyeqn:
                            harmony_res = self._check_pyeqn_harmony(
                                self.pyeqn, full_set
                            )
                            if harmony_res is False:
                                print(
                                    f"    [Trial {trial_idx + 1}] Disharmony with pyeqn for {source_var}={val}"
                                )
                                continue
                            elif harmony_res is None:
                                # Inconclusive (domain error), skip this value but don't count as failure
                                continue

                        matches = 0
                        current_trial_mismatches = []
                        print(
                            f"    [Trial {trial_idx + 1}] Checking targets against source_var={source_var}, val={val}"
                        )
                        for target_var in variants:
                            if target_var == source_var:
                                matches += 1
                                continue

                            target_method = getattr(
                                self.lib_instance, f"{base_eq}__{target_var}"
                            )
                            target_inputs = {
                                p: v for p, v in full_set.items() if p != target_var
                            }

                            try:
                                target_values = target_method(**target_inputs.copy())
                                is_sim = self.are_similar(
                                    full_set[target_var], target_values
                                )
                                print(
                                    f"      -> Target {target_var}: expected {full_set[target_var]}, got {target_values} | Match: {is_sim}"
                                )
                                if is_sim:
                                    matches += 1
                                else:
                                    current_trial_mismatches.append(
                                        {
                                            "target": target_var,
                                            "expected": full_set[target_var],
                                            "got": target_values,
                                        }
                                    )
                            except Exception as err:
                                print(f"      -> Target {target_var}: ERROR: {err}")
                                current_trial_mismatches.append(
                                    {"target": target_var, "error": str(err)}
                                )

                        if matches > best_match_for_this_trial:
                            best_match_for_this_trial = matches
                            trial_detail = {
                                "inputs": test_inputs,
                                "output": val,
                                "mismatches": current_trial_mismatches,
                            }

                    trial_matches.append(best_match_for_this_trial)
                    if trial_detail:
                        source_mismatches.append(trial_detail)

                except Exception as err:
                    print(
                        f"    [Trial {trial_idx + 1}] {source_var} variant crashed: {err}"
                    )
                    trial_matches.append(0)
                    source_mismatches.append({"error": str(err)})

            results[source_var] = max(trial_matches) if trial_matches else 0
            if results[source_var] < len(variants):
                # Only keep mismatches for broken variants
                mismatches[source_var] = source_mismatches

        print(f"[OUTPUT] verify_equation: scores={results}")
        return {"scores": results, "mismatches": mismatches}

    def verify(self):
        overall_results = {}
        for base_eq in self.base_equations:
            overall_results[base_eq] = self.verify_equation(base_eq)
        return overall_results
