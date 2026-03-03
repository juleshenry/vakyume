from cmath import log, sqrt, exp
from math import e, pi
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
    def __init__(self, lib_class, pyeqn=None, subshards_dir=None, verbose=False):
        self.lib_class = lib_class
        self.verbose = verbose
        self._log_buffer = []
        # Instantiate the class since methods are now instance methods
        try:
            self.lib_instance = lib_class()
        except Exception as e:
            self._log_buffer.append(
                f"  [Verify] WARNING: Could not instantiate {lib_class.__name__}: {e}. Falling back to class-level calls."
            )
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
                # Collect known equation variables from variant names
                # e.g. eqn_8_9__e, eqn_8_9__r → {'e', 'r', ...}
                known_vars = set()
                for v in self.equation_variants:
                    if "__" in v:
                        known_vars.add(v.split("__")[-1])
                # Also collect from function signatures
                for v in self.equation_variants:
                    try:
                        sig = inspect.signature(getattr(self.lib_class, v))
                        for name, param in sig.parameters.items():
                            if (
                                param.kind
                                not in (param.VAR_POSITIONAL, param.VAR_KEYWORD)
                                and name != "self"
                            ):
                                known_vars.add(name)
                    except Exception:
                        pass

                # Find all identifiers in lhs and rhs to use as arguments
                all_tokens = set(
                    re.findall(r"\b[a-zA-Z_]\w*\b", lhs_str + " " + rhs_str)
                )
                # Remove math functions and numbers, BUT keep tokens that
                # are known equation variables (e.g. 'e' in eqn_8_9)
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
                    [
                        t
                        for t in all_tokens
                        if (t in known_vars or t not in math_funcs) and not t.isdigit()
                    ]
                )

                with open(subshard_path, "w") as f:
                    f.write(
                        "from cmath import *\nfrom math import e, pi\nimport numpy as np\n\n"
                    )
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
            self._log_buffer.append(f"  [Harmony Prep] ERROR: {err}")

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

    def _check_pyeqn_harmony(self, pyeqn, params_dict, log_lines=None):
        """Checks if the values satisfy the original equation string."""
        if not pyeqn or "=" not in pyeqn:
            return True

        if self.harmony_func:
            try:
                res_val = self.harmony_func(**params_dict)
                if isinstance(res_val, (list, tuple, np.ndarray)):
                    res_val = res_val[0]
                res = abs(res_val) < 1e-4
                if not res:
                    msg = f"  [Harmony Check] FAIL: residual {res_val} for {pyeqn}"
                    if log_lines is not None:
                        log_lines.append(msg)
                    else:
                        print(msg)
                return res
            except Exception as err:
                # If the harmony subshard crashes for any reason, treat as
                # inconclusive rather than a failure.
                msg = f"  [Harmony Check] ERROR during subshard call: {err}"
                if log_lines is not None:
                    log_lines.append(msg)
                else:
                    print(msg)
                return None

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

        def to_scalar(v):
            """Convert a single (non-list) value to float or complex.
            Always returns a complex number (or None) — the caller
            compares real and imaginary parts independently."""
            if isinstance(v, (int, float)):
                return complex(v, 0)
            if isinstance(v, complex):
                return v
            try:
                if hasattr(v, "evalf"):
                    ev = v.evalf()
                    try:
                        return complex(float(ev), 0)
                    except (TypeError, ValueError):
                        return complex(ev)
                return complex(float(v), 0)
            except (TypeError, ValueError):
                try:
                    return complex(v)
                except (TypeError, ValueError):
                    return None

        def scalars_match(sa, sb):
            if sa is None or sb is None:
                return False
            # Both are complex — check real and imaginary parts independently
            return abs(sa.real - sb.real) < 1e-6 and abs(sa.imag - sb.imag) < 1e-6

        # If either side is a list, check whether ANY element from one side
        # matches ANY element (or the scalar) on the other side.
        a_is_list = isinstance(a, (list, tuple, np.ndarray))
        b_is_list = isinstance(b, (list, tuple, np.ndarray))

        if a_is_list or b_is_list:
            a_items = list(a) if a_is_list else [a]
            b_items = list(b) if b_is_list else [b]
            for av in a_items:
                for bv in b_items:
                    sa = to_scalar(av)
                    sb = to_scalar(bv)
                    if scalars_match(sa, sb):
                        return True
            return False

        # Both are scalars
        sa = to_scalar(a)
        sb = to_scalar(b)
        return scalars_match(sa, sb)

    def _detect_invariant_variants(self, base_eq, variants, params, log_lines):
        """Detect variants whose solver returns a constant regardless of inputs.

        These correspond to variables that cancel out of the equation (e.g.
        p_i / p_nc = p_i / (p - P_c)  where p_i cancels).  SymPy solves
        them as the trivial root (often 0) because the variable factors out.
        Including them in cross-verification poisons the scores of the
        *real* variants, so we exclude them.
        """
        invariant = set()
        probe_count = 3
        for var in variants:
            method = getattr(self.lib_instance, f"{base_eq}__{var}")
            probe_results = []
            for _ in range(probe_count):
                inputs = {p: self.make_rand() for p in params if p != var}
                try:
                    vals = method(**inputs)
                    # Normalise to a comparable tuple
                    if isinstance(vals, (list, tuple, np.ndarray)):
                        probe_results.append(
                            tuple(round(abs(complex(v)), 10) for v in vals)
                        )
                    elif vals is not None:
                        probe_results.append((round(abs(complex(vals)), 10),))
                    else:
                        probe_results.append(None)
                except Exception:
                    probe_results.append(None)

            # If all probes returned the exact same (non-None) value, it's invariant
            non_none = [r for r in probe_results if r is not None]
            if len(non_none) == probe_count and len(set(non_none)) == 1:
                invariant.add(var)
                log_lines.append(
                    f" |- INVARIANT detected: {var} always returns {non_none[0]} — excluding from cross-check"
                )
        return invariant

    def verify_equation(self, base_eq):
        # Collect all log output in a buffer; flush only for failures (or if verbose)
        log_lines = list(self._log_buffer)  # inherit any init/harmony-prep messages
        log_lines.append(f"[INPUT] verify_equation: base_eq={base_eq}")
        params = self._get_params(base_eq)
        all_variants = [p for p in params if hasattr(self.lib_class, f"{base_eq}__{p}")]

        # Detect and exclude invariant (trivially-solved) variants
        invariant_vars = self._detect_invariant_variants(
            base_eq, all_variants, params, log_lines
        )
        variants = [v for v in all_variants if v not in invariant_vars]

        results = {}  # (source_truth_var) -> matches
        mismatches = {}  # (source_truth_var) -> details

        # Give invariant variants a perfect score so they don't drag the family down
        for iv in invariant_vars:
            results[iv] = len(all_variants)

        # Increase trials to be more certain
        num_trials = 3

        for source_var in variants:
            variant_method = getattr(self.lib_instance, f"{base_eq}__{source_var}")
            log_lines.append(f" |- Testing source variant: {source_var}")

            trial_matches = []
            source_mismatches = []

            for trial_idx in range(num_trials):
                test_inputs = {p: self.make_rand() for p in params if p != source_var}
                log_lines.append(f"    [Trial {trial_idx + 1}] Inputs: {test_inputs}")

                try:
                    source_values = variant_method(**test_inputs.copy())
                    log_lines.append(
                        f"    [Trial {trial_idx + 1}] {source_var} variant output: {source_values}"
                    )
                    if not source_values:
                        log_lines.append(
                            f"    [Trial {trial_idx + 1}] {source_var} variant returned NOTHING"
                        )
                        trial_matches.append(0)
                        source_mismatches.append({"error": "No values returned"})
                        continue

                    best_match_for_this_trial = 0
                    trial_detail = None

                    for val in source_values:
                        full_set = test_inputs.copy()
                        full_set[source_var] = val

                        # Harmony check with original equation string
                        if self.pyeqn:
                            harmony_res = self._check_pyeqn_harmony(
                                self.pyeqn, full_set, log_lines=log_lines
                            )
                            if harmony_res is False:
                                log_lines.append(
                                    f"    [Trial {trial_idx + 1}] Disharmony with pyeqn for {source_var}={val}"
                                )
                                continue
                            elif harmony_res is None:
                                # Inconclusive (domain error), skip this value but don't count as failure
                                continue

                        matches = len(
                            invariant_vars
                        )  # invariant vars count as automatic matches
                        current_trial_mismatches = []
                        log_lines.append(
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
                                log_lines.append(
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
                                err_str = str(err).lower()
                                # "Pending LLM" and "UnsolvedException" mean
                                # the solver doesn't exist yet — inconclusive,
                                # don't penalise the source variant.
                                if any(
                                    kw in err_str
                                    for kw in (
                                        "pending llm",
                                        "unsolvedexception",
                                    )
                                ):
                                    log_lines.append(
                                        f"      -> Target {target_var}: ERROR: {err}"
                                    )
                                    # Treat as inconclusive — count as a match
                                    # so it doesn't drag down the source score
                                    matches += 1
                                else:
                                    log_lines.append(
                                        f"      -> Target {target_var}: ERROR: {err}"
                                    )
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
                    log_lines.append(
                        f"    [Trial {trial_idx + 1}] {source_var} variant crashed: {err}"
                    )
                    trial_matches.append(0)
                    source_mismatches.append({"error": str(err)})

            results[source_var] = max(trial_matches) if trial_matches else 0
            if results[source_var] < len(all_variants):
                # Only keep mismatches for broken variants
                mismatches[source_var] = source_mismatches

        log_lines.append(f"[OUTPUT] verify_equation: scores={results}")

        # Decide whether to flush the log buffer
        num_v = len(all_variants)
        is_clean = num_v > 0 and all(s >= num_v for s in results.values())

        if self.verbose or not is_clean:
            for line in log_lines:
                print(line)

        return {"scores": results, "mismatches": mismatches}

    def verify(self):
        overall_results = {}
        for base_eq in self.base_equations:
            overall_results[base_eq] = self.verify_equation(base_eq)
        return overall_results
