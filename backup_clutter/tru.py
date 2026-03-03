from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, Basic
from scipy.optimize import newton
from itertools import product
import random
import inspect
import numpy as np

try:
    from tqdm import tqdm
except Exception:
    tqdm = None

from kwasak import kwasak
from config import *


class Verify:
    def __init__(self, lib_class):
        self.lib_class = lib_class
        self.base_equations = [
            name
            for name in dir(lib_class)
            if name.startswith("eqn") and "__" not in name
        ]
        self.equation_variants = [
            name for name in dir(lib_class) if name.startswith("eqn") and "__" in name
        ]

    def _get_params(self, base_eq):
        variants = [v for v in self.equation_variants if v.startswith(base_eq + "__")]
        params = set()
        for v in variants:
            sig = inspect.signature(getattr(self.lib_class, v))
            params.update(sig.parameters.keys())
            params.add(v.split("__")[-1])
        return sorted(list(params))

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
                return abs(va - vb) < 1e-5
            return abs(float(va) - float(vb)) < 1e-5

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
        params = self._get_params(base_eq)
        variants = [p for p in params if hasattr(self.lib_class, f"{base_eq}__{p}")]

        results = {}  # (source_truth_var) -> matches

        # Increase trials to be more certain
        num_trials = 3

        for source_var in variants:
            variant_method = getattr(self.lib_class, f"{base_eq}__{source_var}")

            trial_matches = []
            for trial_idx in range(num_trials):
                test_inputs = {p: self.make_rand() for p in params if p != source_var}

                try:
                    source_values = variant_method(**test_inputs)
                    if not source_values:
                        print(
                            f"  [OOO] {base_eq}__{source_var} returned no values with inputs {test_inputs}"
                        )
                        trial_matches.append(0)
                        continue

                    best_match_for_this_trial = 0
                    for val in source_values:
                        if isinstance(val, complex) and abs(val.imag) > 1e-5:
                            continue

                        full_set = test_inputs.copy()
                        full_set[source_var] = val

                        matches = 0
                        mismatches = []
                        for target_var in variants:
                            if target_var == source_var:
                                matches += 1
                                continue

                            target_method = getattr(
                                self.lib_class, f"{base_eq}__{target_var}"
                            )
                            target_inputs = {
                                p: v for p, v in full_set.items() if p != target_var
                            }

                            try:
                                target_values = target_method(**target_inputs)
                                if self.are_similar(
                                    full_set[target_var], target_values
                                ):
                                    matches += 1
                                else:
                                    mismatches.append((target_var, target_values))
                            except Exception as e:
                                mismatches.append((target_var, f"Error: {e}"))

                        if matches < len(variants):
                            print(
                                f"  [OOO] Trial {trial_idx}: {base_eq}__{source_var} output {val} is inconsistent."
                            )
                            for m_var, m_val in mismatches:
                                print(
                                    f"    |- {base_eq}__{m_var} expected {full_set[m_var]} but got {m_val}"
                                )

                        best_match_for_this_trial = max(
                            best_match_for_this_trial, matches
                        )
                    trial_matches.append(best_match_for_this_trial)
                except Exception as e:
                    print(f"  [OOO] {base_eq}__{source_var} crashed: {e}")
                    trial_matches.append(0)

            results[source_var] = max(trial_matches) if trial_matches else 0

        return results

    def verify(self):
        overall_results = {}
        for base_eq in self.base_equations:
            overall_results[base_eq] = self.verify_equation(base_eq)
        return overall_results
