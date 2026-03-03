import os
import sys
import shutil
import inspect
import importlib
import importlib.util
import json
import re
import timeout_decorator
import argparse
import time
import ast
import py_compile

try:
    from tqdm import tqdm
except ImportError:
    tqdm = None
from sympy import Symbol, solve, Eq, sympify, symbols
from .verifier import Verify
from .config import (
    TAB,
    IntractableSolution,
    OllamaOffline,
    UnsolvedException,
    COOLDOWN_SECONDS,
    MAX_COMP_TIME_SECONDS,
    FUNKTORZ,
)
from .llm import escribir_codigo, repair_codigo, extract_code

# Paths
DEFAULT_MAX_ROUNDS = 10


class VakyumeEncoder(json.JSONEncoder):
    def default(self, obj):
        if hasattr(obj, "evalf"):
            try:
                if hasattr(obj, "is_number") and obj.is_number:
                    if (
                        hasattr(obj, "is_complex")
                        and obj.is_complex
                        and not obj.is_real
                    ):
                        return str(obj)
                    return float(obj.evalf())
                return str(obj)
            except:
                return str(obj)
        if isinstance(obj, (set, tuple)):
            return list(obj)
        if type(obj).__module__ == "numpy":
            if hasattr(obj, "tolist"):
                return obj.tolist()
            try:
                return float(obj)
            except:
                try:
                    return int(obj)
                except:
                    return str(obj)
        return str(obj)


def case_safe_name(name: str) -> str:
    """Returns a case-safe version of a name for filesystems that are case-insensitive.
    Appends _cap to any name that contains uppercase letters.
    """
    if any(c.isupper() for c in name):
        return f"{name}_cap"
    return name


class PipelineContext:
    def __init__(self, project_dir):
        self.project_dir = os.path.abspath(project_dir)
        self.notes_dir = os.path.join(self.project_dir, "notes")
        self.shards_dir = os.path.join(self.project_dir, "shards")
        # Every method itself a shard in a family
        self.reports_dir = os.path.join(self.project_dir, "reports")
        self.repair_prompts_dir = os.path.join(self.project_dir, "repair_prompts")
        self.certified_file = os.path.join(self.project_dir, "vakyume_certified.py")

        # Ensure directories exist
        for d in [
            self.shards_dir,
            self.reports_dir,
            self.repair_prompts_dir,
            self.notes_dir,
        ]:
            if not os.path.exists(d):
                os.makedirs(d)

    def clear_shards(self):
        if os.path.exists(self.shards_dir):
            shutil.rmtree(self.shards_dir)
        os.makedirs(self.shards_dir)
        package_dir = os.path.dirname(__file__)
        kwasak_src = os.path.join(package_dir, "kwasak.py")
        if os.path.exists(kwasak_src):
            shutil.copy(kwasak_src, os.path.join(self.shards_dir, "kwasak.py"))


class Solver:
    def __init__(self):
        self.funktors = FUNKTORZ

    @timeout_decorator.timeout(
        MAX_COMP_TIME_SECONDS,
        timeout_exception=timeout_decorator.timeout_decorator.TimeoutError,
    )
    def get_solns_vanilla_nf(self, nf: str, symb: Symbol, tokens: list = None):
        try:
            # Use explicit symbols to avoid conflicts (like 'Q' being assumptions)
            local_dict = {t: Symbol(t) for t in tokens} if tokens else {}
            expr = sympify(nf, locals=local_dict)
            solns = solve(expr, symb)
            if not solns:
                raise UnsolvedException("Sympy solve returned empty")
            return solns
        except Exception as e:
            raise UnsolvedException(f"Sympy solve failed: {e}")

    def tokenize(self, eqn):
        malos = {"ln", "log"}
        for m in malos:
            eqn = f"{m}".join([o.strip() for o in eqn.split(m)])

        # Dilate functors
        eqn = eqn.split("#")[0]
        dilated = ""
        for i, f in enumerate(eqn):
            if (
                0 < i < len(eqn) - 1
                and eqn[i + 1] != "*"
                and eqn[i - 1] != "*"
                and f in self.funktors
            ):
                dilated += f" {f} "
            else:
                dilated += f
        return dilated.split()

    def get_tokes(self, eqn):
        tokes = set()
        malos = {"ln", "log"}
        for t in self.tokenize(eqn):
            clean = t.strip().replace("(", "").replace(")", "").split("**")[0].strip()
            if clean.isidentifier() and clean not in malos:
                tokes.add(clean)
        return sorted(list(tokes))

    def make_normal_form(self, eqn):
        parts = eqn.split("=")
        if len(parts) != 2:
            return None
        nf = f"({parts[1].split('#')[0].strip()}) - ({parts[0].strip()})"
        return nf.replace("ln(", "log(").replace("ln (", "log(")

    def sympy_failover(self):
        code = [
            f"{TAB}# Placeholder for numerical solver",
            f'{TAB}raise UnsolvedException("Pending LLM/Manual Repair")',
        ]
        return "\n".join(code)


def shard_from_chapters(ctx: PipelineContext):
    print(f"[INPUT] shard_from_chapters: ctx.notes_dir={ctx.notes_dir}")
    solver = Solver()

    if not os.path.exists(ctx.notes_dir):
        print(f"Notes directory not found: {ctx.notes_dir}")
        return

    chapter_files = sorted(os.listdir(ctx.notes_dir))
    created_count = 0

    to_process = []
    for chapter_file in chapter_files:
        if chapter_file.endswith(".py") and not chapter_file.startswith("__"):
            to_process.append(chapter_file)

    if not to_process:
        print(f"[OUTPUT] shard_from_chapters: No chapters to process")
        return

    import_header = (
        "from cmath import log, sqrt, exp\n"
        "from math import e, pi\n"
        "from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp\n"
        "from scipy.optimize import newton\n"
        "import numpy as np\n"
        "from vakyume.config import UnsolvedException\n\n"
    )

    for chapter_file in to_process:
        chapter_path = os.path.join(ctx.notes_dir, chapter_file)
        base_name = os.path.splitext(chapter_file)[0]
        parts = base_name.split("_")
        if parts[0].isdigit():
            parts = parts[1:]
        class_name = "".join(p.capitalize() for p in parts)

        with open(chapter_path, "r") as f:
            lines = f.readlines()

        eqn_number = ""
        for line in lines:
            if x := re.findall(r"\d{1,2}-\d{1,2}\w*", line):
                eqn_number = x[0].replace("-", "_")

            if (
                " = " in line
                and not line.strip().startswith("#")
                and not line.strip().startswith('"""')
            ):
                tokes = solver.get_tokes(line)
                nf = solver.make_normal_form(line)
                if not nf:
                    continue

                family_name = f"{class_name}_eqn_{eqn_number}"
                family_dir = os.path.join(ctx.shards_dir, family_name)
                if not os.path.exists(family_dir):
                    os.makedirs(family_dir)

                # Individual Shards for each variable (Functions, no class)
                for token in tokes:
                    method_name = f"eqn_{eqn_number}__{token}"
                    safe_token = case_safe_name(token)
                    shard_name = f"eqn_{eqn_number}__{safe_token}.py"
                    shard_path = os.path.join(family_dir, shard_name)
                    if os.path.exists(shard_path):
                        continue

                    created_count += 1
                    shard_content = import_header
                    other_args = [t for t in tokes if t != token]
                    args_str = ", ".join(f"{t}: float" for t in other_args)
                    if args_str:
                        args_str = f", {args_str}"
                    header = f"def {method_name}(self{args_str}, **kwargs):"
                    shard_content += header + "\n"
                    shard_content += f"{TAB}# [.pyeqn] {line.strip()}\n"

                    try:
                        solns = solver.get_solns_vanilla_nf(
                            nf, Symbol(token), tokens=tokes
                        )
                        if solns:
                            shard_content += f"{TAB}result = []\n"
                            for soln in solns:
                                shard_content += f"{TAB}{token} = {soln}\n"
                                shard_content += f"{TAB}result.append({token})\n"
                            shard_content += f"{TAB}return result\n"
                        else:
                            shard_content += solver.sympy_failover() + "\n"
                    except Exception:
                        shard_content += solver.sympy_failover() + "\n"

                    with open(shard_path, "w") as sf:
                        sf.write(shard_content)

                    # py_compile check
                    try:
                        py_compile.compile(shard_path, doraise=True)
                    except py_compile.PyCompileError as e:
                        print(f"  [py_compile] FAILED for {shard_path}: {e}")

                # Main shard for the equation (Class and kwasak entry point)
                main_shard_path = os.path.join(family_dir, f"eqn_{eqn_number}.py")
                if not os.path.exists(main_shard_path):
                    with open(main_shard_path, "w") as sf:
                        sf.write(import_header)
                        sf.write("from vakyume.kwasak import kwasak_static\n")
                        for token in tokes:
                            method_name = f"eqn_{eqn_number}__{token}"
                            safe_token = case_safe_name(token)
                            file_name_no_ext = f"eqn_{eqn_number}__{safe_token}"
                            sf.write(f"from .{file_name_no_ext} import {method_name}\n")
                        sf.write(f"\nclass {class_name}:\n")
                        for token in tokes:
                            method_name = f"eqn_{eqn_number}__{token}"
                            sf.write(f"{TAB}{method_name} = {method_name}\n")
                        sf.write(f"\n{TAB}@kwasak_static\n")
                        tokes_str = ", ".join(f"{t}=None" for t in tokes)
                        sf.write(f"{TAB}def eqn_{eqn_number}(self, {tokes_str}):\n")
                        sf.write(f"{TAB * 2}return\n")

    print(f"Scraped {created_count} new shards.")


def verify_family(ctx: PipelineContext, family_name: str, verbose=False):
    """Verifies all shards in an equation family using cross-verification (golden tuple)."""
    family_dir = os.path.join(ctx.shards_dir, family_name)
    if not os.path.exists(family_dir):
        return {"error": f"Family directory {family_name} not found"}

    # We load the main equation file which imports all variants
    class_name = family_name.split("_")[0]
    eqn_num_match = re.search(r"_eqn_(.*)", family_name)
    if not eqn_num_match:
        return {"error": f"Invalid family name {family_name}"}
    eqn_number = eqn_num_match.group(1)

    main_shard_path = os.path.join(family_dir, f"eqn_{eqn_number}.py")
    if not os.path.exists(main_shard_path):
        return {"error": f"Main shard {main_shard_path} not found"}

    # Ensure shards_dir is in sys.path
    if ctx.shards_dir not in sys.path:
        sys.path.append(ctx.shards_dir)

    # We need the parent directory of family_dir in sys.path to import it as a package
    parent_dir = os.path.dirname(family_dir)
    if parent_dir not in sys.path:
        sys.path.append(parent_dir)

    module_name = f"{family_name}.eqn_{eqn_number}"

    try:
        # Import the main shard
        # Because we use relative imports in the main shard (from .eqn_... import ...),
        # we must import it as part of a package.
        module = importlib.import_module(module_name)
        cls = getattr(module, class_name, None)
        if not cls:
            return {"error": f"Class {class_name} not found in {module_name}"}

        # Find pyeqn from one of the variants
        pyeqn = None
        for name, attr in inspect.getmembers(cls):
            if name.startswith(f"eqn_{eqn_number}__"):
                try:
                    source = inspect.getsource(attr)
                    match = re.search(r"# \[\.pyeqn\] (.*)", source)
                    if match:
                        pyeqn = match.group(1).strip()
                        break
                except:
                    pass

        subshards_dir = os.path.join(ctx.shards_dir, "subshards")
        if not os.path.exists(subshards_dir):
            os.makedirs(subshards_dir)

        v = Verify(cls, pyeqn=pyeqn, subshards_dir=subshards_dir, verbose=verbose)
        raw_results = v.verify()

        verification_results = {}
        mismatches = {}

        for base_eq, data in raw_results.items():
            scores = data["scores"]
            mismatches.update(data["mismatches"])
            verification_results[base_eq] = scores

        return {
            "results": verification_results,
            "mismatches": mismatches,
            "shard_errors": {},
        }
    except Exception as e:
        import traceback

        return {"error": f"Verification failed: {e}\n{traceback.format_exc()}"}


def verify_all_shards(ctx: PipelineContext, verbose=False, skip_families=None):
    if skip_families:
        print(
            f"Verifying families... (skipping {len(skip_families)} previously solved)"
        )
    else:
        print("Verifying families...")
    all_results = {}
    family_names = sorted(
        [
            d
            for d in os.listdir(ctx.shards_dir)
            if os.path.isdir(os.path.join(ctx.shards_dir, d)) and d != "subshards"
        ]
    )

    to_verify = [f for f in family_names if not skip_families or f not in skip_families]

    iterator = tqdm(to_verify, desc="Verifying", unit="family") if tqdm else to_verify
    for family_name in iterator:
        all_results[family_name] = verify_family(ctx, family_name, verbose=verbose)

    # For skipped families, insert a synthetic passing result so analyze_results
    # categorizes them as solved without re-running verification.
    if skip_families:
        for family_name in family_names:
            if family_name in skip_families and family_name not in all_results:
                eqn_num_match = re.search(r"_eqn_(.*)", family_name)
                if eqn_num_match:
                    eqn_number = eqn_num_match.group(1)
                    # Build a synthetic result where every variant has a perfect score.
                    # We discover variants from the shard files on disk.
                    family_dir = os.path.join(ctx.shards_dir, family_name)
                    variant_files = [
                        f
                        for f in os.listdir(family_dir)
                        if f.endswith(".py") and "__" in f
                    ]
                    variant_names = []
                    for vf in variant_files:
                        # e.g. eqn_1_1__x.py -> x (strip _cap suffix if present)
                        raw = vf.replace(".py", "").split("__")[-1]
                        if raw.endswith("_cap"):
                            raw = raw[:-4]
                        variant_names.append(raw)
                    num_v = len(variant_names)
                    scores = {v: num_v for v in variant_names}
                    all_results[family_name] = {
                        "results": {f"eqn_{eqn_number}": scores},
                        "mismatches": {},
                        "shard_errors": {},
                    }

    with open(os.path.join(ctx.reports_dir, "verification_results.json"), "w") as rf:
        json.dump(all_results, rf, indent=4, cls=VakyumeEncoder)
    return all_results


def analyze_results(ctx: PipelineContext, all_results):
    analysis = {"solved": [], "inconsistent": {}, "failed": []}

    for family_name, family_res in all_results.items():
        if (
            not family_res
            or not isinstance(family_res, dict)
            or "results" not in family_res
        ):
            analysis["failed"].append(
                {
                    "file": family_name,
                    "error": family_res.get("error") if family_res else "Unknown error",
                }
            )
            continue

        eqn_results_data = family_res.get("results")
        mismatches_data = family_res.get("mismatches", {})
        all_solved = True

        for base_eq, variants_inner in eqn_results_data.items():
            if not variants_inner:
                continue
            best_score = max(variants_inner.values())
            num_v = len(variants_inner)

            trusted = [
                v
                for v, score in variants_inner.items()
                if score == best_score and score >= num_v
            ]
            broken = [v for v, score in variants_inner.items() if score < num_v]

            # Detect placeholder shards that raise UnsolvedException("Pending
            # LLM/Manual Repair") — these get perfect scores from the verifier
            # (inconclusive treated as agreement) but are NOT actually solved.
            eqn_suffix = (
                family_name.split("_eqn_")[1] if "_eqn_" in family_name else None
            )
            if eqn_suffix:
                family_dir = os.path.join(ctx.shards_dir, family_name)
                for v in list(variants_inner.keys()):
                    if v in broken:
                        continue
                    safe_v = case_safe_name(v)
                    shard_path = os.path.join(
                        family_dir, f"eqn_{eqn_suffix}__{safe_v}.py"
                    )
                    if os.path.exists(shard_path):
                        with open(shard_path, "r") as _sf:
                            content = _sf.read()
                        if 'raise UnsolvedException("Pending' in content:
                            if v not in broken:
                                broken.append(v)
                            if v in trusted:
                                trusted.remove(v)

            if broken:
                all_solved = False
                analysis["inconsistent"][family_name] = {
                    "broken": broken,
                    "trusted": trusted,
                    "scores": variants_inner,
                    "mismatches": {
                        v: mismatches_data.get(v)
                        for v in broken
                        if v in mismatches_data
                    },
                }
                break

        if all_solved and eqn_results_data:
            analysis["solved"].append(family_name)

    with open(os.path.join(ctx.reports_dir, "analysis.json"), "w") as af:
        json.dump(analysis, af, indent=4, cls=VakyumeEncoder)
    return analysis


def _generate_derivation_scaffold(pyeqn: str, target_var: str, header: str) -> str:
    """Generate a code scaffold with CONCRETE algebraic derivation for the LLM.

    The scaffold contains step-by-step comments using actual variable names and
    coefficients from the equation, ending with the complete answer expression.
    Phi-3 only needs to copy the final expression (fill-in-the-blank).

    For transcendental / multi-root cases, produces a complete brentq solver
    body that just needs a closing ``return``.
    """
    if not pyeqn or not target_var or not header:
        return ""

    import re as _re

    eq = " ".join(pyeqn.split())
    if " = " not in eq:
        return ""
    lhs, rhs = eq.split(" = ", 1)

    pattern = _re.compile(r"\b" + _re.escape(target_var) + r"\b")
    total = len(pattern.findall(eq))
    if total == 0:
        return ""

    # ── Helper: extract the outer coefficient multiplying the (**exp) term ──
    def _extract_outer_coeff(expr: str, pow_fragment: str) -> str:
        """Given 'A * B * (fragment)**exp', return 'A * B' (everything except
        the power term itself)."""
        # Remove the power sub-expression to get the rest
        rest = expr.replace(pow_fragment, "1").strip()
        # Clean up: 1 * X → X, X * 1 → X
        rest = _re.sub(r"\b1\s*\*\s*", "", rest)
        rest = _re.sub(r"\s*\*\s*1\b", "", rest)
        rest = rest.strip()
        return rest if rest and rest != "1" else ""

    # ── Detect target in exponent → transcendental → complete brentq body ──
    if _re.search(r"\*\*\s*\([^)]*" + _re.escape(target_var) + r"[^)]*\)", eq):
        # Build a complete brentq solver with concrete residual.
        # Replace target_var with target_var + "_val" in the equation.
        rhs_sub = _re.sub(
            r"\b" + _re.escape(target_var) + r"\b",
            target_var + "_val",
            rhs,
        )
        lines = [header]
        lines.append(f"    # [.pyeqn] {pyeqn}")
        lines.append(
            f"    # {target_var} appears in the exponent — use numerical solver"
        )
        lines.append(f"    from scipy.optimize import brentq")
        lines.append(f"    def _res({target_var}_val):")
        lines.append(f"        return {rhs_sub} - {lhs}")
        lines.append(f"    # Scan for sign change in a wide range")
        lines.append(f"    lo, hi = None, None")
        lines.append(f"    prev = _res(1.01)")
        lines.append(f"    for i in range(1, 10000):")
        lines.append(f"        x = 1.0 + i * 0.01")
        lines.append(f"        try:")
        lines.append(f"            cur = _res(x)")
        lines.append(f"        except Exception:")
        lines.append(f"            continue")
        lines.append(f"        if prev * cur < 0:")
        lines.append(f"            lo, hi = x - 0.01, x")
        lines.append(f"            break")
        lines.append(f"        prev = cur")
        lines.append(f"    if lo is None:")
        lines.append(
            f'        raise UnsolvedException("No sign change found for {target_var}")'
        )
        lines.append(f"    {target_var} = brentq(_res, lo, hi)")
        lines.append(f"    return [{target_var}]")
        return "\n".join(lines)

    # ── Pattern: (target / C) ** n  ──
    # e.g. installed_costs = 38000 * (hp / 10) ** 0.45
    # e.g. installation_cost = 16000 * (NS + 2*NC) * (SCON / 1000) ** 0.35
    # e.g. adiabatic_hp = (w/20) * ((P_2 / P_1) ** 0.286 - 1), solve for P_2
    div_pow = _re.search(
        r"\(\s*" + _re.escape(target_var) + r"\s*/\s*([\w.]+)\s*\)"
        r"\s*\*\*\s*([\d.]+)",
        eq,
    )
    if div_pow:
        divisor = div_pow.group(1)
        exp_val = div_pow.group(2)
        inv_exp = f"1.0 / {exp_val}"
        pow_frag = div_pow.group(0)
        # Detect "- 1" pattern: coeff * ((target/C)**n - 1)
        minus1 = _re.search(_re.escape(pow_frag) + r"\s*-\s*1", eq)
        if minus1:
            # Extract outer coeff around the whole "((target/C)**n - 1)" block
            outer_m1 = _extract_outer_coeff(rhs, "(" + pow_frag + " - 1)")
            if not outer_m1:
                outer_m1 = _extract_outer_coeff(rhs, pow_frag + " - 1")
            if outer_m1:
                ratio_expr = f"{lhs} / ({outer_m1}) + 1"
            else:
                ratio_expr = f"{lhs} + 1"
            answer = f"{divisor} * ({ratio_expr}) ** ({inv_exp})"
            lines = [header]
            lines.append(f"    # [.pyeqn] {pyeqn}")
            lines.append(f"    # Solve for {target_var}:")
            lines.append(
                f"    # Step 1: ({target_var} / {divisor}) ** {exp_val}"
                f" - 1 = {lhs} / ({outer_m1 if outer_m1 else '1'})"
            )
            lines.append(
                f"    # Step 2: ({target_var} / {divisor}) ** {exp_val} = {ratio_expr}"
            )
            lines.append(
                f"    # Step 3: {target_var} / {divisor}"
                f" = ({ratio_expr}) ** ({inv_exp})"
            )
            lines.append(f"    # Step 4: {target_var} = {answer}")
            lines.append(f"    {target_var} = {answer}")
            return "\n".join(lines)
        else:
            # No "- 1" — simple (target/C)**n isolation
            outer = _extract_outer_coeff(rhs, pow_frag)
            if outer:
                iso_expr = f"{lhs} / ({outer})"
            else:
                iso_expr = lhs
            answer = f"{divisor} * ({iso_expr}) ** ({inv_exp})"
            lines = [header]
            lines.append(f"    # [.pyeqn] {pyeqn}")
            lines.append(f"    # Solve for {target_var}:")
            lines.append(
                f"    # Step 1: ({target_var} / {divisor}) ** {exp_val} = {iso_expr}"
            )
            lines.append(
                f"    # Step 2: {target_var} / {divisor} = ({iso_expr}) ** ({inv_exp})"
            )
            lines.append(f"    # Step 3: {target_var} = {answer}")
            lines.append(f"    {target_var} = {answer}")
            return "\n".join(lines)

    # ── Pattern: (A / target) ** n — solve for target in denominator of ratio ──
    # e.g. adiabatic_hp = (w/20) * ((P_2 / P_1) ** 0.286 - 1)
    # Solving for P_1:  P_1 = P_2 / (adiabatic_hp * 20 / w + 1) ** (1/0.286)
    # Solving for P_2:  P_2 = P_1 * (adiabatic_hp * 20 / w + 1) ** (1/0.286)
    ratio_pow = _re.search(
        r"\(\s*(\w+)\s*/\s*" + _re.escape(target_var) + r"\s*\)"
        r"\s*\*\*\s*([\d.]+)",
        eq,
    )
    if ratio_pow:
        numerator_var = ratio_pow.group(1)
        exp_val = ratio_pow.group(2)
        inv_exp = f"1.0 / {exp_val}"
        pow_frag = ratio_pow.group(0)
        # Detect "- 1" pattern: coeff * ((A/B)**n - 1)
        minus1 = _re.search(_re.escape(pow_frag) + r"\s*-\s*1", eq)
        if minus1:
            outer_m1 = _extract_outer_coeff(rhs, "(" + pow_frag + " - 1)")
            if not outer_m1:
                outer_m1 = _extract_outer_coeff(rhs, pow_frag + " - 1")
            if outer_m1:
                ratio_expr = f"{lhs} / ({outer_m1}) + 1"
            else:
                ratio_expr = f"{lhs} + 1"
        else:
            outer = _extract_outer_coeff(rhs, pow_frag)
            if outer:
                ratio_expr = f"{lhs} / ({outer})"
            else:
                ratio_expr = lhs
        answer = f"{numerator_var} / ({ratio_expr}) ** ({inv_exp})"
        lines = [header]
        lines.append(f"    # [.pyeqn] {pyeqn}")
        lines.append(f"    # Solve for {target_var}:")
        if minus1:
            outer_label = outer_m1 if outer_m1 else "1"
            lines.append(
                f"    # Step 1: ({numerator_var} / {target_var}) ** {exp_val}"
                f" - 1 = {lhs} / ({outer_label})"
            )
            lines.append(
                f"    # Step 2: ({numerator_var} / {target_var}) ** {exp_val}"
                f" = {ratio_expr}"
            )
        else:
            lines.append(
                f"    # Step 1: ({numerator_var} / {target_var}) ** {exp_val}"
                f" = {ratio_expr}"
            )
        lines.append(f"    # Step 3: {target_var} = {answer}")
        lines.append(f"    {target_var} = {answer}")
        return "\n".join(lines)

    # ── Pattern: solve for NUMERATOR of (A / B) ** n ──
    # e.g. adiabatic_hp = (w/20) * ((P_2 / P_1) ** 0.286 - 1), solve for P_2
    ratio_pow_num = _re.search(
        r"\(\s*" + _re.escape(target_var) + r"\s*/\s*(\w+)\s*\)"
        r"\s*\*\*\s*([\d.]+)",
        eq,
    )
    if ratio_pow_num:
        denom_var = ratio_pow_num.group(1)
        exp_val = ratio_pow_num.group(2)
        inv_exp = f"1.0 / {exp_val}"
        pow_frag = ratio_pow_num.group(0)
        minus1 = _re.search(_re.escape(pow_frag) + r"\s*-\s*1", eq)
        if minus1:
            outer_m1 = _extract_outer_coeff(rhs, pow_frag + " - 1")
            if outer_m1:
                ratio_expr = f"{lhs} / ({outer_m1}) + 1"
            else:
                ratio_expr = f"{lhs} + 1"
        else:
            outer = _extract_outer_coeff(rhs, pow_frag)
            if outer:
                ratio_expr = f"{lhs} / ({outer})"
            else:
                ratio_expr = lhs
        answer = f"{denom_var} * ({ratio_expr}) ** ({inv_exp})"
        lines = [header]
        lines.append(f"    # [.pyeqn] {pyeqn}")
        lines.append(f"    # Solve for {target_var}:")
        if minus1:
            lines.append(
                f"    # Step 1: ({target_var} / {denom_var}) ** {exp_val}"
                f" = {ratio_expr}"
            )
        else:
            lines.append(
                f"    # Step 1: ({target_var} / {denom_var}) ** {exp_val}"
                f" = {ratio_expr}"
            )
        lines.append(
            f"    # Step 2: {target_var} / {denom_var} = ({ratio_expr}) ** ({inv_exp})"
        )
        lines.append(f"    # Step 3: {target_var} = {answer}")
        lines.append(f"    {target_var} = {answer}")
        return "\n".join(lines)

    # ── Pattern: target ** n (bare power, e.g. rho ** 0.84) ──
    # e.g. bhp = bhp_0 * (0.5 + 0.0155 * rho ** 0.84 * mu ** 0.16)
    bare_pow = _re.search(r"\b" + _re.escape(target_var) + r"\s*\*\*\s*([\d.]+)", eq)
    if bare_pow and total == 1:
        exp_val = bare_pow.group(1)
        inv_exp = f"1.0 / {exp_val}"
        # Parse the structure: lhs = outer_coeff * (additive_const + coeff * target**exp * cofactors)
        # We need to produce: target = ((lhs / outer - additive_const) / (coeff * cofactors)) ** inv_exp
        # Try to parse: rhs = A * (B + C * target**exp * D)
        #   where A is outer multiplicand, B is additive constant, C*D are cofactors

        # Find cofactors: everything multiplied with target**exp in the same product
        # The target**exp fragment in the expression
        pow_frag = bare_pow.group(0)  # e.g. "rho ** 0.84"

        # Try to detect if target**exp is inside parentheses with an additive term
        # Pattern: (additive + coeff * target**exp * cofactors)
        paren_match = _re.search(
            r"\(([^()]*?)\b" + _re.escape(pow_frag) + r"([^()]*?)\)",
            eq,
        )
        if paren_match:
            inside_before = paren_match.group(1).strip()
            inside_after = paren_match.group(2).strip()
            # inside_before might be "0.5 + 0.0155 * " → split on +
            # Find the additive constant and the coefficient of target
            parts = _re.split(r"\s*\+\s*", inside_before, maxsplit=1)
            if len(parts) == 2:
                additive_const = parts[0].strip()
                target_coeff = parts[1].rstrip("* ").strip()
            else:
                additive_const = None
                target_coeff = inside_before.rstrip("* ").strip()
            # cofactors after target**exp
            cofactors = inside_after.lstrip("* ").strip()
            # Full parenthesized expression
            full_paren = paren_match.group(0)
            # Outer coefficient: everything outside the parens
            outer = _extract_outer_coeff(rhs, full_paren)

            # Build concrete steps
            lines = [header]
            lines.append(f"    # [.pyeqn] {pyeqn}")
            lines.append(f"    # Solve for {target_var}:")

            step1_rhs = f"{lhs}"
            if outer:
                step1_rhs = f"{lhs} / {outer}"
                lines.append(
                    f"    # Step 1: {step1_rhs} = {inside_before}{pow_frag}{inside_after}"
                )
            else:
                lines.append(
                    f"    # Step 1: {lhs} = {inside_before}{pow_frag}{inside_after}"
                )

            if additive_const:
                iso_base = f"({step1_rhs} - {additive_const})"
                lines.append(
                    f"    # Step 2: {iso_base}"
                    f" = {target_coeff} * {pow_frag}"
                    + (f" * {cofactors}" if cofactors else "")
                )
            else:
                iso_base = step1_rhs

            coeff_divisor = ""
            if target_coeff and target_coeff != "1":
                coeff_divisor += target_coeff
            if cofactors:
                if coeff_divisor:
                    coeff_divisor += f" * {cofactors}"
                else:
                    coeff_divisor = cofactors

            if coeff_divisor:
                isolated = f"({iso_base} / ({coeff_divisor}))"
            else:
                isolated = iso_base

            answer = f"{isolated} ** ({inv_exp})"
            lines.append(f"    # Step 3: {target_var} ** {exp_val} = {isolated}")
            lines.append(f"    # Step 4: {target_var} = {answer}")
            lines.append(f"    {target_var} = {answer}")
            return "\n".join(lines)

        else:
            # target**exp is NOT inside extra parens — simpler structure
            # e.g. lhs = coeff * target ** exp
            outer = _extract_outer_coeff(rhs, pow_frag)
            if outer:
                isolated = f"{lhs} / ({outer})"
            else:
                isolated = lhs
            answer = f"({isolated}) ** ({inv_exp})"
            lines = [header]
            lines.append(f"    # [.pyeqn] {pyeqn}")
            lines.append(f"    # Solve for {target_var}:")
            lines.append(f"    # Step 1: {target_var} ** {exp_val} = {isolated}")
            lines.append(f"    # Step 2: {target_var} = {answer}")
            lines.append(f"    {target_var} = {answer}")
            return "\n".join(lines)

    # ── Pattern: fractional exponent with target in 2+ difference terms ──
    # e.g. S_p = S_Th * ((P-p_s)*(460+T_i) / ((P-p_c)*(460+T_e)))**0.6
    # Target appears in (target - X) AND (target - Y) under one **exp
    # After clearing exponent → linear in target → cross-multiply
    exp_match = _re.search(r"\*\*\s*(0\.\d+)", eq)
    if exp_match and total >= 2:
        exp_val = exp_match.group(1)
        inv_exp_f = round(1.0 / float(exp_val), 6)
        inv_exp = str(inv_exp_f)
        diff_terms = _re.findall(r"\(" + _re.escape(target_var) + r"\s*-\s*(\w+)\)", eq)
        if len(diff_terms) >= 2 and total <= 3:
            # ── Ratio-of-differences pattern (e.g. eqn_10_19 P) ──
            # After clearing exponent we get a ratio of linear terms in target.
            # We need to know which differences are in numerator vs denominator.
            # Parse: big_expr = (...numerator...) / (...denominator...)
            # Inside the **exp, find the full fraction.

            # For eqn_10_19: ((P-p_s)*(460+T_i) / ((P-p_c)*(460+T_e)))**0.6
            # lhs_ratio = S_p / S_Th
            # ratio_cleared = (S_p / S_Th) ** (1/0.6)
            # Then: ratio_cleared = (P-p_s)*(460+T_i) / ((P-p_c)*(460+T_e))
            # R * (P-p_c)*(460+T_e) = (P-p_s)*(460+T_i)

            outer = _extract_outer_coeff(rhs, exp_match.group(0))
            # outer should be like "S_Th * (...)" or just the coeff before **

            # Extract what's inside (...)**exp
            # Find the parenthesized expression before **exp
            paren_depth = 0
            exp_pos = eq.find("**" + exp_val)
            if exp_pos < 0:
                exp_pos = eq.find("** " + exp_val)
            start = exp_pos - 1
            while start >= 0 and eq[start] == " ":
                start -= 1
            if start >= 0 and eq[start] == ")":
                end_p = start
                depth = 1
                start -= 1
                while start >= 0 and depth > 0:
                    if eq[start] == ")":
                        depth += 1
                    elif eq[start] == "(":
                        depth -= 1
                    start -= 1
                inner_expr = eq[start + 2 : end_p]  # content inside outermost parens
            else:
                inner_expr = ""

            if inner_expr:
                # Find the outer multiplier of the **exp term
                # rhs = outer_mult * (inner_expr)**exp
                pow_full = f"({inner_expr})**{exp_val}"
                # Try exact or with spaces
                rhs_stripped = rhs.replace(" ", "")
                pow_stripped = pow_full.replace(" ", "")
                # Get the outer multiplier (everything in rhs except the power term)
                outer_mult = _extract_outer_coeff(rhs, f"({inner_expr})")
                # Actually simpler: the outer multiplier for the whole expression
                # For S_p = S_Th * (...)**0.6, outer_mult = S_Th
                # lhs / outer_mult = (...)**0.6
                # (lhs / outer_mult) ** (1/exp) = inner_expr

                # Identify numerator and denominator of inner_expr
                # inner_expr is like "(P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e))"
                # or "(P-p_0)*(460+T_i)*(P-p_c) / (P * (P-p_s)*(460+T_e))"
                if " / " in inner_expr:
                    num_part, den_part = inner_expr.split(" / ", 1)
                    # Strip outer parens from denominator
                    den_part = den_part.strip()
                    if den_part.startswith("(") and den_part.endswith(")"):
                        # Check if the outer parens are balanced
                        depth = 0
                        for i, c in enumerate(den_part):
                            if c == "(":
                                depth += 1
                            elif c == ")":
                                depth -= 1
                            if depth == 0 and i < len(den_part) - 1:
                                break
                        else:
                            den_part = den_part[1:-1]
                    num_part = num_part.strip()
                else:
                    num_part = inner_expr
                    den_part = ""

                # Determine which diff terms are in numerator vs denominator
                num_diffs = _re.findall(
                    r"\(" + _re.escape(target_var) + r"\s*-\s*(\w+)\)", num_part
                )
                den_diffs = _re.findall(
                    r"\(" + _re.escape(target_var) + r"\s*-\s*(\w+)\)", den_part
                )

                # Non-target factors in num and den
                def _remove_target_diffs(expr, tv):
                    """Remove (target - X) terms, leaving other factors."""
                    cleaned = _re.sub(
                        r"\(\s*" + _re.escape(tv) + r"\s*-\s*\w+\s*\)\s*\*?\s*",
                        "",
                        expr,
                    )
                    cleaned = cleaned.strip().strip("*").strip()
                    return cleaned if cleaned else "1"

                num_other = _remove_target_diffs(num_part, target_var)
                den_other = _remove_target_diffs(den_part, target_var)

                # Also handle bare target_var in denominator (e.g. "P * (P-p_s)...")
                bare_in_den = bool(
                    _re.search(r"\b" + _re.escape(target_var) + r"\b(?!\s*-)", den_part)
                )

                if bare_in_den:
                    # Target appears as a bare multiplier in denominator too
                    # This is the eqn_10_20 P case — needs brentq
                    rhs_sub = _re.sub(
                        r"\b" + _re.escape(target_var) + r"\b",
                        target_var + "_val",
                        rhs,
                    )
                    lines = [header]
                    lines.append(f"    # [.pyeqn] {pyeqn}")
                    lines.append(
                        f"    # {target_var} appears {total} times"
                        f" (including bare) — use numerical solver"
                    )
                    lines.append(f"    from scipy.optimize import brentq")
                    lines.append(f"    def _res({target_var}_val):")
                    lines.append(f"        return {rhs_sub} - {lhs}")
                    lines.append(f"    lo, hi = None, None")
                    lines.append(f"    prev = _res(max(p_0, p_c, p_s) + 1)")
                    lines.append(f"    for i in range(1, 100000):")
                    lines.append(f"        x = max(p_0, p_c, p_s) + 1 + i * 0.1")
                    lines.append(f"        try:")
                    lines.append(f"            cur = _res(x)")
                    lines.append(f"        except Exception:")
                    lines.append(f"            continue")
                    lines.append(f"        if prev * cur < 0:")
                    lines.append(f"            lo, hi = x - 0.1, x")
                    lines.append(f"            break")
                    lines.append(f"        prev = cur")
                    lines.append(f"    if lo is None:")
                    lines.append(
                        f"        raise UnsolvedException("
                        f'"No sign change found for {target_var}")'
                    )
                    lines.append(f"    {target_var} = brentq(_res, lo, hi)")
                    lines.append(f"    return [{target_var}]")
                    return "\n".join(lines)

                # Linear case: only (target - X) terms, no bare target
                # After clearing exponent:
                #   R = (lhs / outer_mult) ** (1/exp)
                #   R = num_part / den_part
                #   R * den_factors * prod(target - den_diffs) = num_factors * prod(target - num_diffs)
                # With exactly 2 diff terms (one num, one den):
                #   R * den_other * (target - den_diff) = num_other * (target - num_diff)
                #   target * (R * den_other - num_other) = R * den_other * den_diff - num_other * num_diff
                #   target = (R * den_other * den_diff - num_other * num_diff) / (R * den_other - num_other)

                # Find the outer coefficient before the **exp term
                # rhs might be "S_Th * (inner)**0.6" → outer_mult = "S_Th"
                # lhs/outer = (inner)**0.6 → R = (lhs/outer)**(1/0.6)
                outer_mult = _extract_outer_coeff(rhs, f"({inner_expr})**{exp_val}")
                if not outer_mult:
                    # Try with spaces around **
                    outer_mult = _extract_outer_coeff(
                        rhs, f"({inner_expr}) ** {exp_val}"
                    )
                # Also try matching without the exp
                if not outer_mult:
                    outer_mult = _extract_outer_coeff(rhs, f"({inner_expr})")

                if outer_mult:
                    R_expr = f"({lhs} / ({outer_mult})) ** ({inv_exp})"
                else:
                    R_expr = f"({lhs}) ** ({inv_exp})"

                lines = [header]
                lines.append(f"    # [.pyeqn] {pyeqn}")
                lines.append(f"    # Solve for {target_var}:")
                lines.append(
                    f"    # Step 1: ({lhs} / {outer_mult if outer_mult else '1'})"
                    f" ** ({inv_exp}) = {inner_expr}"
                )
                lines.append(f"    R = {R_expr}")

                if len(num_diffs) == 1 and len(den_diffs) == 1:
                    nd = num_diffs[0]
                    dd = den_diffs[0]
                    # R * den_other * (target - dd) = num_other * (target - nd)
                    lines.append(
                        f"    # Step 2: R * ({den_other}) * ({target_var} - {dd})"
                        f" = ({num_other}) * ({target_var} - {nd})"
                    )
                    lines.append(
                        f"    # Step 3: {target_var} * (R * ({den_other}) - ({num_other}))"
                        f" = R * ({den_other}) * {dd} - ({num_other}) * {nd}"
                    )
                    answer = (
                        f"(R * ({den_other}) * {dd} - ({num_other}) * {nd})"
                        f" / (R * ({den_other}) - ({num_other}))"
                    )
                    lines.append(f"    # Step 4: {target_var} = {answer}")
                    lines.append(f"    {target_var} = {answer}")
                    return "\n".join(lines)

                elif len(num_diffs) == 2 and len(den_diffs) == 1:
                    # e.g. eqn_10_20 with P in num twice and den once
                    # but P also bare in den → already handled by brentq above
                    # This case shouldn't occur for our equations, fallthrough
                    pass
                elif len(num_diffs) >= 1 and len(den_diffs) >= 1:
                    # General: produce brentq as safe fallback
                    pass

    # ── Pattern: single-occurrence target inside (target - X) with **exp ──
    # e.g. eqn_10_20 solve for p_0: (P - p_0) appears once
    # Clear exponent, isolate (P - p_0), solve for p_0
    if exp_match and total == 1:
        exp_val = exp_match.group(1)
        inv_exp_f = round(1.0 / float(exp_val), 6)
        inv_exp = str(inv_exp_f)

        # Parse the inner expression of (**exp) to find target's role
        exp_pos = eq.find("**" + exp_val)
        if exp_pos < 0:
            exp_pos = eq.find("** " + exp_val)
        if exp_pos >= 0:
            start = exp_pos - 1
            while start >= 0 and eq[start] == " ":
                start -= 1
            if start >= 0 and eq[start] == ")":
                end_p = start
                depth = 1
                start -= 1
                while start >= 0 and depth > 0:
                    if eq[start] == ")":
                        depth += 1
                    elif eq[start] == "(":
                        depth -= 1
                    start -= 1
                inner_expr = eq[start + 2 : end_p]
            else:
                inner_expr = ""
        else:
            inner_expr = ""

        if inner_expr:
            # Get outer multiplier
            outer_mult = _extract_outer_coeff(rhs, f"({inner_expr})**{exp_val}")
            if not outer_mult:
                outer_mult = _extract_outer_coeff(rhs, f"({inner_expr}) ** {exp_val}")
            if not outer_mult:
                outer_mult = _extract_outer_coeff(rhs, f"({inner_expr})")

            if outer_mult:
                R_expr = f"({lhs} / ({outer_mult})) ** ({inv_exp})"
            else:
                R_expr = f"({lhs}) ** ({inv_exp})"

            # inner_expr after clearing exponent = R
            # Find target in inner_expr
            if " / " in inner_expr:
                num_part, den_part = inner_expr.split(" / ", 1)
                den_part = den_part.strip()
                if den_part.startswith("(") and den_part.endswith(")"):
                    depth = 0
                    for i, c in enumerate(den_part):
                        if c == "(":
                            depth += 1
                        elif c == ")":
                            depth -= 1
                        if depth == 0 and i < len(den_part) - 1:
                            break
                    else:
                        den_part = den_part[1:-1]
            else:
                num_part = inner_expr
                den_part = ""

            target_in_num = bool(
                _re.search(r"\b" + _re.escape(target_var) + r"\b", num_part)
            )
            target_in_den = bool(
                _re.search(r"\b" + _re.escape(target_var) + r"\b", den_part)
            )

            # Helper to build factor_expr and answer for a diff term in numerator
            def _handle_num_diff(match_obj, sign, other_var, num_part, den_part):
                """Build scaffold lines for target found in numerator diff term.
                sign='+' means (target - other), sign='-' means (other - target)."""
                num_remaining = num_part.replace(match_obj.group(0), "1")
                num_remaining = _re.sub(r"\b1\s*\*\s*", "", num_remaining)
                num_remaining = _re.sub(r"\s*\*\s*1\b", "", num_remaining)
                num_remaining = num_remaining.strip().strip("*").strip()
                if not num_remaining or num_remaining == "1":
                    num_remaining = ""

                if den_part and num_remaining:
                    factor_expr = f"R * ({den_part}) / ({num_remaining})"
                elif den_part:
                    factor_expr = f"R * ({den_part})"
                elif num_remaining:
                    factor_expr = f"R / ({num_remaining})"
                else:
                    factor_expr = "R"

                if sign == "+":
                    diff_str = f"({target_var} - {other_var})"
                    answer = f"{other_var} + {factor_expr}"
                else:
                    diff_str = f"({other_var} - {target_var})"
                    answer = f"{other_var} - {factor_expr}"

                lines = [header]
                lines.append(f"    # [.pyeqn] {pyeqn}")
                lines.append(f"    # Solve for {target_var}:")
                lines.append(f"    R = {R_expr}")
                lines.append(f"    # After clearing **{exp_val}: R = {inner_expr}")
                lines.append(f"    # {diff_str} = {factor_expr}")
                lines.append(f"    # {target_var} = {answer}")
                lines.append(f"    {target_var} = {answer}")
                return "\n".join(lines)

            # Case: (X - target) in numerator
            rev_diff_num = _re.search(
                r"\(\s*(\w+)\s*-\s*" + _re.escape(target_var) + r"\s*\)", num_part
            )
            if rev_diff_num:
                return _handle_num_diff(
                    rev_diff_num, "-", rev_diff_num.group(1), num_part, den_part
                )

            # Case: (target - X) in numerator
            fwd_diff_num = _re.search(
                r"\(\s*" + _re.escape(target_var) + r"\s*-\s*(\w+)\s*\)", num_part
            )
            if fwd_diff_num:
                return _handle_num_diff(
                    fwd_diff_num, "+", fwd_diff_num.group(1), num_part, den_part
                )

            # Case: (target - X) in denominator
            diff_den = _re.search(
                r"\(\s*" + _re.escape(target_var) + r"\s*-\s*(\w+)\s*\)", den_part
            )
            if diff_den:
                other_var = diff_den.group(1)
                den_remaining = den_part.replace(diff_den.group(0), "1")
                den_remaining = _re.sub(r"\b1\s*\*\s*", "", den_remaining)
                den_remaining = _re.sub(r"\s*\*\s*1\b", "", den_remaining)
                den_remaining = den_remaining.strip().strip("*").strip()
                if not den_remaining or den_remaining == "1":
                    den_remaining = ""

                # (target - other) = num_part / (R * den_remaining)
                if num_part and den_remaining:
                    factor_expr = f"({num_part}) / (R * ({den_remaining}))"
                elif num_part:
                    factor_expr = f"({num_part}) / R"
                elif den_remaining:
                    factor_expr = f"1.0 / (R * ({den_remaining}))"
                else:
                    factor_expr = f"1.0 / R"

                answer = f"{other_var} + {factor_expr}"

                lines = [header]
                lines.append(f"    # [.pyeqn] {pyeqn}")
                lines.append(f"    # Solve for {target_var}:")
                lines.append(f"    R = {R_expr}")
                lines.append(f"    # After clearing **{exp_val}: R = {inner_expr}")
                lines.append(f"    # ({target_var} - {other_var}) = {factor_expr}")
                lines.append(f"    # {target_var} = {answer}")
                lines.append(f"    {target_var} = {answer}")
                return "\n".join(lines)

            # Case: (X - target) in denominator
            rev_diff_den = _re.search(
                r"\(\s*(\w+)\s*-\s*" + _re.escape(target_var) + r"\s*\)", den_part
            )
            if rev_diff_den:
                other_var = rev_diff_den.group(1)
                den_remaining = den_part.replace(rev_diff_den.group(0), "1")
                den_remaining = _re.sub(r"\b1\s*\*\s*", "", den_remaining)
                den_remaining = _re.sub(r"\s*\*\s*1\b", "", den_remaining)
                den_remaining = den_remaining.strip().strip("*").strip()
                if not den_remaining or den_remaining == "1":
                    den_remaining = ""

                # R * den = num, so (other - target) = num / (R * den_remaining)
                if num_part and den_remaining:
                    factor_expr = f"({num_part}) / (R * ({den_remaining}))"
                elif num_part:
                    factor_expr = f"({num_part}) / R"
                elif den_remaining:
                    factor_expr = f"1.0 / (R * ({den_remaining}))"
                else:
                    factor_expr = f"1.0 / R"

                answer = f"{other_var} - {factor_expr}"

                lines = [header]
                lines.append(f"    # [.pyeqn] {pyeqn}")
                lines.append(f"    # Solve for {target_var}:")
                lines.append(f"    R = {R_expr}")
                lines.append(f"    # After clearing **{exp_val}: R = {inner_expr}")
                lines.append(f"    # ({other_var} - {target_var}) = {factor_expr}")
                lines.append(f"    # {target_var} = {answer}")
                lines.append(f"    {target_var} = {answer}")
                return "\n".join(lines)

            # Case: (C + target) in numerator or denominator
            add_num = _re.search(
                r"\((\d+)\s*\+\s*" + _re.escape(target_var) + r"\)", num_part
            )
            add_den = _re.search(
                r"\((\d+)\s*\+\s*" + _re.escape(target_var) + r"\)", den_part
            )
            if add_num:
                const = add_num.group(1)
                num_remaining = num_part.replace(add_num.group(0), "1")
                num_remaining = _re.sub(r"\b1\s*\*\s*", "", num_remaining)
                num_remaining = _re.sub(r"\s*\*\s*1\b", "", num_remaining)
                num_remaining = num_remaining.strip().strip("*").strip()
                if not num_remaining or num_remaining == "1":
                    num_remaining = ""

                if den_part and num_remaining:
                    factor_expr = f"R * ({den_part}) / ({num_remaining})"
                elif den_part:
                    factor_expr = f"R * ({den_part})"
                elif num_remaining:
                    factor_expr = f"R / ({num_remaining})"
                else:
                    factor_expr = "R"

                answer = f"{factor_expr} - {const}"

                lines = [header]
                lines.append(f"    # [.pyeqn] {pyeqn}")
                lines.append(f"    # Solve for {target_var}:")
                lines.append(f"    R = {R_expr}")
                lines.append(f"    # ({const} + {target_var}) = {factor_expr}")
                lines.append(f"    # {target_var} = {answer}")
                lines.append(f"    {target_var} = {answer}")
                return "\n".join(lines)
            elif add_den:
                const = add_den.group(1)
                den_remaining = den_part.replace(add_den.group(0), "1")
                den_remaining = _re.sub(r"\b1\s*\*\s*", "", den_remaining)
                den_remaining = _re.sub(r"\s*\*\s*1\b", "", den_remaining)
                den_remaining = den_remaining.strip().strip("*").strip()
                if not den_remaining or den_remaining == "1":
                    den_remaining = ""

                if num_part and den_remaining:
                    factor_expr = f"({num_part}) / (R * ({den_remaining}))"
                elif num_part:
                    factor_expr = f"({num_part}) / R"
                else:
                    factor_expr = f"1.0 / R"

                answer = f"{factor_expr} - {const}"

                lines = [header]
                lines.append(f"    # [.pyeqn] {pyeqn}")
                lines.append(f"    # Solve for {target_var}:")
                lines.append(f"    R = {R_expr}")
                lines.append(f"    # ({const} + {target_var}) = {factor_expr}")
                lines.append(f"    # {target_var} = {answer}")
                lines.append(f"    {target_var} = {answer}")
                return "\n".join(lines)

    # ── Fallback: general single-occurrence → simple rearrangement ──
    if total == 1:
        lines = [header]
        lines.append(f"    # [.pyeqn] {pyeqn}")
        lines.append(f"    # Solve for {target_var} by rearranging the equation")
        lines.append(f"    {target_var} =")
        return "\n".join(lines)

    # ── Fallback: multiple occurrences → brentq numerical ──
    if total >= 2:
        rhs_sub = _re.sub(
            r"\b" + _re.escape(target_var) + r"\b",
            target_var + "_val",
            rhs,
        )
        lines = [header]
        lines.append(f"    # [.pyeqn] {pyeqn}")
        lines.append(f"    # {target_var} appears {total} times — use numerical solver")
        lines.append(f"    from scipy.optimize import brentq")
        lines.append(f"    def _res({target_var}_val):")
        lines.append(f"        return {rhs_sub} - {lhs}")
        lines.append(f"    lo, hi = None, None")
        lines.append(f"    prev = _res(0.01)")
        lines.append(f"    for i in range(1, 100000):")
        lines.append(f"        x = i * 0.01")
        lines.append(f"        try:")
        lines.append(f"            cur = _res(x)")
        lines.append(f"        except Exception:")
        lines.append(f"            continue")
        lines.append(f"        if prev * cur < 0:")
        lines.append(f"            lo, hi = x - 0.01, x")
        lines.append(f"            break")
        lines.append(f"        prev = cur")
        lines.append(f"    if lo is None:")
        lines.append(
            f'        raise UnsolvedException("No sign change found for {target_var}")'
        )
        lines.append(f"    {target_var} = brentq(_res, lo, hi)")
        lines.append(f"    return [{target_var}]")
        return "\n".join(lines)

    return ""


def attempt_repair_shard(
    ctx: PipelineContext,
    family_name: str,
    broken_variant: str,
    trusted_variants: list,
    scores: dict,
    mismatches: dict,
    error: str = None,
):
    """Repairs a single shard file by showing the LLM a working shard as example."""
    family_dir = os.path.join(ctx.shards_dir, family_name)
    eqn_suffix = family_name.split("_eqn_")[1]
    safe_variant = case_safe_name(broken_variant)
    shard_file = f"eqn_{eqn_suffix}__{safe_variant}.py"
    shard_path = os.path.join(family_dir, shard_file)

    if not os.path.exists(shard_path):
        print(f"  [Repair] Shard {shard_file} not found in {family_name}")
        return {"updated": False}

    with open(shard_path, "r") as f:
        shard_code = f.read()

    pyeqn_match = re.search(r"# \[\.pyeqn\] (.*)", shard_code)
    pyeqn = pyeqn_match.group(1).strip() if pyeqn_match else ""

    # Find a working shard from the same family to use as example
    example_code = None
    for tv in trusted_variants:
        safe_tv = case_safe_name(tv)
        tv_path = os.path.join(family_dir, f"eqn_{eqn_suffix}__{safe_tv}.py")
        if os.path.exists(tv_path):
            with open(tv_path, "r") as f:
                candidate = f.read()
            if "raise UnsolvedException" not in candidate:
                # Strip import lines so Phi-3 doesn't fall back to SymPy
                cleaned_lines = []
                for line in candidate.splitlines():
                    s = line.strip()
                    if s.startswith("import ") or s.startswith("from "):
                        continue
                    cleaned_lines.append(line)
                example_code = "\n".join(cleaned_lines).strip()
                break

    # Extract the function header from the broken shard
    header_match = re.search(
        r"(def eqn_.*?\*\*kwargs\s*,?\s*\)):", shard_code, re.DOTALL
    )
    if header_match:
        # Normalize multi-line headers to single line
        header = re.sub(r"\s+", " ", header_match.group(1)).strip() + ":"
    else:
        header = ""

    # Generate algebraic derivation scaffold (fill-in-the-blank for Phi-3)
    scaffold = _generate_derivation_scaffold(pyeqn, broken_variant, header)

    print(f"\n |- Repairing shard: {shard_file} in family {family_name}")

    # ── Check if scaffold is already a complete solution (bypass LLM) ──
    if scaffold and not scaffold.strip().endswith(f"{broken_variant} ="):
        # Scaffold has a complete expression — either has return[] or needs one added
        has_return = "return [" in scaffold
        has_assignment = f"{broken_variant} =" in scaffold

        if has_assignment or has_return:
            print(f"  [Repair] Scaffold provides complete solution — bypassing LLM")
            code_text = scaffold
            if not has_return:
                code_text += f"\n    return [{broken_variant}]"

            # Build the complete shard file
            import_header = (
                "from cmath import log, sqrt, exp\n"
                "from math import e, pi\n"
                "from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp\n"
                "from scipy.optimize import newton\n"
                "import numpy as np\n"
                "from vakyume.config import UnsolvedException\n\n"
            )

            # Strip only top-level import lines from the scaffold
            # (keep indented imports like 'from scipy.optimize import brentq'
            # that appear inside function bodies)
            body_lines = []
            for line in code_text.splitlines():
                stripped = line.lstrip()
                if line == stripped and (
                    stripped.startswith("import ") or stripped.startswith("from ")
                ):
                    continue
                body_lines.append(line)
            code_text = "\n".join(body_lines)

            # Ensure the pyeqn comment is present
            if "# [.pyeqn]" not in code_text and pyeqn:
                code_text = re.sub(
                    r"(def\s+\w+\(.*?\):)",
                    rf"\1\n    # [.pyeqn] {pyeqn}",
                    code_text,
                    count=1,
                )

            full_code = import_header + code_text + "\n"

            try:
                ast.parse(full_code)
                with open(shard_path, "w") as f:
                    f.write(full_code)
                print(
                    f"  [Repair] Successfully wrote scaffold-derived shard: {shard_file}"
                )
                return {"updated": True}
            except SyntaxError as e:
                print(f"  [Repair] Scaffold syntax error: {e}")
                print(f"  [Repair] Falling through to LLM repair...")

    if scaffold:
        # --- Scaffold approach: give Phi-3 the partial function to complete ---
        system_prompt = (
            "Complete the Python function. Output ONLY the function. "
            "No text. No markdown. Do NOT use sympy. Return as return [value]."
        )

        user_prompt = ""
        if example_code:
            user_prompt += f"Working example (solves for a different variable):\n{example_code}\n\n"

        user_prompt += f"Complete this function:\n{scaffold}\n"

        # Include mismatch diagnostics if available (for retry rounds)
        if mismatches and broken_variant in mismatches:
            variant_mismatches = mismatches[broken_variant]
            if variant_mismatches:
                user_prompt += "\nPrevious attempt failed with these errors:\n"
                for i, trial in enumerate(variant_mismatches[:2]):
                    if "error" in trial:
                        user_prompt += f"  Error: {trial['error']}\n"
                    elif "inputs" in trial:
                        user_prompt += f"  Inputs: {trial['inputs']}\n"
                        user_prompt += f"  Got: {trial.get('output')}\n"
                        for m in trial.get("mismatches", []):
                            if "expected" in m:
                                user_prompt += f"  Expected {m['target']}={m['expected']}, got {m.get('got')}\n"
    else:
        # --- Fallback: free-form prompt when no scaffold pattern matched ---
        system_prompt = (
            "You write Python functions that solve equations. Output ONLY code. No text.\n"
            "Do NOT use sympy. Use direct algebraic rearrangement.\n"
            "Return result as: return [value]\n"
            "Use cmath.sqrt/cmath.log/cmath.exp for complex-safe math.\n"
            "If the variable is in the exponent (transcendental), use "
            "scipy.optimize.brentq with a sign-change scan."
        )

        user_prompt = f"Equation: {pyeqn}\n"
        if example_code:
            user_prompt += f"\nWorking example (solves for a different variable):\n{example_code}\n"

        # Include mismatch diagnostics if available
        if mismatches and broken_variant in mismatches:
            variant_mismatches = mismatches[broken_variant]
            if variant_mismatches:
                user_prompt += "\nPrevious attempt failed with these errors:\n"
                for i, trial in enumerate(variant_mismatches[:2]):
                    if "error" in trial:
                        user_prompt += f"  Error: {trial['error']}\n"
                    elif "inputs" in trial:
                        user_prompt += f"  Inputs: {trial['inputs']}\n"
                        user_prompt += f"  Got: {trial.get('output')}\n"
                        for m in trial.get("mismatches", []):
                            if "expected" in m:
                                user_prompt += f"  Expected {m['target']}={m['expected']}, got {m.get('got')}\n"

        user_prompt += f"\nNow solve for: {broken_variant}\n"
        user_prompt += f"Use this exact header:\n{header}\n"

    from .llm import ask_llm

    raw = ask_llm(system_prompt, user_prompt, stream=False)

    if not raw:
        print(f"  [Repair] LLM returned empty response")
        return {"updated": False}

    method_name = f"eqn_{eqn_suffix}__{broken_variant}"
    code_text = extract_code(raw, target_name=method_name)
    if not code_text.strip():
        # Fallback: try to use whatever came back without target filtering
        code_text = extract_code(raw)
    if not code_text.strip():
        print(f"  [Repair] extract_code returned empty")
        return {"updated": False}

    # --- Normalize the function name ---
    # phi3 may output a different function name; replace it with the correct one.
    any_def = re.search(r"def\s+(\w+)\s*\((.*?)\):", code_text, re.DOTALL)
    if any_def:
        llm_func_name = any_def.group(1)
        params_str = any_def.group(2).strip()

        # Ensure 'self' is the first parameter
        if not params_str.startswith("self"):
            params_str = "self, " + params_str if params_str else "self"

        # Ensure **kwargs is present
        if "**kwargs" not in params_str:
            params_str = params_str.rstrip(", ") + ", **kwargs"

        # Replace the def line with the correct function name and params
        code_text = code_text.replace(
            any_def.group(0), f"def {method_name}({params_str}):"
        )

    # --- Build the complete shard file ---
    import_header = (
        "from cmath import log, sqrt, exp\n"
        "from math import e, pi\n"
        "from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp\n"
        "from scipy.optimize import newton\n"
        "import numpy as np\n"
        "from vakyume.config import UnsolvedException\n\n"
    )

    # Strip any imports the LLM added (we provide our own)
    body_lines = []
    for line in code_text.splitlines():
        s = line.strip()
        if s.startswith("import ") or s.startswith("from "):
            continue
        body_lines.append(line)
    code_text = "\n".join(body_lines)

    # Ensure the pyeqn comment is present
    if "# [.pyeqn]" not in code_text and pyeqn:
        # Insert right after the def line
        code_text = re.sub(
            r"(def\s+\w+\(.*?\):)",
            rf"\1\n    # [.pyeqn] {pyeqn}",
            code_text,
            count=1,
        )

    # Ensure return value is wrapped in a list
    # Phi-3 often writes "return value" instead of "return [value]"
    code_text = re.sub(
        r"return\s+(?!\[)(\w[\w.]*(?:\s*\*\*\s*[\w.()]+)?)\s*$",
        r"return [\1]",
        code_text,
        flags=re.MULTILINE,
    )

    full_code = import_header + code_text + "\n"

    try:
        ast.parse(full_code)
        with open(shard_path, "w") as f:
            f.write(full_code)
        print(f"  [Repair] Successfully wrote repaired shard: {shard_file}")
        return {"updated": True}
    except SyntaxError as e:
        print(f"  [Repair] LLM produced invalid syntax: {e}")
        print(f"  [Repair] Code was:\n{full_code[:500]}")
        return {"updated": False, "error": str(e)}


def run_pipeline(
    project_dir=".",
    max_rounds=DEFAULT_MAX_ROUNDS,
    auto_repair=True,
    overwrite=False,
    verbose=False,
    repair_only=False,
):
    ctx = PipelineContext(project_dir)

    # --repair-only and --overwrite are mutually exclusive in practice
    if repair_only and overwrite:
        print(
            "Warning: --repair-only is ignored when --overwrite is used (shards are being wiped)."
        )
        repair_only = False

    # Load previously solved families when in repair-only mode
    previously_solved = set()
    if repair_only:
        analysis_path = os.path.join(ctx.reports_dir, "analysis.json")
        if os.path.exists(analysis_path):
            with open(analysis_path, "r") as af:
                prior_analysis = json.load(af)
            previously_solved = set(prior_analysis.get("solved", []))

            # Remove families that still have placeholder shards — they need
            # re-verification and repair even if they were "solved" before.
            placeholder_families = set()
            for fam in list(previously_solved):
                if "_eqn_" not in fam:
                    continue
                fam_dir = os.path.join(ctx.shards_dir, fam)
                if not os.path.isdir(fam_dir):
                    continue
                for fname in os.listdir(fam_dir):
                    if not fname.endswith(".py") or "__" not in fname:
                        continue
                    fpath = os.path.join(fam_dir, fname)
                    with open(fpath, "r") as _pf:
                        if 'raise UnsolvedException("Pending' in _pf.read():
                            placeholder_families.add(fam)
                            break
            if placeholder_families:
                previously_solved -= placeholder_families
                print(
                    f"  {len(placeholder_families)} 'solved' families have placeholder shards and will be re-verified."
                )

            num_broken = len(prior_analysis.get("inconsistent", {})) + len(
                prior_analysis.get("failed", [])
            )
            print(
                f"Repair-only mode: {len(previously_solved)} families previously solved, {num_broken} to repair."
            )
        else:
            print(
                "Warning: No prior analysis.json found. Falling back to full pipeline."
            )
            repair_only = False

    if not repair_only:
        if overwrite:
            ctx.clear_shards()
        elif not os.path.exists(ctx.shards_dir) or not os.listdir(ctx.shards_dir):
            ctx.clear_shards()
        shard_from_chapters(ctx)

    for round_idx in range(1, max_rounds + 1):
        print(f"\n--- Round {round_idx}/{max_rounds} ---")
        skip = previously_solved if repair_only else None
        results = verify_all_shards(ctx, verbose=verbose, skip_families=skip)
        analysis = analyze_results(ctx, results)

        print(
            f"Status: {len(analysis['solved'])} families solved, "
            f"{len(analysis['inconsistent'])} inconsistent, "
            f"{len(analysis['failed'])} failed."
        )

        assemble_certified_library(ctx, analysis)

        if not analysis.get("inconsistent") and not analysis.get("failed"):
            print("All equations aligned. Stopping early.")
            break

        # Focus on ONE family at a time for repair
        target_family = None
        if analysis["inconsistent"]:
            target_family = sorted(analysis["inconsistent"].keys())[0]
            repair_info = analysis["inconsistent"][target_family]
            repair_type = "inconsistent"
        elif analysis["failed"]:
            target_family = analysis["failed"][0]["file"]
            repair_info = {
                "broken": [],  # Will be populated
                "trusted": [],
                "scores": {},
                "mismatches": {},
                "error": analysis["failed"][0].get("error"),
            }
            # Populate variants for failed family
            family_dir = os.path.join(ctx.shards_dir, target_family)
            if os.path.exists(family_dir):
                repair_info["broken"] = [
                    f.split("__")[1].replace(".py", "")
                    for f in os.listdir(family_dir)
                    if "__" in f
                ]
            repair_type = "failed"

        if not target_family:
            break

        attempt = 1
        max_attempts = 3
        while attempt <= max_attempts:
            print(
                f" |- [Repair] Family: {target_family} | Attempt: {attempt}/{max_attempts}"
            )

            for broken_variant in repair_info["broken"]:
                attempt_repair_shard(
                    ctx,
                    target_family,
                    broken_variant,
                    repair_info.get("trusted", []),
                    repair_info.get("scores", {}),
                    repair_info.get("mismatches", {}),
                    error=repair_info.get("error"),
                )

            # Re-verify immediately
            print(f" |- Re-verifying family {target_family}...")
            new_res = verify_family(ctx, target_family, verbose=verbose)

            if not new_res.get("error"):
                eqn_name = f"eqn_{target_family.split('_eqn_')[1]}"
                scores = new_res["results"].get(eqn_name, {})
                num_v = len(scores)

                # Check for remaining placeholder shards even if scores look perfect
                still_placeholder = []
                eqn_suffix = (
                    target_family.split("_eqn_")[1]
                    if "_eqn_" in target_family
                    else None
                )
                if eqn_suffix:
                    fam_dir = os.path.join(ctx.shards_dir, target_family)
                    for v in scores:
                        sv = case_safe_name(v)
                        sp = os.path.join(fam_dir, f"eqn_{eqn_suffix}__{sv}.py")
                        if os.path.exists(sp):
                            with open(sp, "r") as _sf:
                                _content = _sf.read()
                            if 'raise UnsolvedException("Pending' in _content:
                                still_placeholder.append(v)

                all_real = (
                    num_v > 0
                    and all(s == num_v for s in scores.values())
                    and not still_placeholder
                )
                if all_real:
                    print(f" |- SUCCESS: Family {target_family} is now certified.")
                    # Track newly solved family so it's skipped in subsequent rounds
                    if repair_only:
                        previously_solved.add(target_family)
                    break
                else:
                    score_broken = [v for v, s in scores.items() if s < num_v]
                    # Merge in any placeholders the verifier missed
                    all_broken = list(set(score_broken + still_placeholder))
                    if still_placeholder:
                        print(f" |- Still has placeholders: {still_placeholder}")
                    else:
                        print(f" |- Still inconsistent: {scores}")
                    # Update repair_info for next attempt
                    repair_info["scores"] = scores
                    repair_info["broken"] = all_broken
                    repair_info["trusted"] = [
                        v
                        for v, s in scores.items()
                        if s == num_v and v not in still_placeholder
                    ]
                    repair_info["mismatches"] = new_res.get("mismatches", {})
            else:
                print(f" |- Still failing: {new_res.get('error')}")
                repair_info["error"] = new_res.get("error")

            attempt += 1
            if COOLDOWN_SECONDS > 0:
                time.sleep(COOLDOWN_SECONDS)

    return analysis


def assemble_certified_library(ctx: PipelineContext, analysis):
    """Combines solved families into a single library file."""
    with open(ctx.certified_file, "w") as out:
        out.write(
            "from cmath import log, sqrt, exp\n"
            "from math import e, pi\n"
            "from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp\n"
            "from scipy.optimize import newton\n"
            "from vakyume.kwasak import kwasak_static\n"
            "from vakyume.config import UnsolvedException\n"
            "import numpy as np\n\n"
        )
        class_groups = {}

        for family_name in sorted(analysis["solved"]):
            class_name = family_name.split("_")[0]
            family_dir = os.path.join(ctx.shards_dir, family_name)

            if class_name not in class_groups:
                class_groups[class_name] = []

            seen_methods = set()

            # Read all methods from the family
            for sf in sorted(os.listdir(family_dir)):
                if not sf.endswith(".py"):
                    continue
                shard_path = os.path.join(family_dir, sf)

                with open(shard_path, "r") as f:
                    content = f.read()
                    try:
                        tree = ast.parse(content)
                    except:
                        continue

                    for node in tree.body:
                        if isinstance(node, ast.FunctionDef):
                            if node.name not in seen_methods:
                                seen_methods.add(node.name)
                                method_source = get_standalone_method_source(
                                    shard_path, node.name
                                )
                                if method_source:
                                    # Indent it for the class
                                    indented = "\n".join(
                                        [
                                            f"{TAB}{line}"
                                            for line in method_source.splitlines()
                                        ]
                                    )
                                    class_groups[class_name].append(indented)
                        elif isinstance(node, ast.ClassDef):
                            # This is the main shard's class
                            for item in node.body:
                                if isinstance(item, ast.FunctionDef):
                                    if item.name not in seen_methods:
                                        seen_methods.add(item.name)
                                        method_source = get_method_source_from_class(
                                            shard_path, class_name, item.name
                                        )
                                        if method_source:
                                            class_groups[class_name].append(
                                                method_source
                                            )

        for class_name, bodies in class_groups.items():
            out.write(f"class {class_name}:\n")
            for body in bodies:
                out.write(body)
                out.write("\n")


def get_standalone_method_source(file_path, method_name):
    """Extracts a top-level function source."""
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
                # Check if it's a new top level thing
                stripped = lines[end_idx].lstrip()
                if stripped.startswith(("def ", "class ", "@")):
                    break
        end_idx += 1
    return "".join(lines[start_idx:end_idx])


def get_method_source_from_class(file_path, class_name, method_name):
    """Helper to extract method source from class in file."""
    with open(file_path, "r") as f:
        lines = f.readlines()

    # Find the class first
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
        # If we hit another top-level class or def, we're done with this class
        if line.strip() and not line.startswith(" "):
            break

    if target_idx == -1:
        return ""

    # Go up to find decorators
    start_idx = target_idx
    while start_idx > class_idx + 1 and lines[start_idx - 1].strip().startswith("@"):
        start_idx -= 1

    # Find end of method by indentation
    line_with_indent = lines[target_idx]
    indent_len = len(line_with_indent) - len(line_with_indent.lstrip())

    end_idx = target_idx + 1
    while end_idx < len(lines):
        line = lines[end_idx]
        if line.strip():
            curr_indent = len(line) - len(line.lstrip())
            if curr_indent <= indent_len:
                # Check if it's a new method or class
                stripped = line.lstrip()
                if stripped.startswith(("def ", "class ", "@")):
                    break
                # Also if indent is 0, it's definitely over
                if curr_indent == 0:
                    break
        end_idx += 1

    return "".join(lines[start_idx:end_idx])


def extract_code_block(text):
    fence = "```"
    if fence not in text:
        return text
    in_code = False
    code_lines = []
    for line in text.splitlines():
        if line.strip().startswith(fence):
            if in_code:
                break
            in_code = True
            continue
