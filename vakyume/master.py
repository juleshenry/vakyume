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
from sympy import Symbol, solve, Eq
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


class PipelineContext:
    def __init__(self, project_dir):
        self.project_dir = os.path.abspath(project_dir)
        self.notes_dir = os.path.join(self.project_dir, "notes")
        self.shards_dir = os.path.join(self.project_dir, "shards")
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
    def get_solns_vanilla_nf(self, nf: str, symb: Symbol):
        try:
            solns = solve(nf, symb)
            if not solns:
                raise UnsolvedException("Sympy solve returned empty")
            return solns
        except Exception:
            raise UnsolvedException("Sympy solve failed")

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
        return f"({parts[1].split('#')[0].strip()}) - ({parts[0].strip()})"

    def sympy_failover(self, eqn_header, normal_form, token, pyeqn=None):
        code = [
            f"{TAB * 2}def func({token}):",
        ]
        if pyeqn:
            code.append(f"{TAB * 3}# [.pyeqn] {pyeqn}")
        code.extend(
            [
                f"{TAB * 3}# Placeholder for numerical solver",
                f"{TAB * 3}pass",
                f'{TAB * 2}raise UnsolvedException("Pending LLM/Manual Repair")',
            ]
        )
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
        "from math import log, sqrt, exp, pow, e\n"
        "from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp\n"
        "from scipy.optimize import newton\n"
        "import numpy as np\n"
        "from kwasak import kwasak_static\n\n"
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

                # Main shard for the equation (kwasak entry point)
                main_shard_path = os.path.join(family_dir, f"eqn_{eqn_number}.py")
                if not os.path.exists(main_shard_path):
                    with open(main_shard_path, "w") as sf:
                        sf.write(import_header)
                        sf.write(f"class {class_name}:\n")
                        sf.write(f"{TAB}@kwasak_static\n")
                        sf.write(
                            f"{TAB}def eqn_{eqn_number}({', '.join(f'{t}=None' for t in tokes)}, **kwargs):\n"
                        )
                        sf.write(f"{TAB * 2}return\n")

                # Shards for each variable
                for token in tokes:
                    shard_name = f"eqn_{eqn_number}__{token}.py"
                    shard_path = os.path.join(family_dir, shard_name)
                    if os.path.exists(shard_path):
                        continue

                    created_count += 1
                    with open(shard_path, "w") as sf:
                        sf.write(import_header)
                        sf.write(f"class {class_name}:\n")
                        other_args = [t for t in tokes if t != token]
                        header = f"{TAB}@staticmethod\n{TAB}def eqn_{eqn_number}__{token}({', '.join(f'{t}: float' for t in other_args)}, **kwargs):"
                        sf.write(header + "\n")
                        sf.write(f"{TAB * 2}# [.pyeqn] {line.strip()}\n")

                        try:
                            solns = solver.get_solns_vanilla_nf(nf, Symbol(token))
                            if solns:
                                sf.write(f"{TAB * 2}result = []\n")
                                for soln in solns:
                                    sf.write(f"{TAB * 2}{token} = {soln}\n")
                                    sf.write(f"{TAB * 2}result.append({token})\n")
                                sf.write(f"{TAB * 2}return result\n")
                            else:
                                sf.write(
                                    solver.sympy_failover(
                                        header, nf, token, pyeqn=line.strip()
                                    )
                                    + "\n"
                                )
                        except Exception:
                            sf.write(
                                solver.sympy_failover(
                                    header, nf, token, pyeqn=line.strip()
                                )
                                + "\n"
                            )

                    # py_compile one by one each shard to verify it at least is valid python
                    try:
                        py_compile.compile(shard_path, doraise=True)
                    except py_compile.PyCompileError as e:
                        print(f"  [py_compile] FAILED for {shard_path}: {e}")

    print(f"Scraped {created_count} new shards.")


def verify_family(ctx: PipelineContext, family_name: str):
    """Verifies all shards in an equation family using cross-verification (golden tuple)."""
    family_dir = os.path.join(ctx.shards_dir, family_name)
    if not os.path.exists(family_dir):
        return {"error": f"Family directory {family_name} not found"}

    shard_files = [f for f in os.listdir(family_dir) if f.endswith(".py")]

    # We need to assemble a MockLib from all shards in the family
    class MockLib:
        pass

    pyeqn = None
    evaled_methods = {}
    errors = {}

    for sf in shard_files:
        shard_path = os.path.join(family_dir, sf)
        # pycompile check
        try:
            py_compile.compile(shard_path, doraise=True)
        except py_compile.PyCompileError as e:
            errors[sf] = f"Syntax error: {e}"
            continue

        module_name = f"shard_{family_name}_{sf[:-3]}"
        try:
            spec = importlib.util.spec_from_file_location(module_name, shard_path)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)

            class_name = family_name.split("_")[0]
            cls = getattr(module, class_name, None)
            if cls:
                for name, attr in inspect.getmembers(cls):
                    if name.startswith("eqn_"):
                        setattr(MockLib, name, attr)
                        evaled_methods[name] = attr
                        if not pyeqn:
                            try:
                                source = inspect.getsource(attr)
                                match = re.search(r"# \[\.pyeqn\] (.*)", source)
                                if match:
                                    pyeqn = match.group(1).strip()
                            except:
                                pass
            else:
                errors[sf] = f"Class {class_name} not found"
        except Exception as e:
            errors[sf] = str(e)

    if not evaled_methods:
        return {"error": "No methods found to verify", "shard_errors": errors}

    try:
        subshards_dir = os.path.join(ctx.shards_dir, "subshards")
        if not os.path.exists(subshards_dir):
            os.makedirs(subshards_dir)

        v = Verify(MockLib, pyeqn=pyeqn, subshards_dir=subshards_dir)
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
            "shard_errors": errors,
        }
    except Exception as err:
        return {"error": f"Verification failed: {err}"}


def verify_all_shards(ctx: PipelineContext):
    print("Verifying families...")
    all_results = {}
    family_names = sorted(
        [
            d
            for d in os.listdir(ctx.shards_dir)
            if os.path.isdir(os.path.join(ctx.shards_dir, d)) and d != "subshards"
        ]
    )

    iterator = (
        tqdm(family_names, desc="Verifying", unit="family") if tqdm else family_names
    )
    for family_name in iterator:
        all_results[family_name] = verify_family(ctx, family_name)

    with open(os.path.join(ctx.reports_dir, "verification_results.json"), "w") as rf:
        json.dump(all_results, rf, indent=4)
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
                if score == best_score and score >= 1
            ]
            broken = [
                v
                for v, score in variants_inner.items()
                if score < best_score or score == 0
            ]

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
        json.dump(analysis, af, indent=4)
    return analysis


def attempt_repair_shard(
    ctx: PipelineContext,
    family_name: str,
    broken_variant: str,
    trusted_variants: list,
    scores: dict,
    mismatches: dict,
):
    """Repairs a single shard file."""
    family_dir = os.path.join(ctx.shards_dir, family_name)
    shard_file = None
    for f in os.listdir(family_dir):
        if f.endswith(f"__{broken_variant}.py"):
            shard_file = f
            break

    if not shard_file:
        print(f"  [Repair] Shard for {broken_variant} not found in {family_name}")
        return {"updated": False}

    shard_path = os.path.join(family_dir, shard_file)
    with open(shard_path, "r") as f:
        shard_code = f.read()

    pyeqn_match = re.search(r"# \[\.pyeqn\] (.*)", shard_code)
    pyeqn = pyeqn_match.group(1).strip() if pyeqn_match else ""

    print(f"\n |- Repairing shard: {shard_file} in family {family_name}")
    raw_iterator = repair_codigo(
        shard_file=shard_file,
        shard_code=shard_code,
        error=None,
        broken_variants=[broken_variant],
        trusted_variants=trusted_variants,
        scores=scores,
        mismatches=mismatches,
        pyeqn=pyeqn,
        stream=True,
    )

    if raw_iterator is None:
        return {"updated": False}

    raw = ""
    for chunk in raw_iterator:
        msg = (
            chunk.get("message")
            if isinstance(chunk, dict)
            else getattr(chunk, "message", None)
        )
        content = (
            msg.get("content")
            if isinstance(msg, dict)
            else getattr(msg, "content", "")
        )
        raw += str(content or "")
        print(content, end="", flush=True)

    code_text = extract_code(raw)
    if not code_text.strip():
        return {"updated": False}

    try:
        ast.parse(code_text)
        with open(shard_path, "w") as f:
            f.write(code_text)
        return {"updated": True}
    except Exception as e:
        print(f"  [Repair] LLM produced invalid syntax: {e}")
        return {"updated": False, "error": str(e)}


def run_pipeline(
    project_dir=".", max_rounds=DEFAULT_MAX_ROUNDS, auto_repair=True, overwrite=False
):
    ctx = PipelineContext(project_dir)
    if overwrite:
        ctx.clear_shards()
    elif not os.path.exists(ctx.shards_dir) or not os.listdir(ctx.shards_dir):
        ctx.clear_shards()

    shard_from_chapters(ctx)

    for round_idx in range(1, max_rounds + 1):
        print(f"\n--- Round {round_idx}/{max_rounds} ---")
        results = verify_all_shards(ctx)
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

        if analysis["inconsistent"]:
            family_name = sorted(analysis["inconsistent"].keys())[0]
            info = analysis["inconsistent"][family_name]
            broken_variant = info["broken"][0]
            attempt_repair_shard(
                ctx,
                family_name,
                broken_variant,
                info["trusted"],
                info["scores"],
                info["mismatches"],
            )
        elif analysis["failed"]:
            family_name = analysis["failed"][0]["file"]
            print(f"Family {family_name} failed. Re-scraping...")
            shutil.rmtree(os.path.join(ctx.shards_dir, family_name), ignore_errors=True)
            shard_from_chapters(ctx)

    return analysis


def assemble_certified_library(ctx: PipelineContext, analysis):
    """Combines solved families into a single library file."""
    with open(ctx.certified_file, "w") as out:
        out.write(
            "from math import log, sqrt, exp, pow, e\n"
            "from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp\n"
            "from scipy.optimize import newton\n"
            "from kwasak import kwasak_static\n"
            "import numpy as np\n\n"
        )
        class_groups = {}

        for family_name in sorted(analysis["solved"]):
            class_name = family_name.split("_")[0]
            family_dir = os.path.join(ctx.shards_dir, family_name)

            if class_name not in class_groups:
                class_groups[class_name] = []

            for sf in sorted(os.listdir(family_dir)):
                if not sf.endswith(".py"):
                    continue
                with open(os.path.join(family_dir, sf), "r") as f:
                    lines = f.readlines()
                    in_class = False
                    for line in lines:
                        if line.startswith(f"class {class_name}:"):
                            in_class = True
                            continue
                        if in_class:
                            # Skip decorators if they are already there or handled by assemble
                            # Actually just take everything inside the class
                            class_groups[class_name].append(line)

        for class_name, bodies in class_groups.items():
            out.write(f"class {class_name}:\n")
            # We need to make sure we don't duplicate methods if multiple shards have them
            # But the new structure has one method per shard, so it should be fine.
            # We might need to handle indentation if it's inconsistent.
            out.writelines(bodies)
            out.write("\n")


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
        if in_code:
            code_lines.append(line)
    return "\n".join(code_lines).strip()
