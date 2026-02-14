import os
import shutil
import inspect
import importlib.util
import json
import re
import timeout_decorator
import argparse
import time

try:
    from tqdm import tqdm
except ImportError:
    tqdm = None
from sympy import Symbol, solve, Eq
from tru import Verify

# Import from llm if available, else mock
from typing import Callable, Optional

try:
    from llm import escribir_codigo as _escribir_codigo
    from llm import repair_codigo as _repair_codigo
    from llm import extract_code as _extract_code
except Exception:
    _escribir_codigo = None
    _repair_codigo = None
    _extract_code = None


def escribir_codigo(*args, **kwargs):
    if _escribir_codigo is None:
        return "# LLM Mock: No Ollama"
    return _escribir_codigo(*args, **kwargs)


def repair_codigo(*args, **kwargs):
    if _repair_codigo is None:
        return "# LLM Mock: No Ollama"
    return _repair_codigo(*args, **kwargs)


def extract_code(text):
    if _extract_code is None:
        return ""
    return _extract_code(text)


from config import *

# Paths
ROOT = os.getcwd()
SHARDS_DIR = os.path.join(ROOT, "shards")
SOURCE_FILE = os.path.join(ROOT, "vakyume.py")
REPORTS_DIR = os.path.join(ROOT, "reports")
CHAPTERS_DIR = os.path.join(ROOT, "chapters")
REPAIR_PROMPTS_DIR = os.path.join(ROOT, "repair_prompts")
DEFAULT_MAX_ROUNDS = 10


def clear_shards():
    if os.path.exists(SHARDS_DIR):
        shutil.rmtree(SHARDS_DIR)
    os.makedirs(SHARDS_DIR)
    if os.path.exists(os.path.join(ROOT, "kwasak.py")):
        shutil.copy(
            os.path.join(ROOT, "kwasak.py"), os.path.join(SHARDS_DIR, "kwasak.py")
        )
    if not os.path.exists(REPORTS_DIR):
        os.makedirs(REPORTS_DIR)
    if not os.path.exists(REPAIR_PROMPTS_DIR):
        os.makedirs(REPAIR_PROMPTS_DIR)


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

    def sympy_failover(self, eqn_header, normal_form, token):
        # In a real scenario, this would call LLM.
        # For now, we leave a placeholder or a numerical hint.
        code = [
            f"{TAB * 2}# [Sympy Failover Placeholder for {token}]",
            f"{TAB * 2}def func({token}):",
            f"{TAB * 3}# Numerical fallback needed for: {normal_form}",
            f"{TAB * 3}return eval(\"{normal_form.replace(token, 'x')}\".replace('x', str({token})))",
            f"{TAB * 2}# result = [newton(func, 1.0)]",
            f'{TAB * 2}raise UnsolvedException("Pending LLM/Manual Repair")',
        ]
        return "\n".join(code)


def shard_from_chapters():
    solver = Solver()

    if not os.path.exists(CHAPTERS_DIR):
        print(f"Chapters directory not found: {CHAPTERS_DIR}")
        return

    chapter_files = sorted(os.listdir(CHAPTERS_DIR))
    created_count = 0

    # Pre-check if anything needs to be done to avoid noisy tqdm
    to_process = []
    for chapter_file in chapter_files:
        if chapter_file.endswith(".py") and not chapter_file.startswith("__"):
            to_process.append(chapter_file)

    if not to_process:
        return

    for chapter_file in to_process:
        chapter_path = os.path.join(CHAPTERS_DIR, chapter_file)
        # Extract class name from filename: 01_vacuum_theory.py -> VacuumTheory
        base_name = os.path.splitext(chapter_file)[0]
        parts = base_name.split("_")[1:]
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

                shard_name = f"{class_name}_eqn_{eqn_number}.py"
                shard_path = os.path.join(SHARDS_DIR, shard_name)
                if os.path.exists(shard_path):
                    continue

                if created_count == 0:
                    print("Scraping new equations from chapters...")

                created_count += 1
                with open(shard_path, "w") as sf:
                    sf.write("from math import log, sqrt, exp, pow, e\n")
                    sf.write(
                        "from sympy import I, Piecewise, LambertW, Eq, symbols, solve\n"
                    )
                    sf.write("from scipy.optimize import newton\n")
                    sf.write("from kwasak import kwasak_static\n")
                    sf.write("import numpy as np\n\n")
                    sf.write(f"class {class_name}:\n")
                    sf.write(f"{TAB}@kwasak_static\n")
                    sf.write(
                        f"{TAB}def eqn_{eqn_number}({', '.join(f'{t}=None' for t in tokes)}, **kwargs):\n"
                    )
                    sf.write(f"{TAB * 2}return\n\n")

                    for token in tokes:
                        other_args = [t for t in tokes if t != token]
                        header = f"{TAB}@staticmethod\n{TAB}def eqn_{eqn_number}__{token}({', '.join(f'{t}: float' for t in other_args)}):"
                        sf.write(header + "\n")
                        sf.write(f"{TAB * 2}# [.pyeqn] {line.strip()}\n")

                        try:
                            solns = solver.get_solns_vanilla_nf(nf, Symbol(token))
                            if solns:
                                sf.write(f"{TAB * 2}result = []\n")
                                for soln in solns:
                                    sf.write(f"{TAB * 2}{token} = {soln}\n")
                                    sf.write(f"{TAB * 2}result.append({token})\n")
                                sf.write(f"{TAB * 2}return result\n\n")
                            else:
                                sf.write(
                                    solver.sympy_failover(header, nf, token) + "\n\n"
                                )
                        except timeout_decorator.timeout_decorator.TimeoutError:
                            sf.write(f"{TAB * 2}# Timeout during Sympy solve\n")
                            sf.write(solver.sympy_failover(header, nf, token) + "\n\n")
                        except Exception as e:
                            sf.write(f"{TAB * 2}# Error during Sympy solve: {e}\n")
                            sf.write(solver.sympy_failover(header, nf, token) + "\n\n")

    if created_count > 0:
        print(f"Scraped {created_count} new equations.")


def verify_all_shards():
    print("Verifying shards...")
    all_results = {}
    shard_files = sorted(
        [f for f in os.listdir(SHARDS_DIR) if f.endswith(".py") and f != "kwasak.py"]
    )

    iterator = (
        tqdm(shard_files, desc="Verifying", unit="shard") if tqdm else shard_files
    )
    for f in iterator:
        shard_path = os.path.join(SHARDS_DIR, f)
        module_name = f[:-3]

        try:
            spec = importlib.util.spec_from_file_location(module_name, shard_path)
            if spec and spec.loader:
                module = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(module)

                for name, obj in inspect.getmembers(module):
                    if (
                        inspect.isclass(obj)
                        and name != "Verify"
                        and obj.__module__ == module_name
                    ):
                        v = Verify(obj)
                        shard_results = v.verify()
                        all_results[f] = shard_results
        except Exception as e:
            # Silence noisy verification errors, they are captured in all_results
            all_results[f] = {"error": str(e)}

    with open(os.path.join(REPORTS_DIR, "verification_results.json"), "w") as rf:
        json.dump(all_results, rf, indent=4)
    return all_results


def analyze_results(all_results):
    print("Analyzing alignment for certification...")
    analysis = {"solved": [], "inconsistent": {}, "failed": []}

    for shard_file, eqn_results in all_results.items():
        if not eqn_results or "error" in eqn_results:
            analysis["failed"].append(
                {
                    "file": shard_file,
                    "error": eqn_results.get("error")
                    if eqn_results
                    else "Unknown error",
                }
            )
            continue

        for eqn_name, variants in eqn_results.items():
            num_variants = len(variants)
            valid_scores = [s for s in variants.values() if s is not None]

            if valid_scores and all(s == num_variants for s in valid_scores):
                analysis["solved"].append(shard_file)
            else:
                inconsistent = [
                    v for v, s in variants.items() if s is None or s < num_variants
                ]
                analysis["inconsistent"][shard_file] = {
                    "eqn": eqn_name,
                    "num_variants": num_variants,
                    "inconsistent_variants": inconsistent,
                    "scores": variants,
                }

    with open(os.path.join(REPORTS_DIR, "analysis.json"), "w") as af:
        json.dump(analysis, af, indent=4)
    return analysis


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


def extract_method_blocks(code_text):
    lines = code_text.splitlines()
    blocks = {}
    i = 0
    while i < len(lines):
        line = lines[i]
        stripped = line.lstrip()
        if stripped.startswith("def eqn_"):
            name = stripped.split("def ", 1)[1].split("(", 1)[0].strip()
            start = i
            if i > 0 and lines[i - 1].lstrip().startswith("@"):
                start = i - 1
            end = i + 1
            while end < len(lines):
                nxt = lines[end]
                nxt_stripped = nxt.lstrip()
                if nxt_stripped.startswith("def eqn_") or nxt_stripped.startswith(
                    "@staticmethod"
                ):
                    break
                end += 1
            block_lines = lines[start:end]
            blocks[name] = block_lines
            i = end
            continue
        i += 1
    return blocks


def normalize_block(block_lines):
    non_empty = [l for l in block_lines if l.strip()]
    if not non_empty:
        return []
    min_indent = min(len(l) - len(l.lstrip()) for l in non_empty)
    normalized = []
    for line in block_lines:
        if line.strip():
            normalized.append(" " * 4 + line[min_indent:].rstrip())
        else:
            normalized.append("")
    return [l + "\n" for l in normalized]


def replace_method_in_file(file_lines, method_name, new_block_lines):
    target = f"def {method_name}("
    for idx, line in enumerate(file_lines):
        if target in line:
            start = idx
            while start > 0 and file_lines[start - 1].lstrip().startswith("@"):
                start -= 1
            end = idx + 1
            while end < len(file_lines):
                nxt = file_lines[end]
                if nxt.startswith("    @") or nxt.startswith("    def "):
                    break
                end += 1
            file_lines[start:end] = new_block_lines + ["\n"]
            return True
    return False


def attempt_repair_shard(
    shard_file, shard_code, error=None, inconsistent_variants=None, scores=None
):
    shard_path = os.path.join(SHARDS_DIR, shard_file)

    # Using the optimized repair prompt
    raw = repair_codigo(
        shard_file=shard_file,
        shard_code=shard_code,
        error=error,
        inconsistent_variants=inconsistent_variants,
        scores=scores,
    )

    code_text = extract_code(raw)
    if not code_text.strip():
        code_text = extract_code_block(raw)

    if not code_text.strip():
        return False

    blocks = extract_method_blocks(code_text)
    if not blocks:
        return False

    with open(shard_path, "r") as sf:
        file_lines = sf.readlines()

    updated = False
    for candidate in blocks.keys():
        # If we have specific variants, check if this candidate matches one
        if inconsistent_variants:
            match = False
            for v in inconsistent_variants:
                if candidate.endswith(f"__{v}"):
                    match = True
                    break
            if not match:
                continue

        # If match or if we are repairing everything
        new_block = normalize_block(blocks[candidate])
        if new_block and replace_method_in_file(file_lines, candidate, new_block):
            updated = True

    if updated:
        with open(shard_path, "w") as sf:
            sf.writelines(file_lines)
    return updated


def run_pipeline(max_rounds=DEFAULT_MAX_ROUNDS, auto_repair=True, overwrite=False):
    analysis = {"solved": [], "inconsistent": {}, "failed": []}

    if not os.path.exists(REPORTS_DIR):
        os.makedirs(REPORTS_DIR)
    if not os.path.exists(REPAIR_PROMPTS_DIR):
        os.makedirs(REPAIR_PROMPTS_DIR)

    if overwrite:
        clear_shards()
    elif not os.path.exists(SHARDS_DIR):
        clear_shards()
    shard_from_chapters()

    for round_idx in range(1, max_rounds + 1):
        if round_idx == 1:
            print("\n--- Initial OOO Status Scan ---")
        else:
            print(f"\n--- Verification Round {round_idx}/{max_rounds} ---")

        results = verify_all_shards()
        analysis = analyze_results(results)

        print(
            f"OOO Status: {len(analysis['solved'])} families fully functional / solved, "
            f"{len(analysis['inconsistent'])} inconsistent, "
            f"{len(analysis['failed'])} failed."
        )

        assemble_certified_library(analysis)

        if not analysis.get("inconsistent") and not analysis.get("failed"):
            print("All equations aligned and no failures. Stopping early.")
            return analysis

        print(f"Starting auto-repair round...")
        any_updates = False

        # Combine inconsistent and failed for repair
        to_repair = []
        for sf, info in analysis.get("inconsistent", {}).items():
            to_repair.append(
                {
                    "file": sf,
                    "variants": info.get("inconsistent_variants", []),
                    "scores": info.get("scores", {}),
                    "error": None,
                }
            )
        for fail_info in analysis.get("failed", []):
            to_repair.append(
                {
                    "file": fail_info["file"],
                    "variants": [],
                    "scores": {},
                    "error": fail_info.get("error"),
                }
            )

        iterator = (
            tqdm(to_repair, desc="Repairing Shards", unit="shard")
            if tqdm
            else to_repair
        )
        for item in iterator:
            shard_file = item["file"]
            shard_path = os.path.join(SHARDS_DIR, shard_file)

            if not os.path.exists(shard_path):
                continue

            # Update progress display
            if not isinstance(iterator, list):
                iterator.set_postfix(shard=shard_file)
            else:
                print(f"Attempting auto-repair for {shard_file}...")

            with open(shard_path, "r") as sf:
                shard_code = sf.read().rstrip()

            updated = attempt_repair_shard(
                shard_file=shard_file,
                shard_code=shard_code,
                error=item["error"],
                inconsistent_variants=item["variants"],
                scores=item["scores"],
            )

            if updated and isinstance(iterator, list):
                print(f"Auto-repair succeeded for {shard_file}.")
            any_updates = any_updates or updated

            # Cool down to prevent hardware overheating during heavy LLM inference
            if COOLDOWN_SECONDS > 0:
                time.sleep(COOLDOWN_SECONDS)

        if not any_updates:
            print("No shard updates produced by auto-repair. Stopping.")
            return analysis

    print("Reached maximum repair rounds.")
    return analysis


def assemble_certified_library(analysis):
    print(
        f"Assembling certified library from {len(analysis['solved'])} aligned equations..."
    )
    certified_file = os.path.join(ROOT, "vakyume_certified.py")

    with open(certified_file, "w") as out:
        out.write("from math import log, sqrt, exp, pow, e\n")
        out.write("from sympy import I, Piecewise, LambertW, Eq, symbols, solve\n")
        out.write("from scipy.optimize import newton\n")
        out.write("from kwasak import kwasak_static\n")
        out.write("import numpy as np\n\n")

        class_groups = {}
        for shard_file in sorted(analysis["solved"]):
            class_name = shard_file.split("_")[0]
            with open(os.path.join(SHARDS_DIR, shard_file), "r") as sf:
                lines = sf.readlines()
                body = []
                in_class = False
                for line in lines:
                    if line.startswith(f"class {class_name}:"):
                        in_class = True
                        continue
                    if in_class:
                        body.append(line)

                if class_name not in class_groups:
                    class_groups[class_name] = []
                class_groups[class_name].extend(body)

        for class_name, bodies in class_groups.items():
            out.write(f"class {class_name}:\n")
            out.writelines(bodies)
            out.write("\n")
    print(f"Certified library created: {certified_file}")


def _parse_bool(value: str) -> bool:
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y"}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--overwrite",
        nargs="?",
        const="true",
        default="false",
        help="Overwrite shards directory (true/false)",
    )
    parser.add_argument(
        "--max-rounds",
        type=int,
        default=DEFAULT_MAX_ROUNDS,
        help="Maximum verification rounds",
    )
    args = parser.parse_args()

    final_analysis = {"solved": [], "inconsistent": {}, "failed": []}
    final_analysis = run_pipeline(
        max_rounds=args.max_rounds,
        auto_repair=True,
        overwrite=_parse_bool(args.overwrite),
    )

    print("\n" + "=" * 40)
    print(f"SOLVED (Aligned):      {len(final_analysis['solved'])}")
    print(f"INCONSISTENT:          {len(final_analysis['inconsistent'])}")
    print(f"UNSOLVABLE/FAILED:     {len(final_analysis['failed'])}")
    print("=" * 40)
