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

try:
    from tqdm import tqdm
except ImportError:
    tqdm = None
from sympy import Symbol, solve, Eq
from .verifier import Verify
from .config import *
from .llm import escribir_codigo, repair_codigo, extract_code

# Paths
DEFAULT_MAX_ROUNDS = 10


class PipelineContext:
    def __init__(self, project_dir):
        self.project_dir = os.path.abspath(project_dir)
        self.notes_dir = os.path.join(self.project_dir, "notes")
        self.shards_dir = os.path.join(self.project_dir, "shards")
        self.subshards_dir = os.path.join(self.shards_dir, "subshards")
        self.reports_dir = os.path.join(self.project_dir, "reports")
        self.repair_prompts_dir = os.path.join(self.project_dir, "repair_prompts")
        self.certified_file = os.path.join(self.project_dir, "vakyume_certified.py")

        # Ensure directories exist
        for d in [
            self.shards_dir,
            self.subshards_dir,
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
        if not os.path.exists(self.subshards_dir):
            os.makedirs(self.subshards_dir)
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

                shard_name = f"{class_name}_eqn_{eqn_number}.py"
                shard_path = os.path.join(ctx.shards_dir, shard_name)
                if os.path.exists(shard_path):
                    continue

                if created_count == 0:
                    print("Scraping new equations from notes...")

                created_count += 1
                prompt_cache_path = os.path.join(
                    ctx.repair_prompts_dir, f"{shard_name}.pyeqn"
                )
                with open(prompt_cache_path, "w") as pf:
                    pf.write(line.strip())

                with open(shard_path, "w") as sf:
                    sf.write("from math import log, sqrt, exp, pow, e\n")
                    sf.write(
                        "from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp\n"
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
                                sf.write(f"{TAB * 2}return result\n\n")
                            else:
                                sf.write(
                                    solver.sympy_failover(
                                        header, nf, token, pyeqn=line.strip()
                                    )
                                    + "\n\n"
                                )
                        except timeout_decorator.timeout_decorator.TimeoutError:
                            sf.write(f"{TAB * 2}# Timeout during Sympy solve\n")
                            sf.write(
                                solver.sympy_failover(
                                    header, nf, token, pyeqn=line.strip()
                                )
                                + "\n\n"
                            )
                        except Exception:
                            sf.write(
                                solver.sympy_failover(
                                    header, nf, token, pyeqn=line.strip()
                                )
                                + "\n\n"
                            )
    print(f"Scraped {created_count} new equations.")


def verify_single_shard(ctx: PipelineContext, shard_file: str):
    print(f"[INPUT] verify_single_shard: shard_file={shard_file}")
    shard_path = os.path.join(ctx.shards_dir, shard_file)
    try:
        with open(shard_path, "r") as f:
            code = f.read()
    except Exception as err:
        return {"error": str(err)}

    # Extract methods
    blocks = extract_method_blocks(code)
    if not blocks:
        return {"error": "No methods found in shard"}

    # Extract pyeqn from comments
    pyeqn = None
    pyeqn_match = re.search(r"# \[\.pyeqn\] (.*)", code)
    if pyeqn_match:
        pyeqn = pyeqn_match.group(1).strip()

    # Inject a robust environment for verification
    import_header = (
        "import math\n"
        "from math import *\n"
        "import numpy as np\n"
        "import sympy\n"
        "from sympy import I, Piecewise, LambertW, Eq, symbols, solve\n"
        "from scipy.optimize import newton\n"
        "class UnsolvedException(Exception): pass\n"
    )

    # Also extract any specific imports from the shard
    shard_imports = []
    for line in code.splitlines():
        if line.strip().startswith(("import ", "from ")) and "kwasak" not in line:
            shard_imports.append(line)

    if shard_imports:
        import_header += "\n".join(shard_imports) + "\n"

    evaled_methods = {}
    errors = {}

    # Prepare shared globals for all methods in the same shard
    global_env = {"__builtins__": __builtins__}
    try:
        exec(import_header, global_env)
    except Exception as err:
        return {"error": f"Import header execution failed: {err}"}

    # Use subshards instead of exec
    if ctx.subshards_dir not in sys.path:
        sys.path.append(ctx.subshards_dir)

    for name, block_lines in blocks.items():
        # Prepare subshard
        block_text = "\n".join(block_lines)
        lines = block_text.splitlines()
        non_empty = [l for l in lines if l.strip()]
        if not non_empty:
            continue
        min_indent = min(len(l) - len(l.lstrip()) for l in non_empty)
        deindented = "\n".join(l[min_indent:] for l in lines)
        deindented = re.sub(r"@staticmethod|@kwasak_static", "", deindented)

        subshard_module_name = f"subshard_{shard_file[:-3]}_{name}"
        subshard_path = os.path.join(ctx.subshards_dir, f"{subshard_module_name}.py")

        try:
            with open(subshard_path, "w") as f:
                f.write(import_header)
                f.write("\n\n")
                f.write(deindented)

            # Import or reload
            importlib.invalidate_caches()
            if subshard_module_name in sys.modules:
                module = importlib.reload(sys.modules[subshard_module_name])
            else:
                module = importlib.import_module(subshard_module_name)

            func = getattr(module, name, None)
            if func:
                evaled_methods[name] = func
                global_env[name] = func
            else:
                errors[name] = f"Method '{name}' not found in subshard module."
        except Exception as err:
            errors[name] = f"{type(err).__name__}: {err}"

    # Create MockLib
    class MockLib:
        pass

    for name, func in evaled_methods.items():
        setattr(MockLib, name, func)

    # Ensure base equations are present for Verify to find them
    base_eqs = set()
    for name in blocks.keys():
        if "__" in name:
            base_eqs.add(name.split("__")[0])
    for base_eq in base_eqs:
        if not hasattr(MockLib, base_eq):
            setattr(MockLib, base_eq, lambda: None)  # Dummy for Verify to find

    try:
        v = Verify(MockLib, pyeqn=pyeqn, subshards_dir=ctx.subshards_dir)
        raw_results = v.verify()

        verification_results = {}
        mismatches = {}

        # The new Verify.verify() returns { base_eq: { "scores": ..., "mismatches": ... } }
        # We need to flatten/adapt this for the rest of the pipeline
        for base_eq, data in raw_results.items():
            scores = data["scores"]
            mismatches.update(data["mismatches"])

            if base_eq not in verification_results:
                verification_results[base_eq] = {}

            for token, score in scores.items():
                verification_results[base_eq][token] = score

        # Merge errors and ensure all variants from blocks are represented
        for name in blocks.keys():
            if "__" in name:
                base_eq, token = name.split("__")
                if base_eq not in verification_results:
                    verification_results[base_eq] = {}
                if name in errors:
                    verification_results[base_eq][token] = -1  # Mark as error
                elif token not in verification_results[base_eq]:
                    pass

        return {"results": verification_results, "mismatches": mismatches}
    except Exception as err:
        print(f"[OUTPUT] verify_single_shard: error={err}")
        return {"error": f"Verification failed: {err}"}


def verify_all_shards(ctx: PipelineContext):
    print("Verifying shards...")
    all_results = {}
    shard_files = sorted(
        [
            f
            for f in os.listdir(ctx.shards_dir)
            if f.endswith(".py") and f != "kwasak.py"
        ]
    )

    iterator = (
        tqdm(shard_files, desc="Verifying", unit="shard") if tqdm else shard_files
    )
    for f in iterator:
        all_results[f] = verify_single_shard(ctx, f)

    with open(os.path.join(ctx.reports_dir, "verification_results.json"), "w") as rf:
        json.dump(all_results, rf, indent=4)
    return all_results


def analyze_results(ctx: PipelineContext, all_results):
    print(f"[INPUT] analyze_results: {len(all_results)} shards in all_results")
    print("Analyzing alignment for certification...")
    analysis = {"solved": [], "inconsistent": {}, "failed": []}

    for shard_file, shard_res in all_results.items():
        if (
            not shard_res
            or not isinstance(shard_res, dict)
            or "results" not in shard_res
        ):
            analysis["failed"].append(
                {
                    "file": shard_file,
                    "error": shard_res.get("error")
                    if (shard_res and isinstance(shard_res, dict))
                    else "Unknown error",
                }
            )
            continue

        eqn_results_data = shard_res.get("results")
        if not isinstance(eqn_results_data, dict):
            analysis["failed"].append(
                {"file": shard_file, "error": "Invalid results format"}
            )
            continue

        mismatches_data = shard_res.get("mismatches", {})
        all_solved = True

        for base_eq, variants_inner in eqn_results_data.items():
            if not isinstance(variants_inner, dict):
                all_solved = False
                continue

            # Improvement: Use relative majority and harmony check
            num_v = len(variants_inner)
            best_score = max(variants_inner.values()) if variants_inner else 0

            # Trusted are those that passed harmony (score >= 1) and have the best score
            trusted = [
                v
                for v, score in variants_inner.items()
                if score == best_score and score >= 1
            ]
            # Broken are those that failed harmony (score == 0) or are inconsistent with the majority
            broken = [
                v
                for v, score in variants_inner.items()
                if score < best_score or score == 0
            ]

            if broken:
                all_solved = False
                analysis["inconsistent"][shard_file] = {
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
            analysis["solved"].append(shard_file)

    with open(os.path.join(ctx.reports_dir, "analysis.json"), "w") as af:
        json.dump(analysis, af, indent=4)
    print(
        f"[OUTPUT] analyze_results: solved={len(analysis['solved'])}, inconsistent={len(analysis['inconsistent'])}, failed={len(analysis['failed'])}"
    )
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
        # Look for def eqn_...
        if stripped.startswith("def "):
            match = re.search(r"def\s+(eqn_\w+)\s*\(", stripped)
            if match:
                name = match.group(1)
                start = i
                # Include preceding decorators
                while start > 0:
                    prev_line = lines[start - 1].lstrip()
                    if prev_line.startswith("@") or not prev_line:
                        start -= 1
                    else:
                        break

                # Find end of method by indentation
                line_indent = len(line) - len(line.lstrip())
                end = i + 1
                while end < len(lines):
                    if lines[end].strip():
                        curr_indent = len(lines[end]) - len(lines[end].lstrip())
                        # If indent is less than or equal to method indent, and not a continuation line
                        if curr_indent <= line_indent:
                            # Check if it's a new method or class or decorator
                            test_line = lines[end].lstrip()
                            if (
                                test_line.startswith("def ")
                                or test_line.startswith("@")
                                or test_line.startswith("class ")
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
    # More robust matching
    target_regex = re.compile(rf"def\s+{method_name}\s*\(")
    for idx, line in enumerate(file_lines):
        if target_regex.search(line):
            start = idx
            # Go up to find decorators
            while start > 0 and file_lines[start - 1].lstrip().startswith("@"):
                start -= 1

            # Find end of existing method
            indent = len(line) - len(line.lstrip())
            end = idx + 1
            while end < len(file_lines):
                if file_lines[end].strip():
                    curr_indent = len(file_lines[end]) - len(file_lines[end].lstrip())
                    if curr_indent <= indent:
                        if (
                            file_lines[end].lstrip().startswith("def ")
                            or file_lines[end].lstrip().startswith("@")
                            or file_lines[end].lstrip().startswith("class ")
                        ):
                            break
                end += 1

            # Replace
            file_lines[start:end] = new_block_lines + ["\n"]
            return True
    print(f"  [Warning] Could not find method {method_name} in file for replacement.")
    return False


def attempt_repair_shard(
    ctx: PipelineContext,
    shard_file,
    shard_code,
    error=None,
    broken_variants=None,
    trusted_variants=None,
    scores=None,
    mismatches=None,
):
    print(
        f"[INPUT] attempt_repair_shard: shard_file={shard_file}, broken={broken_variants}, trusted={trusted_variants}"
    )
    shard_path = os.path.join(ctx.shards_dir, shard_file)
    prompt_cache_path = os.path.join(ctx.repair_prompts_dir, f"{shard_file}.pyeqn")

    # Try to get pyeqn from cache first, then scrape from shard_code, then cache it
    pyeqn = None
    if os.path.exists(prompt_cache_path):
        with open(prompt_cache_path, "r") as f:
            pyeqn = f.read().strip()

    if not pyeqn:
        pyeqn_match = re.search(r"# \[\.pyeqn\] (.*)", shard_code)
        if pyeqn_match:
            pyeqn = pyeqn_match.group(1).strip()
            # Cache it for next time
            with open(prompt_cache_path, "w") as f:
                f.write(pyeqn)

    # Focus on broken variants
    target_broken = broken_variants
    if not broken_variants:
        target_broken = []

    print(f"\n |- Repairing variant(s): {target_broken}")
    raw_iterator = repair_codigo(
        shard_file=shard_file,
        shard_code=shard_code,
        error=error,
        broken_variants=target_broken,
        trusted_variants=trusted_variants,
        scores=scores,
        mismatches=mismatches,
        pyeqn=pyeqn,
        stream=True,
    )

    if raw_iterator is None:
        print("  [Error] LLM returned no response.")
        return {"updated": False}

    print(f"\n--- Code for {shard_file} ---")
    raw = ""
    try:
        for chunk in raw_iterator:
            # Use a safe way to get content that handles both dicts and objects
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
            content = str(content or "")
            raw += content
            print(content, end="", flush=True)
    except Exception as err:
        print(f"\n  [Error] Error during streaming: {err}")
        return {"updated": False}

    print("\n--- End of Code ---")

    code_text = extract_code(raw)
    if not code_text.strip():
        code_text = extract_code_block(raw)
    if not code_text.strip():
        return {"updated": False}

    blocks = extract_method_blocks(code_text)
    if not blocks:
        return {"updated": False}

    with open(shard_path, "r") as sf:
        file_lines = sf.readlines()

    # Work on a copy to validate
    new_file_lines = list(file_lines)
    updated = False
    candidates = sorted(blocks.keys())

    for candidate in candidates:
        # Don't let it overwrite the main kwasak method unless explicitly broken?
        # Actually, let's just stick to broken_variants if they exist.
        if broken_variants:
            match = False
            for v in broken_variants:
                if candidate.endswith(f"__{v}"):
                    match = True
                    break
            if not match:
                continue

        new_block = normalize_block(blocks[candidate])
        if new_block and replace_method_in_file(new_file_lines, candidate, new_block):
            updated = True

    if updated:
        # Validate syntax
        file_content = "".join(new_file_lines)
        try:
            ast.parse(file_content)
        except Exception as err:
            print(f" |- [Warning] LLM produced invalid syntax for {shard_file}: {err}")
            # We still save it so the next iteration can show the LLM its mistake

        with open(shard_path, "w") as sf:
            sf.write(file_content)
        print(f"[OUTPUT] attempt_repair_shard: updated=True for {shard_file}")
        return {"updated": True}
    print(f"[OUTPUT] attempt_repair_shard: updated=False for {shard_file}")
    return {"updated": False}


def run_pipeline(
    project_dir=".", max_rounds=DEFAULT_MAX_ROUNDS, auto_repair=True, overwrite=False
):
    ctx = PipelineContext(project_dir)
    if overwrite:
        ctx.clear_shards()
    elif not os.path.exists(ctx.shards_dir):
        ctx.clear_shards()
    shard_from_chapters(ctx)

    analysis = {"solved": [], "inconsistent": {}, "failed": []}
    for round_idx in range(1, max_rounds + 1):
        print(f"\n--- Round {round_idx}/{max_rounds} ---")
        results = verify_all_shards(ctx)
        analysis = analyze_results(ctx, results)

        print(
            f"Status: {len(analysis['solved'])} solved, "
            f"{len(analysis['inconsistent'])} inconsistent, "
            f"{len(analysis['failed'])} failed."
        )

        assemble_certified_library(ctx, analysis)

        if not analysis.get("inconsistent") and not analysis.get("failed"):
            print("All equations aligned and no failures. Stopping early.")
            return analysis

        to_repair = []
        # Sort keys to ensure consistent order
        inconsistent_files = sorted(analysis.get("inconsistent", {}).keys())
        for sf in inconsistent_files:
            info = analysis["inconsistent"][sf]
            to_repair.append(
                {
                    "file": sf,
                    "broken": info.get("broken", []),
                    "trusted": info.get("trusted", []),
                    "scores": info.get("scores", {}),
                    "mismatches": info.get("mismatches", {}),
                    "error": None,
                }
            )

        failed_shards = sorted(analysis.get("failed", []), key=lambda x: x["file"])
        for fail_info in failed_shards:
            to_repair.append(
                {
                    "file": fail_info["file"],
                    "broken": [],
                    "trusted": [],
                    "scores": {},
                    "error": fail_info.get("error"),
                }
            )

        if not to_repair:
            print("Nothing to repair. Stopping.")
            break

        # Focus on ONE shard at a time
        item = to_repair[0]
        shard_file = item["file"]
        shard_path = os.path.join(ctx.shards_dir, shard_file)
        if not os.path.exists(shard_path):
            print(f"Shard path {shard_path} not found. Skipping.")
            continue

        attempt = 1
        while True:
            sys.stdout.write(f"\r[Repair] Shard: {shard_file} | Attempt: {attempt}")
            sys.stdout.flush()

            if attempt > 5:
                print(
                    f"\n |- Attempt limit reached for {shard_file}. Re-scraping from notes..."
                )
                if os.path.exists(shard_path):
                    os.remove(shard_path)
                shard_from_chapters(ctx)
                # After re-scrape, re-verify this specific shard to get fresh analysis
                res_wrap = verify_single_shard(ctx, shard_file)
                if res_wrap and isinstance(res_wrap, dict) and "results" in res_wrap:
                    is_solved = True
                    results_data = res_wrap.get("results", {})
                    mismatches_data = res_wrap.get("mismatches", {})
                    if not isinstance(results_data, dict):
                        is_solved = False
                    else:
                        for eqn_name, variants_data in results_data.items():
                            if isinstance(variants_data, dict):
                                num_v = len(variants_data)
                                if not all(s == num_v for s in variants_data.values()):
                                    is_solved = False
                                    max_score = max(variants_data.values())
                                    item["trusted"] = [
                                        v
                                        for v, s in variants_data.items()
                                        if s == max_score and s > 1
                                    ]
                                    item["broken"] = [
                                        v
                                        for v, s in variants_data.items()
                                        if s < max_score or s <= 1
                                    ]
                                    item["scores"] = variants_data
                                    item["mismatches"] = {
                                        v: mismatches_data.get(v)
                                        for v in item["broken"]
                                        if isinstance(mismatches_data, dict)
                                        and v in mismatches_data
                                    }
                                    item["error"] = None
                                    break
                    if is_solved:
                        print(f" |- Result: SUCCESS after re-scrape for {shard_file}")
                        break

                    if is_solved:
                        print(f" |- Result: SUCCESS after re-scrape for {shard_file}")
                        break

                    if is_solved:
                        print(f" |- Result: SUCCESS after re-scrape for {shard_file}")
                        break

                    if is_solved:
                        print(f" |- Result: SUCCESS after re-scrape for {shard_file}")
                        break
                attempt = 1  # Reset attempt counter after re-scrape
                continue

            with open(shard_path, "r") as sf:
                shard_code = sf.read().rstrip()

            repair_res = attempt_repair_shard(
                ctx,
                shard_file=shard_file,
                shard_code=shard_code,
                broken_variants=item.get("broken"),
                trusted_variants=item.get("trusted"),
                scores=item.get("scores"),
                mismatches=item.get("mismatches"),
            )

            updated = repair_res.get("updated", False)
            if not updated and "syntax_error" in repair_res:
                item["error"] = repair_res["syntax_error"]

            is_solved = False
            new_results = None

            if updated:
                res_wrap = verify_single_shard(ctx, shard_file)
                if res_wrap and isinstance(res_wrap, dict) and "results" in res_wrap:
                    is_solved = True
                    results_data = res_wrap.get("results", {})
                    mismatches_data = res_wrap.get("mismatches", {})
                    if not isinstance(results_data, dict):
                        is_solved = False
                    else:
                        for eqn_name, variants_data in results_data.items():
                            if isinstance(variants_data, dict):
                                num_v = len(variants_data)
                                if not all(s == num_v for s in variants_data.values()):
                                    is_solved = False
                                    max_score = max(variants_data.values())
                                    item["trusted"] = [
                                        v
                                        for v, s in variants_data.items()
                                        if s == max_score and s > 1
                                    ]
                                    item["broken"] = [
                                        v
                                        for v, s in variants_data.items()
                                        if s < max_score or s <= 1
                                    ]
                                    item["scores"] = variants_data
                                    item["mismatches"] = {
                                        v: mismatches_data.get(v)
                                        for v in item["broken"]
                                        if isinstance(mismatches_data, dict)
                                        and v in mismatches_data
                                    }
                                    item["error"] = None
                                    break
                            else:
                                is_solved = False
                                break
                    if is_solved:
                        print(f"\n |- Result: SUCCESS for {shard_file}")
                        break
                    else:
                        print(f"\n |- Result: FAILED (Still Inconsistent). Retrying...")
                else:
                    err = (
                        res_wrap.get("error")
                        if (res_wrap and isinstance(res_wrap, dict))
                        else "Unknown error"
                    )
                    item["error"] = err
                    print(f"\n |- Result: FAILED (Execution Error: {err}). Retrying...")

            else:
                print(f"\n |- Result: NO CHANGE (LLM produced no update). Retrying...")

            attempt += 1
            if COOLDOWN_SECONDS > 0:
                time.sleep(COOLDOWN_SECONDS)

        # After one shard is fixed, we go back to the top of the round loop (or start next round)
        # Actually, if we fixed it, we should probably re-verify all to get updated analysis
        continue

    return analysis


def assemble_certified_library(ctx: PipelineContext, analysis):
    with open(ctx.certified_file, "w") as out:
        out.write(
            "from math import log, sqrt, exp, pow, e\nfrom sympy import I, Piecewise, LambertW, Eq, symbols, solve\nfrom scipy.optimize import newton\nfrom kwasak import kwasak_static\nimport numpy as np\n\n"
        )
        class_groups = {}
        for shard_file in sorted(analysis["solved"]):
            class_name = shard_file.split("_")[0]
            with open(os.path.join(ctx.shards_dir, shard_file), "r") as sf:
                lines = sf.readlines()
                in_class = False
                for line in lines:
                    if line.startswith(f"class {class_name}:"):
                        in_class = True
                        continue
                    if in_class:
                        if class_name not in class_groups:
                            class_groups[class_name] = []
                        class_groups[class_name].append(line)
        for class_name, bodies in class_groups.items():
            out.write(f"class {class_name}:\n")
            out.writelines(bodies)
            out.write("\n")
