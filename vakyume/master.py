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
                f"{TAB * 3}return eval(\"{normal_form.replace(token, 'x')}\".replace('x', str({token})))",
                f'{TAB * 2}raise UnsolvedException("Pending LLM/Manual Repair")',
            ]
        )
        return "\n".join(code)


def shard_from_chapters(ctx: PipelineContext):
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
                        except Exception as e:
                            sf.write(f"{TAB * 2}# Error during Sympy solve: {e}\n")
                            sf.write(
                                solver.sympy_failover(
                                    header, nf, token, pyeqn=line.strip()
                                )
                                + "\n\n"
                            )

    if created_count > 0:
        print(f"Scraped {created_count} new equations.")


def verify_single_shard(ctx: PipelineContext, shard_file: str):
    shard_path = os.path.join(ctx.shards_dir, shard_file)
    try:
        with open(shard_path, "r") as f:
            code = f.read()
    except Exception as e:
        return {"error": str(e)}

    # Extract methods
    blocks = extract_method_blocks(code)
    if not blocks:
        return {"error": "No methods found in shard"}

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

    for name, block_lines in blocks.items():
        # Prepare for exec
        block_text = "\n".join(block_lines)
        lines = block_text.splitlines()
        non_empty = [l for l in lines if l.strip()]
        if not non_empty:
            continue
        min_indent = min(len(l) - len(l.lstrip()) for l in non_empty)
        deindented = "\n".join(l[min_indent:] for l in lines)
        # Remove decorators for direct exec
        deindented = re.sub(r"@staticmethod|@kwasak_static", "", deindented)

        full_code = import_header + deindented
        try:
            # We use a shared globals for all methods in the same shard to mimic class scope
            # but keep it isolated from the pipeline's globals
            global_env = {"__builtins__": __builtins__}
            exec(import_header, global_env)  # Pre-load imports into env

            locs = {}
            exec(deindented, global_env, locs)
            func = locs.get(name)
            if func:
                evaled_methods[name] = func
            else:
                errors[name] = f"Method '{name}' not found in locs after exec."
        except Exception as e:
            errors[name] = f"{type(e).__name__}: {e}"

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
        v = Verify(MockLib)
        verification_results = v.verify()

        # Merge errors and ensure all variants from blocks are represented
        for name in blocks.keys():
            if "__" in name:
                base_eq, token = name.split("__")
                if base_eq not in verification_results:
                    verification_results[base_eq] = {}
                if name in errors:
                    verification_results[base_eq][token] = -1  # Mark as error
                elif token not in verification_results[base_eq]:
                    # This could happen if it loaded but had no score?
                    # Verify usually gives a score to everything it sees in MockLib
                    pass
        return verification_results
    except Exception as e:
        return {"error": f"Verification failed: {e}"}


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

        all_solved = True
        for eqn_name, variants in eqn_results.items():
            num_variants = len(variants)
            if num_variants == 0:
                continue

            max_score = max(variants.values()) if variants else 0
            if max_score != num_variants:
                all_solved = False
                trusted = [v for v, s in variants.items() if s == max_score and s > 1]
                broken = [v for v, s in variants.items() if s < max_score or s <= 1]

                # If everything is broken, pick the first one as broken to start somewhere
                if not trusted and not broken and variants:
                    broken = list(variants.keys())

                analysis["inconsistent"][shard_file] = {
                    "eqn": eqn_name,
                    "num_variants": num_variants,
                    "trusted_variants": trusted,
                    "broken_variants": broken,
                    "scores": variants,
                }
                # For now, one inconsistent equation in a shard is enough to flag the shard
                break

        if all_solved and eqn_results:
            analysis["solved"].append(shard_file)

    with open(os.path.join(ctx.reports_dir, "analysis.json"), "w") as af:
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
):
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

    # Focus on ONE broken variant if multiple exist
    target_broken = broken_variants
    if broken_variants and len(broken_variants) > 1:
        # Prioritize error variants (-1)
        if scores:
            errors = [v for v in broken_variants if scores.get(v) == -1]
            if errors:
                target_broken = [errors[0]]
            else:
                target_broken = [broken_variants[0]]
        else:
            target_broken = [broken_variants[0]]

    print(f"\n |- Repairing variant(s): {target_broken}")
    raw_iterator = repair_codigo(
        shard_file=shard_file,
        shard_code=shard_code,
        error=error,
        broken_variants=target_broken,
        trusted_variants=trusted_variants,
        scores=scores,
        pyeqn=pyeqn,
        stream=True,
    )

    if raw_iterator is None:
        print("  [Error] LLM returned no response.")
        return False

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
    except Exception as e:
        print(f"\n  [Error] Error during streaming: {e}")
        return False

    print("\n--- End of Code ---")

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
        try:
            ast.parse("".join(new_file_lines))
            with open(shard_path, "w") as sf:
                sf.writelines(new_file_lines)
            return True
        except Exception as e:
            print(f" |- [Warning] LLM produced invalid syntax for {shard_file}: {e}")
            return False
    return False


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
                    "broken": info.get("broken_variants", []),
                    "trusted": info.get("trusted_variants", []),
                    "scores": info.get("scores", {}),
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
                new_results = verify_single_shard(ctx, shard_file)
                if new_results and "error" not in new_results:
                    # Check if it's miraculously solved (unlikely but possible if sympy worked this time)
                    is_solved = True
                    for eqn_name, variants in new_results.items():
                        if isinstance(variants, dict):
                            num_v = len(variants)
                            if not all(s == num_v for s in variants.values()):
                                is_solved = False
                                # Update item for next attempt
                                max_score = max(variants.values())
                                item["trusted"] = [
                                    v
                                    for v, s in variants.items()
                                    if s == max_score and s > 1
                                ]
                                item["broken"] = [
                                    v
                                    for v, s in variants.items()
                                    if s < max_score or s <= 1
                                ]
                                item["scores"] = variants
                                item["error"] = None
                                break
                    if is_solved:
                        print(f" |- Result: SUCCESS after re-scrape for {shard_file}")
                        break
                attempt = 1  # Reset attempt counter after re-scrape
                continue

            with open(shard_path, "r") as sf:
                shard_code = sf.read().rstrip()

            updated = attempt_repair_shard(
                ctx=ctx,
                shard_file=shard_file,
                shard_code=shard_code,
                error=item["error"],
                broken_variants=item["broken"],
                trusted_variants=item["trusted"],
                scores=item["scores"],
            )

            is_solved = False
            new_results = None

            if updated:
                new_results = verify_single_shard(ctx, shard_file)
                if new_results and "error" not in new_results:
                    is_solved = True
                    for eqn_name, variants in new_results.items():
                        if isinstance(variants, dict):
                            num_v = len(variants)
                            if not all(s == num_v for s in variants.values()):
                                is_solved = False
                                # Update item for next attempt
                                max_score = max(variants.values())
                                item["trusted"] = [
                                    v
                                    for v, s in variants.items()
                                    if s == max_score and s > 1
                                ]
                                item["broken"] = [
                                    v
                                    for v, s in variants.items()
                                    if s < max_score or s <= 1
                                ]
                                item["scores"] = variants
                                item["error"] = None
                                break
                        else:
                            is_solved = False
                            break

                if is_solved:
                    print(f"\n |- Result: SUCCESS for {shard_file}")
                    break
                else:
                    # Log the specific failure reason
                    if new_results and "error" in new_results:
                        print(
                            f"\n |- Result: FAILED (Execution Error: {new_results['error']}). Retrying..."
                        )
                    else:
                        print(f"\n |- Result: FAILED (Still Inconsistent). Retrying...")
            else:
                print(f"\n |- Result: NO CHANGE (LLM produced no update). Retrying...")

            attempt += 1
            if COOLDOWN_SECONDS > 0:
                time.sleep(COOLDOWN_SECONDS)

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
