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
        "from math import log, sqrt, exp, pow, e\n"
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
                    shard_name = f"{method_name}.py"
                    shard_path = os.path.join(family_dir, shard_name)
                    if os.path.exists(shard_path):
                        continue

                    created_count += 1
                    shard_content = import_header
                    other_args = [t for t in tokes if t != token]
                    header = f"def {method_name}({', '.join(f'{t}: float' for t in other_args)}, **kwargs):"
                    shard_content += header + "\n"
                    shard_content += f"{TAB}# [.pyeqn] {line.strip()}\n"

                    try:
                        solns = solver.get_solns_vanilla_nf(nf, Symbol(token))
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

                    # One-shot LLM repair if SymPy failed
                    if "UnsolvedException" in shard_content:
                        # Use a stream to collect response but don't log to console
                        raw_iterator = repair_codigo(
                            shard_file=shard_name,
                            shard_code=shard_content,
                            pyeqn=line.strip(),
                            broken_variants=[token],
                            is_subshard=True,
                            stream=True
                        )
                        raw = ""
                        for chunk in raw_iterator:
                            msg = chunk.get("message") if isinstance(chunk, dict) else getattr(chunk, "message", None)
                            content = msg.get("content") if isinstance(msg, dict) else getattr(msg, "content", "")
                            raw += str(content or "")
                        
                        code_text = extract_code(raw, target_name=method_name)
                        if code_text.strip():
                            if "import " not in code_text:
                                shard_content = import_header + code_text
                            else:
                                shard_content = code_text

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
                            sf.write(f"from .{method_name} import {method_name}\n")
                        sf.write(f"\nclass {class_name}:\n")
                        for token in tokes:
                            method_name = f"eqn_{eqn_number}__{token}"
                            sf.write(f"{TAB}{method_name} = staticmethod({method_name})\n")
                        sf.write(f"\n{TAB}@kwasak_static\n")
                        sf.write(
                            f"{TAB}def eqn_{eqn_number}({', '.join(f'{t}=None' for t in tokes)}, **kwargs):\n"
                        )
                        sf.write(f"{TAB * 2}return\n")

    print(f"Scraped {created_count} new shards.")


def verify_family(ctx: PipelineContext, family_name: str):
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

        v = Verify(cls, pyeqn=pyeqn, subshards_dir=subshards_dir)
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
                if score == best_score and score >= num_v
            ]
            broken = [
                v
                for v, score in variants_inner.items()
                if score < num_v
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
    error: str = None
):
    """Repairs a single shard file."""
    family_dir = os.path.join(ctx.shards_dir, family_name)
    shard_file = f"eqn_{family_name.split('_eqn_')[1]}__{broken_variant}.py"
    shard_path = os.path.join(family_dir, shard_file)

    if not os.path.exists(shard_path):
        print(f"  [Repair] Shard {shard_file} not found in {family_name}")
        return {"updated": False}

    with open(shard_path, "r") as f:
        shard_code = f.read()

    pyeqn_match = re.search(r"# \[\.pyeqn\] (.*)", shard_code)
    pyeqn = pyeqn_match.group(1).strip() if pyeqn_match else ""

    print(f"\n |- Repairing shard: {shard_file} in family {family_name}")
    raw_iterator = repair_codigo(
        shard_file=shard_file,
        shard_code=shard_code,
        error=error,
        broken_variants=[broken_variant],
        trusted_variants=trusted_variants,
        scores=scores,
        mismatches=mismatches,
        pyeqn=pyeqn,
        stream=True,
        is_subshard=True,
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

    method_name = f"eqn_{family_name.split('_eqn_')[1]}__{broken_variant}"
    code_text = extract_code(raw, target_name=method_name)
    if not code_text.strip():
        return {"updated": False}

    try:
        ast.parse(code_text)
        # Ensure it has the correct header
        import_header = (
            "from math import log, sqrt, exp, pow, e\n"
            "from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp\n"
            "from scipy.optimize import newton\n"
            "import numpy as np\n"
            "from vakyume.config import UnsolvedException\n\n"
        )
        if "import " not in code_text:
            code_text = import_header + code_text
            
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

        # Focus on ONE family at a time for repair
        target_family = None
        if analysis["inconsistent"]:
            target_family = sorted(analysis["inconsistent"].keys())[0]
            repair_info = analysis["inconsistent"][target_family]
            repair_type = "inconsistent"
        elif analysis["failed"]:
            target_family = analysis["failed"][0]["file"]
            repair_info = {
                "broken": [], # Will be populated
                "trusted": [],
                "scores": {},
                "mismatches": {},
                "error": analysis["failed"][0].get("error")
            }
            # Populate variants for failed family
            family_dir = os.path.join(ctx.shards_dir, target_family)
            if os.path.exists(family_dir):
                repair_info["broken"] = [f.split("__")[1].replace(".py", "") for f in os.listdir(family_dir) if "__" in f]
            repair_type = "failed"

        if not target_family:
            break

        attempt = 1
        max_attempts = 3
        while attempt <= max_attempts:
            print(f" |- [Repair] Family: {target_family} | Attempt: {attempt}/{max_attempts}")
            
            for broken_variant in repair_info["broken"]:
                attempt_repair_shard(
                    ctx,
                    target_family,
                    broken_variant,
                    repair_info.get("trusted", []),
                    repair_info.get("scores", {}),
                    repair_info.get("mismatches", {}),
                    error=repair_info.get("error")
                )
            
            # Re-verify immediately
            print(f" |- Re-verifying family {target_family}...")
            new_res = verify_family(ctx, target_family)
            
            if not new_res.get("error"):
                eqn_name = f"eqn_{target_family.split('_eqn_')[1]}"
                scores = new_res["results"].get(eqn_name, {})
                num_v = len(scores)
                if num_v > 0 and all(s == num_v for s in scores.values()):
                    print(f" |- SUCCESS: Family {target_family} is now certified.")
                    break
                else:
                    print(f" |- Still inconsistent: {scores}")
                    # Update repair_info for next attempt
                    repair_info["scores"] = scores
                    repair_info["broken"] = [v for v, s in scores.items() if s < num_v]
                    repair_info["trusted"] = [v for v, s in scores.items() if s == num_v]
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
            "from math import log, sqrt, exp, pow, e\n"
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
                                method_source = get_standalone_method_source(shard_path, node.name)
                                if method_source:
                                    # Indent it for the class
                                    indented = "\n".join([f"{TAB}{line}" for line in method_source.splitlines()])
                                    class_groups[class_name].append(indented)
                        elif isinstance(node, ast.ClassDef):
                            # This is the main shard's class
                            for item in node.body:
                                if isinstance(item, ast.FunctionDef):
                                    if item.name not in seen_methods:
                                        seen_methods.add(item.name)
                                        method_source = get_method_source_from_class(shard_path, class_name, item.name)
                                        if method_source:
                                            class_groups[class_name].append(method_source)

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
        if line.startswith(f"class {class_name}:") or line.startswith(f"class {class_name}("):
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
        if in_code:
            code_lines.append(line)
    return "\n".join(code_lines).strip()
