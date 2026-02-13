import os
import shutil
import inspect
import importlib.util
import json
import re
import timeout_decorator
from sympy import Symbol, solve, Eq
from tru import Verify

# Import from process.llm if available, else mock
try:
    from process.llm import escribir_codigo, extract_code
except ImportError:
    def escribir_codigo(*args, **kwargs): return "# LLM Mock: No Ollama"
    def extract_code(text): return ""

from config import *

# Paths
ROOT = os.getcwd()
SHARDS_DIR = os.path.join(ROOT, "shards")
SOURCE_FILE = os.path.join(ROOT, "vakyume.py")
REPORTS_DIR = os.path.join(ROOT, "reports")
CHAPTERS_DIR = os.path.join(ROOT, "chapters")

def clear_shards():
    if os.path.exists(SHARDS_DIR):
        shutil.rmtree(SHARDS_DIR)
    os.makedirs(SHARDS_DIR)
    if os.path.exists(os.path.join(ROOT, "kwasak.py")):
        shutil.copy(os.path.join(ROOT, "kwasak.py"), os.path.join(SHARDS_DIR, "kwasak.py"))
    if not os.path.exists(REPORTS_DIR):
        os.makedirs(REPORTS_DIR)

class Solver:
    def __init__(self):
        self.funktors = FUNKTORZ

    @timeout_decorator.timeout(MAX_COMP_TIME_SECONDS, timeout_exception=StopIteration)
    def get_solns_vanilla_nf(self, nf: str, symb: Symbol):
        try:
            return solve(nf, symb)
        except Exception:
            return []

    def tokenize(self, eqn):
        malos = {"ln", "log"}
        for m in malos:
            eqn = f"{m}".join([o.strip() for o in eqn.split(m)])
        
        # Dilate functors
        eqn = eqn.split("#")[0]
        dilated = ""
        for i, f in enumerate(eqn):
            if 0 < i < len(eqn) - 1 and eqn[i + 1] != "*" and eqn[i - 1] != "*" and f in self.funktors:
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
        if len(parts) != 2: return None
        return f"({parts[1].split('#')[0].strip()}) - ({parts[0].strip()})"

    def sympy_failover(self, eqn_header, normal_form, token):
        # In a real scenario, this would call LLM. 
        # For now, we leave a placeholder or a numerical hint.
        code = [
            f"{TAB*2}# [Sympy Failover Placeholder for {token}]",
            f"{TAB*2}def func({token}):",
            f"{TAB*3}# Numerical fallback needed for: {normal_form}",
            f"{TAB*3}return eval(\"{normal_form.replace(token, 'x')}\".replace('x', str({token})))",
            f"{TAB*2}# result = [newton(func, 1.0)]",
            f"{TAB*2}return [] # Pending LLM/Manual Repair"
        ]
        return "\n".join(code)

def shard_from_chapters():
    print("Scraping equations from chapters...")
    solver = Solver()
    
    if not os.path.exists(CHAPTERS_DIR):
        print(f"Chapters directory not found: {CHAPTERS_DIR}")
        return

    for chapter_file in sorted(os.listdir(CHAPTERS_DIR)):
        if not chapter_file.endswith(".py") or chapter_file.startswith("__"):
            continue
            
        chapter_path = os.path.join(CHAPTERS_DIR, chapter_file)
        # Extract class name from filename: 01_vacuum_theory.py -> VacuumTheory
        parts = chapter_file.split("_")[1:]
        class_name = "".join(p.capitalize() for p in parts if p != "py")
        if class_name.endswith("Py"): class_name = class_name[:-2]

        with open(chapter_path, "r") as f:
            lines = f.readlines()

        eqn_number = ""
        for line in lines:
            if x := re.findall(r"\d{1,2}-\d{1,2}\w*", line):
                eqn_number = x[0].replace("-", "_")
            
            if " = " in line and not line.strip().startswith("#") and not line.strip().startswith('"""'):
                tokes = solver.get_tokes(line)
                nf = solver.make_normal_form(line)
                if not nf: continue
                
                shard_name = f"{class_name}_eqn_{eqn_number}.py"
                shard_path = os.path.join(SHARDS_DIR, shard_name)
                
                with open(shard_path, "w") as sf:
                    sf.write("from math import log, sqrt, exp, pow, e\n")
                    sf.write("from sympy import I, Piecewise, LambertW, Eq, symbols, solve\n")
                    sf.write("from scipy.optimize import newton\n")
                    sf.write("from kwasak import kwasak_static\n")
                    sf.write("import numpy as np\n\n")
                    sf.write(f"class {class_name}:\n")
                    sf.write(f"{TAB}@kwasak_static\n")
                    sf.write(f"{TAB}def eqn_{eqn_number}({', '.join(f'{t}=None' for t in tokes)}, **kwargs):\n")
                    sf.write(f"{TAB*2}return\n\n")
                    
                    for token in tokes:
                        other_args = [t for t in tokes if t != token]
                        header = f"{TAB}@staticmethod\n{TAB}def eqn_{eqn_number}__{token}({', '.join(f'{t}: float' for t in other_args)}):"
                        sf.write(header + "\n")
                        sf.write(f"{TAB*2}# [.pyeqn] {line.strip()}\n")
                        
                        try:
                            solns = solver.get_solns_vanilla_nf(nf, Symbol(token))
                            if not solns:
                                sf.write(solver.sympy_failover(header, nf, token) + "\n\n")
                            else:
                                sf.write(f"{TAB*2}result = []\n")
                                for soln in solns:
                                    sf.write(f"{TAB*2}{token} = {soln}\n")
                                    sf.write(f"{TAB*2}result.append({token})\n")
                                sf.write(f"{TAB*2}return result\n\n")
                        except StopIteration:
                            sf.write(f"{TAB*2}# Timeout during Sympy solve\n")
                            sf.write(solver.sympy_failover(header, nf, token) + "\n\n")
                        except Exception as e:
                            sf.write(f"{TAB*2}# Error during Sympy solve: {e}\n")
                            sf.write(solver.sympy_failover(header, nf, token) + "\n\n")

def verify_all_shards():
    print("Verifying shards...")
    all_results = {}
    shard_files = [f for f in os.listdir(SHARDS_DIR) if f.endswith(".py") and f != "kwasak.py"]
    
    for f in shard_files:
        shard_path = os.path.join(SHARDS_DIR, f)
        module_name = f[:-3]
        
        try:
            spec = importlib.util.spec_from_file_location(module_name, shard_path)
            if spec and spec.loader:
                module = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(module)
                
                for name, obj in inspect.getmembers(module):
                    if inspect.isclass(obj) and name != "Verify" and obj.__module__ == module_name:
                        v = Verify(obj)
                        shard_results = v.verify()
                        all_results[f] = shard_results
        except Exception as e:
            print(f"Error verifying {f}: {e}")
            
    with open(os.path.join(REPORTS_DIR, "verification_results.json"), "w") as rf:
        json.dump(all_results, rf, indent=4)
    return all_results

def analyze_results(all_results):
    print("Analyzing alignment for certification...")
    analysis = {"solved": [], "inconsistent": {}, "failed": []}
    
    for shard_file, eqn_results in all_results.items():
        if not eqn_results:
            analysis["failed"].append(shard_file)
            continue
            
        for eqn_name, variants in eqn_results.items():
            num_variants = len(variants)
            valid_scores = [s for s in variants.values() if s is not None]
            
            if valid_scores and all(s == num_variants for s in valid_scores):
                analysis["solved"].append(shard_file)
            else:
                inconsistent = [v for v, s in variants.items() if s is None or s < num_variants]
                analysis["inconsistent"][shard_file] = {
                    "eqn": eqn_name,
                    "num_variants": num_variants,
                    "inconsistent_variants": inconsistent,
                    "scores": variants
                }
    
    with open(os.path.join(REPORTS_DIR, "analysis.json"), "w") as af:
        json.dump(analysis, af, indent=4)
    return analysis

def assemble_certified_library(analysis):
    print(f"Assembling certified library from {len(analysis['solved'])} aligned equations...")
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
                        in_class = True; continue
                    if in_class: body.append(line)
                
                if class_name not in class_groups: class_groups[class_name] = []
                class_groups[class_name].extend(body)

        for class_name, bodies in class_groups.items():
            out.write(f"class {class_name}:\n")
            out.writelines(bodies)
            out.write("\n")
    print(f"Certified library created: {certified_file}")

if __name__ == "__main__":
    # If the user wants to re-shard from chapters, call shard_from_chapters()
    # Otherwise, just verify existing shards.
    # For now, let's keep existing shards if they exist, but allow re-sharding if requested.
    
    # clear_shards()
    # shard_from_chapters()
    
    results = verify_all_shards()
    analysis = analyze_results(results)
    assemble_certified_library(analysis)
    
    print("\n" + "="*40)
    print(f"SOLVED (Aligned):      {len(analysis['solved'])}")
    print(f"INCONSISTENT:          {len(analysis['inconsistent'])}")
    print(f"UNSOLVABLE/FAILED:     {len(analysis['failed'])}")
    print("="*40)
