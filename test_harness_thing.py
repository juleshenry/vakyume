import sys
import importlib
import inspect
import re
import random
import traceback
import math
import requests
import numpy as np
from typing import Dict, Any, List

import os
sys.path.insert(
    0, os.path.join(os.path.dirname(__file__), "projects/VacuumTheory/shards")
)
TARGET_MODULE = "LiquidRing_eqn_10_19"  
# --- Configuration ---
TARGET_MODULE = "liquid_ring"
OLLAMA_API_URL = "http://localhost:11434/api/generate"
MODEL_NAME = "phi3" # or 'mistral', 'llama3'

# --- The "Paranoid" Harness ---
class ParanoidHarness:
    def __init__(self, module_name):
        self.module_name = module_name
        self.module = None
    
    def load_module(self):
        """Force reloads the module from disk."""
        if self.module_name in sys.modules:
            self.module = importlib.reload(sys.modules[self.module_name])
        else:
            self.module = importlib.import_module(self.module_name)
        return self.module

    def find_method_recursively(self, cls, method_name):
        """
        Hunts for a method even if it's hidden in nested classes 
        (like your EquationSolver example).
        """
        # Check current class
        if hasattr(cls, method_name):
            return getattr(cls, method_name)
            
        # Check nested classes
        for name, obj in inspect.getmembers(cls):
            if inspect.isclass(obj):
                # Prevent infinite recursion if generic
                if obj is not cls: 
                    found = self.find_method_recursively(obj, method_name)
                    if found: return found
        return None

    def get_source_location(self, method_name):
        """Finds where the method is defined in the file (Class or Nested Class)."""
        # This is tricky with nested classes. We scan the text file for the def.
        # A robust implementation would use AST, but regex is faster for repairs.
        with open(self.module.__file__, 'r') as f:
            content = f.read()
            
        # Regex to capture the function signature block
        # matches: def eqn_10_19__P(...)
        pattern = re.compile(rf"def {method_name}\s*\(.*?\):", re.DOTALL)
        match = pattern.search(content)
        return match

    def generate_repair_prompt(self, method_name, source_code, error_trace, inputs):
        """
        Constructs a prompt that forces the LLM to look at the Traceback.
        """
        return f"""
        You are a Python Expert fixing broken scientific code.
        
        The method `{method_name}` failed during execution.
        
        CODE:
        {source_code}
        
        ERROR TRACEBACK:
        {error_trace}
        
        CONTEXT:
        The inputs causing failure were: {inputs}
        
        TASK:
        1. Fix the error. If it is a NameError (like P_s vs p_s), correct the variable name.
        2. If it is a convergence error, improve the initial guess.
        3. Return ONLY the corrected python function code.
        """

    def apply_patch(self, method_name, old_code, new_code):
        """
        Replaces the broken function in the file with the LLM's fix.
        """
        file_path = self.module.__file__
        with open(file_path, 'r') as f:
            content = f.read()

        # We need to be careful to replace the correct block.
        # We assume the LLM returns the full "def method(...): ..." block.
        
        # 1. Normalize indentation of new code to match the file's context
        # Find indentation of the old definition
        match = re.search(rf"(\s*)def {method_name}", content)
        if not match:
            print(f"❌ Critical: Could not locate {method_name} in source file to patch.")
            return False
            
        indent = match.group(1)
        indented_new_code = "\n".join([indent + line if line.strip() else line for line in new_code.splitlines()])
        
        # 2. Naive replace: finding the start of def and assuming it ends at next def or unindent
        # (For this demo, we'll rely on the fact that your broken code allows us to locate it via regex)
        # A safer way is to ask LLM to provide the whole class, but that's token-heavy.
        
        # Let's try to locate the specific block of the *broken* code to replace it.
        # Since we don't have the exact broken text in `old_code` (inspect might normalize it),
        # we construct a regex that matches the start of the function.
        
        start_idx = match.start()
        
        # Find the next function definition or class to mark the end
        # This is a heuristic.
        next_def = re.search(r"\n\s*(def|class|@)", content[start_idx + len(match.group(0)):])
        
        if next_def:
            end_idx = start_idx + len(match.group(0)) + next_def.start()
        else:
            end_idx = len(content) # End of file

        new_file_content = content[:start_idx] + indented_new_code + content[end_idx:]
        
        with open(file_path, 'w') as f:
            f.write(new_file_content)
        
        print(f"🩹 Patched {method_name} in {file_path}")
        return True

    def verify_cycle(self, family_name="eqn_10_19"):
        print(f"\n🔬 Verifying Family: {family_name}")
        
        # 1. Generate Random Inputs
        inputs = {
            'P': random.uniform(20, 100), 'S_Th': random.uniform(1, 5),
            'S_p': random.uniform(1, 5), 'T_e': random.uniform(100, 200),
            'T_i': random.uniform(50, 150), 'p_c': random.uniform(1, 5),
            'p_s': random.uniform(1, 5)
        }
        
        # 2. Anchor Step: Calculate S_Th using the (assumed working) S_Th method
        # We need to find the method first
        cls = getattr(self.module, 'LiquidRing')
        s_th_func = self.find_method_recursively(cls, f"{family_name}__S_Th")
        
        try:
            # Generate a Consistent State (Truth Tuple)
            calc_s_th = s_th_func(**inputs)[0]
            inputs['S_Th'] = calc_s_th.real # Flatten complex if necessary
            print(f"⚓ Anchor established. S_Th calculated as {inputs['S_Th']:.4f}")
        except Exception as e:
            print(f"💥 Anchor Failed: {e}")
            return False

        # 3. Test The Problematic Method (eqn_10_19__P)
        target_method_name = f"{family_name}__P"
        target_func = self.find_method_recursively(cls, target_method_name)
        
        if not target_func:
            print(f"❌ Method {target_method_name} not found in class or nested classes!")
            # Should trigger creation/repair, but for now we report missing.
            return False

        print(f"⚔️  Testing {target_method_name}...", end="")
        
        # Prepare inputs (hide P)
        test_inputs = inputs.copy()
        expected_P = test_inputs.pop('P')

        try:
            # --- EXECUTION ---
            result = target_func(**test_inputs)
            
            # --- VALIDATION ---
            # Check if result is empty or garbage
            if not result: raise ValueError("Method returned empty result")
            calculated_P = result[0]
            
            # Math Check
            if isinstance(calculated_P, complex): calculated_P = calculated_P.real
            
            error_margin = abs(calculated_P - expected_P)
            if error_margin > 1.0: # Allow some slip for numerical solvers
                print(f" ❌ MATH FAILURE. Expected {expected_P:.2f}, Got {calculated_P:.2f}")
                raise ValueError(f"Math verification failed. Error margin: {error_margin}")
            
            print(f" ✅ PASS (Error: {error_margin:.4f})")
            return True

        except Exception as e:
            print(f"\n🔥 CRASH DETECTED: {e}")
            
            # Capture the traceback for the LLM
            tb = traceback.format_exc()
            
            # Get source code
            try:
                src = inspect.getsource(target_func)
            except:
                src = "# Could not retrieve source (maybe dynamic?)"

            # Call Repair
            repaired_code = self.call_ollama(target_method_name, src, tb, test_inputs)
            if repaired_code:
                self.apply_patch(target_method_name, src, repaired_code)
                return False # Signal to retry
            return False

    def call_ollama(self, name, code, error, inputs):
        print(f"🚑 Calling {MODEL_NAME} to fix {name}...")
        prompt = self.generate_repair_prompt(name, code, error, inputs)
        
        try:
            res = requests.post(OLLAMA_API_URL, json={
                "model": MODEL_NAME, "prompt": prompt, "stream": False
            })
            ans = res.json()['response']
            # Clean markdown
            ans = ans.replace("```python", "").replace("```", "").strip()
            return ans
        except Exception as e:
            print(f"LLM Error: {e}")
            return None

if __name__ == "__main__":
    harness = ParanoidHarness(TARGET_MODULE)
    
    # Run loop
    for _ in range(3): # Max 3 repair attempts
        harness.load_module()
        success = harness.verify_cycle()
        if success:
            print("\n✨ System Stabilized.")
            break
        print("\n🔄 Reloading for retry...")