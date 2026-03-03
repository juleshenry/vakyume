import os
import sys
from vakyume.verifier import Verify
from vakyume.master import PipelineContext, analyze_results


class Basic:
    @staticmethod
    def eqn_1_1(x=None, y=None, z=None, **kwargs):
        return

    @staticmethod
    def eqn_1_1__z(x: float, y: float, **kwargs):
        return [x - y]  # MESSED UP: should be x + y

    @staticmethod
    def eqn_1_1__x(y: float, z: float, **kwargs):
        return [z - y]  # CORRECT

    @staticmethod
    def eqn_1_1__y(x: float, z: float, **kwargs):
        return [z - x]  # CORRECT


# Mocking the context for subshards
ctx = PipelineContext("dummy_project")
pyeqn = "z = x + y"

v = Verify(Basic, pyeqn=pyeqn, subshards_dir=ctx.subshards_dir)
results = v.verify()
print("Verification results:", results)

# Flatten results like verify_single_shard does
verification_results = {}
mismatches = {}
for base_eq, data in results.items():
    verification_results[base_eq] = data["scores"]
    mismatches.update(data["mismatches"])

# Mocking the input for analyze_results
all_results = {
    "Basic_eqn_1_1.py": {"results": verification_results, "mismatches": mismatches}
}

analysis = analyze_results(ctx, all_results)
import json

print("Analysis:", json.dumps(analysis, indent=2))
