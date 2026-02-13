import importlib.util
import os
import inspect
from tru import Verify

SHARDS_DIR = 'shards'
f = 'FluidFlowVacuumLines_eqn_2_08.py'
shard_path = os.path.join(SHARDS_DIR, f)
module_name = f[:-3]

spec = importlib.util.spec_from_file_location(module_name, shard_path)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

for name, obj in inspect.getmembers(module):
    if inspect.isclass(obj) and name != "Verify":
        print(f"Verifying {name}...")
        v = Verify(obj)
        shard_results = v.verify()
        print(f"Results: {shard_results}")
