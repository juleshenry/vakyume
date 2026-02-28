import os
import subprocess
import shutil

# Paths
PROJECT_DIR = "dummy_project"
NOTES_DIR = os.path.join(PROJECT_DIR, "notes")
SHARDS_DIR = os.path.join(PROJECT_DIR, "shards")

# --- Dummy 1: Basic Addition (Messed Up) ---
BASIC_NOTE = """# Equation 1-1
z = x + y
"""

BASIC_SHARD_MESSED = """from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
from vakyume.kwasak import kwasak_static
import numpy as np

class Basic:
    @kwasak_static
    def eqn_1_1(x=None, y=None, z=None, **kwargs):
        return

    @staticmethod
    def eqn_1_1__z(x: float, y: float, **kwargs):
        # [.pyeqn] z = x + y
        result = []
        z = x - y  # MESSED UP: should be x + y
        result.append(z)
        return result

    @staticmethod
    def eqn_1_1__x(y: float, z: float, **kwargs):
        # [.pyeqn] z = x + y
        result = []
        x = z - y
        result.append(x)
        return result

    @staticmethod
    def eqn_1_1__y(x: float, z: float, **kwargs):
        # [.pyeqn] z = x + y
        result = []
        y = z - x
        result.append(y)
        return result
"""

# --- Dummy 2: Rotary Piston (Failing in VacuumTheory) ---
ROTARY_NOTE = """# Chapter 11 : Rotary Piston and Rotary Vane Pumps
# 11-2 Sizing for Evacuation
t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
"""

# We'll provide a partially working shard.
# We use the correct class name that vakyume would generate.
ROTARY_SHARD_MESSED = """from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from vakyume.kwasak import kwasak_static
import numpy as np

class Rotary:
    @kwasak_static
    def eqn_11_2(Q=None, Q_0=None, Q_external_gas_throughput=None, SP_1=None, SP_2=None, S_vol_pump_speed=None, V=None, t=None, **kwargs):
        return

    @staticmethod
    def eqn_11_2__t(Q: float, Q_0: float, Q_external_gas_throughput: float, SP_1: float, SP_2: float, S_vol_pump_speed: float, V: float, **kwargs):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        result = []
        # Hand-coded error to trigger repair: using + instead of -
        t = V / S_vol_pump_speed * log( (SP_1 + (Q_external_gas_throughput + Q_0))/ (SP_2 + (Q + Q_0)))
        result.append(t)
        return result

    @staticmethod
    def eqn_11_2__V(Q: float, Q_0: float, Q_external_gas_throughput: float, SP_1: float, SP_2: float, S_vol_pump_speed: float, t: float, **kwargs):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        raise Exception("Pending Repair")
"""


def reset_project():
    print("Resetting dummy_project...")
    if os.path.exists(PROJECT_DIR):
        shutil.rmtree(PROJECT_DIR)

    os.makedirs(NOTES_DIR)
    os.makedirs(SHARDS_DIR)
    os.makedirs(os.path.join(PROJECT_DIR, "reports"))
    os.makedirs(os.path.join(PROJECT_DIR, "repair_prompts"))

    with open(os.path.join(NOTES_DIR, "01_basic.py"), "w") as f:
        f.write(BASIC_NOTE)
    with open(os.path.join(NOTES_DIR, "11_rotary.py"), "w") as f:
        f.write(ROTARY_NOTE)

    with open(os.path.join(SHARDS_DIR, "Basic_eqn_1_1.py"), "w") as f:
        f.write(BASIC_SHARD_MESSED)
    with open(os.path.join(SHARDS_DIR, "Rotary_eqn_11_2.py"), "w") as f:
        f.write(ROTARY_SHARD_MESSED)


def run_vakyume():
    print("Running vakyume (NO overwrite to test repair)...")
    # We DO NOT use --overwrite here because we want to keep the messed up shards and fix them.
    # However, shard_from_chapters skips existing files, so it will proceed to verification.
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    vakyume_script = os.path.join(base_dir, "vakyume.py")
    cmd = ["python3", vakyume_script, PROJECT_DIR, "--max-rounds", "3"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(result.stdout)
    if result.stderr:
        print("Errors:")
        print(result.stderr)


def verify_fix():
    print("\n--- Verification ---")

    # Check Basic
    basic_path = os.path.join(SHARDS_DIR, "Basic_eqn_1_1.py")
    if os.path.exists(basic_path):
        with open(basic_path, "r") as f:
            content = f.read()
        if "z = x + y" in content and "z = x - y" not in content:
            print("SUCCESS: Basic equation fixed.")
        else:
            print("FAIL: Basic equation still messed up or not updated.")
    else:
        print("FAIL: Basic shard missing.")

    # Check Rotary
    rotary_path = os.path.join(SHARDS_DIR, "Rotary_eqn_11_2.py")
    if os.path.exists(rotary_path):
        with open(rotary_path, "r") as f:
            content = f.read()
        if "SP_1 - (" in content and "SP_1 + (" not in content:
            print("SUCCESS: Rotary equation fixed (signs corrected).")
        else:
            print("FAIL: Rotary equation still has errors.")
            # Print relevant line
            for line in content.splitlines():
                if "t =" in line and "def" not in line:
                    print(f"Found: {line.strip()}")
    else:
        print("FAIL: Rotary shard missing.")


if __name__ == "__main__":
    reset_project()
    run_vakyume()
    verify_fix()
