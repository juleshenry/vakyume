"""
Shared fixtures and markers for the Vakyume test suite.

Markers
-------
    @pytest.mark.llm   -- test requires a running Ollama instance (integration)
    @pytest.mark.slow   -- test takes a long time to run
"""

import os
import shutil
import tempfile

import pytest


def pytest_configure(config):
    config.addinivalue_line("markers", "llm: requires a running Ollama LLM instance")
    config.addinivalue_line("markers", "slow: long-running test")


# ── Project fixtures ────────────────────────────────────────────────────────

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VACUUM_THEORY_PROJECT = os.path.join(REPO_ROOT, "projects", "VacuumTheory")


@pytest.fixture
def repo_root():
    """Absolute path to the repository root."""
    return REPO_ROOT


@pytest.fixture
def vacuum_theory_project():
    """Absolute path to the VacuumTheory example project (may not exist)."""
    return VACUUM_THEORY_PROJECT


@pytest.fixture
def tmp_project(tmp_path):
    """Create a disposable project directory with the standard sub-folders.

    Yields the path; cleaned up automatically by pytest's tmp_path.
    """
    for sub in ("notes", "shards", "reports", "repair_prompts", "harmony_checks"):
        (tmp_path / sub).mkdir()
    return tmp_path


@pytest.fixture
def dummy_project(tmp_project):
    """A tmp_project pre-populated with the 'Basic addition' and 'Rotary' dummy shards.

    Useful for verifier and repair tests.
    """
    basic_note = "# Equation 1-1\nz = x + y\n"
    rotary_note = (
        "# Chapter 11 : Rotary Piston and Rotary Vane Pumps\n"
        "# 11-2 Sizing for Evacuation\n"
        "t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))\n"
    )

    (tmp_project / "notes" / "01_basic.py").write_text(basic_note)
    (tmp_project / "notes" / "11_rotary.py").write_text(rotary_note)

    basic_shard = (
        "from math import log, sqrt, exp, pow, e\n"
        "from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp\n"
        "from scipy.optimize import newton\n"
        "from vakyume.kwasak import kwasak\n"
        "import numpy as np\n\n"
        "class Basic:\n"
        "    @kwasak\n"
        "    def eqn_1_1(x=None, y=None, z=None, **kwargs):\n"
        "        return\n\n"
        "    @staticmethod\n"
        "    def eqn_1_1__z(x: float, y: float, **kwargs):\n"
        "        # [.pyeqn] z = x + y\n"
        "        z = x - y  # deliberate error\n"
        "        return [z]\n\n"
        "    @staticmethod\n"
        "    def eqn_1_1__x(y: float, z: float, **kwargs):\n"
        "        x = z - y\n"
        "        return [x]\n\n"
        "    @staticmethod\n"
        "    def eqn_1_1__y(x: float, z: float, **kwargs):\n"
        "        y = z - x\n"
        "        return [y]\n"
    )

    shard_dir = tmp_project / "shards" / "Basic_eqn_1_1"
    shard_dir.mkdir(parents=True, exist_ok=True)
    (shard_dir / "Basic_eqn_1_1.py").write_text(basic_shard)

    return tmp_project
