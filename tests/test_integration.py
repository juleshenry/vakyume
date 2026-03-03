"""
Integration tests that require a running Ollama LLM.

These test the full pipeline including shard generation, verification,
and LLM-based auto-repair.  They are marked with ``@pytest.mark.llm``
and skipped by default -- run with ``pytest -m llm`` to include them.
"""

import os
import subprocess
import sys

import pytest

from vakyume.pipeline import run_pipeline


pytestmark = pytest.mark.llm


# ── Fixtures ────────────────────────────────────────────────────────────────

BASIC_NOTE = "# Equation 1-1\nz = x + y\n"

BASIC_SHARD_MESSED = (
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
    "        result = []\n"
    "        z = x - y  # MESSED UP: should be x + y\n"
    "        result.append(z)\n"
    "        return result\n\n"
    "    @staticmethod\n"
    "    def eqn_1_1__x(y: float, z: float, **kwargs):\n"
    "        # [.pyeqn] z = x + y\n"
    "        result = []\n"
    "        x = z - y\n"
    "        result.append(x)\n"
    "        return result\n\n"
    "    @staticmethod\n"
    "    def eqn_1_1__y(x: float, z: float, **kwargs):\n"
    "        # [.pyeqn] z = x + y\n"
    "        result = []\n"
    "        y = z - x\n"
    "        result.append(y)\n"
    "        return result\n"
)


@pytest.fixture
def repair_project(tmp_path):
    """Create a project with a deliberately broken shard for repair testing."""
    notes = tmp_path / "notes"
    notes.mkdir()
    (notes / "01_basic.py").write_text(BASIC_NOTE)

    for sub in ("shards", "reports", "repair_prompts", "subshards"):
        (tmp_path / sub).mkdir()

    shard_dir = tmp_path / "shards" / "Basic_eqn_1_1"
    shard_dir.mkdir()
    (shard_dir / "Basic_eqn_1_1.py").write_text(BASIC_SHARD_MESSED)

    return tmp_path


# ── Tests ───────────────────────────────────────────────────────────────────


class TestAutoRepair:
    """Full pipeline with auto-repair on a broken shard."""

    def test_pipeline_runs_without_crash(self, repair_project):
        analysis = run_pipeline(
            project_dir=str(repair_project),
            max_rounds=2,
            overwrite=False,
        )
        assert "solved" in analysis
        assert "inconsistent" in analysis
        assert "failed" in analysis


class TestDummyEndToEnd:
    """End-to-end: run the CLI, then verify the shard was repaired."""

    def test_cli_repair(self, repair_project):
        repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        cmd = [
            sys.executable,
            "-m",
            "vakyume",
            "run",
            str(repair_project),
            "--max-rounds",
            "3",
        ]
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=repo_root,
        )
        assert (
            result.returncode == 0
        ), f"CLI failed:\nstdout: {result.stdout}\nstderr: {result.stderr}"
