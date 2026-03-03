"""
Unit tests for the Vakyume verifier.

These exercise the ``Verify`` class on small inline equation classes
without requiring Ollama or any project on disk.
"""

import json
import pytest

from vakyume.verifier import Verify
from vakyume.master import PipelineContext, analyze_results


# ── Inline equation classes used by the tests ───────────────────────────────


class _BasicCorrect:
    """z = x + y  (all solvers correct)."""

    @staticmethod
    def eqn_1_1(x=None, y=None, z=None, **kwargs):
        return

    @staticmethod
    def eqn_1_1__z(x: float, y: float, **kwargs):
        return [x + y]

    @staticmethod
    def eqn_1_1__x(y: float, z: float, **kwargs):
        return [z - y]

    @staticmethod
    def eqn_1_1__y(x: float, z: float, **kwargs):
        return [z - x]


class _BasicBroken:
    """z = x + y  (``__z`` solver deliberately wrong: x - y)."""

    @staticmethod
    def eqn_1_1(x=None, y=None, z=None, **kwargs):
        return

    @staticmethod
    def eqn_1_1__z(x: float, y: float, **kwargs):
        return [x - y]  # deliberate error

    @staticmethod
    def eqn_1_1__x(y: float, z: float, **kwargs):
        return [z - y]

    @staticmethod
    def eqn_1_1__y(x: float, z: float, **kwargs):
        return [z - x]


class _RotaryBroken:
    """Rotary piston equation with wrong signs and a pending-repair variant."""

    @staticmethod
    def eqn_11_2(
        Q=None,
        Q_0=None,
        Q_external_gas_throughput=None,
        SP_1=None,
        SP_2=None,
        S_vol_pump_speed=None,
        V=None,
        t=None,
        **kwargs,
    ):
        return

    @staticmethod
    def eqn_11_2__t(
        Q: float,
        Q_0: float,
        Q_external_gas_throughput: float,
        SP_1: float,
        SP_2: float,
        S_vol_pump_speed: float,
        V: float,
        **kwargs,
    ):
        from math import log

        t = (
            V
            / S_vol_pump_speed
            * log((SP_1 + (Q_external_gas_throughput + Q_0)) / (SP_2 + (Q + Q_0)))
        )
        return [t]

    @staticmethod
    def eqn_11_2__V(
        Q: float,
        Q_0: float,
        Q_external_gas_throughput: float,
        SP_1: float,
        SP_2: float,
        S_vol_pump_speed: float,
        t: float,
        **kwargs,
    ):
        raise Exception("Pending Repair")


# ── Tests ───────────────────────────────────────────────────────────────────


class TestVerifyCorrect:
    """Verify passes on a fully correct equation class."""

    def test_all_variants_consistent(self, tmp_path):
        v = Verify(_BasicCorrect, pyeqn="z = x + y", subshards_dir=str(tmp_path))
        results = v.verify()
        for eqn, data in results.items():
            scores = data["scores"]
            valid = [s for s in scores.values() if s is not None]
            assert len(valid) > 0, f"No valid scores for {eqn}"
            assert all(s == valid[0] for s in valid), (
                f"Scores should all match for correct equation: {scores}"
            )


class TestVerifyBroken:
    """Verify catches a broken solver variant."""

    def test_detects_mismatch(self, tmp_path):
        v = Verify(_BasicBroken, pyeqn="z = x + y", subshards_dir=str(tmp_path))
        results = v.verify()
        for eqn, data in results.items():
            mismatches = data.get("mismatches", {})
            # The broken __z variant should be flagged
            assert len(mismatches) > 0, (
                "Expected at least one mismatch for the broken equation"
            )


class TestVerifyPendingRepair:
    """Verify gracefully handles a variant that raises an exception."""

    def test_exception_variant_returns_none_score(self, tmp_path):
        pyeqn = (
            "t = V / S_vol_pump_speed * ln( "
            "(SP_1 - (Q_external_gas_throughput + Q_0))"
            "/ (SP_2 - (Q + Q_0)))"
        )
        v = Verify(_RotaryBroken, pyeqn=pyeqn, subshards_dir=str(tmp_path))
        results = v.verify()
        # At least one variant should exist; the one that raises should
        # produce a None score or be excluded.
        for eqn, data in results.items():
            scores = data["scores"]
            # eqn_11_2__V raises, so its score should be None or absent
            v_score = scores.get("eqn_11_2__V")
            assert v_score is None or v_score == 0, (
                f"Expected None/0 for raising variant, got {v_score}"
            )


class TestAnalyzeResults:
    """analyze_results correctly categorises solved / inconsistent / failed."""

    def test_broken_is_inconsistent(self, tmp_path):
        ctx = PipelineContext(str(tmp_path))
        v = Verify(_BasicBroken, pyeqn="z = x + y", subshards_dir=str(tmp_path))
        raw = v.verify()

        verification_results = {}
        mismatches = {}
        for base_eq, data in raw.items():
            verification_results[base_eq] = data["scores"]
            mismatches.update(data["mismatches"])

        all_results = {
            "BasicBroken_eqn_1_1.py": {
                "results": verification_results,
                "mismatches": mismatches,
            }
        }
        analysis = analyze_results(ctx, all_results)
        assert len(analysis["inconsistent"]) > 0, (
            "Broken equation should be classified as inconsistent"
        )
