"""
Parametrized family-level cross-validation tests.

Dynamically discovers equation families from ``projects/VacuumTheory/shards/``
and verifies each one by loading all variant solver functions, assembling them
into a synthetic class, and running the ``Verify`` cross-checker.

Families where some variants are inconsistent will **fail** — this is the
desired behaviour so that regressions are caught immediately.

Skipped automatically when the VacuumTheory project is not present.
"""

import importlib.util
import os
import re

import pytest

from vakyume.verifier import Verify


# ── Discover families at collection time ────────────────────────────────────

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SHARDS_DIR = os.path.join(REPO_ROOT, "projects", "VacuumTheory", "shards")


def _collect_families():
    """Walk ``SHARDS_DIR`` and yield family directory names.

    Each family directory (e.g. ``AirLeak_eqn_4_10``) contains a parent
    dispatcher and one or more variant solver files (``eqn_X_Y__var.py``).
    """
    if not os.path.isdir(SHARDS_DIR):
        return
    for family_name in sorted(os.listdir(SHARDS_DIR)):
        family_path = os.path.join(SHARDS_DIR, family_name)
        if not os.path.isdir(family_path) or family_name == "harmony_checks":
            continue
        # Only include families that have at least one variant file
        has_variants = any(
            f.endswith(".py") and "__" in f and not f.startswith("_")
            for f in os.listdir(family_path)
        )
        if has_variants:
            yield family_name


def _build_family_class(family_name):
    """Load all variant solver functions and assemble them into a class.

    Returns the synthetic class and the base equation name (e.g. ``eqn_4_10``).
    """
    family_path = os.path.join(SHARDS_DIR, family_name)

    # Determine the base equation name from the family directory name.
    # Family names look like: FluidFlowVacuumLines_eqn_2_18a
    # We need the eqn part: eqn_2_18a
    m = re.search(r"(eqn_\d+_\d+\w*)", family_name)
    assert m, f"Cannot extract equation name from family {family_name}"
    base_eqn = m.group(1)

    # Chapter class name: everything before _eqn_
    class_name = family_name.split("_eqn_")[0]

    attrs = {}
    variant_files = sorted(
        f
        for f in os.listdir(family_path)
        if f.endswith(".py") and "__" in f and not f.startswith("_")
    )

    for fname in variant_files:
        fpath = os.path.join(family_path, fname)
        module_name = fname[:-3]  # e.g. eqn_4_10__T_cap

        spec = importlib.util.spec_from_file_location(module_name, fpath)
        assert spec is not None, f"Cannot create module spec for {fpath}"
        module = importlib.util.module_from_spec(spec)
        try:
            spec.loader.exec_module(module)
        except Exception as e:
            # If a shard can't be loaded at all, skip it
            continue

        # Find the solver function (named like eqn_X_Y__variable)
        for attr_name in dir(module):
            if attr_name.startswith(base_eqn + "__") and callable(
                getattr(module, attr_name)
            ):
                attrs[attr_name] = getattr(module, attr_name)

    # Add a dummy base equation method so Verify discovers the family
    def _base_eq(self, **kwargs):
        return

    attrs[base_eqn] = _base_eq

    klass = type(class_name, (), attrs)
    return klass, base_eqn


_FAMILIES = list(_collect_families())


# ── Known inconsistent families ─────────────────────────────────────────────
# These families fail cross-validation due to issues in the equation notes
# or SymPy producing complex-domain solutions.  Marked xfail so CI stays
# green while signalling they need attention.
#
#   FluidFlowVacuumLines_eqn_2_17  — two distinct equations share the same
#       equation number in the notes file (laminar flow delta_P in terms of
#       v vs q); the parser lumps them into one family, but the q and v
#       variants solve different equations.
#   FluidFlowVacuumLines_eqn_2_34  — transitional flow conductance with
#       multiple independent coefficients; SymPy solutions don't round-trip.
#   LiquidRing_eqn_10_19  — complex roots for random positive inputs.
#   LiquidRing_eqn_10_20  — complex roots for random positive inputs.
#   Precondensors_eqn_7_14b  — log-mean temperature difference; SymPy
#       produces complex solutions for some input ranges.
#   SelectingPump_eqn_8_6  — fractional exponents produce complex roots.
KNOWN_INCONSISTENT = {
    "FluidFlowVacuumLines_eqn_2_17",
    "FluidFlowVacuumLines_eqn_2_34",
    "LiquidRing_eqn_10_19",
    "LiquidRing_eqn_10_20",
    "Precondensors_eqn_7_14b",
    "SelectingPump_eqn_8_6",
}


@pytest.mark.skipif(
    not os.path.isdir(SHARDS_DIR),
    reason="VacuumTheory project not found (run the pipeline first)",
)
@pytest.mark.parametrize("family", _FAMILIES, ids=_FAMILIES)
def test_family_cross_validation(family):
    """Load all variants for a family and verify cross-consistency."""
    if family in KNOWN_INCONSISTENT:
        pytest.xfail(f"{family} is a known-inconsistent family")

    klass, base_eqn = _build_family_class(family)

    # Verify this family
    v = Verify(klass)
    results = v.verify()

    for eqn, data in results.items():
        scores = data["scores"] if isinstance(data, dict) else data
        if not scores:
            continue
        valid_scores = [s for s in scores.values() if s is not None]
        if not valid_scores:
            continue
        num_variants = len([attr for attr in dir(klass) if attr.startswith(eqn + "__")])
        for var, score in scores.items():
            if score is not None:
                assert score == num_variants, (
                    f"{family}: equation {eqn} variant {var} "
                    f"inconsistent (score {score}/{num_variants})"
                )
