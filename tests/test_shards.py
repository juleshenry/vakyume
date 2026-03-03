"""
Parametrized shard-consistency tests.

Dynamically discovers shard files from ``projects/VacuumTheory/shards/``
(subdirectory layout) and verifies each one using the ``Verify`` class.

Skipped automatically when the VacuumTheory project is not present.
"""

import importlib.util
import os

import pytest

from vakyume.verifier import Verify


# ── Discover shards at collection time ──────────────────────────────────────

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SHARDS_DIR = os.path.join(REPO_ROOT, "projects", "VacuumTheory", "shards")


def _collect_shard_files():
    """Walk ``SHARDS_DIR`` and yield (family_dir, shard_file) tuples."""
    if not os.path.isdir(SHARDS_DIR):
        return
    for family_name in sorted(os.listdir(SHARDS_DIR)):
        family_path = os.path.join(SHARDS_DIR, family_name)
        if not os.path.isdir(family_path):
            continue
        for fname in sorted(os.listdir(family_path)):
            if fname.endswith(".py") and not fname.startswith("_"):
                yield family_name, fname


_SHARD_PARAMS = list(_collect_shard_files())


@pytest.mark.skipif(
    not os.path.isdir(SHARDS_DIR),
    reason="VacuumTheory project not found (run the pipeline first)",
)
@pytest.mark.parametrize(
    "family,shard_file",
    _SHARD_PARAMS,
    ids=[f"{f}/{s}" for f, s in _SHARD_PARAMS],
)
def test_shard_consistency(family, shard_file):
    """Load a shard module and verify all variant solvers agree."""
    shard_path = os.path.join(SHARDS_DIR, family, shard_file)
    module_name = shard_file[:-3]

    spec = importlib.util.spec_from_file_location(module_name, shard_path)
    assert spec is not None, f"Cannot create module spec for {shard_path}"
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    for name, obj in module.__dict__.items():
        if isinstance(obj, type) and name != "Verify":
            v = Verify(obj)
            res = v.verify()
            for eqn, data in res.items():
                scores = data["scores"] if isinstance(data, dict) else data
                if not scores:
                    continue
                valid_scores = [s for s in scores.values() if s is not None]
                if not valid_scores:
                    continue
                max_score = max(valid_scores)
                for var, score in scores.items():
                    if score is not None:
                        assert score == max_score, (
                            f"{module_name}: equation {eqn} variant {var} "
                            f"inconsistent (score {score}/{max_score})"
                        )
