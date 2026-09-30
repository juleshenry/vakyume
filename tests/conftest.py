import os
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
PROJECTS = REPO_ROOT / "projects"

# The old shard/verifier suite needs the `legacy` extra; opt in explicitly.
collect_ignore_glob = [] if os.environ.get("VAKYUME_LEGACY_TESTS") else ["legacy/*"]


@pytest.fixture(autouse=True, scope="session")
def _solver_cache(tmp_path_factory):
    # keep sympy.solve results out of the user's ~/.cache during tests
    os.environ.setdefault("VAKYUME_CACHE_DIR", str(tmp_path_factory.mktemp("solve-cache")))
