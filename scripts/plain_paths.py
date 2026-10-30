"""Portable locations for archived scripts; defaults contain no study data."""

import os
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def plain_path(root, relative=""):
    if root == "PYTHON":
        return os.environ.get("PLAIN_PYTHON", sys.executable)
    defaults = {
        "REPO": REPO,
        "PROJECT": REPO / "private_workspace/project",
        "HEATMAP": REPO / "private_workspace/heatmap",
        "LEGACY_WORK": REPO / "private_workspace/legacy_work",
        "WORKSPACE": REPO / "private_workspace/workspace",
        "DESKTOP_OUTPUTS": REPO / "private_workspace/outputs",
        "R_LIBRARY": REPO / "private_workspace/r_library",
    }
    if root not in defaults:
        raise ValueError("Unknown path root: " + root)
    base = Path(os.environ.get("PLAIN_" + root + "_ROOT", str(defaults[root]))).expanduser()
    result = str(base / relative)
    # Some original R/Python code concatenates the next filename fragment.
    if relative.endswith("/") and not result.endswith(os.sep):
        result += os.sep
    return result
