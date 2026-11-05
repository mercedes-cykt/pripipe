"""Small helpers shared by stages."""

from __future__ import annotations

import os
import sys
from pathlib import Path


def env_with_env_bin() -> dict[str, str]:
    """PATH with the current interpreter's bin dir first.

    Lets `python -m pripipe` find gffread/kallisto/... even when the conda
    environment was not `activate`d (only its interpreter was used).
    """
    env = os.environ.copy()
    bin_dir = str(Path(sys.executable).parent)
    parts = env.get("PATH", "").split(os.pathsep)
    if bin_dir not in parts:
        env["PATH"] = bin_dir + os.pathsep + env.get("PATH", "")
    return env


def missing_tool_message(tool: str) -> str:
    return (f"external tool not found: {tool} (searched PATH and "
            f"{Path(sys.executable).parent}); run `conda activate pri-miRNA-pipe`")
