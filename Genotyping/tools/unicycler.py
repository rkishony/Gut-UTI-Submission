"""Unicycler wrapper for genome assembly."""

import os
import subprocess
from pathlib import Path

from tools.run_utils import get_tool, run_cmd


def run_unicycler(
    r1: Path,
    r2: Path,
    out_dir: Path,
    threads: int = 4,
) -> subprocess.CompletedProcess:
    """Run Unicycler on a pair of FASTQ files."""
    unicycler = get_tool("unicycler")
    cmd = [
        unicycler,
        "-1", str(r1),
        "-2", str(r2),
        "-o", str(out_dir),
        "-t", str(threads),
    ]

    # Unicycler's conda env has dependencies (spades, racon, etc.)
    # that must be on PATH.
    env = os.environ.copy()
    env["PATH"] = str(unicycler.parent) + os.pathsep + env.get("PATH", "")

    return run_cmd(cmd, check=False, env=env)
