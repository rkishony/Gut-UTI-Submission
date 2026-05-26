"""ResFinder / PointFinder wrapper.

Uses Marta's older resfinder installation at
/zdata/user-data/marta/marta/pathometagen_pipeline/resfindertool/
with her own python env and KMA binary.
"""

from pathlib import Path

from config import RESFINDER_DB, POINTFINDER_DB, TOOLS
from tools.run_utils import run_cmd

RESFINDER_SCRIPT = TOOLS["resfinder"]
RESFINDER_PYTHON = TOOLS["resfinder-python"]
RESFINDER_KMA = TOOLS["resfinder-kma"]


def run_resfinder(
    r1: Path,
    r2: Path,
    output_dir: Path,
    *,
    min_cov: float = 0.99,
) -> Path:
    """Run ResFinder for acquired resistance genes (KMA on paired fastq)."""
    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        RESFINDER_PYTHON, RESFINDER_SCRIPT,
        "-ifq", str(r1), str(r2),
        "-o", str(output_dir),
        "-db_res", str(RESFINDER_DB),
        "-db_res_kma", str(RESFINDER_DB),
        "-k", RESFINDER_KMA,
        "-acq",
        "-l", str(min_cov),
    ]
    run_cmd(cmd)
    return output_dir


def run_pointfinder(
    r1: Path,
    r2: Path,
    output_dir: Path,
    species: str = "Escherichia coli",
) -> Path:
    """Run PointFinder for chromosomal point mutations (KMA on paired fastq)."""
    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        RESFINDER_PYTHON, RESFINDER_SCRIPT,
        "-ifq", str(r1), str(r2),
        "-o", str(output_dir),
        "-db_res", str(RESFINDER_DB),
        "-db_res_kma", str(RESFINDER_DB),
        "-db_point", str(POINTFINDER_DB),
        "-db_point_kma", str(POINTFINDER_DB),
        "-k", RESFINDER_KMA,
        "-c",
        "-s", species,
    ]
    run_cmd(cmd)
    return output_dir
