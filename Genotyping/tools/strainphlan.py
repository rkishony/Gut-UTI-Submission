"""StrainPhlAn wrapper using Marta's installation and database.

Uses her manuscript_general conda env for strainphlan,
her pre-extracted DB markers (t__SGB10068), and her reference genomes.
"""

import os
from pathlib import Path

from config import (
    STRAINPHLAN_DB_MARKERS,
    STRAINPHLAN_REF_GENOMES_DIR,
    STRAINPHLAN_DEFAULT_CLADE,
)
from tools.metaphlan import METAPHLAN_DB_FOLDER
from tools.run_utils import get_tool, run_cmd


def _strainphlan_env() -> dict:
    env = os.environ.copy()
    env["DEFAULT_DB_FOLDER"] = str(METAPHLAN_DB_FOLDER)
    manuscript_bin = "/zdata/user-data/marta/.conda/envs/manuscript_general/bin"
    env["PATH"] = manuscript_bin + ":" + env.get("PATH", "")
    return env


def run_strainphlan(
    consensus_markers: list[Path],
    output_dir: Path,
    clade: str = STRAINPHLAN_DEFAULT_CLADE,
    db_markers: Path = STRAINPHLAN_DB_MARKERS,
    ref_genomes_dir: Path = STRAINPHLAN_REF_GENOMES_DIR,
    nproc: int = 8,
    mutation_rates: bool = True,
) -> Path:
    """Run strainphlan on consensus marker .pkl files.

    Returns the output tree path.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    ref_fnas = sorted(ref_genomes_dir.glob("*.fna"))

    cmd = [
        get_tool("strainphlan"),
        "-s", *[str(p) for p in consensus_markers],
        "-m", str(db_markers),
        "-o", str(output_dir),
        "-n", str(nproc),
        "-c", clade,
    ]
    for ref in ref_fnas:
        cmd += ["-r", str(ref)]
    if mutation_rates:
        cmd.append("--mutation_rates")

    run_cmd(cmd, env=_strainphlan_env(), stream=True, prefix="[strainphlan] ")

    return output_dir / f"RAxML_bestTree.{clade}.StrainPhlAn4.tre"
