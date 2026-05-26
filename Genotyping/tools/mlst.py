"""Genome-level MLST wrapper (tseemann/mlst).

Uses Marta's installation at
/zdata/user-data/marta/Manuscript_final/Manuscript/analyses/MLST_isolates/mlst/
with her manuscript_mlst_iso conda env for Perl dependencies.
"""

import os
from pathlib import Path

from tools.run_utils import run_cmd

MLST_BIN = Path(
    "/zdata/user-data/marta/Manuscript_final/Manuscript/analyses/MLST_isolates/mlst/bin/mlst"
)
MLST_CONDA_ENV = Path("/zdata/user-data/marta/.conda/envs/manuscript_mlst_iso")


def _mlst_env() -> dict:
    """Return env dict that puts the conda env's bin/lib on PATH/PERL5LIB."""
    env = os.environ.copy()
    env_bin = str(MLST_CONDA_ENV / "bin")
    env["PATH"] = env_bin + ":" + env.get("PATH", "")
    perl_lib = str(MLST_CONDA_ENV / "lib" / "perl5")
    env["PERL5LIB"] = perl_lib + ":" + env.get("PERL5LIB", "")
    return env


def run_mlst(
    assembly: Path,
    output_tsv: Path,
    *,
    scheme: str | None = None,
    nopath: bool = True,
    minid: float | None = None,
    mincov: float | None = None,
    threads: int | None = None,
) -> Path:
    """Run MLST on an assembled genome FASTA.

    Returns *output_tsv* path.  Output columns (TSV):
        FILE  SCHEME  ST  allele1  allele2  ...
    """
    output_tsv.parent.mkdir(parents=True, exist_ok=True)

    cmd: list[str] = [str(MLST_BIN)]
    if scheme:
        cmd += ["--scheme", scheme]
    if nopath:
        cmd.append("--nopath")
    if minid is not None:
        cmd += ["--minid", str(minid)]
    if mincov is not None:
        cmd += ["--mincov", str(mincov)]
    if threads is not None:
        cmd += ["--threads", str(threads)]
    cmd.append(str(assembly))
    result = run_cmd(cmd, env=_mlst_env())

    output_tsv.write_text(result.stdout)
    return output_tsv
