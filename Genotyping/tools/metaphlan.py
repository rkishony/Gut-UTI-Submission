"""MetaPhlAn profiling and sample2markers wrappers."""

import os
from pathlib import Path
from tools.run_utils import get_tool, run_cmd

METAPHLAN_DB_FOLDER = Path(
    "/zdata/user-data/marta/.conda/envs/manuscript_general/lib/python3.7/site-packages/metaphlan/metaphlan_databases"
)


_MANUSCRIPT_GENERAL_BIN = str(Path("/zdata/user-data/marta/.conda/envs/manuscript_general/bin"))


def _metaphlan_env() -> dict:
    """Return env dict with DEFAULT_DB_FOLDER set and manuscript_general/bin on PATH."""
    env = os.environ.copy()
    env["DEFAULT_DB_FOLDER"] = str(METAPHLAN_DB_FOLDER)
    env["PATH"] = _MANUSCRIPT_GENERAL_BIN + os.pathsep + env.get("PATH", "")
    return env


def metaphlan_profile(
    r1: Path,
    r2: Path,
    output_profile: Path,
    sam_out: Path,
    bowtie2_out: Path,
    nproc: int = 4,
) -> None:
    """Run MetaPhlAn on paired-end reads.

    Produces a taxonomic profile TSV, a SAM file (bz2), and bowtie2 output (bz2).
    """
    for p in (output_profile, sam_out, bowtie2_out):
        p.parent.mkdir(parents=True, exist_ok=True)

    run_cmd(
        [get_tool("metaphlan"),
         f"{r1},{r2}",
         "--input_type", "fastq",
         "-s", str(sam_out),
         "--bowtie2out", str(bowtie2_out),
         "--nproc", str(nproc),
         "-o", str(output_profile)],
        env=_metaphlan_env(),
        stream=True, prefix="[metaphlan] ",
    )


_SAMPLE2MARKERS_WRAPPER = (
    "import multiprocessing; multiprocessing.set_start_method('forkserver'); "
    "from metaphlan.utils.sample2markers import main; import sys; sys.exit(main())"
)

_MANUSCRIPT_GENERAL_PYTHON = str(
    Path("/zdata/user-data/marta/.conda/envs/manuscript_general/bin/python")
)


def sample2markers(
    sam_bz2: Path,
    output_dir: Path,
    nproc: int = 8,
) -> None:
    """Extract strain-level markers from a MetaPhlAn SAM file.

    Uses 'forkserver' start method to avoid pysam/htslib fork-safety bug
    in the old cmseq multiprocessing (corrupts file handles after fork).
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    run_cmd(
        [_MANUSCRIPT_GENERAL_PYTHON, "-c", _SAMPLE2MARKERS_WRAPPER,
         "-i", str(sam_bz2),
         "-o", str(output_dir),
         "-n", str(nproc)],
        env=_metaphlan_env(),
        stream=True, prefix="[sample2markers] ",
    )
