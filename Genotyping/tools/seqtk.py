"""Seqtk wrapper for FASTQ downsampling."""

from pathlib import Path

from tools.run_utils import get_tool, run_cmd


def downsample_fastq(
    input_path: Path,
    output_path: Path,
    n_reads: int,
    seed: int = 42,
) -> None:
    """Downsample a FASTQ file to *n_reads* using seqtk sample."""
    seqtk = get_tool("seqtk")
    cmd = f"{seqtk} sample -s {seed} {input_path} {n_reads} | gzip > {output_path}"
    run_cmd(cmd, shell=True)
