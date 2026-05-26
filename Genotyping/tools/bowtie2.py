"""Bowtie2 / samtools / bedtools wrappers for paired-end read alignment."""

from collections import defaultdict
from pathlib import Path
from typing import Optional

import numpy as np

from tools.run_utils import get_tool, run_cmd


def bowtie2_build(ref_fasta: Path, index_prefix: Path, threads: int = 4) -> None:
    """Build bowtie2 index from a FASTA reference."""
    index_prefix.parent.mkdir(parents=True, exist_ok=True)
    run_cmd([get_tool("bowtie2-build"), "--quiet", "--threads", str(threads),
             str(ref_fasta), str(index_prefix)])


def bowtie2_align_paired(
    index_prefix: Path,
    r1: Path,
    r2: Path,
    output_sam: Path,
    threads: int = 6,
    met_file: Optional[Path] = None,
) -> None:
    """Align paired-end reads with bowtie2 --very-fast."""
    output_sam.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        get_tool("bowtie2"), "--very-fast",
        "-p", str(threads),
        "-x", str(index_prefix),
        "-1", str(r1),
        "-2", str(r2),
        "-S", str(output_sam),
    ]
    if met_file is not None:
        met_file.parent.mkdir(parents=True, exist_ok=True)
        cmd += ["--met-file", str(met_file)]

    result = run_cmd(cmd, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"bowtie2 failed: {result.stderr.strip()}")
    return result


def sam_to_sorted_bam(
    sam_path: Path,
    sorted_bam_path: Path,
    threads: int = 6,
) -> None:
    """Convert SAM → sorted BAM and index it."""
    samtools = get_tool("samtools")
    run_cmd([samtools, "sort", "-@", str(threads),
             "-o", str(sorted_bam_path), str(sam_path)])
    run_cmd([samtools, "index", str(sorted_bam_path)])


def bedtools_genomecov(
    sorted_bam: Path,
    output_path: Path,
) -> None:
    """Run bedtools genomecov -d and save as compressed .npz.

    The .npz contains:
      coverage       - flat int16 array (all contigs concatenated)
      contig_names   - string array of contig names, in order
      contig_offsets - int array of start positions into coverage
    """
    cmd = [get_tool("bedtools"), "genomecov", "-d", "-ibam", str(sorted_bam)]
    result = run_cmd(cmd, check=True)

    # Parse per-base depths, grouped by contig (dict preserves insertion order)
    contig_depths: dict[str, list[int]] = defaultdict(list)
    for line in result.stdout.splitlines():
        contig, _pos, depth = line.split("\t")
        contig_depths[contig].append(int(depth))

    # Build flat array + offsets
    contig_names = list(contig_depths.keys())
    coverage = np.concatenate(list(contig_depths.values()))
    assert coverage.max() <= np.iinfo(np.int32).max, f"Depth {coverage.max()} exceeds int32 range"
    coverage = coverage.astype(np.int32)
    offsets = np.cumsum([0] + [len(a) for a in list(contig_depths.values())[:-1]])

    np.savez_compressed(
        output_path,
        coverage=coverage,
        contig_names=np.array(contig_names, dtype=str),
        contig_offsets=offsets,
    )
