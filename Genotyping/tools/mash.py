"""Mash kmer sketch wrapper."""

from pathlib import Path

from tools.run_utils import get_tool, run_cmd


def mash_sketch(
    input_file: Path,
    output_prefix: Path,
    sketch_size: int = 10000,
) -> Path:
    """Run `mash sketch` on a FASTA/FASTQ file.

    Args:
        input_file: Input FASTA or FASTQ (may be gzipped).
        output_prefix: Output prefix (mash appends '.msh').
        sketch_size: Sketch size (-s).

    Returns:
        Path to the generated .msh sketch file.
    """
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        get_tool("mash"), "sketch",
        "-s", str(sketch_size),
        "-o", str(output_prefix),
        str(input_file),
    ]

    run_cmd(cmd)

    return output_prefix.with_suffix(".msh")


def mash_paste(output_prefix: Path, sketch_files: list[Path]) -> Path:
    """Combine multiple .msh sketches into one."""
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    cmd = [get_tool("mash"), "paste", str(output_prefix)] + [str(f) for f in sketch_files]
    run_cmd(cmd)
    return output_prefix.with_suffix(".msh")


def mash_dist(
    reference: Path,
    query: Path,
    output_file: Path,
    name_map: dict[str, str] | None = None,
) -> Path:
    """Run `mash dist` and save tab-delimited output to *output_file*.

    Args:
        name_map: If provided, replaces ref/query paths with mapped names.
    """
    output_file.parent.mkdir(parents=True, exist_ok=True)
    cmd = [get_tool("mash"), "dist", str(reference), str(query)]
    result = run_cmd(cmd)
    if name_map:
        lines = []
        for line in result.stdout.strip().splitlines():
            ref, qry, *rest = line.split("\t")
            ref = name_map.get(ref, ref)
            qry = name_map.get(qry, qry)
            lines.append("\t".join([ref, qry] + rest))
        output_file.write_text("\n".join(lines) + "\n")
    else:
        output_file.write_text(result.stdout)
    return output_file
