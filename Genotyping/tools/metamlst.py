"""MetaMLST wrapper: align reads to MLST loci and run MetaMLST typing.

Reproduces the original bash workflow:
  bowtie2 --very-sensitive-local -a --no-unal -x index -U R1.fq.gz R2.fq.gz \
      | samtools view -bS - > sample.bam
  metamlst.py sample.bam
  metamlst-merge.py ./out
"""

import os
from pathlib import Path

from config import METAMLST_BT2_INDEX
from tools.run_utils import get_tool, run_cmd


def _metamlst_env():
    """Return env dict with metamlst-samtools on PATH."""
    samtools_dir = str(get_tool("metamlst-samtools").parent)
    env = os.environ.copy()
    env["PATH"] = samtools_dir + ":" + env.get("PATH", "")
    return env


def metamlst_align(
    r1: Path,
    r2: Path,
    output_bam: Path,
) -> Path:
    """Align reads to MetaMLST bowtie2 index → BAM.

    Uses -U (unpaired) with both R1 and R2 files, matching the original workflow.
    """
    output_bam.parent.mkdir(parents=True, exist_ok=True)

    bowtie2 = get_tool("metamlst-bowtie2")
    samtools = get_tool("metamlst-samtools")

    cmd = (
        f"{bowtie2} --very-sensitive-local -a --no-unal "
        f"-x {METAMLST_BT2_INDEX} "
        f"-U {r1} {r2} "
        f"| {samtools} view -bS - > {output_bam}"
    )
    run_cmd(cmd, shell=True)

    return output_bam


def metamlst_type(
    bam_file: Path,
    output_dir: Path,
) -> Path:
    """Run metamlst.py on a BAM to reconstruct MLST profiles."""
    output_dir.mkdir(parents=True, exist_ok=True)
    run_cmd(
        [get_tool("metamlst-python"), get_tool("metamlst"),
         str(bam_file),
         "-o", str(output_dir),
         "--log"],
        env=_metamlst_env(),
    )
    return output_dir


def metamlst_merge(output_dir: Path) -> Path:
    """Run metamlst-merge.py to combine all per-sample results.

    metamlst-merge.py expects .nfo files directly inside the given folder,
    but our pipeline stores them under <run>/typing/.  We symlink them into a
    flat staging directory so the merge script can find them.
    """
    staging = output_dir / "_merge_staging"
    staging.mkdir(parents=True, exist_ok=True)

    for nfo in output_dir.glob("*/typing/*.nfo"):
        sample_name = nfo.parents[1].name
        dest = staging / f"{sample_name}.nfo"
        text = nfo.read_text()
        text = text.replace("\taligned\t", f"\t{sample_name}\t")
        dest.write_text(text)

    run_cmd([get_tool("metamlst-python"), get_tool("metamlst-merge"), str(staging)],
            env=_metamlst_env())

    merged_src = staging / "merged"
    merged_dst = output_dir / "merged"
    if merged_dst.exists():
        import shutil
        shutil.rmtree(merged_dst)
    if merged_src.exists():
        merged_src.rename(merged_dst)

    import shutil
    shutil.rmtree(staging, ignore_errors=True)

    return output_dir
