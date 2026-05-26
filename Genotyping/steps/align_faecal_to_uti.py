#!/usr/bin/env python3
"""Step: Align faecal metagenomic reads to UTI assemblies."""

import threading

import pandas as pd
from pathlib import Path

from config_utils import get_output_dir
from steps.per_run_step import PerRunStep
from steps.step import run_cli
from steps.assemble import AssembleStep
from steps.create_pairs import CreatePairsStep
from steps.extract_run_names import ExtractRunNamesStep
from utils.metagenome_io import get_metagenome_reads
from tools.bowtie2 import bowtie2_build, bowtie2_align_paired, sam_to_sorted_bam, bedtools_genomecov
from utils.file_utils import create_if_missing


def align_paired_and_genomecov(
    ref_fasta: Path,
    r1: Path,
    r2: Path,
    output_dir: Path,
    index_prefix: Path,
    threads: int = 6,
    keep_bam: bool = False,
) -> Path:
    """Full pipeline: align → sort → genomecov → cleanup.

    Returns path to coverage.npz output.
    """
    create_if_missing(output_dir)
    sam_path = output_dir / "aligned.sam"
    sorted_bam = output_dir / "aligned.sorted.bam"
    genomecov_path = output_dir / "coverage.npz"
    metrics_path = output_dir / "bowtie2.metrics.txt"
    log_path = output_dir / "bowtie2.log"

    result = bowtie2_align_paired(
        index_prefix, r1, r2, sam_path,
        threads=threads, met_file=metrics_path
    )
    log_path.write_text(result.stderr)
    sam_to_sorted_bam(sam_path, sorted_bam, threads=threads)
    bedtools_genomecov(sorted_bam, genomecov_path)

    # Cleanup intermediates
    sam_path.unlink(missing_ok=True)
    if not keep_bam:
        sorted_bam.unlink(missing_ok=True)
        bai = Path(str(sorted_bam) + ".bai")
        bai.unlink(missing_ok=True)
    return genomecov_path


class AlignFaecalToUTIStep(PerRunStep):
    """Align faecal reads to UTI assemblies for each pair."""

    output_dir_name = "alignment"
    cores_per_run = 4
    output_paths = {"coverage": "coverage.npz", "metrics": "bowtie2.metrics.txt", "log": "bowtie2.log"}
    sort_by = ["match_num", "uti_run_name"]
    ascending = [True, True]

    def __init__(self, csv_file=None, **kwargs):
        super().__init__(csv_file=csv_file or CreatePairsStep.output_path(), **kwargs)
        self.threads = self.threads_per_worker
        self._runs_df = None
        self._index_lock = threading.Lock()

    def _should_process(self, row: pd.Series) -> bool:
        """Use uti/faecal run names for filtering instead of compound pair name."""
        uti_row = self._get_run_row(row["uti_run_name"])
        faecal_row = self._get_run_row(row["faecal_run_name"])
        return super()._should_process(uti_row) and super()._should_process(faecal_row)

    def setup(self) -> None:
        alignment_dir = get_output_dir("alignment")
        create_if_missing(alignment_dir)
        create_if_missing(alignment_dir / "indices")
        self._runs_df = pd.read_csv(ExtractRunNamesStep.output_path())

    def _get_run_row(self, run_name: str) -> pd.Series:
        return self._runs_df[self._runs_df["run_name"] == run_name].iloc[0]

    def _create_index_if_missing(self, uti_row: pd.Series) -> Path:
        """Build bowtie2 index for a UTI assembly (thread-safe, once per UTI)."""
        uti_run = uti_row["run_name"]
        index_prefix = get_output_dir("alignment") / "indices" / uti_run
        with self._index_lock:
            if Path(f"{index_prefix}.1.bt2").exists():
                return index_prefix
            assembly_fasta = AssembleStep.run_output_path(uti_row)["assembly"]
            if not assembly_fasta.exists():
                raise FileNotFoundError(
                    f"Assembly not found for {uti_run}: {assembly_fasta}"
                )
            self.print(f"Building index for {uti_run}")
            bowtie2_build(assembly_fasta, index_prefix, threads=self.threads)
        return index_prefix

    def run_input_paths(self, row: pd.Series) -> list:
        uti_row = self._get_run_row(row["uti_run_name"])
        faecal_row = self._get_run_row(row["faecal_run_name"])
        assembly = AssembleStep.run_output_path(uti_row)["assembly"]
        r1, r2 = get_metagenome_reads(faecal_row)
        return [assembly, r1, r2]

    def process_run(self, row: pd.Series) -> None:

        # UTI:
        uti_run = row["uti_run_name"]
        uti_row = self._get_run_row(uti_run)
        index_prefix = self._create_index_if_missing(uti_row)
        assembly_fasta = AssembleStep.run_output_path(uti_row)["assembly"]

        # Faecal:
        faecal_run = row["faecal_run_name"]
        faecal_row = self._get_run_row(faecal_run)
        r1, r2 = get_metagenome_reads(faecal_row)

        pair_dir = self.run_folder(row)
        align_paired_and_genomecov(
            ref_fasta=assembly_fasta,
            r1=r1,
            r2=r2,
            output_dir=pair_dir,
            index_prefix=index_prefix,
            threads=self.threads,
        )


if __name__ == "__main__":
    run_cli(AlignFaecalToUTIStep)
