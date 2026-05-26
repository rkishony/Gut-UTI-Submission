#!/usr/bin/env python3
"""Step: Downsample metagenomes if read count exceeds target based on coverage_factor."""

import pandas as pd
from pathlib import Path

from config import BASE_READS_PER_COVERAGE_THRESHOLD, BASE_READS_PER_COVERAGE_TARGET
from steps.extract_run_names import ExtractRunNamesStep
from tools.cutadapt import load_cutadapt_report
from steps.per_run_step import PerRunStep
from steps.step import run_cli
from steps.trim_sequences import TrimSequencesStep
from tools.seqtk import downsample_fastq


class DownSampleStep(PerRunStep):
    output_dir_name = "downsampled"
    sort_by = None

    def __init__(self, csv_file=None, **kwargs):
        super().__init__(csv_file=csv_file or ExtractRunNamesStep.output_path(), **kwargs)

    output_paths = {
        "R1": "R1.fastq.gz",
        "R2": "R2.fastq.gz",
        "not_downsampled": "not_downsampled.txt",
    }

    def run_input_paths(self, row: pd.Series) -> dict:
        return TrimSequencesStep.run_output_path(row)

    def filter_df(self, df: pd.DataFrame) -> pd.DataFrame:
        return df[(~df["is_isolate"]) & (df["is_faecal"])]

    def _run_output_paths_satisfied(self, row: pd.Series) -> bool:
        out = self.run_output_path(row)
        return (out["R1"].exists() and out["R2"].exists()) or out["not_downsampled"].exists()

    def _write_not_downsampled(self, out: dict, msg: str) -> None:
        out["not_downsampled"].parent.mkdir(parents=True, exist_ok=True)
        out["not_downsampled"].write_text(msg + "\n")

    def process_run(self, row: pd.Series) -> None:
        trim = self.run_input_paths(row)
        out = self.run_output_path(row)
        report_path = trim["report"]

        report = load_cutadapt_report(report_path)
        assert report is not None
        num_reads = report["read_counts"]["output"]

        cov = row.get("coverage_factor")
        if pd.isna(cov) or not isinstance(cov, (int, float)) or cov <= 0:
            self._write_not_downsampled(out, "no coverage_factor")
            return

        if BASE_READS_PER_COVERAGE_THRESHOLD is None or BASE_READS_PER_COVERAGE_TARGET is None:
            self._write_not_downsampled(out, "downsampling disabled (threshold/target is None)")
            return

        threshold_reads = int(BASE_READS_PER_COVERAGE_THRESHOLD * float(cov))
        if num_reads <= threshold_reads:
            self._write_not_downsampled(out, f"num_reads ({num_reads:,}) <= threshold ({threshold_reads:,})")
            return

        target_reads = int(BASE_READS_PER_COVERAGE_TARGET * float(cov))
        self.print(f"downsampling {num_reads:,} -> {target_reads:,} read pairs (threshold={threshold_reads:,})", row)
        out["R1"].parent.mkdir(parents=True, exist_ok=True)
        downsample_fastq(trim["R1"], out["R1"], target_reads)
        downsample_fastq(trim["R2"], out["R2"], target_reads)


if __name__ == "__main__":
    run_cli(DownSampleStep)
