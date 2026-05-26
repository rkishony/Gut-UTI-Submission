#!/usr/bin/env python3
"""Step 2: Trim sequencing reads using cutadapt."""

import pandas as pd

from get_raw_files import get_raw_files
from steps.extract_run_names import ExtractRunNamesStep
from steps.per_run_step import PerRunStep
from steps.step import run_cli
from tools.cutadapt import clip_and_trim_paired


class TrimSequencesStep(PerRunStep):
    output_dir_name = "trimmed"
    sort_by = None

    def __init__(self, csv_file=None, **kwargs):
        super().__init__(csv_file=csv_file or ExtractRunNamesStep.output_path(), **kwargs)
    output_paths = {"R1": "R1.fastq.gz", "R2": "R2.fastq.gz", "report": "cutadapt_report.json"}

    def run_input_paths(self, row: pd.Series) -> tuple:
        return get_raw_files(row["raw_path"], row["batch"], raise_error=False)

    def process_run(self, row: pd.Series) -> None:
        r1_file, r2_file = self.run_input_paths(row)
        out = self.run_output_path(row)
        clip_and_trim_paired(r1_file, r2_file, out["R1"], out["R2"], report_path=out["report"])


if __name__ == "__main__":
    run_cli(TrimSequencesStep)
