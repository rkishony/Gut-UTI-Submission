"""Step: Run MetaPhlAn profiling on metagenome (non-isolate) samples."""

import pandas as pd
from pathlib import Path

from steps.per_run_step import PerRunStep
from utils.metagenome_io import get_metagenome_reads
from steps.step import run_cli
from tools.metaphlan import metaphlan_profile


class MetaphlanStep(PerRunStep):
    cores_per_run = 4

    output_paths = "profile.tsv"

    def run_input_paths(self, row: pd.Series) -> tuple:
        return get_metagenome_reads(row)

    def filter_df(self, df: pd.DataFrame) -> pd.DataFrame:
        return df[~df["is_isolate"]]

    def process_run(self, row: pd.Series) -> None:
        r1, r2 = self.run_input_paths(row)
        run_dir = self.run_folder(row)
        profile = self.run_output_path(row)
        sam = run_dir / "sam.bz2"
        bowtie2 = run_dir / "bowtie2.bz2"

        metaphlan_profile(r1, r2, profile, sam, bowtie2,
                          nproc=self.threads_per_worker)
        sam.unlink()
        bowtie2.unlink()


if __name__ == "__main__":
    run_cli(MetaphlanStep)
