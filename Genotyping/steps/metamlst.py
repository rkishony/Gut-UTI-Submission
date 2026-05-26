"""Step: Run MetaMLST typing on metagenome (non-isolate) samples."""

import pandas as pd
from pathlib import Path

from config_utils import get_output_dir
from steps.per_run_step import PerRunStep
from utils.metagenome_io import get_metagenome_reads
from steps.step import Step, run_cli
from tools.metamlst import metamlst_align, metamlst_type, metamlst_merge


class MetaMLSTStep(PerRunStep):
    output_dir_name = "metamlst"
    output_paths = "typing"

    def run_input_paths(self, row: pd.Series) -> tuple:
        return get_metagenome_reads(row)

    def filter_df(self, df: pd.DataFrame) -> pd.DataFrame:
        return df[~df["is_isolate"]]

    def process_run(self, row: pd.Series) -> None:
        r1, r2 = self.run_input_paths(row)
        run_dir = self.run_folder(row)
        bam = run_dir / "aligned.bam"
        typing_dir = self.run_output_path(row)

        metamlst_align(r1, r2, bam)
        metamlst_type(bam, typing_dir)
        bam.unlink(missing_ok=True)
        bam.with_suffix(".bam.bai").unlink(missing_ok=True)


class MetaMLSTMergeStep(Step):
    """Merge all per-sample MetaMLST typing results."""
    allow_failure = True

    @classmethod
    def output_path(cls) -> Path:
        return get_output_dir("metamlst") / "merged"

    def _run(self) -> None:
        metamlst_dir = get_output_dir("metamlst")

        typing_dirs = sorted(d for d in metamlst_dir.glob("*/typing") if d.is_dir())
        if not typing_dirs:
            raise RuntimeError("No per-sample typing directories found")

        metamlst_merge(metamlst_dir)
        return f"{len(typing_dirs)} merged samples"


if __name__ == "__main__":
    run_cli(MetaMLSTStep, MetaMLSTMergeStep)
