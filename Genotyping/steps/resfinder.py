"""Step: Run ResFinder + PointFinder on metagenome (non-isolate) samples."""

import shutil
import pandas as pd
from pathlib import Path

from steps.per_run_step import PerRunStep
from utils.metagenome_io import get_metagenome_reads
from steps.step import run_cli
from tools.resfinder import run_resfinder, run_pointfinder


class ResfinderStep(PerRunStep):
    output_dir_name = "resfinder"
    output_paths = ["acq", "pointEcoli"]

    def run_input_paths(self, row: pd.Series) -> tuple:
        return get_metagenome_reads(row)

    def filter_df(self, df: pd.DataFrame) -> pd.DataFrame:
        return df[~df["is_isolate"]]

    def process_run(self, row: pd.Series) -> None:
        r1, r2 = self.run_input_paths(row)
        acq_dir, point_dir = self.run_output_path(row)

        run_resfinder(r1, r2, acq_dir)
        run_pointfinder(r1, r2, point_dir)

        for kma_dir in (acq_dir / "resfinder_kma", point_dir / "pointfinder_kma"):
            if kma_dir.is_dir():
                for f in kma_dir.iterdir():
                    if f.is_file() and f.suffix != ".res":
                        f.unlink(missing_ok=True)


if __name__ == "__main__":
    run_cli(ResfinderStep)
