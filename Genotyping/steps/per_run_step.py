"""PerRunStep: MultiStep defaulting to the trim-stats CSV."""

import pandas as pd
from pathlib import Path

from config_utils import get_output_dir
from steps.step import MultiStep


class PerRunStep(MultiStep):
    output_dir_name: str  # subclasses must set this
    sort_by: str | list[str] | None = "num_reads"
    ascending: bool | list[bool] | None = False

    def __init__(self, csv_file: str | Path | None = None, **kwargs):
        if csv_file is None:
            from steps.trim_stats import TrimStatsStep
            csv_file = TrimStatsStep.output_path()[0]
        super().__init__(csv_file=csv_file, **kwargs)

    @classmethod
    def run_folder(cls, row: pd.Series) -> Path:
        return get_output_dir(cls.output_dir_name) / row["run_name"]
