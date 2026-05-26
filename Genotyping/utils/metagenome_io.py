"""Helper to resolve metagenome read paths (downsampled or trimmed)."""

from pathlib import Path
from typing import Tuple

import pandas as pd

from steps.downsample import DownSampleStep
from steps.trim_sequences import TrimSequencesStep


def get_metagenome_reads(row: pd.Series) -> Tuple[Path, Path]:
    """Return (R1, R2): downsampled if available, else trimmed."""
    down = DownSampleStep.run_output_path(row)
    if down["R1"].exists() and down["R2"].exists():
        return down["R1"], down["R2"]
    trim = TrimSequencesStep.run_output_path(row)
    return trim["R1"], trim["R2"]
