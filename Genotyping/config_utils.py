"""Helpers that read/mutate config values at runtime."""
from pathlib import Path

import config


ALL_FOLDERS = (
    "", "trimmed", "downsampled", "assembly", "alignment", "summary",
    "metaphlan", "mash", "metamlst", "resfinder", "mlst", "strainphlan",
    "sra_upload", "logs",
)


def get_output_dir(name: str) -> Path:
    """Return OUTPUT_DIR / name, always reflecting current OUTPUT_DIR."""
    assert name in ALL_FOLDERS, f"Unknown output directory: {name}"
    return config.OUTPUT_DIR / name
