from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
COLORS_TXT   = PROJECT_ROOT / "colors.txt"
META_CSV     = PROJECT_ROOT / "output" / "summary" / "trim_stats.csv"
MASH_CSV     = PROJECT_ROOT / "output" / "mash" / "mash_distance_matrix.csv"
ALIGN_CSV    = PROJECT_ROOT / "output" / "summary" / "alignment_stats_coalesced.csv"
UNCOV_THRESH = 1e-4


COLORS: dict[str, str] = {}
for _line in COLORS_TXT.read_text().splitlines():
    _line = _line.strip()
    if _line:
        _name, _hex = _line.split(",", 1)
        COLORS[_name.strip()] = _hex.strip()


def load_color(name: str) -> str:
    """Return '#RRGGBB' hex color by name from colors.txt."""
    return COLORS[name]


def get_patient_meta(fid: int) -> pd.DataFrame:
    meta = pd.read_csv(META_CSV)
    meta = meta[(meta["fid"] == fid) & (meta["is_isolate"] == True)].copy()
    return meta


def load_mash_distances() -> pd.DataFrame:
    return pd.read_csv(MASH_CSV, index_col=0)


def load_uti_coverage() -> pd.Series:
    """Return minimum uncov_frac per UTI run (matched pairs only)."""
    align = pd.read_csv(ALIGN_CSV)
    matched = align[align["pair_type"] == "matched"].copy()
    matched["uncov_frac"] = matched["uncov_len"] / matched["genome_len"].replace(0, np.nan)
    return matched.groupby("uti_run_name")["uncov_frac"].min()


def classify_uti_coverage(run_name: str, uncov_by_uti: pd.Series) -> bool | None:
    """True = covered, False = uncovered, None = no alignment data."""
    if run_name not in uncov_by_uti.index:
        return None
    frac = uncov_by_uti[run_name]
    if np.isnan(frac):
        return None
    return frac < UNCOV_THRESH


def get_all_isolate_meta() -> pd.DataFrame:
    meta = pd.read_csv(META_CSV)
    return meta[meta["is_isolate"] == True].copy()
