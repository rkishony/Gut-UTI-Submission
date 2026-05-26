#!/usr/bin/env python3
"""Step: Collect alignment statistics from bowtie2 log files."""

import re

import numpy as np
import pandas as pd
from pathlib import Path

from config_utils import get_output_dir
from steps.step import Step, run_cli
from steps.create_pairs import CreatePairsStep
from steps.align_faecal_to_uti import AlignFaecalToUTIStep


FILENAME = "alignment_stats.csv"
COALESCED_FILENAME = "alignment_stats_coalesced.csv"

SUM_FIELDS = ["reads", "conc_1", "conc_gt1", "conc_0"]

STAT_COLS = [
    "reads",
    "conc_1", "conc_1_pct",
    "conc_gt1", "conc_gt1_pct",
    "conc_0", "conc_0_pct",
    "genome_len", "uncov_len",
]

_LOG_PATTERNS = [
    ("reads",      r"(\d+) reads; of these:",                                     int,   0),
    ("align_rate", r"([\d.]+)% overall alignment rate",                            float, 0.0),
    ("conc_1",     r"(\d+) \(([\d.]+)%\) aligned concordantly exactly 1 time",     int,   0),
    ("conc_1_pct", r"(\d+) \(([\d.]+)%\) aligned concordantly exactly 1 time",     float, 0.0),
    ("conc_gt1",   r"(\d+) \(([\d.]+)%\) aligned concordantly >1 times",           int,   0),
    ("conc_gt1_pct", r"(\d+) \(([\d.]+)%\) aligned concordantly >1 times",         float, 0.0),
    ("conc_0",     r"(\d+) \(([\d.]+)%\) aligned concordantly 0 times",            int,   0),
    ("conc_0_pct", r"(\d+) \(([\d.]+)%\) aligned concordantly 0 times",            float, 0.0),
]


def _parse_log(log_text: str) -> dict:
    """Parse a bowtie2 log file into a stats dict."""
    stats = {}
    for name, pattern, typ, default in _LOG_PATTERNS:
        m = re.search(pattern, log_text)
        group_idx = 2 if name.endswith("_pct") and m and m.lastindex >= 2 else 1
        stats[name] = typ(m.group(group_idx)) if m else default
    return stats


class AlignmentStatsStep(Step):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    @classmethod
    def output_path(cls) -> Path:
        return get_output_dir("summary") / FILENAME

    @classmethod
    def coalesced_output_path(cls) -> Path:
        return get_output_dir("summary") / COALESCED_FILENAME

    def _run(self) -> None:
        pairs_df = pd.read_csv(CreatePairsStep.output_path())
        rows = []
        coalesced_rows = []
        prev_key = None
        prev_cov = None

        for _, row in pairs_df.iterrows():
            align = AlignFaecalToUTIStep.run_output_path(row)

            stats = _parse_log(align["log"].read_text() if align["log"].exists() else "")
            stats["status"] = "ok" if align["log"].exists() else "no_log"

            npz_path = align["coverage"]
            cov = None
            if npz_path.exists():
                cov = np.load(npz_path)["coverage"]
                stats["genome_len"] = len(cov)
                stats["uncov_len"] = int(np.sum(cov == 0))
            else:
                stats["genome_len"] = stats["uncov_len"] = 0

            stats["run_name"] = row["run_name"]
            stats["uti_run_name"] = row["uti_run_name"]
            stats["faecal_run_name"] = row["faecal_run_name"]
            rows.append(stats)

            # --- coalesced: accumulate bio_reps sharing (uti, pair_type, match_num) ---
            uncov_rep = stats["uncov_len"] if cov is not None else np.nan
            key = (row["uti_run_name"], row["pair_type"], row["match_num"])
            if key == prev_key:
                c = coalesced_rows[-1]
                c["n_bio_reps"] += 1
                c["faecal_run_names"] += f",{row['faecal_run_name']}"
                for f in SUM_FIELDS:
                    c[f] += stats[f]
                c["uncov_len2"] = uncov_rep
                prev_cov = (prev_cov + cov if prev_cov is not None else cov) if cov is not None else prev_cov
                c["uncov_len"] = int(np.sum(prev_cov == 0)) if prev_cov is not None else np.nan
                c["status"] = "ok" if c["status"] == "ok" and stats["status"] == "ok" else "partial"
            else:
                coalesced_rows.append({
                    "uti_run_name": row["uti_run_name"],
                    "pair_type": row["pair_type"],
                    "match_num": row["match_num"],
                    "faecal_run_names": row["faecal_run_name"],
                    "n_bio_reps": 1,
                    "status": stats["status"],
                    **{f: stats[f] for f in SUM_FIELDS},
                    "genome_len": stats["genome_len"],
                    "uncov_len": uncov_rep,
                    "uncov_len1": uncov_rep,
                    "uncov_len2": np.nan,
                })
                prev_cov = cov
            prev_key = key

        # --- per-run table ---
        out = pd.DataFrame(rows)
        out = out.merge(pairs_df[["run_name", "pair_type", "match_num", "bio_rep"]], on="run_name", how="left")
        out = out[["run_name", "uti_run_name", "faecal_run_name", "pair_type", "match_num", "bio_rep",
                    "status", "align_rate"] + STAT_COLS]
        out.to_csv(self.output_path(), index=False)

        # --- coalesced table ---
        coal = pd.DataFrame(coalesced_rows)
        for fld in ["conc_1", "conc_gt1", "conc_0"]:
            coal[f"{fld}_pct"] = np.where(coal["reads"] > 0, coal[fld] / coal["reads"] * 100, 0.0).round(2)
        coal = coal[["uti_run_name", "faecal_run_names", "pair_type", "match_num", "n_bio_reps", "status"]
                     + STAT_COLS + ["uncov_len1", "uncov_len2"]]
        coal.to_csv(self.coalesced_output_path(), index=False)

        return f"{len(out)} per-run, {len(coal)} coalesced"


if __name__ == "__main__":
    run_cli(AlignmentStatsStep)
