#!/usr/bin/env python3
"""Step: Collect trimming statistics from cutadapt JSON reports."""

import pandas as pd
from pathlib import Path

import config
from config_utils import get_output_dir
from steps.downsample import DownSampleStep
from steps.extract_run_names import ExtractRunNamesStep
from steps.step import Step, run_cli
from steps.trim_sequences import TrimSequencesStep
from tools.cutadapt import extract_trim_stats, load_cutadapt_report


FILENAME = "trim_stats.csv"
SORTED_FILENAME = "trim_stats_sorted.csv"

def _get_sample_name(row: pd.Series) -> str:
    fid = row['fid']
    if row['is_isolate'] and not row['is_faecal']:
        return f"P{fid}.U{int(row['uid'])}.{row['isolate']}"
    elif row['is_isolate'] and row['is_faecal']:
        return f"P{fid}.F.{row['isolate']}"
    else:
        bio = row['bio_rep']
        suffix = chr(ord('a') + bio - 2) if bio > 1 else ''
        return f"P{fid}.F{suffix}"


class TrimStatsStep(Step):
    COLUMNS_FROM_RUN_INFO = ["run_name", "fid", "uid", "is_isolate", "is_faecal", "isolate", "replicate", "batch", "CFU", "coverage_factor"]
    def __init__(self, csv_file: str | Path = None, **kwargs):
        super().__init__(**kwargs)
        self.csv_file = csv_file or ExtractRunNamesStep.output_path()

    @classmethod
    def output_path(cls) -> list[Path]:
        summary_dir = get_output_dir("summary")
        return [summary_dir / FILENAME, summary_dir / SORTED_FILENAME]

    def _run(self) -> None:
        df = pd.read_csv(self.csv_file)
        rows = []
        for _, row in df.iterrows():
            trim = TrimSequencesStep.run_output_path(row)
            report = load_cutadapt_report(trim["report"])

            has_two_fasta = trim["R1"].exists() and trim["R2"].exists()
            stats = extract_trim_stats(report)
            stats["run_name"] = row["run_name"]
            stats["has_two_fasta"] = has_two_fasta

            # Downsampling info
            down = DownSampleStep.run_output_path(row)
            downsampled = down["R1"].exists() and down["R2"].exists()
            stats["downsampled"] = downsampled
            if downsampled and config.BASE_READS_PER_COVERAGE_TARGET is not None:
                cov = row.get("coverage_factor")
                target = int(config.BASE_READS_PER_COVERAGE_TARGET * float(cov)) if pd.notna(cov) else pd.NA
                stats["downsampled_num_reads"] = target
            else:
                stats["downsampled_num_reads"] = pd.NA
            stats["effective_num_reads"] = (
                stats["downsampled_num_reads"] if downsampled else stats["num_reads"]
            )
            eff = stats["effective_num_reads"]
            stats["effective_total_bases"] = (
                int(eff * stats["read_len"] * 2) if pd.notna(eff) and eff else 0
            )

            for col in self.COLUMNS_FROM_RUN_INFO:
                stats[col] = row[col]
            rows.append(stats)

        out = pd.DataFrame(rows)

        # Mark replicate rows (same fid x uid)
        num_reads = out["effective_num_reads"]
        out["is_rep"] = True
        for group_idx in out.groupby(["fid", "uid", "isolate"], dropna=False).groups.values():
            group_num_reads = num_reads.loc[group_idx]
            rep_idx = group_num_reads.idxmax()
            out.loc[rep_idx, "is_rep"] = False

        # Mark biological replicates: hs vs ns are different platings (bio reps),
        # hs1 vs hs2 are same plating (tech reps, not bio reps).
        out["_plating"] = out["batch"].str.replace(r"\d+$", "", regex=True)
        bio_group_cols = [c for c in self.COLUMNS_FROM_RUN_INFO if c not in ("batch", "run_name", "replicate")]
        out["bio_rep"] = 0
        for _, grp in out.groupby(bio_group_cols, dropna=False):
            if grp["_plating"].nunique() > 1:
                reps = grp[~grp["is_rep"]]
                assert len(reps) == 1, f"Expected exactly 1 representative, got {len(reps)}: {grp['run_name'].tolist()}"
                rep_plating = reps["_plating"].iloc[0]
                platings = [rep_plating] + [p for p in grp["_plating"].unique() if p != rep_plating]
                plating_to_bio = {p: i for i, p in enumerate(platings, start=1)}
                for plating, sub in grp.groupby("_plating"):
                    best = sub["effective_num_reads"].idxmax()
                    out.loc[best, "bio_rep"] = plating_to_bio[plating]
        out.drop(columns=["_plating"], inplace=True)

        out['sample_name'] = out.apply(_get_sample_name, axis=1)

        # reorder columns
        cols = self.COLUMNS_FROM_RUN_INFO + [
            "sample_name", "is_rep", "bio_rep",
            "has_two_fasta", "raw_num_reads", "num_reads", "percentage_reads_trimmed",
            "raw_read_len", "read_len", "total_bases", "percentage_with_adapter",
            "downsampled", "downsampled_num_reads", "effective_num_reads", "effective_total_bases"]
        out = out[cols]
        out["uid"] = out["uid"].astype("Int64")
        out.to_csv(self.output_path()[0], index=False)

        # sort by raw_num_reads descending
        out = out.sort_values(["is_faecal", "is_isolate", "raw_num_reads"], ascending=[True, True, False])
        out.to_csv(self.output_path()[1], index=False)
        return f"{len(out)} runs"


if __name__ == "__main__":
    run_cli(TrimStatsStep)
