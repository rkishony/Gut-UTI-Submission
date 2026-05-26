"""Step: Run Mash sketch on assembled genomes and compute distance matrix."""

import sys

import numpy as np
import pandas as pd
from pathlib import Path

from config_utils import get_output_dir
from steps.assemble import AssembleStep
from steps.extract_run_names import ExtractRunNamesStep
from steps.per_run_step import PerRunStep
from steps.step import Step, run_cli
from tools.mash import mash_sketch, mash_paste, mash_dist


class MashStep(PerRunStep):
    output_dir_name = "mash"
    output_paths = "sketch.msh"

    def run_input_paths(self, row: pd.Series) -> Path:
        return AssembleStep.run_output_path(row)["assembly"]

    def filter_df(self, df: pd.DataFrame) -> pd.DataFrame:
        return df[df["is_isolate"]]

    def process_run(self, row: pd.Series) -> None:
        assembly_fasta = self.run_input_paths(row)
        output_prefix = self.run_folder(row) / "sketch"
        mash_sketch(assembly_fasta, output_prefix)


class MashDistStep(Step):
    """Collect all per-run Mash sketches and compute a pairwise distance matrix."""

    def __init__(self, csv_file: str | Path = None, **kwargs):
        super().__init__(**kwargs)
        self.csv_file = csv_file or ExtractRunNamesStep.output_path()

    @classmethod
    def output_path(cls) -> Path:
        return get_output_dir("mash") / "mash_distance_matrix.csv"

    @classmethod
    def summary_path(cls) -> Path:
        return get_output_dir("summary") / "mash_stats.csv"

    def _run(self) -> None:
        df = pd.read_csv(self.csv_file)
        df = df[df["is_isolate"]]

        sketches = []
        run_names = []
        name_map = {}
        skipped = []
        for _, row in df.iterrows():
            msh = MashStep.run_output_path(row)
            if not msh.exists():
                skipped.append(row['run_name'])
                continue
            sketches.append(msh)
            run_names.append(row["run_name"])
            name_map[str(AssembleStep.run_output_path(row)["assembly"])] = row["run_name"]
        if skipped:
            print(f"[MashDist] Skipping {len(skipped)} runs with missing sketches: {skipped}")

        combined = get_output_dir("mash") / "all_sketches"
        combined.with_suffix(".msh").unlink(missing_ok=True)
        mash_paste(combined, sketches)
        combined_msh = combined.with_suffix(".msh")

        raw_path = get_output_dir("mash") / "mash_dist_raw.tsv"
        mash_dist(combined_msh, combined_msh, raw_path, name_map=name_map)

        name_to_idx = {n: i for i, n in enumerate(run_names)}
        n = len(run_names)
        dist_matrix = np.zeros((n, n))
        for line in raw_path.read_text().strip().splitlines():
            ref, query, distance, p_value, shared = line.split("\t")
            if ref in name_to_idx and query in name_to_idx:
                i, j = name_to_idx[ref], name_to_idx[query]
                dist_matrix[i, j] = float(distance)
                dist_matrix[j, i] = float(distance)

        dist_df = pd.DataFrame(dist_matrix, index=run_names, columns=run_names)
        dist_df.index.name = "run_name"
        dist_df.to_csv(self.output_path())

        self._write_summary(df, dist_matrix, run_names, skipped)
        return f"{n}x{n} matrix"

    def _write_summary(self, df: pd.DataFrame, dist_matrix: np.ndarray,
                       run_names: list[str], skipped: list[str]) -> None:
        included = set(run_names)
        rows = [{"run_name": row["run_name"], "has_sketch": row["run_name"] in included}
                for _, row in df.iterrows()]
        pd.DataFrame(rows).to_csv(self.summary_path(), index=False)


if __name__ == "__main__":
    run_cli(MashStep)
