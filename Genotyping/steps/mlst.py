"""Steps: Run MLST on isolate assemblies + collect summary table."""

import pandas as pd
from pathlib import Path

from config_utils import get_output_dir
from steps.assemble import AssembleStep
from steps.per_run_step import PerRunStep
from steps.extract_run_names import ExtractRunNamesStep
from steps.step import Step, run_cli
from tools.mlst import run_mlst


class MLSTStep(PerRunStep):
    output_dir_name = "mlst"
    cores_per_run = 2

    output_paths = "mlst.tsv"

    def run_input_paths(self, row: pd.Series) -> Path:
        return AssembleStep.run_output_path(row)["assembly"]

    def filter_df(self, df: pd.DataFrame) -> pd.DataFrame:
        return df[df["is_isolate"]]

    def process_run(self, row: pd.Series) -> None:
        assembly = self.run_input_paths(row)
        output_tsv = self.run_output_path(row)
        run_mlst(assembly, output_tsv, threads=self.threads_per_worker)


class CollectMLSTStep(Step):
    """Collect all per-isolate MLST outputs into one summary table."""
    FIELDS_TO_COPY = ["fid", "uid", "is_isolate", "is_faecal", "isolate", "replicate"]

    def __init__(self, csv_file: str | Path | None = None, **kwargs):
        super().__init__(**kwargs)
        self.csv_file = csv_file or ExtractRunNamesStep.output_path()

    @classmethod
    def output_path(cls) -> Path:
        return get_output_dir("summary") / "mlst_isolates.csv"

    def _run(self) -> str:
        df = pd.read_csv(self.csv_file)
        df = df[df["is_isolate"]]

        rows = []
        for _, row in df.iterrows():
            run_name = row["run_name"]
            mlst_path = MLSTStep.run_output_path(row)

            summary = {field: row[field] for field in self.FIELDS_TO_COPY}
            summary.update({
                "run_name": run_name,
                "species": pd.NA,
                "ST": pd.NA,
                "status": "missing_file",
            })

            if mlst_path.exists():
                text = mlst_path.read_text().strip()
                if not text:
                    summary["status"] = "empty"
                else:
                    parts = text.split("\t")
                    if len(parts) >= 3:
                        species_raw = parts[1].strip()
                        st_raw = parts[2].strip()
                        summary["species"] = species_raw if species_raw != "-" else pd.NA
                        summary["ST"] = pd.to_numeric(st_raw, errors="coerce")
                        summary["status"] = "no_call" if (species_raw == "-" and st_raw == "-") else "ok"
                    else:
                        summary["status"] = "malformed"

            rows.append(summary)

        out = pd.DataFrame(rows)
        out["uid"] = pd.to_numeric(out["uid"], errors="coerce").astype("Int64")
        out["ST"] = out["ST"].astype("Int64")
        out.to_csv(self.output_path(), index=False)
        return f"{len(out)} isolates"


if __name__ == "__main__":
    run_cli(MLSTStep, CollectMLSTStep)
