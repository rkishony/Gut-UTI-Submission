#!/usr/bin/env python3
"""Step: Assemble trimmed reads using Unicycler."""

import pandas as pd
from pathlib import Path

from config import MAX_READS_FOR_ASSEMBLY
from steps.per_run_step import PerRunStep
from steps.step import run_cli
from steps.trim_sequences import TrimSequencesStep
from tools.seqtk import downsample_fastq
from tools.unicycler import run_unicycler


REASSONS_TO_ERROR_MESSAGES = {
    "spades_failed": "Error: SPAdes failed to produce assemblies",
}


def get_error_reason_from_log(log_text: str) -> str | None:
    for reason, error_message in REASSONS_TO_ERROR_MESSAGES.items():
        if error_message in log_text:
            return reason
    return None


class AssembleStep(PerRunStep):
    output_dir_name = "assembly"
    cores_per_run = 4

    output_paths = {"assembly": "assembly.fasta", "log": "unicycler.log"}

    def run_input_paths(self, row: pd.Series) -> dict:
        trim = TrimSequencesStep.run_output_path(row)
        return {"R1": trim["R1"], "R2": trim["R2"]}

    def filter_df(self, df: pd.DataFrame) -> pd.DataFrame:
        return df[df["is_isolate"]]

    def _should_process(self, row: pd.Series) -> bool:
        log_file = self.run_output_path(row)["log"]
        if not log_file.exists():
            return True
        with open(log_file, "r") as f:
            log_text = f.read()
        reason = get_error_reason_from_log(log_text)
        if reason is not None:
            self.print(f"{reason}, skipping", row)
            return False
        return True

    @staticmethod
    def _needs_downsampling(n_reads: int) -> bool:
        return MAX_READS_FOR_ASSEMBLY is not None and n_reads > MAX_READS_FOR_ASSEMBLY

    def _run_output_exists(self, row: pd.Series) -> bool:
        if not super()._run_output_exists(row):
            return False
        # Re-run if downsampling is needed but wasn't done or was done with a different target
        if self._needs_downsampling(row["num_reads"]):
            downsampled_flag = self.run_folder(row) / "downsampled.txt"
            if not downsampled_flag.exists():
                self.print("needs downsampling, re-running", row)
                return False
            if f"to {MAX_READS_FOR_ASSEMBLY:,}" not in downsampled_flag.read_text():
                self.print("downsampled with different target, re-running", row)
                return False
        return True

    def process_run(self, row: pd.Series) -> None:
        run_name = row["run_name"]
        inputs = self.run_input_paths(row)
        r1, r2 = inputs["R1"], inputs["R2"]
        out_dir = self.run_folder(row)

        n_reads = row["num_reads"]
        downsample = self._needs_downsampling(n_reads)
        if downsample:
            self.print(f"downsampling {n_reads:,} -> {MAX_READS_FOR_ASSEMBLY:,} read pairs", row)
            sub_r1 = out_dir / "sub.R1.fastq.gz"
            sub_r2 = out_dir / "sub.R2.fastq.gz"
            downsample_fastq(r1, sub_r1, MAX_READS_FOR_ASSEMBLY)
            downsample_fastq(r2, sub_r2, MAX_READS_FOR_ASSEMBLY)
            r1, r2 = sub_r1, sub_r2

            (out_dir / "downsampled.txt").write_text(
                f"Downsampled from {n_reads:,} to {MAX_READS_FOR_ASSEMBLY:,} read pairs"
            )

        result = run_unicycler(r1, r2, out_dir, self.threads_per_worker)

        for gfa in out_dir.glob("*.gfa"):
            gfa.unlink()

        # delete the downsampled files if created
        if downsample:
            sub_r1.unlink()
            sub_r2.unlink()

        if result.returncode != 0:
            raise RuntimeError(
                f"Unicycler failed for {run_name} (exit {result.returncode}):\n{result.stderr}"
            )


if __name__ == "__main__":
    run_cli(AssembleStep)
