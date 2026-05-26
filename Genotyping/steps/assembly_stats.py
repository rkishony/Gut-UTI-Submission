#!/usr/bin/env python3
"""Step: Collect assembly statistics from Unicycler log files."""

import re
import pandas as pd
from pathlib import Path

from config_utils import get_output_dir
from steps.extract_run_names import ExtractRunNamesStep
from steps.step import Step, run_cli
from steps.assemble import AssembleStep, get_error_reason_from_log


FILENAME = "assembly_stats.csv"


def _parse_stats_parts(parts: list[str]) -> dict:
    def _parse_int(s: str) -> int:
        return int(s.replace(",", ""))

    return {
        "num_segments": _parse_int(parts[1]),
        "num_links": _parse_int(parts[2]),
        "total_length": _parse_int(parts[3]),
        "N50": _parse_int(parts[4]),
        "longest_segment": _parse_int(parts[5]),
    }


class AssemblyStatsStep(Step):
    def __init__(self, csv_file: str | Path = None, **kwargs):
        super().__init__(**kwargs)
        self.csv_file = csv_file or ExtractRunNamesStep.output_path()

    @classmethod
    def output_path(cls) -> Path:
        return get_output_dir("summary") / FILENAME

    @staticmethod
    def _parse_log(log_text: str) -> dict:
        """Parse the component table from the 'Bridged assembly graph' section."""
        # Find the component table
        pattern = (
            r"Component\s+Segments\s+Links\s+Length\s+N50\s+Longest segment\s+Status\s*\n"
            r"(.*?)(?:\n\n|\Z)"
        )
        m = re.search(pattern, log_text, re.DOTALL)
        if not m:
            return {
                "num_segments": 0, "num_links": 0, "total_length": 0,
                "N50": 0, "longest_segment": 0,
                "num_components": 0, "num_complete": 0, "num_incomplete": 0,
            }

        lines = [ln for ln in m.group(1).strip().splitlines() if ln.strip()]
        total_stats = {}
        num_components = 0
        num_complete = 0
        num_incomplete = 0

        for line in lines:
            parts = line.split()
            if not parts:
                continue
            if parts[0] == "total":
                total_stats = _parse_stats_parts(parts)
            else:
                num_components += 1
                status = parts[-1] if parts[-1] in ("complete", "incomplete") else ""
                if status == "complete":
                    num_complete += 1
                elif status == "incomplete":
                    num_incomplete += 1

        # When there's only one component, unicycler omits the "total" row
        if not total_stats and num_components == 1:
            parts = [ln for ln in lines if ln.split()[0] != "total"][0].split()
            total_stats = _parse_stats_parts(parts)

        total_stats.update({
            "num_components": num_components,
            "num_complete": num_complete,
            "num_incomplete": num_incomplete,
        })
        return total_stats

    def _run(self) -> None:
        df = pd.read_csv(self.csv_file)
        df = df[df["is_isolate"]]
        rows = []
        for _, row in df.iterrows():
            run_name = row["run_name"]

            asm = AssembleStep.run_output_path(row)
            fasta_path = asm["assembly"]
            has_fasta = fasta_path.exists()

            log_path = asm["log"]
            if log_path.exists():
                log_text = log_path.read_text()
                stats = self._parse_log(log_text)
                error_reason = get_error_reason_from_log(log_text)
            else:
                stats = self._parse_log("")
                error_reason = "no_log"

            stats["run_name"] = run_name
            stats["has_fasta"] = has_fasta
            stats["status"] = "ok" if error_reason is None else error_reason
            rows.append(stats)

        out = pd.DataFrame(rows)
        cols = ["run_name", "has_fasta", "status", "num_segments", "num_links", "total_length",
                "N50", "longest_segment", "num_components", "num_complete", "num_incomplete"]
        out = out[cols]
        out.to_csv(self.output_path(), index=False)
        return f"{len(out)} runs"


if __name__ == "__main__":
    run_cli(AssemblyStatsStep)
