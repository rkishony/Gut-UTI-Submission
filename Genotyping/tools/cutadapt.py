"""Cutadapt wrapper for adapter/quality trimming."""

import json
from pathlib import Path
from typing import Optional

from tools.run_utils import get_tool, run_cmd


def load_cutadapt_report(path: Path) -> dict | None:
    """Load cutadapt JSON report. Returns None if path does not exist."""
    if not path.exists():
        return None
    with open(path) as f:
        return json.load(f)


def extract_trim_stats(report: dict | None) -> dict:
    """Extract trimming statistics from cutadapt JSON report."""
    if report is None:
        return {
            "raw_num_reads": 0,
            "num_reads": 0,
            "percentage_reads_trimmed": 0,
            "raw_read_len": 0,
            "read_len": 0,
            "total_bases": 0,
            "percentage_with_adapter": 0,
        }
    rc = report["read_counts"]
    bc = report["basepair_counts"]
    raw_reads = rc["input"]
    out_reads = rc["output"]
    r1_adapter = rc.get("read1_with_adapter", 0) or 0
    r2_adapter = rc.get("read2_with_adapter", 0) or 0

    return {
        "raw_num_reads": raw_reads,
        "num_reads": out_reads,
        "percentage_reads_trimmed": round((raw_reads - out_reads) / raw_reads * 100, 2) if raw_reads else 0,
        "raw_read_len": round(bc["input"] / (raw_reads * 2), 1) if raw_reads else 0,
        "read_len": round(bc["output"] / (out_reads * 2), 1) if out_reads else 0,
        "total_bases": bc["output"],
        "percentage_with_adapter": round((r1_adapter + r2_adapter) / (raw_reads * 2) * 100, 2) if raw_reads else 0,
    }


def clip_and_trim_paired(
    r1_in: Path,
    r2_in: Path,
    r1_out: Path,
    r2_out: Path,
    report_path: Optional[Path] = None,
    quality_threshold: int = 20,
    min_length: int = 35,
    adapter: str = "CTGTCTCTTATA",
) -> None:
    """Trim R1 and R2 together with cutadapt (paired-end)."""
    cmd = [
        get_tool("cutadapt"),
        "-a", adapter,
        "-A", adapter,
        "-q", str(quality_threshold),
        "-m", str(min_length),
        "-o", str(r1_out),
        "-p", str(r2_out),
    ]
    if report_path:
        cmd += ["--json", str(report_path)]
    cmd += [str(r1_in), str(r2_in)]
    run_cmd(cmd)
