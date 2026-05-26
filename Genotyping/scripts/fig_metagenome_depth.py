import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

PROJECT_ROOT = Path(__file__).resolve().parent.parent

FONTTYPE = "Arial"
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = [FONTTYPE]
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"
sys.path.insert(0, str(PROJECT_ROOT))
from config import BASE_READS_PER_COVERAGE_THRESHOLD, BASE_READS_PER_COVERAGE_TARGET

CSV_PATH = PROJECT_ROOT / "output/summary/trim_stats_sorted.csv"
OUT_PATH = PROJECT_ROOT / "figures/metagenome_depth.svg"
NORM_OUT_PATH = PROJECT_ROOT / "figures/metagenome_depth_normalized.svg"

RUN_TYPE_COLORS = {"hs1": "steelblue", "hs2": "mediumseagreen", "ns": "darkorange"}
LEGEND_HANDLES = [Patch(facecolor=c, label=rt) for rt, c in RUN_TYPE_COLORS.items()]


def _plot_panel(ax, data: pd.DataFrame, values: pd.Series, *,
                title: str, ylabel: str, log: bool,
                ymin: float = None, ymax: float = None,
                threshold: float = None, show_legend: bool = False):
    colors = [RUN_TYPE_COLORS.get(rt, "gray") for rt in data["batch"]]
    ax.bar(range(len(data)), values, color=colors, edgecolor="none")
    if threshold is not None:
        ax.axhline(y=threshold, color="red", linestyle="--", alpha=0.7)
    ax.set_xlabel("Sample (ranked)")
    ax.set_xlim(-0.5, len(data) - 0.5)
    if ymin is not None or ymax is not None:
        ax.set_ylim(ymin, ymax)
    if log:
        ax.set_yscale("log")
    scale = "log" if log else "linear"
    ax.set_title(f"{title} (n={len(data)}) — {scale}")
    ax.set_ylabel(ylabel)
    if show_legend:
        ax.legend(handles=LEGEND_HANDLES)


def main():
    df = pd.read_csv(CSV_PATH)
    meta = df[~df["is_isolate"]].sort_values("num_reads", ascending=False).reset_index(drop=True)
    iso = df[df["is_isolate"]].sort_values("num_reads", ascending=False).reset_index(drop=True)

    # --- Raw depth figure ---
    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
    for data, label, row, threshold, ymin, ymax, ymin_log, ymax_log in [
        (meta, "Metagenome", 0, 30, 0, 100, 1,   300),
        (iso,  "Isolate",    1, 1,  0, 5,   0.1, 10),
    ]:
        vals = data["num_reads"] / 1e6
        _plot_panel(axes[row, 0], data, vals, title=label, ylabel="M reads",
                    log=False, ymin=ymin, ymax=ymax, threshold=threshold, show_legend=(row == 0))
        _plot_panel(axes[row, 1], data, vals, title=label, ylabel="M reads (log)",
                    log=True, ymin=ymin_log, ymax=ymax_log, threshold=threshold)
    fig.tight_layout()
    fig.savefig(OUT_PATH)
    print(f"Saved to {OUT_PATH}")

    # --- Normalized figure (metagenome only) ---
    fig2, axes2 = plt.subplots(2, 2, figsize=(14, 9))
    for row, norm_col in enumerate(["CFU", "coverage_factor"]):
        sub = meta[meta[norm_col].notna() & (meta[norm_col] > 0)].copy()
        sub = sub.sort_values("num_reads", ascending=False).reset_index(drop=True)
        vals = sub["num_reads"] / 1e6 / sub[norm_col]
        _plot_panel(axes2[row, 0], sub, vals, title=f"reads / {norm_col}",
                    ylabel=f"M reads / {norm_col}", log=False, show_legend=(row == 0))
        _plot_panel(axes2[row, 1], sub, vals, title=f"reads / {norm_col}",
                    ylabel=f"M reads / {norm_col} (log)", log=True)
    fig2.tight_layout()
    fig2.savefig(NORM_OUT_PATH)
    print(f"Saved to {NORM_OUT_PATH}")

    # --- Scatter: num_reads vs CFU / coverage_factor ---
    scatter_out = PROJECT_ROOT / "figures" / "metagenome_depth_scatter.svg"
    fig3, axes3 = plt.subplots(1, 2, figsize=(14, 5))
    for ax, col in zip(axes3, ["CFU", "coverage_factor"]):
        sub = meta[meta[col].notna() & (meta[col] > 0)].copy()
        colors = [RUN_TYPE_COLORS.get(rt, "gray") for rt in sub["batch"]]
        ax.scatter(sub[col], sub["num_reads"] / 1e6, c=colors, edgecolor="none", alpha=0.7)
        ax.set_xlabel(col)
        ax.set_ylabel("M reads")
        ax.set_title(f"num_reads vs {col} (n={len(sub)})")
        ax.set_xscale("log")
        ax.set_yscale("log")
    ax_cov = axes3[1]  # coverage_factor panel
    xlims = ax_cov.get_xlim()
    xs = np.logspace(np.log10(xlims[0]), np.log10(xlims[1]), 200)
    if BASE_READS_PER_COVERAGE_THRESHOLD is not None:
        ax_cov.plot(xs, BASE_READS_PER_COVERAGE_THRESHOLD * xs / 1e6,
                    color="red", linestyle="--", alpha=0.7,
                    label=f"threshold ({BASE_READS_PER_COVERAGE_THRESHOLD:,} × cov)")
    if BASE_READS_PER_COVERAGE_TARGET is not None:
        ax_cov.plot(xs, BASE_READS_PER_COVERAGE_TARGET * xs / 1e6,
                    color="blue", linestyle="--", alpha=0.7,
                    label=f"target ({BASE_READS_PER_COVERAGE_TARGET:,} × cov)")
    ax_cov.set_xlim(xlims)
    ax_cov.legend()
    axes3[0].legend(handles=LEGEND_HANDLES)
    fig3.tight_layout()
    fig3.savefig(scatter_out)
    print(f"Saved to {scatter_out}")


if __name__ == "__main__":
    main()
