import sys
import argparse
from collections import Counter
from pathlib import Path
from typing import Any, Callable

import io

from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import Phylo

from pycirclize import Circos
from scipy.cluster.hierarchy import linkage, to_tree
from scipy.spatial.distance import squareform

PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR      = PROJECT_ROOT / "figures" / "patient_trees"
MAX_DIST       = 0.1
SCALE_BAR_VAL  = 0.02

FONTTYPE = "Arial"
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = [FONTTYPE]
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_utils import (
    load_color, get_patient_meta, load_mash_distances, load_uti_coverage,
    classify_uti_coverage, get_all_isolate_meta, UNCOV_THRESH,
)
FECAL_COLOR         = load_color("FAECAL_same")
UTI_COVERED_COLOR   = load_color("UTI_covered")
UTI_UNCOVERED_COLOR = load_color("UTI_covered")
UTI_NOT_TESTED_COLOR = (0.7, 0.7, 0.7)


def plot_tree_legend(out_svg: Path, is_small: bool = False) -> None:
    """Standalone legend for the circular patient trees."""
    scale = 0.43 if is_small else 1.0
    shape_scale = 0.4 if is_small else 1.0
    fig, ax = plt.subplots(figsize=(3.2 * scale, 1.1 * scale))
    fig.subplots_adjust(left=0.05, right=1, top=1, bottom=0)
    ax.set_axis_off()

    entries = [
        ("o",  60 * shape_scale, FECAL_COLOR,         True,  "Faecal isolate"),
        ("s", 200 * shape_scale, UTI_COVERED_COLOR,   True,  "Urine isolate, microbiome-covered  (\u03B7\u2080<1e-4)"),
        ("s", 200 * shape_scale, UTI_UNCOVERED_COLOR, False, "Urine isolate, microbiome-uncovered  (\u03B7\u2080>1e-4)"),
    ]

    lw = 1.4 if is_small else 2
    fs = 6 if is_small else 8
    y_positions = [0.75, 0.45, 0.15]
    for (marker, ms, color, filled, label), y in zip(entries, y_positions):
        fc = color if filled else "none"
        ax.scatter([0.04], [y], marker=marker, s=ms, facecolors=fc,
                   edgecolors=color, lw=lw, clip_on=False,
                   transform=ax.transAxes, zorder=5)
        ax.text(0.10, y, label, va="center", ha="left", fontsize=fs,
                transform=ax.transAxes)

    out_svg.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_svg, facecolor="none" if is_small else "white",
                bbox_inches="tight", pad_inches=0.01)
    plt.close(fig)
    print(f"Saved legend: {out_svg}")


def linkage_to_newick(root, labels: list[str]) -> str:
    """Convert scipy linkage ClusterNode to newick string with branch lengths."""
    def recurse(node, parent_dist: float) -> str:
        branch = parent_dist - node.dist
        if node.is_leaf():
            return f"{labels[node.id]}:{branch:.8f}"
        left  = recurse(node.left,  node.dist)
        right = recurse(node.right, node.dist)
        return f"({left},{right}):{branch:.8f}"
    return recurse(root, root.dist) + ";"


def plot_circular_tree(
    nwk: str,
    leaf_to_symbol: dict[str, object],
    leaf_to_color: dict[str, tuple[float, float, float]],
    out_svg: Path,
    leaf_to_filled: dict[str, bool] | None = None,
    tree_line_kws: dict[str, Any] | None = None,
    tree_label_formatter: Callable[[str], str] | None = None,
    bottom_left_label: str | None = None,
    scale_bar_value: float | None = None,
    max_radial_depth: float | None = None,
    is_small: bool = False,
) -> None:
    min_fig_size = 2.6
    max_fig_size = 8.0
    size_per_leaf = 0.14

    tree = Phylo.read(io.StringIO(nwk), "newick")
    if tree_line_kws is None:
        tree_line_kws = {"color": "k", "lw": 0.8}
    if tree_label_formatter is None:
        tree_label_formatter = lambda t: t.replace("_", " ")

    total_depth = max(tree.depths().values())
    if max_radial_depth and total_depth < max_radial_depth:
        r_inner = 100 * (1 - total_depth / max_radial_depth)
    else:
        r_inner = 0

    circos, tv = Circos.initialize_from_tree(
        tree,
        r_lim=(r_inner, 100),
        leaf_label_size=0,
        line_kws=tree_line_kws,
        label_formatter=tree_label_formatter,
    )

    # Draw leaf-end ticks + labels first, then place symbols outside labels.
    leaf_labels = tv.leaf_labels
    leaf_x = {label: tv.name2xr[label][0] for label in leaf_labels}
    leaf_x_vals = [leaf_x[label] for label in leaf_labels]

    tv.track.xticks(
        leaf_x_vals,
        labels=leaf_labels,
        tick_length=1.5,
        outer=True,
        show_bottom_line=False,
        label_size=6 * (0.85 if is_small else 1),
        label_margin=4.0,
        label_orientation="vertical",
        line_kws=dict(color="#737373", lw=0.8),
        text_kws=dict(color="k"),
    )

    fig = circos.plotfig()
    fig_size = float(np.clip(size_per_leaf * len(leaf_labels), min_fig_size, max_fig_size))
    if is_small:
        fig_size *= 0.5
    fig.set_size_inches(fig_size, fig_size)
    ax = fig.axes[0]

    # Place symbols on an outer circle beyond the labels.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    max_r = 0
    for txt in ax.texts:
        bbox = txt.get_window_extent(renderer)
        corners_data = ax.transData.inverted().transform(bbox.corners())
        max_r = max(max_r, corners_data[:, 1].max())
    marker_r = max_r + 32 / fig_size
    for leaf_name in leaf_labels:
        if leaf_name not in leaf_to_symbol or leaf_name not in leaf_to_color:
            continue
        rad = tv.track.x_to_rad(leaf_x[leaf_name])
        marker = leaf_to_symbol[leaf_name]
        shape_scale = 0.4 if is_small else 1
        marker_size = 60 * shape_scale
        if marker == "s":
            marker = (4, 0, 45 - np.degrees(rad))
            marker_size = 200 * shape_scale
        filled = leaf_to_filled.get(leaf_name, True) if leaf_to_filled else True
        color = leaf_to_color[leaf_name]
        ax.scatter(
            [rad],
            [marker_r],
            facecolors=color if filled else "none",
            edgecolors=color,
            marker=marker,
            s=marker_size,
            lw=1.4 if is_small else 2,
            clip_on=False,
            zorder=6,
        )

    # ── scale bar ──────────────────────────────────────────────────────────
    bar_val = scale_bar_value
    effective_depth = max_radial_depth or total_depth
    if effective_depth > 0 and bar_val:
        bar_r = (bar_val / effective_depth) * 100  # data r-units
        disp_0 = ax.transData.transform((0, 50))
        disp_1 = ax.transData.transform((0, 50 + bar_r))
        fig_inv = fig.transFigure.inverted()
        f0, f1 = fig_inv.transform(disp_0), fig_inv.transform(disp_1)
        bar_fig_len = np.hypot(f1[0] - f0[0], f1[1] - f0[1])
        x0, y0 = 0.06, 0.04
        fig.add_artist(Line2D(
            [x0, x0 + bar_fig_len], [y0, y0],
            transform=fig.transFigure, color="k", lw=1.5,
            solid_capstyle="butt", zorder=10,
        ))
        label = f"{bar_val:g}"
        fig.text(x0 + bar_fig_len / 2, y0 + 0.012, label,
                 ha="center", va="bottom", fontsize=5,
                 transform=fig.transFigure)

    out_svg.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_svg, facecolor="none" if is_small else "white")
    plt.close(fig)
    print(f"Saved: {out_svg}")


def plot_patient_tree(fid: int, dist_df: pd.DataFrame | None = None) -> bool:
    # ── load metadata ────────────────────────────────────────────────────────
    meta = get_patient_meta(fid)
    if meta.empty:
        print(f"Patient {fid}: no isolate metadata, skipping.")
        return False

    faecal_types = set(meta["is_faecal"].dropna().astype(bool).unique().tolist())
    if not {True, False}.issubset(faecal_types):
        print(f"Patient {fid}: missing UTI or faecal isolates, skipping.")
        return False

    # ── load distance matrix ─────────────────────────────────────────────────
    if dist_df is None:
        dist_df = load_mash_distances()
    run_names = [r for r in meta["run_name"] if r in dist_df.index]
    meta      = meta[meta["run_name"].isin(run_names)].set_index("run_name")
    sample_name_by_run = meta["sample_name"].to_dict()
    base_leaf_names = [sample_name_by_run[r] for r in run_names]
    duplicated_labels = {
        label for label, count in Counter(base_leaf_names).items() if count > 1
    }
    leaf_name_by_run = {
        run_name: (
            run_name if sample_name_by_run[run_name] in duplicated_labels else sample_name_by_run[run_name]
        )
        for run_name in run_names
    }
    leaf_names = [leaf_name_by_run[r] for r in run_names]

    n = len(run_names)
    print(f"Patient {fid}: {n} isolates")
    if n < 2:
        print(f"Patient {fid}: need at least 2 isolates, skipping.")
        return False

    D = dist_df.loc[run_names, run_names].values.astype(float)
    np.fill_diagonal(D, 0.0)
    if MAX_DIST is not None:
        D = np.clip(D, 0, MAX_DIST)

    # ── build linkage tree → newick ──────────────────────────────────────────
    Z    = linkage(squareform(D), method="average")   # UPGMA
    root = to_tree(Z, rd=False)
    nwk  = linkage_to_newick(root, leaf_names)

    uncov_by_uti = load_uti_coverage()

    leaf_to_symbol: dict[str, object] = {}
    leaf_to_color: dict[str, tuple[float, float, float]] = {}
    leaf_to_filled: dict[str, bool] = {}
    for run_name in run_names:
        leaf_name = leaf_name_by_run[run_name]
        is_faecal = run_name in meta.index and bool(meta.loc[run_name, "is_faecal"])
        if is_faecal:
            leaf_to_symbol[leaf_name] = "o"
            leaf_to_color[leaf_name] = FECAL_COLOR
            leaf_to_filled[leaf_name] = True
        else:
            leaf_to_symbol[leaf_name] = "s"
            cov = classify_uti_coverage(run_name, uncov_by_uti)
            if cov is None:
                leaf_to_color[leaf_name] = UTI_NOT_TESTED_COLOR
                leaf_to_filled[leaf_name] = False
            else:
                leaf_to_color[leaf_name] = UTI_COVERED_COLOR if cov else UTI_UNCOVERED_COLOR
                leaf_to_filled[leaf_name] = cov

    out_svg = OUT_DIR / f"circular_tree_P{fid}.svg"
    common_kw = dict(
        leaf_to_filled=leaf_to_filled,
        bottom_left_label=f"Patient: {fid}",
        scale_bar_value=SCALE_BAR_VAL,
        max_radial_depth=MAX_DIST,
    )
    plot_circular_tree(nwk, leaf_to_symbol, leaf_to_color, out_svg, **common_kw)

    if fid == 6:
        out_small = OUT_DIR.parent / f"circular_tree_P{fid}_small.svg"
        plot_circular_tree(nwk, leaf_to_symbol, leaf_to_color, out_small,
                           is_small=True, **common_kw)
    return True


def get_fids_with_uti_and_faecal_isolates() -> list[int]:
    meta = get_all_isolate_meta()
    by_fid = meta.groupby("fid")["is_faecal"].agg(
        lambda s: set(s.dropna().astype(bool))
    )
    return sorted(fid for fid, types in by_fid.items() if {True, False}.issubset(types))


def main(fid: int | None = None, legend_only: bool = False) -> None:
    plot_tree_legend(OUT_DIR / "tree_legend.svg")
    plot_tree_legend(OUT_DIR.parent / "tree_legend_small.svg", is_small=True)

    if legend_only:
        return

    dist_df = load_mash_distances()
    if fid is not None:
        plot_patient_tree(fid, dist_df=dist_df)
        return

    fids = get_fids_with_uti_and_faecal_isolates()
    if not fids:
        print("No patients found with both UTI and faecal isolates.")
        return

    print(f"Plotting circular trees for {len(fids)} patients.")
    for patient_fid in fids:
        plot_patient_tree(patient_fid, dist_df=dist_df)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Circular phylogenetic tree from MASH distances.")
    parser.add_argument("--fid", type=int, default=None, help="Single patient FID. If omitted, run all eligible patients.")
    parser.add_argument("--legend-only", action="store_true", help="Generate only the legend SVG.")
    args = parser.parse_args()
    main(args.fid, legend_only=args.legend_only)
