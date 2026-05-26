"""CDF of MASH distance from each UTI isolate to its closest same-patient
fecal isolate, stratified by metagenome coverage (covered vs uncovered)."""

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_utils import (
    load_color, load_mash_distances, load_uti_coverage,
    classify_uti_coverage, get_all_isolate_meta,
)
from record_value import record_value

PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUT_PNG      = PROJECT_ROOT / "figures" / "uti_distance_to_closest_fecal.svg"

FONTTYPE = "Arial"
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = [FONTTYPE]
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"

FIRST_ISOLATE_ONLY = True

UTI_COVERED_COLOR   = load_color("UTI_covered")
UTI_UNCOVERED_COLOR = load_color("UTI_uncovered")


def _collect_distances(meta, dist_df, uncov_by_uti):
    """Return (covered_dists, uncovered_dists, covered_fids, uncovered_fids)."""
    covered_dists, uncovered_dists = [], []
    covered_fids, uncovered_fids = [], []
    for fid, grp in meta.groupby("fid"):
        uti_runs = [r for r in grp.loc[grp["is_faecal"] == False, "run_name"]
                    if r in dist_df.index]
        fecal_runs = [r for r in grp.loc[grp["is_faecal"] == True, "run_name"]
                      if r in dist_df.index]
        if not uti_runs or not fecal_runs:
            continue
        for u in uti_runs:
            cov = classify_uti_coverage(u, uncov_by_uti)
            if cov is None:
                continue
            min_dist = dist_df.loc[u, fecal_runs].min()
            if cov:
                covered_dists.append(min_dist)
                covered_fids.append(fid)
            else:
                uncovered_dists.append(min_dist)
                uncovered_fids.append(fid)
    return (np.array(covered_dists), np.array(uncovered_dists),
            np.array(covered_fids), np.array(uncovered_fids))


def _filter_first_isolate(meta):
    """Keep only the UTI isolate with the lowest letter per (fid, uid)."""
    uti_mask = meta["is_faecal"] == False
    uti_meta = meta[uti_mask].sort_values("isolate")
    first_idx = uti_meta.groupby(["fid", "uid"]).head(1).index
    return meta[~uti_mask | meta.index.isin(first_idx)]


def compute_min_distances():
    """Return (all_covered, all_uncovered, filt_covered, filt_uncovered,
            all_cov_fids, all_uncov_fids)."""
    meta = get_all_isolate_meta()
    dist_df = load_mash_distances()
    uncov_by_uti = load_uti_coverage()

    all_cov, all_uncov, all_cov_fids, all_uncov_fids = \
        _collect_distances(meta, dist_df, uncov_by_uti)

    if FIRST_ISOLATE_ONLY:
        filt_meta = _filter_first_isolate(meta)
        filt_cov, filt_uncov, _, _ = _collect_distances(filt_meta, dist_df, uncov_by_uti)
    else:
        filt_cov, filt_uncov = all_cov, all_uncov

    return all_cov, all_uncov, filt_cov, filt_uncov, all_cov_fids, all_uncov_fids


def _pval_stars(p: float) -> str:
    if p < 1e-3:
        return "***"
    if p < 1e-2:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def plot_strip(covered: np.ndarray, uncovered: np.ndarray,
               filt_covered: np.ndarray, filt_uncovered: np.ndarray,
               cov_fids: np.ndarray = None, uncov_fids: np.ndarray = None,
               highlight_fid: int = 6) -> None:
    fig, ax = plt.subplots(figsize=(1.8, 2.8))

    rng = np.random.default_rng(42)
    jitter_w = 0.25
    for pos, vals, fids, color, label in [
        (1, covered,   cov_fids,   UTI_COVERED_COLOR,   f"Covered (n={len(covered)})"),
        (2, uncovered, uncov_fids, UTI_UNCOVERED_COLOR, f"Uncovered (n={len(uncovered)})"),
    ]:
        jx = rng.normal(0, jitter_w / 2.5, size=len(vals)).clip(-jitter_w, jitter_w)
        ax.scatter(pos + jx, vals, s=20, alpha=0.7,
                   color=color, edgecolors="none", zorder=3, label=label)
        if fids is not None:
            mask = fids == highlight_fid
            if mask.any():
                ax.scatter(pos + jx[mask], vals[mask], s=20, facecolors="none",
                           edgecolors="black", linewidths=0.5, zorder=5)
        ax.plot([pos - 0.2, pos + 0.2], [np.median(vals)] * 2,
                color="k", lw=1.5, zorder=4)

    ax.set_yscale("log")
    ax.set_xticks([1, 2])
    ax.set_xticklabels(["Covered", "Uncovered"], fontsize=6)
    ax.tick_params(axis="both", labelsize=6)
    ax.set_xlabel("Urine isolate coverage\nby pathometagenome", fontsize=6)
    ax.set_ylabel("MASH distance to closest fecal isolate", fontsize=6)

    # Mann–Whitney U test (using filtered data — one isolate per UTI)
    _, pval = mannwhitneyu(filt_covered, filt_uncovered, alternative="two-sided")
    stars = _pval_stars(pval)
    record_value("UTI dist to fecal, P-val mannwhitney", "uti_dist_to_feacal.pval", pval, "%.1e")

    y_top = max(covered.max(), uncovered.max()) * 5
    bar_y = max(covered.max(), uncovered.max()) * 1.8
    ax.plot([1, 1, 2, 2], [bar_y * 0.9, bar_y, bar_y, bar_y * 0.9],
            lw=1, color="k")
    ax.text(1.5, bar_y * 1.02, stars, ha="center", va="bottom", fontsize=6)
    ax.set_ylim(top=y_top)

    fig.tight_layout()
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, facecolor="white")
    print(f"Saved: {OUT_PNG}")
    plt.close(fig)


def main() -> None:
    all_cov, all_uncov, filt_cov, filt_uncov, cov_fids, uncov_fids = compute_min_distances()
    record_value("UTI dist to fecal", "uti_dist_to_feacal.n_covered", len(all_cov))
    record_value("UTI dist to fecal", "uti_dist_to_feacal.n_uncovered", len(all_uncov))
    record_value("UTI dist to fecal", "uti_dist_to_feacal.median_covered", float(np.median(all_cov)), "%.4f")
    record_value("UTI dist to fecal", "uti_dist_to_feacal.median_uncovered", float(np.median(all_uncov)), "%.4f")
    plot_strip(all_cov, all_uncov, filt_cov, filt_uncov, cov_fids, uncov_fids)


if __name__ == "__main__":
    main()
