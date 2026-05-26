"""
Radial-line MLST pie chart: one line per sample, colored by ST category.

Usage:
    python scripts/fig_mlst_pie.py
"""

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from pathlib import Path

from plot_utils import load_color

FONTTYPE = "Arial"
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = [FONTTYPE]
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"
from fig_strainphlan_tree import load_tree_data, OUTPUT_DIR
from record_value import record_value

MULTI_ALT_COLORS = [load_color('st_multiple1'), load_color('st_multiple2')]
UNRESOLVED_COLOR = load_color('st_unresolved')


def _load_pie_data():
    """Load tree + MLST data and return sorted DataFrame with colors."""
    _tree, remaining, cats, leaf_sts, _sample_name_map = load_tree_data()

    rows = [{'sample': n, 'cat': cats[n], 'st': leaf_sts.get(n)} for n in remaining]
    df = pd.DataFrame(rows)

    df['st'] = df['st'].fillna(-1).astype(int)
    st_counts = df.groupby('st')['sample'].transform('count')
    df['count'] = st_counts.astype(int)

    cat_order = {'st_multiple': 0, 'st_unique_known': 1, 'st_unique_new': 2, 'st_none': 3}
    df['cat_order'] = df['cat'].map(cat_order)
    df = df.sort_values(['cat_order', 'count', 'st'],
                        ascending=[True, False, True]).reset_index(drop=True)
    df.drop(columns='cat_order', inplace=True)

    color_new = load_color('st_unique_new')
    color_known = load_color('st_unique_known')

    colors, prev_st, ci = [], None, 0
    for _, r in df.iterrows():
        if r['cat'] == 'st_none':
            colors.append(UNRESOLVED_COLOR)
        elif r['cat'] == 'st_unique_new':
            colors.append(color_new)
        elif r['cat'] == 'st_unique_known':
            colors.append(color_known)
        else:
            if r['st'] != prev_st:
                ci += 1
                prev_st = r['st']
            colors.append(MULTI_ALT_COLORS[ci % 2])
    df['color'] = colors

    group = []
    for _, r in df.iterrows():
        if r['cat'] == 'st_none':
            group.append('unresolved')
        elif r['cat'] == 'st_multiple':
            group.append(f"ST {int(r['st'])}")
        else:
            group.append(r['cat'])
    df['group'] = group
    return df


def _record_pie_values(df):
    """Record ST group counts and validate against typing_summary.csv."""
    from fig_strainphlan_tree import BASE_DIR, TRIM_STATS

    trim_df = pd.read_csv(TRIM_STATS)
    n_total = len(trim_df[~trim_df['is_isolate']])
    n_replicates = int(trim_df[~trim_df['is_isolate']]['is_rep'].sum())
    n_valid = n_total - n_replicates
    n_not_in_tree = n_valid - len(df)

    desc = 'MLST pie'
    record_value(desc, 'MLST_pie.n_total_metagenome_runs', n_total)
    record_value(desc, 'MLST_pie.n_replicates', n_replicates)
    record_value(desc, 'MLST_pie.n_valid', n_valid)
    record_value(desc, 'MLST_pie.n_not_in_tree', n_not_in_tree)
    record_value(desc, 'MLST_pie.n_in_tree_and_pie', len(df))
    record_value(desc, 'MLST_pie.n_st_unique_known', (df['cat'] == 'st_unique_known').sum())
    record_value(desc, 'MLST_pie.n_st_unique_new', (df['cat'] == 'st_unique_new').sum())
    record_value(desc, 'MLST_pie.n_st_unresolved', (df['cat'] == 'st_none').sum())
    record_value(desc, 'MLST_pie.n_st_multiple', (df['cat'] == 'st_multiple').sum())

    typing_csv = BASE_DIR / 'output/summary/typing_summary.csv'
    ts = pd.read_csv(typing_csv)
    expected = set(ts.loc[ts['in_pie_line_chart'], 'run_name'])
    actual = set(df['sample'])
    assert actual == expected, (
        f"Pie samples mismatch with typing_summary.csv: "
        f"extra={actual - expected}, missing={expected - actual}"
    )



def plot_mlst_pie(output_path=None):
    """Radial-line pie: one line per sample, colored by ST category."""
    df = _load_pie_data()
    _record_pie_values(df)
    n = len(df)

    fig, ax = plt.subplots(figsize=(3.5, 2.9))
    ax.set_aspect('equal')
    ax.set_xlim(-1.6, 1.6)
    ax.set_ylim(-1.6, 1.6)
    ax.axis('off')

    angles = np.pi / 2 - np.linspace(0, 2 * np.pi, n, endpoint=False)

    for i, (_, row) in enumerate(df.iterrows()):
        ax.plot([0, np.cos(angles[i])], [0, np.sin(angles[i])],
                color=row['color'], linewidth=1)

    breaks = [0] + [i for i in range(1, n) if df.iloc[i]['group'] != df.iloc[i - 1]['group']]
    for j, start in enumerate(breaks):
        end = breaks[j + 1] if j + 1 < len(breaks) else n
        mid_idx = (start + end - 1) / 2.0
        mid_angle = np.interp(mid_idx, range(n), angles)

        grp = df.iloc[start]['group']
        if grp == 'st_unique_new':
            label = 'Unique unknown (new) STs'
        elif grp == 'st_unique_known':
            label = 'Unique known STs'
        elif grp == 'unresolved':
            label = 'Unresolved'
        else:
            label = grp

        r_label = 1.15
        cx, cy = np.cos(mid_angle), np.sin(mid_angle)
        ha = 'left' if cx >= 0 else 'right'
        ax.text(r_label * cx, r_label * cy, label,
                fontsize=5, ha=ha, va='center')

    if output_path is None:
        output_path = OUTPUT_DIR / 'fig_mlst_pie_lines.svg'
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, bbox_inches='tight', transparent=True)
    plt.close(fig)
    print(f"Saved: {output_path}")


def main():
    plot_mlst_pie()


if __name__ == '__main__':
    main()
