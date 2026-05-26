"""
Unrooted radial phylogenetic tree from StrainPhlAn via iTOL batch API.
Leaves labeled/colored by MetaMLST ST category.

Samples: metagenome is_rep=False from trim_stats.csv.

Usage:
    python scripts/fig_strainphlan_tree.py                 # upload + export
    python scripts/fig_strainphlan_tree.py --tree-id ID    # re-export existing
"""

import argparse
import re
import tempfile
import zipfile
import time
import requests
import pandas as pd
from pathlib import Path
from PIL import Image
from Bio import Phylo
from io import StringIO
from collections import Counter

from plot_utils import load_color

BASE_DIR = Path(__file__).resolve().parent.parent

TREE_FILE = BASE_DIR / "output/strainphlan/tree/RAxML_bestTree.t__SGB10068.StrainPhlAn4.tre"
METAMLST_REPORT = BASE_DIR / "output/metamlst/merged/escherichia1_report.txt"
METAMLST_REPORT_E2 = BASE_DIR / "output/metamlst/merged/escherichia2_report.txt"
TRIM_STATS = BASE_DIR / "output/summary/trim_stats.csv"
OUTPUT_DIR = BASE_DIR / "figures"

ITOL_UPLOAD_URL = "https://itol.embl.de/batch_uploader.cgi"
ITOL_EXPORT_URL = "https://itol.embl.de/batch_downloader.cgi"

SVG_SIZE_INCHES = (6, 7)
SVG_DPI = 96

NEW_ST_THRESHOLD = 100000
E2_ST_OFFSET = 200000
INCLUDE_REFERENCE = False
REFERENCE_ID = 'NZ_OW967975.1'

ST_CAT_NAMES = ['st_unique_known', 'st_unique_new', 'st_multiple',
                'st_none', 'st_multiple1', 'st_unresolved']
CAT_COLORS = {name: load_color(name) for name in ST_CAT_NAMES}

TREE_LINESTYLE = {
    'style': 'solid',               # 'solid' (branches are always solid in iTOL)
    'line_width': 3.5,
    'color': '#000000',
}

EXTENSION_LINESTYLE = {
    'align_labels': 1,
    'style': 'dashed',              # 'dashed' or 'solid'
    'line_width': 1.5,
    'color': '#888888',
}

FONTTYPE = 'Arial'

SAMPLE_LABEL_TEXT_STYLE = {
    'font_size': 35,
    'font_name': FONTTYPE,
    'font_style': 'bold',
    'color': '#000000',
}

SAMPLE_TYPE_TEXT_STYLE = {
    'position': -1,
    'font_size': 60,
    'font_style': 'bold',
    'font_style_unresolved': 'normal',
}

ITOL_EXPORT_PARAMS = {
    # layout
    'display_mode': 3,              # 1=rect, 2=circular, 3=unrooted
    'arc': 360,                     # circular/unrooted arc (degrees)
    'rotation': 0,                  # circular rotation (degrees)
    'unrooted_rotation': 0,         # unrooted rotation (degrees)
    'inverted': 0,                  # 1=inverted (circular/rect only)
    'ignore_branch_length': 0,      # 1=cladogram
    # scaling
    'horizontal_scale_factor': 1,
    'vertical_shift_factor': 1,
    # branches
    'line_width': TREE_LINESTYLE['line_width'],
    'default_branch_color': TREE_LINESTYLE['color'],
    # labels
    'label_display': 1,             # 0=hide, 1=show
    'align_labels': EXTENSION_LINESTYLE['align_labels'],
    'dashed_lines': 1 if EXTENSION_LINESTYLE['style'] == 'dashed' else 0,
    'current_font_size': SAMPLE_LABEL_TEXT_STYLE['font_size'],
    'current_font_name': SAMPLE_LABEL_TEXT_STYLE['font_name'],
    'current_font_style': SAMPLE_LABEL_TEXT_STYLE['font_style'],
    'default_label_color': SAMPLE_LABEL_TEXT_STYLE['color'],
    # sorting
    'leaf_sorting': 1,              # 1=sort, 2=original order
    # bootstrap / metadata
    'bootstrap_display': 0,         # 0=hide, 1=show
    # branch lengths
    'branchlength_display': 0,      # 0=hide, 1=show on branches
    # internal scale
    'internal_scale': 0,            # 0=off, 1=on
    # datasets
    'datasets_visible': '0',        # show first dataset (ST labels)
}


def get_api_key():
    key_file = BASE_DIR / '.itol_api_key'
    if key_file.exists():
        return key_file.read_text().strip()
    import os
    return os.environ.get('ITOL_API_KEY', '').strip()


# ── Data ──────────────────────────────────────────────────────────────────

def _read_report(path):
    st_map = {}
    with open(path) as f:
        next(f)
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                st_map[parts[2]] = int(parts[0])
    return st_map


def load_metamlst_st():
    st_map = _read_report(METAMLST_REPORT)
    if METAMLST_REPORT_E2.exists():
        for sample, st in _read_report(METAMLST_REPORT_E2).items():
            if sample not in st_map:
                st_map[sample] = st + E2_ST_OFFSET
    return st_map


def get_valid_samples():
    df = pd.read_csv(TRIM_STATS)
    valid = df[~df['is_isolate'] & ~df['is_rep']]
    return valid['run_name'], valid.set_index('run_name')['sample_name'].to_dict()


def classify_sts(st_map, leaves):
    leaf_sts = {l: st_map[l] for l in leaves if l in st_map}
    st_counts = Counter(leaf_sts.values())
    cats = {}
    for leaf in leaves:
        st = leaf_sts.get(leaf)
        if st is None:
            cats[leaf] = 'st_none'
        elif st_counts[st] > 1:
            cats[leaf] = 'st_multiple'
        elif (st % E2_ST_OFFSET) >= NEW_ST_THRESHOLD:
            cats[leaf] = 'st_unique_new'
        else:
            cats[leaf] = 'st_unique_known'
    return cats, leaf_sts


def make_st_label(name, leaf_sts):
    st = leaf_sts.get(name)
    if st is None:
        return None
    if st >= E2_ST_OFFSET:
        raw = st - E2_ST_OFFSET
        return f"ST NEW (P)" if raw >= NEW_ST_THRESHOLD else f"ST {raw} (P)"
    return "ST NEW" if st >= NEW_ST_THRESHOLD else f"ST {st}"


def prune_tree(tree, keep):
    for name in [t.name for t in tree.get_terminals() if t.name not in keep]:
        tree.prune(name)


def load_tree_data():
    """Single entry point: load tree, filter samples, classify STs.

    Returns (tree, remaining, cats, leaf_sts).
    """
    st_map = load_metamlst_st()
    valid, sample_name_map = get_valid_samples()
    tree = Phylo.read(str(TREE_FILE), 'newick')
    keep = {t.name for t in tree.get_terminals()} & set(valid)
    if INCLUDE_REFERENCE:
        keep |= {REFERENCE_ID}
    prune_tree(tree, keep)
    remaining = sorted(t.name for t in tree.get_terminals())
    cats, leaf_sts = classify_sts(st_map, remaining)
    return tree, remaining, cats, leaf_sts, sample_name_map


# ── iTOL annotation files ────────────────────────────────────────────────

def write_tree_file(path, tree):
    buf = StringIO()
    Phylo.write(tree, buf, 'newick')
    path.write_text(buf.getvalue())


def write_labels(path, leaves, sample_name_map):
    """Sample-name labels close to the leaves, in black."""
    lines = ["LABELS", "SEPARATOR COMMA", "DATA"]
    for name in sorted(leaves):
        lines.append(f"{name},{sample_name_map.get(name, name)}")
    path.write_text("\n".join(lines) + "\n")


def write_st_dataset(path, leaves, cats, leaf_sts):
    """ST labels as a colored text dataset (rendered outside the Fxxx labels)."""
    lines = [
        "DATASET_TEXT",
        "SEPARATOR COMMA",
        "DATASET_LABEL,ST types",
        "COLOR,#000000",
        "SIZE_FACTOR,1.0",
        "DATA",
    ]
    pos = SAMPLE_TYPE_TEXT_STYLE['position']
    fs = SAMPLE_TYPE_TEXT_STYLE['font_style']
    fs_unres = SAMPLE_TYPE_TEXT_STYLE['font_style_unresolved']
    sf = SAMPLE_TYPE_TEXT_STYLE['font_size'] / SAMPLE_LABEL_TEXT_STYLE['font_size']
    for name in sorted(leaves):
        st_label = make_st_label(name, leaf_sts)
        if st_label is not None:
            cat = cats.get(name, 'st_none')
            color = CAT_COLORS['st_multiple1'] if cat == 'st_multiple' else CAT_COLORS[cat]
            lines.append(f"{name},{st_label},{pos},{color},{fs},{sf}")
        elif cats.get(name) == 'st_none':
            lines.append(f"{name},unresolved,{pos},{CAT_COLORS['st_unresolved']},{fs_unres},{sf}")
    path.write_text("\n".join(lines) + "\n")


# ── iTOL API ──────────────────────────────────────────────────────────────

def upload_tree(zip_path, api_key, project='Sample project'):
    with open(zip_path, 'rb') as f:
        resp = requests.post(ITOL_UPLOAD_URL, files={
            'zipFile': ('upload.zip', f, 'application/zip'),
        }, data={
            'APIkey': api_key,
            'projectName': project,
            'treeName': 'StrainPhlAn_Ecoli_SGB10068',
        })
    text = resp.text.strip()
    print(f"Upload response: {text}")
    for line in text.split('\n'):
        if line.startswith('SUCCESS:'):
            return line.split('SUCCESS:')[1].strip()
    raise RuntimeError(f"Upload failed: {text}")


def png_composite_white(path: Path) -> None:
    """Flatten iTOL PNG transparency onto white background."""
    img = Image.open(path)
    if img.mode == 'P':
        img = img.convert('RGBA')
    if img.mode == 'RGBA':
        bg = Image.new('RGB', img.size, (255, 255, 255))
        bg.paste(img, mask=img.split()[3])
        bg.save(path, 'PNG')
    elif img.mode == 'LA':
        rgba = img.convert('RGBA')
        bg = Image.new('RGB', rgba.size, (255, 255, 255))
        bg.paste(rgba, mask=rgba.split()[3])
        bg.save(path, 'PNG')


_DASHED_PATH_RE = re.compile(
    r'(<path\s[^>]*stroke-dasharray=")(\d+,\d+)("[^>]*>)'
)


def autofit_svg_viewbox(svg_text, padding=20):
    """Recompute viewBox by rendering to bitmap and finding content bounds."""
    import cairosvg
    from PIL import Image as _Image
    from io import BytesIO

    big_vb = "-2000,-2000,6000,6000"
    big_svg = re.sub(r'width="[^"]*"', 'width="6000"', svg_text, count=1)
    big_svg = re.sub(r'height="[^"]*"', 'height="6000"', big_svg, count=1)
    big_svg = re.sub(r'viewBox="[^"]*"', f'viewBox="{big_vb}"', big_svg, count=1)

    _Image.MAX_IMAGE_PIXELS = None
    png_data = cairosvg.svg2png(bytestring=big_svg.encode())
    img = _Image.open(BytesIO(png_data)).convert('RGBA')
    bbox = img.getbbox()
    if bbox is None:
        return svg_text

    vb_x0, vb_y0, vb_w_total, vb_h_total = -2000, -2000, 6000, 6000
    px_w, px_h = img.size
    x_min = vb_x0 + bbox[0] / px_w * vb_w_total - padding
    y_min = vb_y0 + bbox[1] / px_h * vb_h_total - padding
    x_max = vb_x0 + bbox[2] / px_w * vb_w_total + padding
    y_max = vb_y0 + bbox[3] / px_h * vb_h_total + padding
    vb_w, vb_h = x_max - x_min, y_max - y_min

    w_in, h_in = SVG_SIZE_INCHES
    svg_text = re.sub(r'width="[^"]*"', f'width="{w_in}in"', svg_text, count=1)
    svg_text = re.sub(r'height="[^"]*"', f'height="{h_in}in"', svg_text, count=1)
    svg_text = re.sub(
        r'viewBox="[^"]*"',
        f'viewBox="{x_min:.2f},{y_min:.2f},{vb_w:.2f},{vb_h:.2f}"',
        svg_text, count=1)
    return svg_text


def restyle_extension_lines(svg_text):
    """Patch extension-line attributes in the exported SVG.

    The iTOL batch API only supports dashed_lines=0/1.  Width and color
    of the extension lines are controlled here via EXTENSION_LINESTYLE.
    """
    color = EXTENSION_LINESTYLE['color']
    width = EXTENSION_LINESTYLE['line_width']
    style = EXTENSION_LINESTYLE['style']

    def _patch(m):
        tag = m.group(0)
        tag = re.sub(r'stroke="[^"]*"', f'stroke="{color}"', tag, count=1)
        tag = re.sub(r'stroke-width="[^"]*"', f'stroke-width="{width}"', tag, count=1)
        if style != 'dashed':
            tag = re.sub(r'stroke-dasharray="[^"]*"',
                         'stroke-dasharray="none"', tag, count=1)
        return tag

    return _DASHED_PATH_RE.sub(_patch, svg_text)


def export_tree(tree_id, fmt, output_path, use_saved_view=False, **overrides):
    params = {'tree': tree_id, 'format': fmt}
    if not use_saved_view:
        params.update(ITOL_EXPORT_PARAMS)
    params.update(overrides)
    resp = requests.post(ITOL_EXPORT_URL, data=params)
    if 'text/html' in resp.headers.get('Content-Type', ''):
        raise RuntimeError(f"Export error: {resp.text}")
    output_path.write_bytes(resp.content)
    print(f"Exported: {output_path} ({len(resp.content)} bytes)")


def export_patched_svg(tree_id, output_path, use_saved_view=False, **overrides):
    """Export SVG from iTOL with extension lines restyled."""
    params = {'tree': tree_id, 'format': 'svg'}
    if not use_saved_view:
        params.update(ITOL_EXPORT_PARAMS)
    params.update(overrides)
    resp = requests.post(ITOL_EXPORT_URL, data=params)
    if 'text/html' in resp.headers.get('Content-Type', ''):
        raise RuntimeError(f"Export error: {resp.text}")
    svg_text = restyle_extension_lines(resp.text)
    svg_text = autofit_svg_viewbox(svg_text)
    output_path.write_text(svg_text)
    print(f"Exported: {output_path}")


# ── Main ──────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--tree-id', help='Re-export existing iTOL tree')
    parser.add_argument('--saved-view', action='store_true',
                        help='Use saved default view from iTOL browser (tweak there first)')
    args = parser.parse_args()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    if args.tree_id:
        tree_id = args.tree_id
    else:
        api_key = get_api_key()
        if not api_key:
            raise SystemExit(
                "No API key. Set ITOL_API_KEY env var or create .itol_api_key")

        tree, remaining, cats, leaf_sts, sample_name_map = load_tree_data()
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            tree_path = tmpdir / "tree.tree"
            labels_path = tmpdir / "01_labels.txt"
            st_path = tmpdir / "02_st_labels.txt"

            write_tree_file(tree_path, tree)
            write_labels(labels_path, remaining, sample_name_map)
            write_st_dataset(st_path, remaining, cats, leaf_sts)

            zip_path = tmpdir / "upload.zip"
            with zipfile.ZipFile(zip_path, 'w') as zf:
                for p in [tree_path, labels_path, st_path]:
                    zf.write(p, p.name)

            tree_id = upload_tree(zip_path, api_key)

        print(f"Tree ID: {tree_id}")
        print(f"View: https://itol.embl.de/tree/{tree_id}")
        time.sleep(2)

    export_patched_svg(tree_id, OUTPUT_DIR / "fig_strainphlan_tree.svg",
                       use_saved_view=args.saved_view)

    for fmt, ext in [('pdf', '.pdf'), ('png', '.png')]:
        export_tree(tree_id, fmt, OUTPUT_DIR / f"fig_strainphlan_tree{ext}",
                    use_saved_view=args.saved_view)


if __name__ == '__main__':
    main()
