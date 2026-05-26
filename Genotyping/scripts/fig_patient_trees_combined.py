"""Combine patient tree SVGs into a PDF using maximal-rectangles bin packing."""

from pathlib import Path
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import inch
from reportlab.pdfgen import canvas
from reportlab.graphics import renderPDF
from reportlab.graphics.shapes import Group, String
from svglib.svglib import svg2rlg

PX_TO_PT = 72 / 96
FONTTYPE = "Helvetica"


def _fix_font_sizes(node):
    """svglib reads SVG px font sizes as pt; scale them back to correct size."""
    if isinstance(node, String) and hasattr(node, 'fontSize'):
        node.fontSize *= 1 / PX_TO_PT
    if isinstance(node, Group):
        for child in node.contents:
            _fix_font_sizes(child)

INPUT_DIR = Path(__file__).resolve().parent.parent / "figures" / "patient_trees"
OUTPUT_PDF = INPUT_DIR.parent / "patient_trees_combined.pdf"
PAGE_W, PAGE_H = A4
MARGIN = 0.25 * inch
TITLE_H = 0.18 * inch
GAP = 6
LEGEND_SVG = INPUT_DIR / "tree_legend.svg"
LEGEND_H = 0.45 * inch
USABLE_W = PAGE_W - 2 * MARGIN
USABLE_H = PAGE_H - 2 * MARGIN - LEGEND_H

MAX_IMG_FRACTION = 0.73


class MaxRectsBin:
    """Maximal rectangles bin packing with best-short-side-fit heuristic."""

    def __init__(self, w, h):
        self.w = w
        self.h = h
        self.free = [(0, 0, w, h)]  # (x, y, w, h)  y=0 is top

    def insert(self, w, h):
        best = None
        best_ssf = float('inf')
        for i, (fx, fy, fw, fh) in enumerate(self.free):
            if w <= fw + 0.5 and h <= fh + 0.5:
                ssf = min(fw - w, fh - h)
                if ssf < best_ssf:
                    best_ssf = ssf
                    best = (fx, fy, i)
        if best is None:
            return None
        px, py, _ = best
        self._split_and_prune(px, py, w, h)
        return (px, py)

    def _split_and_prune(self, px, py, w, h):
        new_free = []
        for (fx, fy, fw, fh) in self.free:
            # if no overlap with placed rect, keep as is
            if px >= fx + fw or px + w <= fx or py >= fy + fh or py + h <= fy:
                new_free.append((fx, fy, fw, fh))
                continue
            # split into up to 4 maximal rects around the placed rect
            if px > fx:
                new_free.append((fx, fy, px - fx, fh))
            if px + w < fx + fw:
                new_free.append((px + w, fy, fx + fw - px - w, fh))
            if py > fy:
                new_free.append((fx, fy, fw, py - fy))
            if py + h < fy + fh:
                new_free.append((fx, py + h, fw, fy + fh - py - h))

        # remove rects fully contained in another
        self.free = []
        for i, a in enumerate(new_free):
            contained = False
            for j, b in enumerate(new_free):
                if i != j and _contains(b, a):
                    contained = True
                    break
            if not contained:
                self.free.append(a)


def _contains(outer, inner):
    return (inner[0] >= outer[0] and inner[1] >= outer[1] and
            inner[0] + inner[2] <= outer[0] + outer[2] and
            inner[1] + inner[3] <= outer[1] + outer[3])


def load_sizes(svgs):
    sizes = {}
    for p in svgs:
        drawing = svg2rlg(str(p))
        sizes[p] = (drawing.width, drawing.height)
    return sizes


def main():
    svgs = list(INPUT_DIR.glob("circular_tree_*.svg"))
    sizes = load_sizes(svgs)
    svgs.sort(key=lambda p: sizes[p][0] * sizes[p][1], reverse=True)
    print(f"Found {len(svgs)} images")

    max_w = max(w for w, _ in sizes.values())
    scale = (USABLE_W * MAX_IMG_FRACTION) / max_w
    print(f"Uniform scale: {scale:.4f}")

    # pack into pages
    pages = []  # list of [(svg_path, x, y, img_w, img_h), ...]
    remaining = list(svgs)

    while remaining:
        b = MaxRectsBin(USABLE_W, USABLE_H)
        page_items = []
        not_placed = []

        for img_path in remaining:
            w_px, h_px = sizes[img_path]
            img_w = w_px * scale + GAP
            img_h = h_px * scale + TITLE_H + GAP
            pos = b.insert(img_w, img_h)
            if pos is not None:
                bx, by = pos
                page_items.append((
                    img_path,
                    MARGIN + bx,
                    PAGE_H - MARGIN - by - img_h + GAP / 2,
                    w_px * scale,
                    h_px * scale
                ))
            else:
                not_placed.append(img_path)

        pages.append(page_items)
        if not_placed and len(not_placed) == len(remaining):
            # can't place even one — force it on a new page alone
            img_path = not_placed.pop(0)
            w_px, h_px = sizes[img_path]
            pages.append([(img_path, MARGIN, MARGIN, w_px * scale, h_px * scale)])
        remaining = not_placed

    print(f"Pages: {len(pages)}")

    legend_drawing = svg2rlg(str(LEGEND_SVG)) if LEGEND_SVG.exists() else None
    if legend_drawing:
        _fix_font_sizes(legend_drawing)

    c = canvas.Canvas(str(OUTPUT_PDF), pagesize=A4)
    for i, page_items in enumerate(pages):
        if i > 0:
            c.showPage()
        for svg_path, x, y, img_w, img_h in page_items:
            patient_label = svg_path.stem.replace("circular_tree_", "")
            c.setFont(f"{FONTTYPE}-Bold", 7)
            c.drawCentredString(x + img_w / 2, y + img_h + 3, patient_label)
            drawing = svg2rlg(str(svg_path))
            _fix_font_sizes(drawing)
            sx = img_w / drawing.width
            sy = img_h / drawing.height
            drawing.width = img_w
            drawing.height = img_h
            drawing.scale(sx, sy)
            renderPDF.draw(drawing, c, x, y)

        if legend_drawing:
            ld = svg2rlg(str(LEGEND_SVG))
            _fix_font_sizes(ld)
            ls = LEGEND_H / ld.height
            ld.width *= ls
            ld.height = LEGEND_H
            ld.scale(ls, ls)
            lx = MARGIN
            ly = MARGIN * 0.5 + 0.5 * inch / 2.54
            renderPDF.draw(ld, c, lx, ly)

    c.save()
    print(f"Saved → {OUTPUT_PDF}")


if __name__ == "__main__":
    main()
