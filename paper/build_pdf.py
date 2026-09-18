"""Render manuscript.md to PDF with ReportLab.

Pandoc and a TeX toolchain are not installed here, and installing one to lay out
a markdown file would be a heavier dependency than the job needs. ReportLab is
already present, so this reads the subset of markdown the manuscript actually
uses -- headings, paragraphs, bullets, blockquotes, pipe tables, images with
captions, indented equation blocks, and inline bold/italic/code -- and lays it
out directly.

    python paper/build_pdf.py [manuscript.md] [-o out.pdf] [--single-column]

Two-column is the default because that is the shape a journal expects: the
title, abstract and index terms run full width on page one, the body flows in
two columns after it, and each page carries a running head and page number.

ponytail: no journal class file. Every venue ships its own LaTeX template with
its own margins and font stack, so matching one before the venue is chosen is
work that gets thrown away -- pick the venue, then port.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.platypus import (
    BaseDocTemplate,
    Frame,
    Image,
    NextPageTemplate,
    PageBreak,
    PageTemplate,
    Paragraph,
    Preformatted,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)

HERE = Path(__file__).resolve().parent


RUNNING_HEAD = "Atmosphere-Aware QUBO Formulations for THz Array Design"


def _rule_and_footer(canvas, doc, head: str = "") -> None:
    canvas.saveState()
    canvas.setStrokeColor(colors.HexColor("#1f4e79"))
    canvas.setLineWidth(0.5)
    if head:
        y = A4[1] - 1.1 * cm
        canvas.line(1.5 * cm, y, A4[0] - 1.5 * cm, y)
        canvas.setFont("Helvetica", 7)
        canvas.setFillColor(colors.HexColor("#555555"))
        canvas.drawString(1.5 * cm, y + 0.12 * cm, head)
    else:
        canvas.line(1.5 * cm, 1.25 * cm, A4[0] - 1.5 * cm, 1.25 * cm)
    canvas.setFont("Helvetica", 7)
    canvas.setFillColor(colors.HexColor("#555555"))
    canvas.drawRightString(A4[0] - 1.5 * cm, 0.85 * cm, f"Page {doc.page}")
    canvas.restoreState()


class TwoColumnPaper(BaseDocTemplate):
    """Full-width front matter on page one, two columns thereafter."""

    def __init__(self, filename: str, **kw):
        super().__init__(filename, pagesize=A4, leftMargin=1.5 * cm,
                         rightMargin=1.5 * cm, topMargin=1.5 * cm,
                         bottomMargin=1.6 * cm, **kw)
        gutter = 0.6 * cm
        col = (self.width - gutter) / 2.0
        self.column_width = col
        self.addPageTemplates([
            PageTemplate(id="front",
                         frames=[Frame(self.leftMargin, self.bottomMargin,
                                       self.width, self.height, id="full")],
                         onPage=lambda c, d: _rule_and_footer(c, d)),
            PageTemplate(id="body",
                         frames=[Frame(self.leftMargin, self.bottomMargin, col,
                                       self.height, id="c1"),
                                 Frame(self.leftMargin + col + gutter,
                                       self.bottomMargin, col, self.height, id="c2")],
                         onPage=lambda c, d: _rule_and_footer(c, d, RUNNING_HEAD)),
        ])


def styles():
    base = getSampleStyleSheet()
    s = {
        "title": ParagraphStyle("title", parent=base["Title"], fontSize=17, leading=21,
                                spaceAfter=10),
        "author": ParagraphStyle("author", parent=base["Normal"], fontSize=11,
                                 alignment=TA_CENTER, spaceAfter=14),
        "h1": ParagraphStyle("h1", parent=base["Heading1"], fontSize=13, leading=16,
                             spaceBefore=14, spaceAfter=6),
        "h2": ParagraphStyle("h2", parent=base["Heading2"], fontSize=11.5, leading=14,
                             spaceBefore=10, spaceAfter=4),
        "h3": ParagraphStyle("h3", parent=base["Heading3"], fontSize=10.5, leading=13,
                             spaceBefore=8, spaceAfter=3),
        "body": ParagraphStyle("body", parent=base["BodyText"], fontSize=9.6, leading=13.4,
                               alignment=TA_JUSTIFY, spaceAfter=6),
        "bullet": ParagraphStyle("bullet", parent=base["BodyText"], fontSize=9.6,
                                 leading=13.2, leftIndent=14, bulletIndent=4,
                                 spaceBefore=1, spaceAfter=4),
        "quote": ParagraphStyle("quote", parent=base["BodyText"], fontSize=8.8,
                                leading=12, leftIndent=14, rightIndent=10,
                                textColor=colors.HexColor("#444444"), spaceAfter=8),
        "caption": ParagraphStyle("caption", parent=base["BodyText"], fontSize=8.4,
                                  leading=11, alignment=TA_CENTER,
                                  textColor=colors.HexColor("#333333"), spaceAfter=10),
        "eq": ParagraphStyle("eq", parent=base["Code"], fontSize=8.8, leading=11.6,
                             leftIndent=18, spaceBefore=3, spaceAfter=7),
        "cell": ParagraphStyle("cell", parent=base["BodyText"], fontSize=8.1, leading=10),
        "cellh": ParagraphStyle("cellh", parent=base["BodyText"], fontSize=8.1,
                                leading=10, textColor=colors.white),
    }
    return s


def inline(text: str) -> str:
    """Markdown inline spans to ReportLab markup, escaping XML first."""
    t = (text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))
    t = re.sub(r"`([^`]+)`", r'<font face="Courier">\1</font>', t)
    t = re.sub(r"\*\*([^*]+)\*\*", r"<b>\1</b>", t)
    t = re.sub(r"(?<!\*)\*([^*]+)\*(?!\*)", r"<i>\1</i>", t)
    t = re.sub(r"\[([^\]]+)\]\((https?://[^)]+)\)", r'<link href="\2">\1</link>', t)
    t = re.sub(r"\[([^\]]+)\]\((?!https?://)[^)]+\)", r"\1", t)   # local links: keep text
    return t


def make_table(rows: list, st: dict, width: float) -> Table:
    """Pipe table to a Table, first row treated as the header."""
    head, *body = rows
    data = [[Paragraph(inline(c), st["cellh"]) for c in head]]
    data += [[Paragraph(inline(c), st["cell"]) for c in r] for r in body]
    n = max(len(r) for r in data)
    data = [r + [Paragraph("", st["cell"])] * (n - len(r)) for r in data]

    t = Table(data, colWidths=[width / n] * n, repeatRows=1, hAlign="LEFT")
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#33475b")),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("GRID", (0, 0), (-1, -1), 0.3, colors.HexColor("#b0b8c0")),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1),
         [colors.white, colors.HexColor("#f2f4f7")]),
        ("LEFTPADDING", (0, 0), (-1, -1), 4),
        ("RIGHTPADDING", (0, 0), (-1, -1), 4),
        ("TOPPADDING", (0, 0), (-1, -1), 3),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
    ]))
    return t


def figure(path: Path, caption: str, st: dict, width: float) -> list:
    """Scale an image to the frame width, preserving aspect ratio."""
    from reportlab.lib.utils import ImageReader

    if not path.exists():
        return [Paragraph(f"[missing figure: {path.name}]", st["caption"])]
    iw, ih = ImageReader(str(path)).getSize()
    w = min(width, 15.5 * cm)     # real aspect ratio from the file, never forced square
    return [Spacer(1, 4), Image(str(path), width=w, height=w * ih / iw),
            Paragraph(inline(caption), st["caption"])]


def is_separator(line: str) -> bool:
    return bool(re.fullmatch(r"\|[\s:|-]+\|", line.strip()))


def split_row(line: str) -> list:
    return [c.strip() for c in line.strip().strip("|").split("|")]


def build(md_path: Path, out_path: Path, two_column: bool = True) -> Path:
    st = styles()
    if two_column:
        doc = TwoColumnPaper(str(out_path), title=md_path.stem, author="Manasa Vamsi")
        width = doc.column_width
        # a 7.9 cm column needs smaller type than a 16.6 cm page
        for key, size, lead in (("body", 8.2, 11.0), ("bullet", 8.2, 11.2),
                                ("cell", 6.9, 8.6), ("cellh", 6.9, 8.6),
                                ("caption", 7.2, 9.4), ("eq", 7.6, 10.0),
                                ("quote", 7.6, 10.2)):
            st[key].fontSize, st[key].leading = size, lead
        st["h1"].fontSize, st["h1"].leading = 11.0, 13.0
        st["h2"].fontSize, st["h2"].leading = 9.4, 11.5
        st["h3"].fontSize, st["h3"].leading = 8.6, 10.6
    else:
        doc = SimpleDocTemplate(
            str(out_path), pagesize=A4,
            leftMargin=2.2 * cm, rightMargin=2.2 * cm,
            topMargin=2.0 * cm, bottomMargin=2.0 * cm,
            title=md_path.stem, author="Manasa Vamsi",
        )
        width = doc.width
    flow: list = []
    lines = md_path.read_text(encoding="utf-8").splitlines()

    i, para, quote, bullets, table = 0, [], [], [], []
    pending_break = False

    def flush():
        """Emit whatever block is open."""
        nonlocal para, quote, bullets, table
        if para:
            flow.append(Paragraph(inline(" ".join(para)), st["body"]))
            para = []
        if quote:
            flow.append(Paragraph(inline(" ".join(quote)), st["quote"]))
            quote = []
        if bullets:
            for b in bullets:
                flow.append(Paragraph(inline(b), st["bullet"], bulletText="•"))
            flow.append(Spacer(1, 5))
            bullets = []
        if table:
            flow.append(Spacer(1, 3))
            flow.append(make_table(table, st, width))
            flow.append(Spacer(1, 9))
            table = []

    while i < len(lines):
        raw = lines[i]
        line = raw.strip()

        if not line:
            flush()
            if pending_break:
                flow.extend([NextPageTemplate("body"), PageBreak()])
                pending_break = False
            i += 1
            continue

        # image with caption
        m = re.fullmatch(r"!\[(.*?)\]\((.+?)\)", line)
        if m:
            flush()
            flow.extend(figure((md_path.parent / m.group(2)).resolve(),
                               m.group(1), st, width))
            i += 1
            continue

        # pipe table
        if line.startswith("|") and i + 1 < len(lines) and is_separator(lines[i + 1]):
            flush()
            table = [split_row(line)]
            i += 2
            while i < len(lines) and lines[i].strip().startswith("|"):
                table.append(split_row(lines[i]))
                i += 1
            flush()
            continue

        # headings
        m = re.fullmatch(r"(#{1,4})\s+(.*)", line)
        if m:
            flush()
            level = len(m.group(1))
            key = {1: "title", 2: "h1", 3: "h2", 4: "h3"}[level]
            flow.append(Paragraph(inline(m.group(2)), st[key]))
            i += 1
            continue

        if line.startswith(">"):
            flush() if (para or bullets or table) else None
            quote.append(line.lstrip("> ").strip())
            i += 1
            continue

        if re.match(r"[-*]\s+", line):
            if para or table:
                flush()
            bullets.append(re.sub(r"^[-*]\s+", "", line))
            i += 1
            continue

        # indented block: equations / code
        if raw.startswith("    ") and not raw.strip().startswith("|"):
            flush()
            block = []
            while i < len(lines) and (lines[i].startswith("    ") or not lines[i].strip()):
                if lines[i].strip():
                    block.append(lines[i][4:].rstrip())
                elif block:
                    break
                i += 1
            # Paragraph, not Preformatted: a 7.9 cm column is narrower than
            # several of these formulae, and Preformatted clips rather than wraps.
            # escape only -- markdown emphasis must not run inside a formula,
            # or `2*pi*sigma` silently becomes `2<i>pi</i>sigma`.
            def esc(x: str) -> str:
                return (x.replace("&", "&amp;").replace("<", "&lt;")
                         .replace(">", "&gt;").replace("  ", "&nbsp;&nbsp;"))

            body = "<br/>".join(esc(b) for b in block)
            flow.append(Paragraph(f'<font face="Courier">{body}</font>', st["eq"]))
            continue

        # author line: the single bold line right after the title
        if line.startswith("**") and line.endswith("**") and len(flow) == 1:
            flow.append(Paragraph(inline(line.strip("*")), st["author"]))
            i += 1
            continue

        para.append(line)
        i += 1

        # front matter ends after the index terms, which wrap over more than one
        # source line -- so arm the break here and fire it at the blank line.
        if two_column and line.startswith("**Index Terms:**"):
            pending_break = True

    flush()
    doc.build(flow)
    return out_path


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("source", nargs="?", default=str(HERE / "manuscript.md"))
    ap.add_argument("-o", "--out", default=None)
    ap.add_argument("--single-column", action="store_true",
                    help="one wide column instead of the journal-style two")
    a = ap.parse_args()

    src = Path(a.source).resolve()
    out = Path(a.out).resolve() if a.out else src.with_suffix(".pdf")
    build(src, out, two_column=not a.single_column)
    print(f"wrote {out} ({out.stat().st_size / 1024:.0f} kB)")


if __name__ == "__main__":
    main()
