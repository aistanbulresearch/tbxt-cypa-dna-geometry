"""Build the editable manuscript from Marked's GFM token stream."""

import argparse
import html
import json
import re
from pathlib import Path

from docx import Document
from docx.enum.style import WD_STYLE_TYPE
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.opc.constants import RELATIONSHIP_TYPE
from docx.shared import Inches, Pt, RGBColor
from PIL import Image


CAPTION = re.compile(r"^(?:\*\*)?(Figure|Table)\s+S?\d+[.:]", re.I)
BODY_WIDTH = 6.78


def xml(tag, **attrs):
    element = OxmlElement(tag)
    for key, value in attrs.items():
        element.set(qn(f"w:{key}"), str(value))
    return element


def style(document, name, font, size, *, above=0, below=0, lines=1.2,
          bold=False, keep=False, centered=False, color="000000"):
    styles = document.styles
    spec = styles[name] if name in styles else styles.add_style(name, WD_STYLE_TYPE.PARAGRAPH)
    spec.font.name = font
    spec.font.size = Pt(size)
    spec.font.bold = bold
    spec.font.italic = False
    spec.font.color.rgb = RGBColor.from_string(color)
    fmt = spec.paragraph_format
    fmt.line_spacing = lines
    fmt.space_before = Pt(above)
    fmt.space_after = Pt(below)
    fmt.keep_with_next = keep
    fmt.keep_together = keep
    fmt.widow_control = True
    fmt.alignment = WD_ALIGN_PARAGRAPH.CENTER if centered else WD_ALIGN_PARAGRAPH.LEFT
    for node in spec.element.xpath("./w:pPr/w:pBdr"):
        node.getparent().remove(node)
    return spec


def design(document):
    style(document, "Normal", "Georgia", 11, below=7, lines=1.18, color="171717")
    style(document, "Title", "Arial", 24, below=11, lines=1.08, bold=True, keep=True)
    style(document, "Heading 1", "Arial", 15, above=16, below=6, bold=True, keep=True)
    style(document, "Heading 2", "Arial", 12, above=11, below=5, bold=True, keep=True)
    style(document, "Heading 3", "Arial", 11, above=9, below=4, bold=True, keep=True)
    style(document, "Metadata", "Arial", 10, below=11, keep=True)
    style(document, "Figure", "Arial", 9, above=9, lines=1, keep=True, centered=True)
    style(document, "Caption", "Arial", 9.5, above=5, below=10, color="222222")
    style(document, "Table Caption", "Arial", 9.5, above=8, below=4, keep=True, color="222222")
    style(document, "Table Cell", "Arial", 9, lines=1.11, color="171717")
    style(document, "References", "Georgia", 10, below=3, lines=1.08)
    ref_format = document.styles["References"].paragraph_format
    ref_format.left_indent = Pt(15)
    ref_format.first_line_indent = Pt(-15)
    ref_format.keep_together = True
    style(document, "Footer", "Arial", 9, color="4A4A4A")
    sect = document.sections[0]
    sect.page_width = Inches(8.5)
    sect.page_height = Inches(11)
    sect.left_margin = sect.right_margin = Inches(.86)
    sect.top_margin = sect.bottom_margin = Inches(.76)
    sect.footer_distance = Inches(.39)
    footer = sect.footer.paragraphs[0]
    footer.alignment = WD_ALIGN_PARAGRAPH.RIGHT
    number = xml("w:fldSimple", instr="PAGE")
    run = xml("w:r")
    value = xml("w:t")
    value.text = "1"
    run.append(value)
    number.append(run)
    footer._p.append(number)
    lang = xml("w:lang", val="en-US")
    document.styles["Normal"].element.get_or_add_rPr().append(lang)


def add_link(paragraph, token, bold, italic):
    href = html.unescape(token["href"])
    if not href.startswith(("https://", "http://", "mailto:")):
        add_inline(paragraph, token["tokens"], bold=bold, italic=italic)
        return
    link = OxmlElement("w:hyperlink")
    rel_id = paragraph.part.relate_to(href, RELATIONSHIP_TYPE.HYPERLINK, is_external=True)
    link.set(qn("r:id"), rel_id)
    start = len(paragraph._p)
    add_inline(paragraph, token["tokens"], bold=bold, italic=italic)
    for child in list(paragraph._p)[start:]:
        if child.tag != qn("w:r"):
            continue
        if child.rPr is None or child.rPr.find(qn("w:color")) is None:
            child.get_or_add_rPr().append(xml("w:color", val="234D77"))
        link.append(child)
    paragraph._p.append(link)


def add_inline(paragraph, tokens, *, bold=False, italic=False):
    for token in tokens:
        kind = token["type"]
        if kind in {"strong", "em", "del"}:
            start = len(paragraph.runs)
            add_inline(paragraph, token["tokens"], bold=bold or kind == "strong",
                       italic=italic or kind == "em")
            if kind == "del":
                for run in paragraph.runs[start:]:
                    run.font.strike = True
            continue
        if kind == "link":
            add_link(paragraph, token, bold, italic)
            continue
        if kind == "br":
            paragraph.add_run().add_break()
            continue
        if kind == "text" and token.get("tokens"):
            add_inline(paragraph, token["tokens"], bold=bold, italic=italic)
            continue
        value = html.unescape(token.get("text", token["raw"]))
        if kind == "html" and re.fullmatch(r"<(?:br\s*/?|/?sup|/?sub)>", value, re.I):
            if value.lower().startswith("<br"):
                paragraph.add_run().add_break()
            continue
        run = paragraph.add_run(value.replace("\n", " "))
        if bold:
            run.bold = True
        if italic:
            run.italic = True
        if kind == "codespan":
            run.font.name = "Consolas"


def add_picture(document, token, directory, separate):
    relative = Path(html.unescape(token["href"]))
    path = (directory / relative).resolve()
    if not path.is_relative_to(directory):
        raise ValueError(f"Image must be in manuscript directory: {relative}")
    with Image.open(path) as image:
        width_pixels, height_pixels = image.size
    max_width = 5.80 if relative.stem.startswith("figure1-") else BODY_WIDTH
    width = min(max_width, 5.05 * width_pixels / height_pixels)
    paragraph = document.add_paragraph(style="Figure")
    paragraph.paragraph_format.keep_with_next = bool(token["text"] or separate)
    shape = paragraph.add_run().add_picture(str(path), width=Inches(width))
    alt = html.unescape(token["text"])
    shape._inline.docPr.set("descr", alt or relative.stem.replace("-", " "))
    if alt and not separate:
        cap = document.add_paragraph(style="Caption")
        cap.paragraph_format.keep_together = True
        label = re.match(r"^(Figure\s+S?\d+[.:])\s+", alt, re.I)
        if label and all(item["type"] == "text" for item in token.get("tokens", [])):
            cap.add_run(label[1]).bold = True
            cap.add_run(" " + alt[label.end():])
        else:
            add_inline(cap, token.get("tokens") or [{"type": "text", "raw": alt}])


def table_widths(token):
    weights = [
        max(len(cell["text"]) for cell in (header, *(row[i] for row in token["rows"])))
        for i, header in enumerate(token["header"])
    ]
    proportions = [min(80, max(14, weight)) ** .62 for weight in weights]
    return [Inches(BODY_WIDTH * p / sum(proportions)) for p in proportions]


def add_table(document, token):
    table = document.add_table(rows=0, cols=len(token["header"]))
    table.autofit = False
    tbl = table._tbl.tblPr
    borders = xml("w:tblBorders")
    margins = xml("w:tblCellMar")
    for edge in ("top", "left", "bottom", "right", "insideH", "insideV"):
        borders.append(xml(f"w:{edge}", val="single", sz=5, color="D9D9D9"))
        if edge in {"top", "bottom", "left", "right"}:
            margins.append(xml(f"w:{edge}", w=88 if edge in {"top", "bottom"} else 104, type="dxa"))
    tbl.insert_element_before(borders, "w:tblLayout", "w:tblLook")
    tbl.insert_element_before(margins, "w:tblLook")
    widths = table_widths(token)
    for column, width in zip(table.columns, widths):
        column.width = width
    for row_index, row_tokens in enumerate([token["header"], *token["rows"]]):
        row = table.add_row()
        row_style = row._tr.get_or_add_trPr()
        row_style.append(xml("w:cantSplit"))
        if row_index == 0:
            row_style.append(xml("w:tblHeader", val="true"))
        for i, (cell, cell_token, width) in enumerate(zip(row.cells, row_tokens, widths)):
            cell.width = width
            cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
            paragraph = cell.paragraphs[0]
            paragraph.style = "Table Cell"
            paragraph.paragraph_format.keep_with_next = (
                row_index == 0 or (len(token["rows"]) <= 6 and row_index < len(token["rows"]))
            )
            if i < len(token["align"]) and token["align"][i] in {"right", "center"}:
                paragraph.alignment = (WD_ALIGN_PARAGRAPH.RIGHT if token["align"][i] == "right"
                                       else WD_ALIGN_PARAGRAPH.CENTER)
            if row_index == 0:
                cell._tc.get_or_add_tcPr().append(xml("w:shd", val="clear", fill="EDF1F5"))
            add_inline(paragraph, cell_token["tokens"], bold=row_index == 0)
    spacer = document.add_paragraph()
    spacer.paragraph_format.space_after = Pt(3)
    spacer.paragraph_format.line_spacing = Pt(2)
    spacer.add_run().font.size = Pt(2)


def add_list(document, token, reference):
    start = token.get("start") or 1
    for i, item in enumerate(token["items"], start=start):
        name = "References" if reference else "Normal"
        paragraph = document.add_paragraph(style=name)
        paragraph.paragraph_format.left_indent = Pt(23)
        paragraph.paragraph_format.first_line_indent = Pt(-23)
        paragraph.paragraph_format.tab_stops.add_tab_stop(Pt(23))
        paragraph.add_run((f"{i}." if token.get("ordered") else "•") + "\t")
        for j, child in enumerate(item["tokens"]):
            if child["type"] == "list":
                add_list(document, child, reference)
            elif child["type"] != "space":
                if j:
                    paragraph.add_run().add_break()
                add_inline(paragraph, child.get("tokens") or [child])


def render(tokens, directory):
    document = Document()
    design(document)
    document.core_properties.author = "Ozge A. Cavus, Aysegul Kuskucu"
    document.core_properties.last_modified_by = ""
    document.core_properties.comments = ""
    document.core_properties.language = "en-US"
    reference = False
    for i, token in enumerate(tokens):
        kind = token["type"]
        previous = tokens[i - 1] if i else {}
        following = tokens[i + 1] if i + 1 < len(tokens) else {}
        if kind == "heading":
            reference = token["text"].strip().lower() == "references"
            name = "Title" if token["depth"] == 1 else f"Heading {min(3, token['depth'] - 1)}"
            if name == "Title":
                document.core_properties.title = token["text"]
            add_inline(document.add_paragraph(style=name), token["tokens"])
        elif kind == "paragraph":
            inline = token["tokens"]
            if len(inline) == 1 and inline[0]["type"] == "image":
                match = CAPTION.match(following.get("text", ""))
                add_picture(document, inline[0], directory, bool(match and match[1].lower() == "figure"))
                continue
            name = "Normal"
            if CAPTION.match(token["text"]):
                name = "Table Caption" if following.get("type") == "table" else "Caption"
            elif reference:
                name = "References"
            elif previous.get("type") == "heading" and previous["depth"] == 1:
                name = "Metadata"
            paragraph = document.add_paragraph(style=name)
            if name == "Caption":
                paragraph.paragraph_format.keep_together = True
            add_inline(paragraph, inline)
        elif kind == "table":
            add_table(document, token)
        elif kind == "list":
            add_list(document, token, reference)
        elif kind in {"blockquote", "code"}:
            add_inline(document.add_paragraph(), token.get("tokens") or [token])
        elif kind != "hr":
            raise ValueError(f"Unsupported Markdown block {kind!r}")
    return document


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tokens", type=Path)
    parser.add_argument("manuscript", type=Path)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    items = json.loads(args.tokens.read_text(encoding="utf-8"))
    document = render(
        [item for item in items if item["type"] != "space"], args.manuscript.resolve().parent
    )
    document.save(args.destination)
