#!/usr/bin/env python3
"""Validate the publishable introduction + Chapters 1–27 lecture tree."""

from __future__ import annotations

import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SRC = ROOT / "src"
EXPECTED_NEW = {
    10: "research_design_architecture",
    11: "measurement_digital_traces",
    15: "sensitivity_partial_identification",
    18: "interference_marketplace_experiments",
}
PROHIBITED = {
    "dash": re.compile(r"——|—|–"),
    "Chinese in math text": re.compile(r"\\text\{[^}]*[一-鿿]"),
    "zero-padded prose chapter": re.compile(r"第\s*0[1-9]\s*章"),
    "external script direction": re.compile(r"完整脚本|本节代码|见\s*code/"),
}
LOCAL_LINK = re.compile(r"(?:href|src)=[\"']([^\"'#?]+)")
MD_IMAGE = re.compile(r"!\[[^]]*\]\(([^)\s]+)")
CHAPTER_REF = re.compile(r"Chapter\s+(\d+)")


def fail(errors: list[str], message: str) -> None:
    errors.append(message)


def numbered_sources() -> dict[int, Path]:
    result: dict[int, Path] = {}
    for path in SRC.glob("[0-9][0-9]_*.md"):
        number = int(path.name[:2])
        if number in result:
            raise ValueError(f"duplicate chapter number {number}: {result[number]}, {path}")
        result[number] = path
    return result


def validate_sources(errors: list[str]) -> None:
    intro = SRC / "introduction.md"
    if not intro.exists():
        fail(errors, "missing src/introduction.md")
    if (SRC / "00_overview_roadmap.md").exists() or list(SRC.glob("00_*.md")):
        fail(errors, "Chapter 00 source still exists")

    try:
        chapters = numbered_sources()
    except ValueError as exc:
        fail(errors, str(exc))
        chapters = {}
    missing = sorted(set(range(1, 28)) - set(chapters))
    extra = sorted(set(chapters) - set(range(1, 28)))
    if missing:
        fail(errors, f"missing numbered sources: {missing}")
    if extra:
        fail(errors, f"unexpected numbered sources: {extra}")

    for number, slug in EXPECTED_NEW.items():
        expected = SRC / f"{number:02d}_{slug}.md"
        if not expected.exists():
            fail(errors, f"missing required new chapter: {expected.relative_to(ROOT)}")

    publishable = [p for p in sorted(SRC.glob("*.md")) if p.name != "demo.md"]
    for path in publishable:
        text = path.read_text(encoding="utf-8")
        if re.match(r"\d{2}_", path.name):
            number = int(path.name[:2])
            for field in ("seriesline", "chapterlabel"):
                match = re.search(rf'^{field}:\s*"([^"]+)"', text, re.M)
                if not match or not re.search(rf"Chapter {number}(?!\d)", match.group(1)):
                    fail(errors, f"{path.relative_to(ROOT)}: {field} does not match Chapter {number}")
        for label, pattern in PROHIBITED.items():
            if match := pattern.search(text):
                line = text.count("\n", 0, match.start()) + 1
                fail(errors, f"{path.relative_to(ROOT)}:{line}: prohibited {label}")
        for ref in map(int, CHAPTER_REF.findall(text)):
            if not 1 <= ref <= 27:
                fail(errors, f"{path.relative_to(ROOT)}: invalid chapter reference {ref}")
        for image in MD_IMAGE.findall(text):
            if image.startswith(("http://", "https://", "data:")):
                continue
            target = (ROOT / image).resolve()
            if not target.exists():
                fail(errors, f"{path.relative_to(ROOT)}: missing image {image}")
        if path.name != "introduction.md" and path.name != "27_quick_reference.md":
            if "::: {.apa-refs}" not in text:
                fail(errors, f"{path.relative_to(ROOT)}: missing APA reference block")


def validate_html(errors: list[str]) -> None:
    expected = [ROOT / "introduction.html"]
    expected.extend(ROOT / f"{n:02d}_{numbered_sources().get(n, Path('missing')).stem[3:]}.html" for n in range(1, 28) if n in numbered_sources())
    for path in expected:
        if not path.exists():
            fail(errors, f"missing rendered HTML: {path.name}")

    for path in [ROOT / "index.html", ROOT / "assets" / "template.html"]:
        if not path.exists():
            fail(errors, f"missing shared page: {path.relative_to(ROOT)}")
            continue
        text = path.read_text(encoding="utf-8")
        if "00_overview_roadmap" in text:
            fail(errors, f"{path.relative_to(ROOT)}: stale Chapter 00 link")
        required_links = ["introduction.html"] + [
            f"{n:02d}_{numbered_sources()[n].stem[3:]}.html" for n in range(1, 28)
        ]
        for required in required_links:
            if f'href="{required}"' not in text:
                fail(errors, f"{path.relative_to(ROOT)}: missing navigation link {required}")
        for link in LOCAL_LINK.findall(text):
            if link.startswith(("http://", "https://", "mailto:", "javascript:")):
                continue
            if "$" in link:
                continue
            # Both template and index are rendered/served from the book root.
            target = (ROOT / link).resolve()
            if not target.exists():
                fail(errors, f"{path.relative_to(ROOT)}: broken local link {link}")


def main() -> int:
    errors: list[str] = []
    validate_sources(errors)
    validate_html(errors)
    if errors:
        print(f"VALIDATION FAILED ({len(errors)} errors)")
        for error in errors:
            print(f"- {error}")
        return 1
    print("VALIDATION OK: introduction + Chapters 1–27")
    return 0


if __name__ == "__main__":
    sys.exit(main())
