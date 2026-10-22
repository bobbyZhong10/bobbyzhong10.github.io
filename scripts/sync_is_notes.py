#!/usr/bin/env python3
"""Replay the website's presentation layer on top of the English IS-economics notes.

The notes are written and rendered in iCloud (notes_v3_en). The site copy under
notes/is-economics-methods/ adds a presentation layer: site YAML fields, a template with the theme
script, versioned stylesheet, topbar and notes.js, a restyled cover page, the site's own style.css,
a notes.js whose chapter order follows the book, and static redirect pages for the v2 file names.
This script rebuilds that layer from the notes instead of hand-editing html.

  build    copy a notes workspace into --out, apply the site layer, render every page with the
           book's pandoc command, then write the redirect pages
  install  mirror the publish set of a build into the collection folder, deleting everything else
  check    link-check the collection inside a Jekyll _site; exits 1 when it finds errors

Run from the repository root, for example:
  python3 scripts/sync_is_notes.py build --src "$V3E" --site notes/is-economics-methods \\
      --redirects scripts/is_notes_redirects.tsv --version notes-ui-20261021 --out ~/v3site_build
  python3 scripts/sync_is_notes.py install --build ~/v3site_build --dest notes/is-economics-methods
  python3 scripts/sync_is_notes.py check --site-root _site --collection notes/is-economics-methods \\
      --redirects scripts/is_notes_redirects.tsv
"""

from __future__ import annotations

import argparse
import csv
import filecmp
import html
import json
import os
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit


# Site layer, verbatim from the pre-release site (commit 6964e91).
SERIES = "Foundations of Information Systems Economics: Empirical Methods for Information Systems Economics"
SIDEBARTITLE = "Foundations of Information Systems Economics"
SIDEBARSUB = "Theory · Methods · Applications"
INTRO_BODYCLASS = "notes-introduction"
THEME_INIT = "<script>(function(){var t=null;try{t=localStorage.getItem('theme');}catch(e){}if(t!=='dark'&&t!=='light'){t=(window.matchMedia&&window.matchMedia('(prefers-color-scheme: dark)').matches)?'dark':'light';}document.documentElement.setAttribute('data-theme',t);})();</script>"  # noqa: E501
TOPBAR = """<div class="notes-topbar">
  <a class="notes-home" href="https://bobbyzhong10.github.io/">← Yunzhou (Bobby) Zhong</a>
  <div class="notes-topbar__right">
    <a href="https://bobbyzhong10.github.io/notes/">Research Notes</a>
    <button id="theme-toggle" class="notes-theme-toggle" type="button" aria-label="Toggle theme"><span></span></button>
  </div>
</div>"""
BACKLINK = '<a class="site-back-link" href="/notes/"><span aria-hidden="true">←</span> Back to Research Notes</a>'

# The cover page and the redirect pages carry the topbar on one line, as the site's index.html does.
TOPBAR_INLINE = re.sub(r">\s+<", "><", TOPBAR)
STYLESHEET = '  <link rel="stylesheet" href="assets/style.css" />'
PANDOC_ARGS = ["-f", "markdown+autolink_bare_uris", "--template=assets/template.html", "--katex=assets/katex/",
               "--lua-filter=assets/wrap-tables.lua", "--toc", "--toc-depth=2"]
PANDOC_VERSION = "3.12.1"
VERSION = re.compile(r"notes-ui-\d{8}[a-z]?")  # a letter marks a second release on the same day
VERSION_QUERY = re.compile(r"\?v=notes-ui-\d{8}[a-z]?")

BUILD_MARKER = ".sync_is_notes_build"  # "incomplete" while build runs, "complete" once it succeeded
BUILD_COMPLETE = "complete"
PUBLISH_ASSETS = ("fig", "katex", "style.css", "notes.js", "template.html", "wrap-tables.lua")
JUNK = {".DS_Store", "__pycache__"}
V2_CHAPTERS = 27
TOOLING_DIRS = ("src", "code", "tools")
TOOLING_FILES = ("template.html", "wrap-tables.lua")
TOOLING_SUFFIXES = (".py", ".tsv")
PAGE_COUNT = 37  # cover, introduction and 35 chapters

REDIRECT_COLUMNS = ["old_file", "kind", "targets", "note"]
PAGE_NAME = re.compile(r"\d\d_[a-z0-9_]+\.html")
FRONT_MATTER = re.compile(r"---\n(.*?\n)---\n", re.S)
STYLE_BLOCK = re.compile(r"^[ \t]*<style>.*?</style>[ \t]*\n", re.M | re.S)
COVER_SUB_EN = re.compile(r'^[ \t]*<div class="cover-sub-en">.*</div>[ \t]*\n', re.M)
COVER_FOOTER = re.compile(r'<footer class="cover-footer">.*?</footer>', re.S)
CHAPTER_ORDER = re.compile(r"(var CHAPTER_ORDER = \[\n).*?(\n[ \t]*\];)", re.S)
FILE_NAV = re.compile(r'<div class="file-nav">(.*?)</nav>', re.S)
FILE_NAV_LINK = re.compile(r'<a href="(\d\d_[a-z0-9_]+\.html)">([^<]*)</a>')
EXTERNAL = re.compile(r"^(?:[a-zA-Z][a-zA-Z0-9+.-]*:|//)")
SITEMAP_LOC = re.compile(r"<loc>\s*([^<]+?)\s*</loc>")
CANONICAL = re.compile(r'<link rel="canonical" href="([^"]*)"')
NOINDEX = '<meta name="robots" content="noindex"'


class BuildError(Exception):
    """A problem with the inputs or the target folder; reported without a traceback."""


def read(path: Path) -> str:
    with open(path, encoding="utf-8", newline="") as handle:
        return handle.read()


def write(path: Path, text: str) -> None:
    with open(path, "w", encoding="utf-8", newline="") as handle:
        handle.write(text)


def html_escape(text: str) -> str:
    return html.escape(text, quote=False)


def _replace_once(text: str, old: str, new: str, where: str) -> str:
    count = text.count(old)
    if count != 1:
        raise ValueError(f"{where}: expected exactly one {old.strip()!r}, found {count}")
    return text.replace(old, new)


def _sub_once(pattern: re.Pattern, new: str, text: str, where: str) -> str:
    out, count = pattern.subn(lambda _match: new, text)
    if count != 1:
        raise ValueError(f"{where}: expected exactly one match of {pattern.pattern!r}, found {count}")
    return out


def _stylesheet(version: str) -> str:
    return f'  <link rel="stylesheet" href="assets/style.css?v={version}" />'


def _notes_script(version: str) -> str:
    return f'<script defer src="assets/notes.js?v={version}"></script>'


def site_yaml(text: str, intro: bool = False) -> str:
    """Set the three site fields in a chapter's front matter; the introduction also gets a body class."""
    match = FRONT_MATTER.match(text)
    if not match:
        raise ValueError("no YAML front matter at the top of the file")
    head = match.group(1)
    for key, value in (("series", SERIES), ("sidebartitle", SIDEBARTITLE), ("sidebarsub", SIDEBARSUB)):
        head = _sub_once(re.compile(rf"^{key}: .*$", re.M), f'{key}: "{value}"', head, "front matter")
    if intro:
        head = re.sub(r"^bodyclass: .*\n", "", head, flags=re.M) + f'bodyclass: "{INTRO_BODYCLASS}"\n'
    return f"---\n{head}---\n" + text[match.end():]


def site_template(v3_template: str, version: str) -> str:
    """Add the theme script, versioned stylesheet, body class, topbar and notes.js to the book's template."""
    out = _replace_once(v3_template, STYLESHEET + "\n", f"{THEME_INIT}\n{_stylesheet(version)}\n", "template")
    out = _replace_once(out, "<body>\n", '<body$if(bodyclass)$ class="$bodyclass$"$endif$>\n' + TOPBAR + "\n",
                        "template")
    return _replace_once(out, "</body>\n", _notes_script(version) + "\n</body>\n", "template")


def site_index(v3_index: str, version: str) -> str:
    """Restyle the book's cover page for the site (its cover styles live in the site's style.css)."""
    out = _sub_once(STYLE_BLOCK, "", v3_index, "index.html")
    out = _replace_once(out, STYLESHEET + "\n", f"{THEME_INIT}\n{_stylesheet(version)}\n", "index.html")
    out = _replace_once(out, "<body>\n", f"<body>\n{TOPBAR_INLINE}\n", "index.html")
    out = _replace_once(out, '<div class="cover-wrapper">\n', f'<div class="cover-wrapper">\n  {BACKLINK}\n',
                        "index.html")
    out = _sub_once(COVER_SUB_EN, "", out, "index.html")
    out = _sub_once(COVER_FOOTER, f'<footer class="cover-footer">{html_escape(SERIES)}</footer>', out, "index.html")
    return _replace_once(out, "</body>\n", _notes_script(version) + "\n</body>\n", "index.html")


def _js_string(text: str) -> str:
    return "'" + text.replace("\\", "\\\\").replace("'", "\\'") + "'"


def notes_js(old_js: str, toc_rows: list, version: str) -> str:
    """Regenerate CHAPTER_ORDER (cover, introduction, then the toc) and bump any versioned asset URL."""
    entries = [("index.html", "Cover"), ("introduction.html", "Introduction · How to Use")]
    entries += [(f"{row['number']}_{row['slug']}.html", f"{row['number']} {row['en_short']}") for row in toc_rows]
    array = ",\n".join(f"    {{ file: {_js_string(name)}, label: {_js_string(label)} }}" for name, label in entries)
    out, count = CHAPTER_ORDER.subn(lambda match: match.group(1) + array + match.group(2), old_js)
    if count != 1:
        raise ValueError(f"notes.js: expected exactly one CHAPTER_ORDER array, found {count}")
    return VERSION_QUERY.sub(f"?v={version}", out)


def _and_join(items: list) -> str:
    return items[0] if len(items) == 1 else ", ".join(items[:-1]) + " and " + items[-1]


def _link_chapters(text: str, titles: dict) -> str:
    """Escape a note and link each "Chapter N" it names to that chapter's page, when there is one."""
    pages = {int(name[:2]): name for name in titles}
    return re.sub(r"Chapter (\d+)", lambda match: f'<a href="{pages[int(match.group(1))]}">{match.group(0)}</a>'
                  if int(match.group(1)) in pages else match.group(0), html_escape(text))


def stub_html(row: dict, titles: dict, version: str, *, canonical_base: str, descs: dict | None = None,
              edition: int = 35) -> str:
    """The static page left at a v2 file name; its canonical link is canonical_base + the first target.

    A one-to-one row redirects at once and keeps the #anchor. A split row is a landing page
    (decision 3): it lists every new chapter with its description, does not redirect and carries noindex.
    """
    targets = row["targets"]
    first = targets[0]
    old_number = int(row["old_file"][:2])

    def chapter(name: str) -> str:
        return f"Chapter {int(name[:2])} · {titles[name]}"

    head = ['<!DOCTYPE html>', '<html lang="en">', '<head>', '  <meta charset="utf-8" />',
            '  <meta name="viewport" content="width=device-width, initial-scale=1.0" />']
    canonical = f'  <link rel="canonical" href="{canonical_base}{first}" />'
    sub = f"This page was Chapter {old_number} of the {V2_CHAPTERS}-chapter edition."
    note = [f"    <p>{_link_chapters(row['note'], titles)}</p>"] if row["note"] else []
    if row["kind"] == "one":
        head += [f"  <title>{html_escape(chapter(first))} | {SIDEBARTITLE}</title>",
                 canonical,
                 f"  <script>location.replace({json.dumps(first)} + location.hash)</script>",
                 f'  <meta http-equiv="refresh" content="0; url={first}" />']
        heading = "This chapter has moved"
        body = ['  <div class="cover-intro">',
                f'    <p>This chapter is now <a href="{first}">{html_escape(chapter(first))}</a> '
                f"in the {edition}-chapter edition.</p>",
                *note,
                "  </div>"]
    else:
        missing = [name for name in targets if not (descs or {}).get(name)]
        if missing:
            raise ValueError(f"{row['old_file']}: a landing page needs a description for {', '.join(missing)}")
        numbers = _and_join([str(int(name[:2])) for name in targets])
        head += ['  <meta name="robots" content="noindex" />',
                 f"  <title>Chapter {old_number} is now Chapters {numbers} | {SIDEBARTITLE}</title>",
                 canonical]
        heading = "This chapter has been split"
        cards = [f'    <a class="chapter-card" href="{name}"><div class="ch-no">{name[:2]}</div>'
                 f'<div class="ch-title">{html_escape(titles[name])}</div>'
                 f'<div class="ch-desc">{html_escape(descs[name])}</div></a>' for name in targets]
        body = ['  <div class="cover-intro">',
                f"    <p>In the {edition}-chapter edition its material is spread over Chapters {numbers}. "
                "Pick the one you were looking for:</p>",
                "  </div>",
                '  <div class="chapter-grid">', *cards, "  </div>",
                *(['  <div class="cover-intro">', *note, "  </div>"] if note else [])]
    head += [THEME_INIT, _stylesheet(version), "</head>"]
    page = [*head, "<body>", TOPBAR_INLINE, '<div class="cover-wrapper">', '  <header class="cover-header">',
            f"    <h1>{heading}</h1>", f'    <div class="cover-sub">{sub}</div>', "  </header>", *body,
            f'  <footer class="cover-footer"><a href="index.html">{html_escape(SERIES)}</a></footer>', "</div>",
            _notes_script(version), "</body>", "</html>"]
    return "\n".join(page) + "\n"


def load_redirects(path: Path) -> list:
    """Rows of the redirect table: old_file, kind (one | split), targets (first is canonical), note."""
    with open(path, encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t", quoting=csv.QUOTE_NONE)
        if next(reader, None) != REDIRECT_COLUMNS:
            raise ValueError(f"{path}: the header must be {' / '.join(REDIRECT_COLUMNS)}, tab-separated")
        rows, seen = [], set()
        for line, fields in enumerate(reader, start=2):
            if not any(field.strip() for field in fields):
                continue
            if len(fields) not in (3, 4):
                raise ValueError(f"{path}:{line}: expected 3 or 4 tab-separated fields, found {len(fields)}")
            old, kind, targets = (field.strip() for field in fields[:3])
            targets = [name.strip() for name in targets.split(",") if name.strip()]
            note = fields[3].strip() if len(fields) == 4 else ""
            bad = [name for name in [old, *targets] if not PAGE_NAME.fullmatch(name)]
            if bad:
                raise ValueError(f"{path}:{line}: not a chapter page name: {', '.join(bad)}")
            if kind not in ("one", "split"):
                raise ValueError(f"{path}:{line}: kind must be one or split, not {kind!r}")
            if kind == "one" and len(targets) != 1:
                raise ValueError(f"{path}:{line}: a one-to-one row needs exactly one target")
            if kind == "split" and len(targets) < 2:
                raise ValueError(f"{path}:{line}: a split row needs at least two targets")
            if old in seen:
                raise ValueError(f"{path}:{line}: {old} appears twice")
            seen.add(old)
            rows.append({"old_file": old, "kind": kind, "targets": targets, "note": note})
    return rows


def file_nav(template: str) -> list:
    """(file, label) for each chapter link in the template's sidebar."""
    block = FILE_NAV.search(template)
    if not block:
        raise ValueError('template: no <div class="file-nav"> block')
    return [(name, html.unescape(label)) for name, label in FILE_NAV_LINK.findall(block.group(1))]


def load_toc(workspace: Path, template: str) -> list:
    """Chapter rows from toc_v3.tsv, which must agree with the template's sidebar.

    The frozen v2 notes have no toc file; for them the sidebar itself is the chapter list.
    """
    nav = file_nav(template)
    toc = workspace / "toc_v3.tsv"
    if not toc.is_file():
        return [{"number": name[:2], "slug": name[3:-5], "en_short": label[3:], "en_desc": "", "status": "done"}
                for name, label in nav]
    with open(toc, encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t", quoting=csv.QUOTE_NONE))
    unfinished = [row["number"] for row in rows if row["status"] not in ("done", "legacy")]
    if unfinished:
        raise BuildError(f"toc_v3.tsv: chapters not ready to publish: {', '.join(unfinished)}")
    expected = [(f"{row['number']}_{row['slug']}.html", f"{row['number']} {row['en_short']}") for row in rows]
    if nav != expected:
        raise BuildError("the template's file-nav does not match toc_v3.tsv; regenerate the template first")
    return rows


def yaml_field(text: str, key: str) -> str:
    match = FRONT_MATTER.match(text)
    field = match and re.search(rf"^{key}: (.*)$", match.group(1), re.M)
    if not field:
        raise ValueError(f"front matter has no {key}")
    value = field.group(1).strip()
    if value.startswith('"'):
        return json.loads(value)
    if value.startswith("'"):
        return value[1:-1].replace("''", "'")
    return value


def canonical_base(site: Path) -> str:
    """The collection's absolute URL: url and baseurl from the site's _config.yml, then the folder's path."""
    site = site.resolve()
    config = next((folder / "_config.yml" for folder in site.parents if (folder / "_config.yml").is_file()), None)
    if config is None:
        raise BuildError(f"no _config.yml above {site}; the redirect pages need the site's url")
    text = read(config)

    def value(key: str) -> str:
        match = re.search(rf"^{key}[ \t]*:[ \t]*[\"']?([^\"'#\s]*)", text, re.M)
        return match.group(1) if match else ""

    url = value("url").rstrip("/")
    if not url:
        raise BuildError(f"{config} has no top-level url")
    path = "/".join(part for part in (value("baseurl").strip("/"), site.relative_to(config.parent).as_posix()) if part)
    return f"{url}/{path}/"


def _inside(path: Path, root: Path) -> bool:
    return path == root or root in path.parents


def _prepare_out(out: Path) -> None:
    if out.exists():
        if not out.is_dir():
            raise BuildError(f"--out {out} is not a folder")
        if any(out.iterdir()) and not (out / BUILD_MARKER).is_file():
            raise BuildError(f"--out {out} is not empty and is not an earlier build; refusing to delete it")
        shutil.rmtree(out)
    out.mkdir(parents=True)
    write(out / BUILD_MARKER, "incomplete\n")


def _copy_ignore(src: Path):
    def ignore(directory: str, names: list) -> set:
        here = Path(directory)
        skip = {name for name in names if name in JUNK}
        if here in (src, src / "src"):
            skip |= {name for name in names if name.startswith("demo")}
        if here == src:
            skip |= {"RENUMBERED_V3"} & set(names)
        return skip
    return ignore


def _pandoc(out: Path, stem: str) -> list:
    command = ["pandoc", f"src/{stem}.md", *PANDOC_ARGS, "-o", f"{stem}.html"]
    result = subprocess.run(command, cwd=out, capture_output=True, text=True)
    if result.returncode != 0:
        raise BuildError(f"pandoc failed on src/{stem}.md:\n{result.stderr}")
    return [f"{stem}: {line}" for line in result.stderr.splitlines()
            if line.strip() and not line.startswith("[WARNING] Deprecated: --katex")]


def _check_pandoc() -> None:
    if not shutil.which("pandoc"):
        raise BuildError("pandoc not found on PATH")
    result = subprocess.run(["pandoc", "--version"], capture_output=True, text=True)
    if result.returncode != 0:
        raise BuildError(f"pandoc --version failed:\n{result.stderr}")
    version = result.stdout.split("\n")[0]
    if version != f"pandoc {PANDOC_VERSION}":
        print(f"warning: {version}; the golden test was recorded with pandoc {PANDOC_VERSION}", file=sys.stderr)


def render(out: Path, stems: list) -> None:
    """The book's pandoc command for every page, then inline_figures.py, as tools/render_all.sh does."""
    with ThreadPoolExecutor(max_workers=os.cpu_count() or 4) as pool:
        warnings = [line for lines in pool.map(lambda stem: _pandoc(out, stem), stems) for line in lines]
    for line in warnings:
        print(line, file=sys.stderr)
    if subprocess.run([sys.executable, "inline_figures.py"], cwd=out).returncode != 0:
        raise BuildError("inline_figures.py failed")
    print(f"rendered {len(stems)} pages, {len(warnings)} pandoc warning lines besides the --katex notice")


def build(src: Path, site: Path, redirects: Path, version: str, out: Path) -> None:
    src, site, out = src.resolve(), site.resolve(), out.expanduser().resolve()
    if not VERSION.fullmatch(version):
        raise BuildError(f"--version must look like notes-ui-YYYYMMDD, optionally with one letter, not {version!r}")
    for required in (src / "index.html", src / "assets/template.html", src / "src/introduction.md",
                     src / "inline_figures.py", site / "assets/style.css", site / "assets/notes.js"):
        if not required.is_file():
            raise BuildError(f"missing {required}")
    if _inside(out, src) or _inside(src, out) or _inside(out, site) or _inside(site, out):
        raise BuildError("--out must be a separate folder, outside --src and --site")
    rows = load_redirects(redirects)
    base = canonical_base(site) if rows else ""
    _check_pandoc()

    _prepare_out(out)
    shutil.copytree(src, out, ignore=_copy_ignore(src), copy_function=shutil.copyfile, dirs_exist_ok=True)

    for md in sorted((out / "src").glob("*.md")):
        write(md, site_yaml(read(md), intro=md.stem == "introduction"))
    template = read(out / "assets/template.html")
    toc = load_toc(out, template)
    write(out / "assets/template.html", site_template(template, version))

    stems = ["introduction"] + [f"{row['number']}_{row['slug']}" for row in toc]
    render(out, stems)
    write(out / "index.html", site_index(read(out / "index.html"), version))
    write(out / "assets/notes.js", notes_js(read(site / "assets/notes.js"), toc, version))
    shutil.copyfile(site / "assets/style.css", out / "assets/style.css")

    pages = {"index.html", *(f"{stem}.html" for stem in stems)}
    titles = {f"{stem}.html": yaml_field(read(out / "src" / f"{stem}.md"), "title") for stem in stems[1:]}
    descs = {f"{row['number']}_{row['slug']}.html": row["en_desc"] for row in toc}
    for row in rows:
        if row["old_file"] in pages:
            raise BuildError(f"redirect {row['old_file']} would overwrite a page of this edition")
        unknown = [name for name in row["targets"] if name not in titles]
        if unknown:
            raise BuildError(f"redirect {row['old_file']}: no such chapter {', '.join(unknown)}")
        write(out / row["old_file"], stub_html(row, titles, version, canonical_base=base, descs=descs,
                                               edition=len(toc)))

    leftover = sorted({path.name for path in out.glob("*.html")} - pages - {row["old_file"] for row in rows})
    if leftover:
        raise BuildError(f"html files in --src that are neither toc pages nor redirects: {', '.join(leftover)}")
    write(out / BUILD_MARKER, f"{BUILD_COMPLETE}\n")
    print(f"built {out}: {len(pages)} pages and {len(rows)} redirect pages, version {version}")


def _files(root: Path, base: Path) -> list:
    return [path.relative_to(base).as_posix() for path in root.rglob("*")
            if path.is_file() and not JUNK & set(path.relative_to(base).parts)]


def publish_set(build_dir: Path) -> set:
    """Pages, assets, sources, code and inline_figures.py; never tools/, toc_v3.tsv or validate_book.py."""
    wanted = {path.name for path in build_dir.glob("*.html") if not path.name.startswith("demo")}
    wanted |= {f"src/{path.name}" for path in (build_dir / "src").glob("*.md") if not path.name.startswith("demo")}
    for name in PUBLISH_ASSETS:
        path = build_dir / "assets" / name
        wanted |= set(_files(path, build_dir)) if path.is_dir() else {f"assets/{name}"}
    wanted |= set(_files(build_dir / "code", build_dir)) | {"inline_figures.py"}
    missing = sorted(name for name in wanted | {"index.html", "introduction.html"} if not (build_dir / name).is_file())
    if missing:
        raise BuildError(f"the build lacks {', '.join(missing)}")
    return wanted


def install(build_dir: Path, dest: Path) -> None:
    build_dir, dest = build_dir.resolve(), dest.resolve()
    if _inside(dest, build_dir) or _inside(build_dir, dest):
        raise BuildError("--build and --dest must not contain each other")
    marker = build_dir / BUILD_MARKER
    if not marker.is_file():
        raise BuildError(f"{build_dir} is not a folder written by `build` (no {BUILD_MARKER})")
    if read(marker).strip() != BUILD_COMPLETE:
        raise BuildError(f"the build in {build_dir} did not finish; run build again")
    wanted = publish_set(build_dir)
    cover = dest / "index.html"
    if dest.exists() and any(dest.iterdir()) and not (cover.is_file() and f"<h1>{SIDEBARTITLE}</h1>" in read(cover)):
        raise BuildError(f"{dest} is not empty and its index.html is not this book's cover; "
                         "refusing to delete anything in it")
    dest.mkdir(parents=True, exist_ok=True)

    removed = copied = 0
    for path in sorted(dest.rglob("*"), reverse=True):
        if path.is_dir() and not path.is_symlink():
            if not any(path.iterdir()):
                path.rmdir()
        elif path.relative_to(dest).as_posix() not in wanted:
            path.unlink()
            removed += 1
    for name in sorted(wanted):
        source, target = build_dir / name, dest / name
        if target.is_file() and filecmp.cmp(source, target, shallow=False):
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        copied += 1
    print(f"installed {len(wanted)} files into {dest}: {copied} written, {removed} removed")


class _PageParser(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.links = []
        self.ids = set()

    def handle_starttag(self, tag, attrs):
        for name, value in attrs:
            if name == "id" or (tag == "a" and name == "name"):
                self.ids.add(value)
            elif name in ("href", "src") and value:
                self.links.append(value)

    handle_startendtag = handle_starttag


def _parse(path: Path) -> _PageParser:
    parser = _PageParser()
    parser.feed(read(path))
    parser.close()
    return parser


def _tooling(rel: Path) -> bool:
    """True for paths that _config.yml must keep out of _site: sources, code, tools and build files."""
    return rel.parts[0] in TOOLING_DIRS or rel.suffix in TOOLING_SUFFIXES or rel.name in TOOLING_FILES


def check_site(site_root: Path, collection: str, redirects: Path) -> list:
    """Errors in the built collection: broken links and anchors, redirects, published tooling, sitemap."""
    site_root = Path(site_root).resolve()
    collection = collection.strip("/")
    coll = site_root / collection
    if not coll.is_dir():
        return [f"{collection}: no such folder under {site_root}"]
    rows = load_redirects(redirects)
    errors = []

    parsed = {}

    def page(path: Path) -> _PageParser:
        if path not in parsed:
            parsed[path] = _parse(path)
        return parsed[path]

    for path in sorted(coll.rglob("*")):
        rel = path.relative_to(coll)
        if path.is_dir() and len(rel.parts) == 1 and rel.name in TOOLING_DIRS:
            errors.append(f"{collection}/{rel.as_posix()}/: must not be published")
        elif path.is_file() and rel.parts[0] not in TOOLING_DIRS and _tooling(rel):
            errors.append(f"{collection}/{rel.as_posix()}: must not be published")

    for source in sorted(coll.rglob("*.html")):
        if _tooling(source.relative_to(coll)):
            continue
        name = source.relative_to(site_root).as_posix()
        for link in page(source).links:
            if EXTERNAL.match(link):
                continue
            parts = urlsplit(link)
            target = source
            if parts.path:
                path = unquote(parts.path)
                target = Path(os.path.normpath(site_root / path.lstrip("/") if path.startswith("/")
                                               else source.parent / path))
                if target.is_dir():
                    target = target / "index.html"
                if not _inside(target, site_root) or not target.is_file():
                    errors.append(f"{name}: broken link {link}")
                    continue
            fragment = unquote(parts.fragment)
            if fragment and target.suffix == ".html" and fragment not in page(target).ids:
                errors.append(f"{name}: missing anchor {link}")

    sitemap = site_root / "sitemap.xml"
    locs = SITEMAP_LOC.findall(read(sitemap)) if sitemap.is_file() else []
    prefix = f"/{collection}/"
    # The collection's absolute URL, as the site itself publishes it in sitemap.xml.
    bases = {loc[:loc.index(prefix) + len(prefix)] for loc in locs if urlsplit(loc).path.startswith(prefix)}
    base = bases.pop() if len(bases) == 1 else None

    stubs = {row["old_file"] for row in rows}
    for row in rows:
        stub = coll / row["old_file"]
        if not stub.is_file():
            errors.append(f"{collection}/{row['old_file']}: redirect page missing")
            continue
        text, first, name = read(stub), row["targets"][0], f"{collection}/{row['old_file']}"
        canonical = CANONICAL.search(text)
        if not canonical:
            errors.append(f"{name}: not a redirect page (no canonical link)")
        elif (not urlsplit(canonical.group(1)).netloc or not canonical.group(1).endswith(prefix + first)
              or (base and canonical.group(1) != base + first)):
            errors.append(f"{name}: the canonical link must be the absolute URL {base or '...' + prefix}{first}")
        if row["kind"] == "split":
            if NOINDEX not in text:
                errors.append(f"{name}: a landing page needs noindex (decision 3)")
            if "http-equiv" in text or "location.replace" in text:
                errors.append(f"{name}: a landing page must not redirect (decision 3)")
        else:
            if NOINDEX in text:
                errors.append(f"{name}: a one-to-one redirect must not carry noindex")
            if (f"location.replace({json.dumps(first)} + location.hash)" not in text
                    or f'content="0; url={first}"' not in text):
                errors.append(f"{name}: must redirect to {first} at once and keep the anchor (location.hash)")
        for name in row["targets"]:
            if not (coll / name).is_file():
                errors.append(f"{collection}/{row['old_file']}: redirect target {name} missing")

    if not sitemap.is_file():
        errors.append("sitemap.xml missing")
    else:
        listed = [path[len(prefix):] or "index.html" for path in (unquote(urlsplit(loc).path) for loc in locs)
                  if path.startswith(prefix)]
        pages = {path.name for path in coll.glob("*.html")} - stubs
        if len(listed) != PAGE_COUNT:
            errors.append(f"sitemap.xml lists {len(listed)} pages under {collection}, expected {PAGE_COUNT}")
        for name in sorted(set(listed) & stubs):
            errors.append(f"sitemap.xml lists redirect page {name}")
        for name in sorted(set(listed) - pages - stubs):
            errors.append(f"sitemap.xml lists {name}, which is not a page of {collection}")
        for name in sorted(pages - set(listed)):
            errors.append(f"sitemap.xml does not list {collection}/{name}")
    return errors


def main(argv: list | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    commands = parser.add_subparsers(dest="command", required=True)
    build_args = commands.add_parser("build", help="render a notes workspace with the site layer")
    build_args.add_argument("--src", type=Path, required=True, help="notes workspace, e.g. notes_v3_en")
    build_args.add_argument("--site", type=Path, required=True, help="the site's collection folder (style.css, notes.js)")
    build_args.add_argument("--redirects", type=Path, required=True, help="redirect table (TSV)")
    build_args.add_argument("--version", required=True, help="cache-busting version, notes-ui-YYYYMMDD")
    build_args.add_argument("--out", type=Path, required=True, help="build folder; an earlier build is replaced")
    install_args = commands.add_parser("install", help="mirror a build's publish set into the collection folder")
    install_args.add_argument("--build", type=Path, required=True)
    install_args.add_argument("--dest", type=Path, required=True)
    check_args = commands.add_parser("check", help="link-check the collection in a Jekyll _site")
    check_args.add_argument("--site-root", type=Path, required=True)
    check_args.add_argument("--collection", required=True)
    check_args.add_argument("--redirects", type=Path, required=True)
    args = parser.parse_args(argv)

    try:
        if args.command == "build":
            build(args.src, args.site, args.redirects, args.version, args.out)
        elif args.command == "install":
            install(args.build, args.dest)
        else:
            errors = check_site(args.site_root, args.collection, args.redirects)
            for error in errors:
                print(error)
            print(f"check: {len(errors)} error(s) in {args.collection}")
            return 1 if errors else 0
    except (BuildError, ValueError, OSError) as error:
        print(f"sync_is_notes.py {args.command}: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
