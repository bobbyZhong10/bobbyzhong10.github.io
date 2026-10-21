#!/usr/bin/env python3
"""Tests for scripts/sync_is_notes.py.

Run from the repository root:  python3 -m unittest discover -s scripts/tests

The golden tests rebuild the frozen v2 English notes (notes_en) with the site layer and compare the
result with the pre-release site, read from commit 6964e91 with `git show` / `git ls-tree`, never
from the working tree. The v3 tests build notes_v3_en. Both sources live in iCloud; IS_NOTES_EN and
IS_NOTES_V3_EN override their paths. A test whose source or commit is missing is skipped.
"""

import csv
import hashlib
import os
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SCRIPTS = REPO / "scripts"
sys.path.insert(0, str(SCRIPTS))
import sync_is_notes as sync  # noqa: E402


NOTES_ROOT = (Path.home() / "Library/Mobile Documents/com~apple~CloudDocs/Research/Method"
              / "Micro Theory_Microeconometrics")
NOTES_EN = Path(os.environ.get("IS_NOTES_EN", NOTES_ROOT / "notes_en"))
NOTES_V3_EN = Path(os.environ.get("IS_NOTES_V3_EN", NOTES_ROOT / "notes_v3_en"))
GOLDEN_COMMIT = "6964e91"
COLLECTION = "notes/is-economics-methods"
OLD_VERSION = "notes-ui-20260822"
NEW_VERSION = "notes-ui-20261021"
REDIRECTS = SCRIPTS / "is_notes_redirects.tsv"

# Appendix C of the Phase 4 plan: v2 page -> (kind, v3 targets, first target is canonical).
APPENDIX_C = {
    "05_estimation_mle_gmm.html": ("split", ["05_m_estimation_nonlinear_regression.html",
                                             "06_maximum_likelihood.html",
                                             "07_gmm_minimum_distance.html",
                                             "12_binary_multinomial_response.html"]),
    "06_inference.html": ("one", ["08_inference.html"]),
    "07_panel_dynamic_panel.html": ("split", ["10_linear_panel_fe_re.html", "11_dynamic_panel.html"]),
    "08_nonparametrics.html": ("one", ["14_nonparametrics.html"]),
    "09_potential_outcomes.html": ("one", ["15_potential_outcomes.html"]),
    "10_research_design_architecture.html": ("one", ["16_research_design_architecture.html"]),
    "11_measurement_digital_traces.html": ("one", ["17_measurement_digital_traces.html"]),
    "12_iv_late_mte.html": ("one", ["18_iv_late_mte.html"]),
    "13_did_event_studies.html": ("one", ["19_did_event_studies.html"]),
    "14_rdd_synthetic_control.html": ("split", ["21_regression_discontinuity.html",
                                                "22_synthetic_control_factor_models.html"]),
    "15_sensitivity_partial_identification.html": ("one", ["23_sensitivity_partial_identification.html"]),
    "16_causal_ml_dml.html": ("one", ["24_causal_ml_dml.html"]),
    "17_digital_experimentation.html": ("one", ["25_digital_experimentation.html"]),
    "18_interference_marketplace_experiments.html": ("one", ["26_interference_marketplace_experiments.html"]),
    "19_discrete_choice.html": ("one", ["27_discrete_choice.html"]),
    "20_demand_blp.html": ("one", ["28_demand_blp.html"]),
    "21_supply_markups_mergers.html": ("one", ["29_supply_markups_mergers.html"]),
    "22_entry_partial_id.html": ("one", ["30_entry_partial_id.html"]),
    "23_dynamic_structural.html": ("one", ["31_dynamic_structural.html"]),
    "24_platforms_networks.html": ("one", ["32_platforms_networks.html"]),
    "25_search_advertising.html": ("one", ["33_search_advertising.html"]),
    "26_ai_algorithms_data.html": ("one", ["34_ai_algorithms_data.html"]),
    "27_quick_reference.html": ("one", ["35_quick_reference.html"]),
}

# Golden normalizations, applied in this order to both sides before an html page is compared. Each
# entry is (id, what it absorbs, one concrete example, pattern, replacement). Nothing else differs.
NORMALIZATIONS = [
    ("W1", "whitespace that contains a line break, between two tags, is dropped",
     'online topbar is one line `</a><div class="notes-topbar__right">`; the template writes '
     '`</a>\\n  <div class="notes-topbar__right">`',
     r">[ \t]*\n\s*<", "><"),
    ("W2", "every other whitespace run folds to one space",
     'pandoc 3.12.1 wraps the longer site sidebar title at 72 columns: '
     '`Foundations of Information Systems\\nEconomics</div>`',
     r"\s+", " "),
    ("P1", "pandoc 3.12.1 trims the line breaks that open and close a multi-line $$ block",
     '09_potential_outcomes: online `<span class="math display">\\n\\underbrace{...}_{\\text{selection bias}}.\\n'
     '</span>` vs `<span class="math display">\\underbrace{...}_{\\text{selection bias}}.</span>`',
     r'(<span class="math display">) ?([^<]*?) ?</span>', r"\1\2</span>"),
    ("P2", "skylighting re-tokenized R numbers and signs (dv/fl/sc spans; the code text is unchanged)",
     '01_consumer_producer_welfare: online `,.<span class="dv">4</span>` vs `,<span class="fl">.4</span>`; '
     '20_demand_blp: online `<span class="fl">1e-12</span>` vs '
     '`<span class="dv">1</span>e<span class="sc">-</span><span class="dv">12</span>`',
     r'<span class="(?:dv|fl|sc)">([^<]*)</span>', r"\1"),
    ("P3", "pandoc 3.12.1 rounds the total relative width of a multiline table to 99%",
     '05_estimation_mle_gmm: online `<table style="width:100%;">` vs `<table style="width:99%;">`',
     r'<table style="width:(?:99|100)%;">', '<table style="width:100%;">'),
    ("S1", "stale online <title>: 16 pages were rendered with the old sidebartitle and only the visible "
     "sidebar was hand-edited; the replay takes the title from the site's own src YAML",
     '03_information_mechanism_design: online `| IS Economics · Empirical Methods</title>` vs '
     '`| Foundations of Information Systems Economics</title>`',
     r"\| IS Economics · Empirical Methods</title>", "| Foundations of Information Systems Economics</title>"),
]

# The frozen notes_en still ships KaTeX 0.16.22; commit 6964e91 upgraded the site copy to 0.19.0
# (notes_v3_en ships 0.19.0 as well). These three files are the only non-overlay bytes that differ.
KATEX_UPGRADED = {"assets/katex/katex.min.js", "assets/katex/katex.min.css", "assets/katex/auto-render.min.js"}

# Files that the site layer rewrites; everything else in the publish set is copied byte for byte.
OVERLAY_SUFFIXES = (".html", ".md")
OVERLAY_FILES = {"assets/notes.js", "assets/template.html", "assets/style.css"}


def norm(text: str, skip: str = None) -> str:
    for rule, _what, _example, pattern, replacement in NORMALIZATIONS:
        if rule != skip:
            text = re.sub(pattern, replacement, text)
    return text.strip()


def git(*args: str, binary: bool = False):
    result = subprocess.run(["git", "-C", str(REPO), *args], check=True, capture_output=True)
    return result.stdout if binary else result.stdout.decode("utf-8")


def have_golden_commit() -> bool:
    return subprocess.run(["git", "-C", str(REPO), "cat-file", "-e", f"{GOLDEN_COMMIT}^{{commit}}"],
                          capture_output=True).returncode == 0


def online(path: str, binary: bool = False):
    return git("show", f"{GOLDEN_COMMIT}:{COLLECTION}/{path}", binary=binary)


def online_tree() -> dict:
    """Published path -> git blob id, for the collection at the golden commit."""
    tree = {}
    for line in git("ls-tree", "-r", GOLDEN_COMMIT, COLLECTION).splitlines():
        meta, path = line.split("\t", 1)
        tree[path[len(COLLECTION) + 1:]] = meta.split()[2]
    return tree


def blob_id(path: Path) -> str:
    data = path.read_bytes()
    return hashlib.sha1(b"blob %d\0" % len(data) + data).hexdigest()


def files_under(root: Path) -> set:
    return {p.relative_to(root).as_posix() for p in root.rglob("*") if p.is_file()}


def make_online_site(root: Path) -> Path:
    """The two site-owned files `build --site` reads, taken from the golden commit."""
    for name in ("assets/style.css", "assets/notes.js"):
        target = root / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(online(name, binary=True))
    return root


def run_cli(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, str(SCRIPTS / "sync_is_notes.py"), *map(str, args)],
                          capture_output=True, text=True)


def read_v3_toc() -> list:
    with open(NOTES_V3_EN / "toc_v3.tsv", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t", quoting=csv.QUOTE_NONE))


def chapter_order(js: str) -> list:
    block = re.search(r"var CHAPTER_ORDER = \[(.*?)\];", js, re.S).group(1)
    return re.findall(r"\{ file: '([^']+)', label: '((?:[^'\\]|\\.)*)' \}", block)


class TestGoldenV2Overlay(unittest.TestCase):
    """build --src notes_en reproduces the pre-release site at 6964e91."""

    @classmethod
    def setUpClass(cls):
        if not NOTES_EN.is_dir():
            raise unittest.SkipTest(f"frozen v2 source not on this machine: {NOTES_EN}")
        if not have_golden_commit():
            raise unittest.SkipTest(f"commit {GOLDEN_COMMIT} not in this repository")
        cls.tmp = tempfile.TemporaryDirectory()
        tmp = Path(cls.tmp.name)
        site = make_online_site(tmp / "site")
        empty = tmp / "no_redirects.tsv"
        empty.write_text("old_file\tkind\ttargets\tnote\n", encoding="utf-8")
        cls.out = tmp / "build"
        result = run_cli("build", "--src", NOTES_EN, "--site", site, "--redirects", empty,
                         "--version", OLD_VERSION, "--out", cls.out)
        if result.returncode != 0:
            cls.tmp.cleanup()
            raise AssertionError(f"build failed:\n{result.stdout}\n{result.stderr}")
        cls.tree = online_tree()

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_src_markdown_byte_identical(self):
        online_md = sorted(p for p in self.tree if p.startswith("src/"))
        self.assertEqual(online_md, sorted(f"src/{p.name}" for p in (self.out / "src").glob("*.md")))
        for path in online_md:
            with self.subTest(path=path):
                self.assertEqual(online(path, binary=True), (self.out / path).read_bytes())

    def test_notes_js_byte_identical(self):
        self.assertEqual(online("assets/notes.js", binary=True), (self.out / "assets/notes.js").read_bytes())

    def test_style_css_byte_identical(self):
        self.assertEqual(online("assets/style.css", binary=True), (self.out / "assets/style.css").read_bytes())

    def test_pages_equal_online_after_norm(self):
        pages = sorted(p for p in self.tree if "/" not in p and p.endswith(".html") and p != "index.html")
        self.assertEqual(len(pages), 28)
        for page in pages:
            with self.subTest(page=page):
                self.assertEqual(norm(online(page)), norm((self.out / page).read_text(encoding="utf-8")))

    def test_every_normalization_is_needed(self):
        pages = sorted(p for p in self.tree if "/" not in p and p.endswith(".html"))
        texts = [(online(page), (self.out / page).read_text(encoding="utf-8")) for page in pages]
        for rule, *_ in NORMALIZATIONS:
            with self.subTest(rule=rule):
                self.assertTrue(any(norm(a, skip=rule) != norm(b, skip=rule) for a, b in texts),
                                f"{rule} absorbs no difference; drop it from NORMALIZATIONS")

    def test_index_equal_online_after_norm(self):
        self.assertEqual(norm(online("index.html")),
                         norm((self.out / "index.html").read_text(encoding="utf-8")))

    def test_template_differs_only_by_theme_init_line(self):
        built = (self.out / "assets/template.html").read_text(encoding="utf-8").splitlines(keepends=True)
        self.assertEqual(sum(1 for line in built if line == sync.THEME_INIT + "\n"), 1)
        without = "".join(line for line in built if line != sync.THEME_INIT + "\n")
        self.assertEqual(online("assets/template.html"), without)

    def test_install_reproduces_online_file_set(self):
        with tempfile.TemporaryDirectory() as dest:
            result = run_cli("install", "--build", self.out, "--dest", dest)
            self.assertEqual(result.returncode, 0, result.stderr)
            published = files_under(Path(dest))
            # Decision 7: the v2 validator is no longer published.
            self.assertEqual(set(self.tree) - {"validate_book.py"}, published)
            differing = sorted(p for p in published
                               if not p.endswith(OVERLAY_SUFFIXES) and p not in OVERLAY_FILES
                               and blob_id(Path(dest) / p) != self.tree[p])
            self.assertEqual(sorted(KATEX_UPGRADED), differing)


class TestV3Build(unittest.TestCase):
    """build --src notes_v3_en with the real redirect table and a new version."""

    @classmethod
    def setUpClass(cls):
        if not (NOTES_V3_EN / "toc_v3.tsv").is_file():
            raise unittest.SkipTest(f"v3 English workspace not on this machine: {NOTES_V3_EN}")
        if not have_golden_commit():
            raise unittest.SkipTest(f"commit {GOLDEN_COMMIT} not in this repository")
        cls.tmp = tempfile.TemporaryDirectory()
        tmp = Path(cls.tmp.name)
        site = make_online_site(tmp / "site")
        cls.out = tmp / "build"
        result = run_cli("build", "--src", NOTES_V3_EN, "--site", site, "--redirects", REDIRECTS,
                         "--version", NEW_VERSION, "--out", cls.out)
        if result.returncode != 0:
            cls.tmp.cleanup()
            raise AssertionError(f"build failed:\n{result.stdout}\n{result.stderr}")
        cls.toc = read_v3_toc()
        cls.pages = [f"{row['number']}_{row['slug']}.html" for row in cls.toc]

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def page(self, name: str) -> str:
        return (self.out / name).read_text(encoding="utf-8")

    def test_notes_js_has_37_entries_in_toc_order(self):
        entries = chapter_order(self.page("assets/notes.js"))
        self.assertEqual(len(entries), 37)
        self.assertEqual(entries[0], ("index.html", "Cover"))
        self.assertEqual(entries[1], ("introduction.html", "Introduction · How to Use"))
        expected = [(f"{row['number']}_{row['slug']}.html",
                     f"{row['number']} {row['en_short']}".replace("'", "\\'")) for row in self.toc]
        self.assertEqual(expected, entries[2:])
        template = self.page("assets/template.html")
        self.assertEqual(re.findall(r"\?v=([\w-]+)", template), [NEW_VERSION, NEW_VERSION])
        self.assertIn(f'assets/style.css?v={NEW_VERSION}"', template)
        self.assertIn(f'assets/notes.js?v={NEW_VERSION}"', template)

    def test_every_page_carries_new_version(self):
        stubs = [row["old_file"] for row in sync.load_redirects(REDIRECTS)]
        for name in ["index.html", "introduction.html", *self.pages, *stubs]:
            with self.subTest(page=name):
                html = self.page(name)
                self.assertIn(f'href="assets/style.css?v={NEW_VERSION}"', html)
                self.assertIn(f'src="assets/notes.js?v={NEW_VERSION}"', html)
                self.assertEqual(set(re.findall(r"\?v=(notes-ui-\d+)", html)), {NEW_VERSION})
        stale = sorted(p.relative_to(self.out).as_posix() for p in self.out.rglob("*")
                       if p.is_file() and OLD_VERSION.encode() in p.read_bytes())
        self.assertEqual([], stale)

    def test_redirect_targets_are_v3_pages(self):
        for row in sync.load_redirects(REDIRECTS):
            with self.subTest(old=row["old_file"]):
                self.assertTrue(set(row["targets"]) <= set(self.pages), row["targets"])

    def test_split_stub_lists_all_targets(self):
        descs = {f"{row['number']}_{row['slug']}.html": row["en_desc"] for row in self.toc}
        for old, (kind, targets) in APPENDIX_C.items():
            if kind != "split":
                continue
            with self.subTest(old=old):
                html = self.page(old)
                hrefs = re.findall(r'<a class="chapter-card" href="([^"]+)"', html)
                self.assertEqual(targets, hrefs)
                for target in targets:
                    self.assertIn(sync.html_escape(descs[target]), html)
                self.assertIn('<meta name="robots" content="noindex" />', html)
                self.assertIn(f'<link rel="canonical" href="{targets[0]}" />', html)
                self.assertNotIn("http-equiv", html)
                self.assertNotIn("location.replace", html)

    def test_one_to_one_stub_keeps_hash(self):
        html = self.page("06_inference.html")
        self.assertIn('<script>location.replace("08_inference.html" + location.hash)</script>', html)
        self.assertIn('<meta http-equiv="refresh" content="0; url=08_inference.html" />', html)
        self.assertIn('<link rel="canonical" href="08_inference.html" />', html)
        self.assertIn('<meta name="robots" content="noindex" />', html)
        self.assertRegex(html, r'This chapter is now <a href="08_inference.html">Chapter 8 · Inference: '
                               r'[^<]+</a> in the 35-chapter edition\.')
        for old, (kind, targets) in APPENDIX_C.items():
            if kind == "one":
                with self.subTest(old=old):
                    self.assertIn(f'location.replace("{targets[0]}" + location.hash)', self.page(old))

    def test_file_nav_markers_kept(self):
        template = self.page("assets/template.html")
        block = re.search(r"<!-- TOC:FILE-NAV -->\n(.*?)<!-- /TOC:FILE-NAV -->", template, re.S)
        self.assertIsNotNone(block)
        self.assertEqual(re.findall(r'<a href="([^"]+)">', block.group(1)), self.pages)

    def test_rendered_pages_carry_site_layer(self):
        for name in ["introduction.html", *self.pages]:
            with self.subTest(page=name):
                html = self.page(name)
                self.assertIn(sync.THEME_INIT, html)
                self.assertIn('<div class="notes-topbar">', html)
                self.assertIn(f"<div class=\"sidebar-title\">{sync.SIDEBARTITLE}", norm(html))
        self.assertIn(f'<body class="{sync.INTRO_BODYCLASS}">', self.page("introduction.html"))

    def test_index_carries_site_layer(self):
        html = self.page("index.html")
        self.assertNotIn("<style>", html)
        self.assertNotIn("cover-sub-en", html)
        self.assertIn(sync.THEME_INIT, html)
        self.assertIn(sync.BACKLINK, html)
        self.assertIn(f'<footer class="cover-footer">{sync.SERIES}</footer>', html)
        self.assertEqual(re.findall(r'<a class="chapter-card" href="(\d\d_[^"]+)"', html), self.pages)

    def test_install_publishes_only_the_publish_set(self):
        stubs = [row["old_file"] for row in sync.load_redirects(REDIRECTS)]
        with tempfile.TemporaryDirectory() as dest:
            result = run_cli("install", "--build", self.out, "--dest", dest)
            self.assertEqual(result.returncode, 0, result.stderr)
            published = files_under(Path(dest))
        top_html = sorted(p for p in published if "/" not in p and p.endswith(".html"))
        self.assertEqual(sorted(["index.html", "introduction.html", *self.pages, *stubs]), top_html)
        self.assertEqual(len(top_html), 2 + 35 + 23)
        for path in published:
            with self.subTest(path=path):
                self.assertFalse(path.startswith(("tools/", "demo")), path)
                self.assertNotIn(path, {"toc_v3.tsv", "validate_book.py", "RENUMBERED_V3"})
                self.assertTrue(path in top_html or path == "inline_figures.py"
                                or path.startswith(("src/", "code/", "assets/fig/", "assets/katex/"))
                                or path in {"assets/style.css", "assets/notes.js", "assets/template.html",
                                            "assets/wrap-tables.lua"}, path)


class TestRedirectTable(unittest.TestCase):

    def test_redirects_cover_all_v2_names(self):
        rows = sync.load_redirects(REDIRECTS)
        self.assertEqual(len(rows), 23)
        self.assertEqual({row["old_file"]: (row["kind"], row["targets"]) for row in rows}, APPENDIX_C)
        if have_golden_commit():
            names = git("ls-tree", "--name-only", GOLDEN_COMMIT, f"{COLLECTION}/").split()
            v2_pages = {Path(n).name for n in names if re.fullmatch(r"\d\d_.+\.html", Path(n).name)
                        and int(Path(n).name[:2]) >= 5}
            self.assertEqual(v2_pages, {row["old_file"] for row in rows})

    def test_notes_name_the_moved_material(self):
        notes = {row["old_file"]: row["note"] for row in sync.load_redirects(REDIRECTS)}
        self.assertIn("Chapter 20", notes["13_did_event_studies.html"])
        self.assertIn("DID_M", notes["13_did_event_studies.html"])
        self.assertIn("Chapter 12", notes["05_estimation_mle_gmm.html"])
        self.assertIn("Helix", notes["05_estimation_mle_gmm.html"])

    def test_load_rejects_bad_rows(self):
        bad = {
            "unknown kind": "06_inference.html\tmove\t08_inference.html\t\n",
            "one with two targets": "06_inference.html\tone\t08_inference.html,09_x.html\t\n",
            "split with one target": "07_panel.html\tsplit\t10_linear_panel_fe_re.html\t\n",
            "duplicate old file": "06_inference.html\tone\t08_inference.html\t\n" * 2,
        }
        for name, body in bad.items():
            with self.subTest(case=name), tempfile.TemporaryDirectory() as tmp:
                path = Path(tmp) / "r.tsv"
                path.write_text("old_file\tkind\ttargets\tnote\n" + body, encoding="utf-8")
                with self.assertRaises(ValueError):
                    sync.load_redirects(path)


CHAPTER_YAML = """---
title: "Inference: Delta Method, Bootstrap, Clustering & Testing"
subtitle: "Honest Standard Errors"
seriesline: "Foundations of Information Systems Economics · Chapter 8"
series: "Foundations of Information Systems Economics"
sidebartitle: "Foundations of IS Economics"
sidebarsub: "Graduate Empirical Methods"
chapterlabel: "Chapter 8 · Inference: Delta Method, Bootstrap, Clustering & Testing"
---

## Introduction

series: "body text is not front matter"
"""


class TestSiteYaml(unittest.TestCase):

    def test_site_yaml_fields(self):
        out = sync.site_yaml(CHAPTER_YAML)
        self.assertEqual(out, CHAPTER_YAML
                         .replace('series: "Foundations of Information Systems Economics"\n',
                                  f'series: "{sync.SERIES}"\n', 1)
                         .replace('sidebartitle: "Foundations of IS Economics"',
                                  f'sidebartitle: "{sync.SIDEBARTITLE}"')
                         .replace('sidebarsub: "Graduate Empirical Methods"',
                                  f'sidebarsub: "{sync.SIDEBARSUB}"'))
        self.assertNotIn("bodyclass", out)
        self.assertIn('series: "body text is not front matter"', out)

    def test_site_yaml_intro_gets_bodyclass(self):
        out = sync.site_yaml(CHAPTER_YAML, intro=True)
        self.assertIn(f'chapterlabel: "Chapter 8 · Inference: Delta Method, Bootstrap, Clustering & Testing"\n'
                      f'bodyclass: "{sync.INTRO_BODYCLASS}"\n---\n', out)
        self.assertEqual(out, sync.site_yaml(out, intro=True))

    def test_site_yaml_requires_each_field(self):
        for key in ("series", "sidebartitle", "sidebarsub"):
            with self.subTest(key=key):
                text = re.sub(rf"^{key}: .*\n", "", CHAPTER_YAML, count=1, flags=re.M)
                with self.assertRaises(ValueError):
                    sync.site_yaml(text)


MINI_TEMPLATE = """<html>
<head>
  <title>$title$</title>
  <link rel="stylesheet" href="assets/style.css" />
</head>
<body>
<div class="file-nav">
      <!-- TOC:FILE-NAV -->
      <a href="01_a.html">01 A</a>
      <!-- /TOC:FILE-NAV -->
</div>
$body$
</body>
</html>
"""


class TestSiteTemplate(unittest.TestCase):

    def test_site_template_adds_the_five_pieces(self):
        out = sync.site_template(MINI_TEMPLATE, NEW_VERSION)
        self.assertEqual(out, MINI_TEMPLATE
                         .replace('  <link rel="stylesheet" href="assets/style.css" />\n',
                                  f'{sync.THEME_INIT}\n'
                                  f'  <link rel="stylesheet" href="assets/style.css?v={NEW_VERSION}" />\n')
                         .replace("<body>\n", '<body$if(bodyclass)$ class="$bodyclass$"$endif$>\n'
                                  + sync.TOPBAR + "\n")
                         .replace("</body>\n", f'<script defer src="assets/notes.js?v={NEW_VERSION}"></script>\n'
                                  "</body>\n"))

    def test_site_template_rejects_a_missing_anchor(self):
        for anchor in ('  <link rel="stylesheet" href="assets/style.css" />\n', "<body>\n", "</body>\n"):
            with self.subTest(anchor=anchor):
                with self.assertRaises(ValueError):
                    sync.site_template(MINI_TEMPLATE.replace(anchor, ""), NEW_VERSION)
        with self.assertRaises(ValueError):
            sync.site_template(sync.site_template(MINI_TEMPLATE, NEW_VERSION), NEW_VERSION)


OLD_JS = """(function () {
  var CHAPTER_ORDER = [
    { file: 'index.html', label: 'Cover' },
    { file: '01_old.html', label: '01 Old' }
  ];
  var other = 'assets/x.css?v=notes-ui-20260822';
}());
"""


class TestNotesJs(unittest.TestCase):

    def test_notes_js_replaces_only_chapter_order(self):
        rows = [{"number": "01", "slug": "a", "en_short": "Consumer, Producer & Welfare"},
                {"number": "02", "slug": "b", "en_short": "Sant'Anna"}]
        out = sync.notes_js(OLD_JS, rows, NEW_VERSION)
        self.assertEqual(out, """(function () {
  var CHAPTER_ORDER = [
    { file: 'index.html', label: 'Cover' },
    { file: 'introduction.html', label: 'Introduction · How to Use' },
    { file: '01_a.html', label: '01 Consumer, Producer & Welfare' },
    { file: '02_b.html', label: '02 Sant\\'Anna' }
  ];
  var other = 'assets/x.css?v=notes-ui-20261021';
}());
""")

    def test_notes_js_needs_the_array(self):
        with self.assertRaises(ValueError):
            sync.notes_js("(function () {}());\n", [], NEW_VERSION)


TITLES = {"08_inference.html": "Inference: Delta Method, Bootstrap, Clustering & Testing",
          "21_regression_discontinuity.html": "Regression Discontinuity Designs",
          "22_synthetic_control_factor_models.html": "Synthetic Control and Panel Factor Models"}
DESCS = {"21_regression_discontinuity.html": "Continuity-based identification & more",
         "22_synthetic_control_factor_models.html": "Synthetic control <inference>"}


class TestStubHtml(unittest.TestCase):

    def test_one_to_one_stub(self):
        row = {"old_file": "06_inference.html", "kind": "one", "targets": ["08_inference.html"],
               "note": "Tests & notes."}
        html = sync.stub_html(row, TITLES, NEW_VERSION)
        self.assertIn('<meta name="robots" content="noindex" />', html)
        self.assertIn('<link rel="canonical" href="08_inference.html" />', html)
        self.assertIn('<meta http-equiv="refresh" content="0; url=08_inference.html" />', html)
        self.assertIn('<script>location.replace("08_inference.html" + location.hash)</script>', html)
        self.assertIn('This chapter is now <a href="08_inference.html">Chapter 8 · Inference: Delta Method, '
                      'Bootstrap, Clustering &amp; Testing</a> in the 35-chapter edition.', html)
        self.assertIn("Tests &amp; notes.", html)
        self.assertIn(f'href="assets/style.css?v={NEW_VERSION}"', html)
        self.assertIn(sync.THEME_INIT, html)
        self.assertLess(html.index("location.replace"), html.index("http-equiv"))

    def test_split_landing_page(self):
        row = {"old_file": "14_rdd_synthetic_control.html", "kind": "split",
               "targets": ["21_regression_discontinuity.html", "22_synthetic_control_factor_models.html"],
               "note": ""}
        html = sync.stub_html(row, TITLES, NEW_VERSION, descs=DESCS)
        self.assertEqual(re.findall(r'<a class="chapter-card" href="([^"]+)"', html), row["targets"])
        self.assertIn('<link rel="canonical" href="21_regression_discontinuity.html" />', html)
        self.assertIn('<meta name="robots" content="noindex" />', html)
        self.assertIn("Continuity-based identification &amp; more", html)
        self.assertIn("Synthetic control &lt;inference&gt;", html)
        self.assertNotIn("http-equiv", html)
        self.assertNotIn("location", html)
        self.assertIn(f'src="assets/notes.js?v={NEW_VERSION}"', html)

    def test_split_landing_page_needs_descriptions(self):
        row = {"old_file": "14_rdd_synthetic_control.html", "kind": "split",
               "targets": ["21_regression_discontinuity.html", "22_synthetic_control_factor_models.html"],
               "note": ""}
        with self.assertRaises(ValueError):
            sync.stub_html(row, TITLES, NEW_VERSION)


def make_build(root: Path) -> Path:
    """A minimal build tree, with tooling the publish set must leave out."""
    files = {
        "index.html": "<html></html>", "introduction.html": "i", "01_a.html": "a", "06_inference.html": "stub",
        "assets/style.css": "css", "assets/notes.js": "js", "assets/template.html": "t",
        "assets/wrap-tables.lua": "lua", "assets/fig/f.svg": "<svg/>", "assets/katex/katex.min.js": "k",
        "src/01_a.md": "md", "src/introduction.md": "md", "code/01_a/run.R": "R", "inline_figures.py": "py",
        "tools/toc.py": "py", "toc_v3.tsv": "tsv", "validate_book.py": "py", "assets/fig/.DS_Store": "",
        "code/01_a/__pycache__/x.pyc": "", "demo.html": "demo", "src/demo.md": "demo", "RENUMBERED_V3": "",
    }
    for name, body in files.items():
        (root / name).parent.mkdir(parents=True, exist_ok=True)
        (root / name).write_text(body, encoding="utf-8")
    (root / sync.BUILD_MARKER).write_text(sync.BUILD_COMPLETE + "\n", encoding="utf-8")
    return root


class TestInstall(unittest.TestCase):

    def test_publish_set_excludes_tooling(self):
        with tempfile.TemporaryDirectory() as tmp:
            build = make_build(Path(tmp) / "build")
            dest = Path(tmp) / "dest"
            for name, body in {"index.html": f"<h1>{sync.SIDEBARTITLE}</h1>", "validate_book.py": "v2", "27_quick_reference.html": "v2",
                               "src/27_quick_reference.md": "v2", "assets/fig/old.svg": "v2"}.items():
                (dest / name).parent.mkdir(parents=True, exist_ok=True)
                (dest / name).write_text(body, encoding="utf-8")
            result = run_cli("install", "--build", build, "--dest", dest)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(files_under(dest), {
                "index.html", "introduction.html", "01_a.html", "06_inference.html", "assets/style.css",
                "assets/notes.js", "assets/template.html", "assets/wrap-tables.lua", "assets/fig/f.svg",
                "assets/katex/katex.min.js", "src/01_a.md", "src/introduction.md", "code/01_a/run.R",
                "inline_figures.py"})
            self.assertEqual((dest / "index.html").read_text(encoding="utf-8"), "<html></html>")

    def test_install_refuses_a_foreign_dest(self):
        with tempfile.TemporaryDirectory() as tmp:
            build = make_build(Path(tmp) / "build")
            dest = Path(tmp) / "elsewhere"
            dest.mkdir()
            (dest / "precious.txt").write_text("keep", encoding="utf-8")
            result = run_cli("install", "--build", build, "--dest", dest)
            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(files_under(dest), {"precious.txt"})
            # A sibling collection has an index.html too, but not this book's cover.
            (dest / "index.html").write_text("<h1>From Choices to Counterfactuals</h1>", encoding="utf-8")
            result = run_cli("install", "--build", build, "--dest", dest)
            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(files_under(dest), {"precious.txt", "index.html"})

    def test_install_refuses_a_tree_that_is_not_a_build(self):
        with tempfile.TemporaryDirectory() as tmp:
            build = make_build(Path(tmp) / "build")
            (build / sync.BUILD_MARKER).unlink()
            result = run_cli("install", "--build", build, "--dest", Path(tmp) / "dest")
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse((Path(tmp) / "dest").exists())

    def test_install_refuses_an_unfinished_build(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            for name in ("ws/index.html", "ws/assets/template.html", "ws/src/introduction.md", "ws/inline_figures.py",
                         "site/assets/style.css", "site/assets/notes.js"):
                (tmp / name).parent.mkdir(parents=True, exist_ok=True)
                (tmp / name).write_text("no front matter, no file-nav", encoding="utf-8")
            redirects = write_redirects(tmp / "r.tsv")
            result = run_cli("build", "--src", tmp / "ws", "--site", tmp / "site", "--redirects", redirects,
                             "--version", NEW_VERSION, "--out", tmp / "out")
            self.assertEqual(result.returncode, 2, result.stderr)
            self.assertTrue((tmp / "out" / sync.BUILD_MARKER).is_file())
            result = run_cli("install", "--build", tmp / "out", "--dest", tmp / "dest")
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("did not finish", result.stderr)
            self.assertFalse((tmp / "dest").exists())


SITE_URL = "https://bobbyzhong10.github.io"


def make_site(root: Path, stub_target: str = "08_c.html") -> Path:
    """A minimal _site that passes check: 37 pages, one one-to-one stub, a sitemap without it."""
    coll = root / COLLECTION
    (coll / "assets").mkdir(parents=True)
    (coll / "assets/style.css").write_text("css", encoding="utf-8")
    (root / "notes").mkdir(exist_ok=True)
    (root / "notes/index.html").write_text("<html></html>", encoding="utf-8")
    chapters = [f"{n:02d}_c.html" for n in range(1, 36)]
    pages = ["index.html", "introduction.html", *chapters]
    for i, page in enumerate(pages):
        nxt = pages[(i + 1) % len(pages)]
        (coll / page).write_text(
            f'<html><head><link rel="stylesheet" href="assets/style.css?v={NEW_VERSION}" /></head><body>'
            f'<a href="/notes/">back</a><a href="https://example.org/x">ext</a><a href="#top">top</a>'
            f'<h2 id="top">T</h2><a href="{nxt}#top">next</a><img src="data:image/png;base64,AAAA" />'
            f'</body></html>', encoding="utf-8")
    (coll / "06_old.html").write_text(
        f'<html><head><meta name="robots" content="noindex" /><link rel="canonical" href="{stub_target}" />'
        f'</head><body><a href="{stub_target}">moved</a></body></html>', encoding="utf-8")
    urls = [f"{SITE_URL}/{COLLECTION}/"] + [f"{SITE_URL}/{COLLECTION}/{p}" for p in pages[1:]]
    urls.append(f"{SITE_URL}/about/")
    (root / "sitemap.xml").write_text(
        '<?xml version="1.0" encoding="UTF-8"?><urlset>'
        + "".join(f"<url><loc>{u}</loc></url>" for u in urls) + "</urlset>", encoding="utf-8")
    return root


def write_redirects(path: Path, target: str = "08_c.html") -> Path:
    path.write_text(f"old_file\tkind\ttargets\tnote\n06_old.html\tone\t{target}\t\n", encoding="utf-8")
    return path


class TestCheck(unittest.TestCase):

    def check(self, site: Path, redirects: Path) -> subprocess.CompletedProcess:
        return run_cli("check", "--site-root", site, "--collection", COLLECTION, "--redirects", redirects)

    def test_check_passes_a_clean_site(self):
        with tempfile.TemporaryDirectory() as tmp:
            site = make_site(Path(tmp) / "_site")
            result = self.check(site, write_redirects(Path(tmp) / "r.tsv"))
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertEqual(sync.check_site(site, COLLECTION, Path(tmp) / "r.tsv"), [])

    def test_check_flags_broken_link(self):
        with tempfile.TemporaryDirectory() as tmp:
            site = make_site(Path(tmp) / "_site")
            page = site / COLLECTION / "03_c.html"
            page.write_text(page.read_text(encoding="utf-8").replace("</body>", '<a href="missing.html">x</a></body>'),
                            encoding="utf-8")
            result = self.check(site, write_redirects(Path(tmp) / "r.tsv"))
            self.assertEqual(result.returncode, 1)
            self.assertIn("missing.html", result.stdout)
            self.assertIn("03_c.html", result.stdout)

    def test_check_flags_missing_anchor(self):
        with tempfile.TemporaryDirectory() as tmp:
            site = make_site(Path(tmp) / "_site")
            page = site / COLLECTION / "04_c.html"
            page.write_text(page.read_text(encoding="utf-8").replace("</body>", '<a href="#nowhere">x</a>'
                                                                     '<a href="05_c.html#gone">y</a></body>'),
                            encoding="utf-8")
            errors = sync.check_site(site, COLLECTION, write_redirects(Path(tmp) / "r.tsv"))
            self.assertEqual(len(errors), 2, errors)
            self.assertTrue(any("#nowhere" in e for e in errors), errors)
            self.assertTrue(any("05_c.html#gone" in e for e in errors), errors)

    def test_check_flags_published_tooling(self):
        with tempfile.TemporaryDirectory() as tmp:
            site = make_site(Path(tmp) / "_site")
            coll = site / COLLECTION
            for name in ("src/01_c.md", "code/run.R", "tools/toc.py", "toc_v3.tsv", "inline_figures.py",
                         "assets/template.html", "assets/wrap-tables.lua"):
                (coll / name).parent.mkdir(parents=True, exist_ok=True)
                (coll / name).write_text('<link rel="stylesheet" href="$css$" />', encoding="utf-8")
            (coll / "tools/report.html").write_text('<a href="missing.html">x</a>', encoding="utf-8")
            errors = sync.check_site(site, COLLECTION, write_redirects(Path(tmp) / "r.tsv"))
            for name in ("src", "code", "tools", "toc_v3.tsv", "inline_figures.py", "assets/template.html",
                         "assets/wrap-tables.lua"):
                self.assertTrue(any(name in e for e in errors), (name, errors))
            self.assertTrue(all(e.endswith("must not be published") for e in errors), errors)

    def test_check_flags_redirect_without_target(self):
        with tempfile.TemporaryDirectory() as tmp:
            site = make_site(Path(tmp) / "_site", stub_target="99_gone.html")
            errors = sync.check_site(site, COLLECTION, write_redirects(Path(tmp) / "r.tsv", "99_gone.html"))
            self.assertTrue(any("99_gone.html" in e and "06_old.html" in e for e in errors), errors)

    def test_check_flags_redirect_that_is_not_a_stub(self):
        with tempfile.TemporaryDirectory() as tmp:
            site = make_site(Path(tmp) / "_site")
            (site / COLLECTION / "06_old.html").write_text('<html><a href="08_c.html">old page</a></html>',
                                                         encoding="utf-8")
            errors = sync.check_site(site, COLLECTION, write_redirects(Path(tmp) / "r.tsv"))
            self.assertTrue(any("06_old.html" in e for e in errors), errors)

    def test_check_flags_sitemap(self):
        with tempfile.TemporaryDirectory() as tmp:
            site = make_site(Path(tmp) / "_site")
            sitemap = site / "sitemap.xml"
            sitemap.write_text(sitemap.read_text(encoding="utf-8").replace(
                "</urlset>", f"<url><loc>{SITE_URL}/{COLLECTION}/06_old.html</loc></url></urlset>"),
                encoding="utf-8")
            errors = sync.check_site(site, COLLECTION, write_redirects(Path(tmp) / "r.tsv"))
            self.assertTrue(any("sitemap" in e and "06_old.html" in e for e in errors), errors)
            self.assertTrue(any("sitemap" in e and "38" in e for e in errors), errors)


if __name__ == "__main__":
    unittest.main()
