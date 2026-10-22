# Site Architecture

This repository contains a Jekyll-based academic website plus two self-contained lecture-note sites. GitHub Pages builds the repository root and publishes the generated files.

## Build Flow

```text
_pages/*.md
    + _layouts/default.html
    + _includes/*
    + _data/navigation.yml
    + assets/css/main.scss + _sass/*
                    |
                    v
              Jekyll build
                    |
                    v
                  _site/
```

`_site/` is generated output and is ignored by Git. Make permanent changes in the source files, not in `_site/`.

The `notes/` directory is different: its HTML, CSS, JavaScript, and local assets are checked-in source files. Jekyll copies them into `_site/notes/` without applying the main site layout.

## Main Directories

| Path | Purpose |
| --- | --- |
| `_pages/` | Main page content and route front matter |
| `_layouts/default.html` | Shared page shell, body classes, main content, and footer |
| `_includes/` | Header, navigation, home hero, SEO, scripts, and reusable fragments |
| `_data/navigation.yml` | Top navigation labels and URLs |
| `_config.yml` | Site identity, author metadata, plugins, build rules, and exclusions |
| `assets/css/main.scss` | Stylesheet entry point: imports the theme modules, then `_sass/_site.scss` |
| `_sass/_site.scss` | All site-specific styles and the color tokens (light and dark) |
| `_sass/` (other files) | Base theme modules imported by `main.scss` |
| `images/` | Avatar, favicons, manifest, and other site images |
| `notes/` | Standalone lecture-note collections |
| `docs/` | Internal documentation; excluded from the published site |

## Page and Route Map

| URL | Source |
| --- | --- |
| `/` | `_pages/about.md` |
| `/research-projects/` | `_pages/research-projects.md` |
| `/education/` | `_pages/education.md` |
| `/awards-services/` | `_pages/awards-services.md` |
| `/teaching/` | `_pages/teaching.md` |
| `/cv/` | `_pages/cv.md` |

Each page uses YAML front matter for its permalink, title, and layout settings. The default layout adds a route-derived body class such as `page-education`, which can be used for page-specific styles.

## Common Modification Points

- Edit page text and cards in the corresponding file under `_pages/`.
- For a paper whose authors are listed alphabetically, add `{% include alphabetical-order-marker.html %}` immediately after its author list in `_pages/research-projects.md`. Keep the `† Authors listed in alphabetical order.` legend in the same publication section.
- Add, remove, or reorder top navigation items in `_data/navigation.yml`.
- Update name, institution, avatar, profile links, or SEO metadata in `_config.yml`.
- Edit the homepage identity block in `_includes/home-hero.html`.
- Edit the shared navigation shell in `_includes/masthead.html`.
- Edit theme switching and the abstract toggles in `_includes/scripts.html`.
- Add site-wide or page-specific styles to the matching section of `_sass/_site.scss`. Reuse the color tokens at the top of that file instead of hard-coding new colors.
- Change low-level theme behavior in `_sass/` only when the change should affect the whole site.

## Brand Assets

- `images/brand-mark.svg` is the generated vector mark used as the preferred browser favicon; it is intentionally not displayed inside the page layout.
- `images/profile-yunzhou-zhong.png` is the homepage portrait and social-sharing image.
- `scripts/generate_brand_assets.py` regenerates the PNG favicons, Apple/Android icons, Windows tile, and `.ico` file from the same navy, crimson, and white BZ monogram design.
- To change the monogram geometry, typography, or colors, edit `scripts/generate_brand_assets.py`, run it, and rebuild the site.

## Lecture Notes

### Foundations of Information Systems Economics

Location: `notes/is-economics-methods/`

The collection is the introduction plus 35 numbered chapters. Its master copy is the English workspace `notes_v3_en` in iCloud. `scripts/sync_is_notes.py` generates the site copy from it and replays the site's presentation layer on top of the book: the site metadata, the topbar and theme switch, the versioned stylesheet and `notes.js`, and the restyled cover.

- `index.html` is the collection cover and chapter index.
- Root-level chapter `.html` files are the published Pandoc outputs.
- `src/*.md` contains the chapter sources with the site metadata. They are not published.
- `code/` holds each chapter's R code. It is mirrored into the repository for reference but, like `src/`, excluded from the served site.
- `assets/template.html` defines the shared chapter page structure.
- `assets/style.css` controls the complete lecture-note theme. It is the site's own file: the sync script reads it from this folder, so edit it here.
- `assets/notes.js` provides the topbar, theme switch, mobile sidebar, and previous/next chapter links. Its chapter order is generated from the book's table of contents.
- `assets/katex/` contains local math-rendering assets.
- `inline_figures.py` is a post-Pandoc figure-processing utility.
- The 23 numbered file names of the earlier 27-chapter edition, listed in `scripts/is_notes_redirects.tsv`, are redirect pages. Each names its new page as an absolute canonical URL built from `url` and `baseurl` in `_config.yml`. A moved chapter redirects at once and keeps the `#anchor`. A split chapter (old 05, 07, and 14) is a landing page that links its new chapters, does not redirect, carries `noindex`, and names the first new chapter as canonical. Defaults in `_config.yml` keep all 23 pages out of `sitemap.xml`.

The book's validator `validate_book.py`, its `tools/`, and `toc_v3.tsv` stay in the workspace and are not published.

For chapter content changes, edit the book in `notes_v3_en`, not the HTML here, then sync from the repository root:

```bash
python3 scripts/sync_is_notes.py build --src "<path to notes_v3_en>" --site notes/is-economics-methods --redirects scripts/is_notes_redirects.tsv --version notes-ui-YYYYMMDD --out <build folder>
python3 scripts/sync_is_notes.py install --build <build folder> --dest notes/is-economics-methods
python3 scripts/sync_is_notes.py check --site-root _site --collection notes/is-economics-methods --redirects scripts/is_notes_redirects.tsv
```

- `build` renders every page with the book's Pandoc command. A new `--version` date makes returning readers' browsers fetch the new stylesheet and `notes.js`; add a letter (`notes-ui-YYYYMMDDb`) for a second release on the same day.
- `build` expects Pandoc 3.12.1. The script's golden test registers the output differences of that version, and another version prints a warning.
- `install` mirrors the build into this folder and deletes everything else in it.
- Run `check` after a complete Jekyll build. It must report 0 errors: it checks links and anchors, the redirect pages, that no sources or tooling are published, and that `sitemap.xml` lists exactly the 37 pages.
- The script's tests run with `python3 -m unittest discover -s scripts/tests`.

### From Choices to Counterfactuals

Location: `notes/structural-econometrics/`

- `index.html` is the collection cover and chapter index.
- Numbered root-level `.html` files are the published chapters.
- `assets/style.css` is the shared lecture-note stylesheet.
- `assets/notes.js` provides lecture-note interactions.
- `assets/katex/` contains local math-rendering assets.

Keep relative asset and chapter links intact when moving or renaming files in either collection.

## Local Build and Preview

The local toolchain mirrors GitHub Pages: `mise.toml` pins Ruby 3.3.4 and the `Gemfile` pins `github-pages` 232. When GitHub Pages moves to a new release (listed at <https://pages.github.com/versions/>), update both together. With mise activated in the shell:

```bash
mise install
BUNDLE_PATH=vendor/bundle bundle install
BUNDLE_PATH=vendor/bundle bundle exec jekyll build
python3 -m http.server 4000 --bind 127.0.0.1 --directory _site
```

Preview the site at `http://localhost:4000/`. Rebuild after changing Jekyll pages, includes, configuration, or SCSS. Static note files are also recopied during the build.

## Deployment and Safety Rules

GitHub Pages publishes after changes are pushed to the configured GitHub repository. Before pushing:

1. Run a complete Jekyll build.
2. Check the affected routes and note links locally.
3. Review `git status` and `git diff --check`.
4. Do not commit `_site/`, `vendor/`, `.bundle/`, `.DS_Store`, or other generated caches.

When `_config.yml` changes, restart or rebuild the site because configuration is not reloaded incrementally.
