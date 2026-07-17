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
| `assets/css/main.scss` | Main stylesheet entry point and site-specific visual overrides |
| `_sass/` | Base theme modules imported by `main.scss` |
| `assets/js/` | Main-site JavaScript and bundled plugins |
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
- Add, remove, or reorder top navigation items in `_data/navigation.yml`.
- Update name, institution, avatar, profile links, or SEO metadata in `_config.yml`.
- Edit the homepage identity block in `_includes/home-hero.html`.
- Edit the shared navigation shell in `_includes/masthead.html`.
- Edit theme switching and heading-anchor behavior in `_includes/scripts.html`.
- Add site-wide or page-specific styles near the end of `assets/css/main.scss`. Reuse the existing CSS variables instead of hard-coding new colors.
- Change low-level theme behavior in `_sass/` only when the change should affect the whole site.

## Lecture Notes

### Foundations of Information Systems Economics

Location: `notes/is-economics-methods/`

- `index.html` is the collection cover and chapter index.
- `src/*.md` contains the editable chapter sources and metadata.
- Root-level chapter `.html` files are the published Pandoc outputs.
- `assets/template.html` defines the shared chapter page structure.
- `assets/style.css` controls the complete lecture-note theme.
- `assets/katex/` contains local math-rendering assets.
- `validate_book.py` checks the collection structure.
- `inline_figures.py` is a post-Pandoc figure-processing utility.

For chapter content changes, edit the matching file in `src/` and regenerate its HTML. If HTML must be edited directly, keep the corresponding Markdown metadata and content synchronized so a future regeneration does not undo the change.

### From Choices to Counterfactuals

Location: `notes/structural-econometrics/`

- `index.html` is the collection cover and chapter index.
- Numbered root-level `.html` files are the published chapters.
- `assets/style.css` is the shared lecture-note stylesheet.
- `assets/notes.js` provides lecture-note interactions.
- `assets/katex/` contains local math-rendering assets.

Keep relative asset and chapter links intact when moving or renaming files in either collection.

## Local Build and Preview

The dependency lock works with Ruby 3.1. On macOS with Homebrew:

```bash
export PATH="/opt/homebrew/opt/ruby@3.1/bin:$PATH"
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
