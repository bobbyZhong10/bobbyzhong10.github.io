# Yunzhou (Bobby) Zhong — Personal Academic Website

Source code for [bobbyzhong10.github.io](https://bobbyzhong10.github.io), my personal academic website. The site is built with [Jekyll](https://jekyllrb.com/) and published through GitHub Pages.

## Table of Contents

- [Site Map](#site-map)
- [Repository Layout](#repository-layout)
- [Local Development](#local-development)
- [Editing Content](#editing-content)
- [Brand Assets](#brand-assets)
- [Citation Statistics (Optional)](#citation-statistics-optional)
- [Deployment](#deployment)
- [License](#license)
- [Contact](#contact)

## Site Map

| Route | Source file | What it contains |
| --- | --- | --- |
| `/` | `_pages/about.md` | Introduction, PhD application note, and news |
| `/research-projects/` | `_pages/research-projects.md` | Research interests, working papers, work in progress, conference proceedings, presentations, and thesis |
| `/education/` | `_pages/education.md` | Academic background and PhD-level coursework |
| `/awards-services/` | `_pages/awards-services.md` | Honors and awards, academic service, and professional memberships |
| `/teaching/` | `_pages/teaching.md` | Teaching experience |
| `/cv/` | `_pages/cv.md` | Embedded CV with a link to the PDF |
| `/notes/` | `_pages/notes.md` | Notes index |

The top navigation is defined in `_data/navigation.yml`. Every page uses the single `_layouts/default.html` layout, which adds a route-derived body class (for example `page-education`) that can be targeted for page-specific styles.

## Repository Layout

```text
.
├── _config.yml              # Site identity, author links, plugins, build exclusions
├── _data/navigation.yml     # Top navigation labels and URLs
├── _pages/                  # One Markdown file per page (front matter sets the permalink)
├── _layouts/default.html    # Shared page shell: head, masthead, hero, content, footer
├── _includes/               # Masthead, home hero, SEO tags, scripts, reusable fragments
├── _sass/                   # Base theme modules imported by assets/css/main.scss
├── assets/
│   ├── css/main.scss        # Stylesheet entry point plus all site-specific overrides
│   ├── js/                  # Main-site JavaScript and bundled plugins
│   └── fonts/               # Font Awesome and Academicons webfonts
├── images/                  # Portrait, favicons, app icons, web manifest, brand mark
├── notes/                   # Standalone note collections, copied to the site as-is
├── scripts/generate_brand_assets.py   # Regenerates the favicon and icon set
├── google_scholar_crawler/  # Optional citation-count crawler used by the workflow
├── .github/workflows/       # GitHub Actions (citation crawler, manual trigger only)
├── docs/                    # Internal documentation and design plans (not published)
├── Gemfile / Gemfile.lock   # Ruby dependencies pinned to the GitHub Pages gem
└── run_server.sh            # Shortcut for a live-reloading local server
```

`_site/` is generated output and is ignored by Git. Make permanent changes in the source files, never in `_site/`.

## Local Development

### Prerequisites

- Ruby 3.1 (the lockfile was resolved against this version)
- Bundler

On macOS with Homebrew:

```bash
brew install ruby@3.1
export PATH="/opt/homebrew/opt/ruby@3.1/bin:$PATH"
```

### Install dependencies

```bash
BUNDLE_PATH=vendor/bundle bundle install
```

Installing into `vendor/bundle` keeps gems local to the project. Do not commit `vendor/` or `.bundle/`.

### Run a local server

Live-reloading server (this is what `run_server.sh` does):

```bash
bundle exec jekyll liveserve
```

Or a plain Jekyll server:

```bash
bundle exec jekyll serve --host 127.0.0.1 --port 4000
```

Then open `http://127.0.0.1:4000/`. Jekyll rebuilds when pages, includes, or SCSS change. `_config.yml` is not reloaded automatically, so restart the server after editing it.

### Build only

```bash
BUNDLE_PATH=vendor/bundle bundle exec jekyll build
python3 -m http.server 4000 --bind 127.0.0.1 --directory _site
```

## Editing Content

Common changes and where to make them:

| Task | Where |
| --- | --- |
| Edit page text, cards, news items, or paper lists | The matching file under `_pages/` |
| Add, remove, or reorder navigation items | `_data/navigation.yml` |
| Update name, affiliation, email, Scholar, ORCID, or LinkedIn links | `author:` block in `_config.yml` |
| Change the homepage identity block (name, links, portrait) | `_includes/home-hero.html` |
| Change the navigation shell or theme toggle button | `_includes/masthead.html` |
| Adjust theme switching or heading-anchor behavior | `_includes/scripts.html` |
| Add site-wide or page-specific styles | End of `assets/css/main.scss`, reusing the existing CSS variables |
| Change low-level theme behavior for the whole site | `_sass/` |
| Replace the portrait or social-sharing image | `images/profile-yunzhou-zhong.png` |
| Update the CV | Replace the Google Drive link and "Last updated" date in `_pages/cv.md` |

Conventions worth keeping:

- Pages are written in Markdown with inline HTML blocks. Wrap Markdown inside HTML containers with `markdown="1"` so Jekyll still renders it.
- For a paper whose authors are listed alphabetically, add `{% include alphabetical-order-marker.html %}` right after the author list in `_pages/research-projects.md`, and keep the "Authors listed in alphabetical order" legend in the same section.
- The site supports light and dark themes. Reuse the existing CSS variables rather than hard-coding colors so both themes stay consistent.

A fuller description of the build flow and modification points lives in `docs/SITE_ARCHITECTURE.md`.

## Brand Assets

The favicons, Apple and Android icons, Windows tile, and `.ico` file all derive from one BZ monogram design (white initials on a navy rounded square with a crimson rule). `images/brand-mark.svg` is the vector version used as the preferred browser favicon.

To change the monogram's geometry, typography, or colors, edit and rerun the generator, which requires Python 3 with Pillow and the Georgia Bold font available on macOS:

```bash
python3 scripts/generate_brand_assets.py
```

## Citation Statistics (Optional)

`google_scholar_crawler/main.py` fetches citation counts from Google Scholar with the `scholarly` package and writes them as JSON. The GitHub Actions workflow in `.github/workflows/google_scholar_crawler.yaml` runs the crawler and force-pushes the results to a `google-scholar-stats` branch.

This workflow is currently paused. It only runs on manual trigger (`workflow_dispatch`), and the include that reads the statistics is commented out in `_includes/scripts.html`. To re-enable it:

1. Add a `GOOGLE_SCHOLAR_ID` repository secret with the Scholar profile ID.
2. Run the workflow manually from the Actions tab so the `google-scholar-stats` branch exists.
3. Uncomment the `fetch_google_scholar_stats.html` include in `_includes/scripts.html`.

## Deployment

GitHub Pages builds the `main` branch from the repository root and publishes the result at [bobbyzhong10.github.io](https://bobbyzhong10.github.io). There is no separate deploy step.

Before pushing:

1. Run a complete Jekyll build locally.
2. Check the affected routes and any note links.
3. Review `git status` and `git diff --check`.
4. Do not commit `_site/`, `vendor/`, `.bundle/`, `.DS_Store`, or other generated caches.

## License

The site's source code is available under the MIT License. See [LICENSE](LICENSE). Written content and personal images are my own work; please get in touch before reusing them.

## Contact

- Email: [bobbyzyz@sas.upenn.edu](mailto:bobbyzyz@sas.upenn.edu)
- Google Scholar: [profile](https://scholar.google.com/citations?user=6YaV4FYAAAAJ&hl)
- ORCID: [0009-0001-8029-7531](https://orcid.org/0009-0001-8029-7531)
- LinkedIn: [yunzhou-zhong](https://www.linkedin.com/in/yunzhou-zhong-659647229)
