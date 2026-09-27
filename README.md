# Yunzhou (Bobby) Zhong — Academic Website

Source for [bobbyzhong10.github.io](https://bobbyzhong10.github.io), a static site built with [Jekyll](https://jekyllrb.com/) and deployed by GitHub Pages.

## Local Preview

Requires Ruby and Bundler.

```bash
bundle install
bundle exec jekyll serve
```

Then open <http://127.0.0.1:4000/>. `./run_server.sh` starts the same server with live reload. Restart the server after editing `_config.yml`.

If the build stops with `Invalid US-ASCII character`, run it with a UTF-8 locale: `LANG=en_US.UTF-8 bundle exec jekyll serve`.

## Where Things Live

| Path | Contents |
| --- | --- |
| `_pages/` | Page content: home, research, education, awards and service, teaching, CV, notes |
| `_config.yml` | Site settings, author profile and links, Google Analytics ID, build exclusions |
| `_data/navigation.yml` | Top navigation |
| `_layouts/`, `_includes/` | Page shell, header, home hero, SEO tags, and page scripts |
| `_sass/_site.scss` | All site styles, including the light and dark color tokens |
| `_sass/` (other files) | Base theme modules |
| `assets/` | Stylesheet entry point, icon CSS, and icon fonts |
| `images/` | Portrait, favicons, and brand mark |
| `notes/` | Research-note collections, published as pre-built HTML |
| `scripts/` | Generator for the favicon and brand-mark images |
| `docs/` | Internal documentation (not published) |

## Common Edits

- **Page text:** edit the matching file in `_pages/`.
- **Name, affiliation, and profile links:** the `author` block in `_config.yml`.
- **Styles and colors:** `_sass/_site.scss`. Reuse the color tokens at the top of the file rather than adding new colors.
- **Abstracts and paper links:** `_pages/research-projects.md`. Each abstract and link button sits under its paper; copy an existing entry's markup.
- **Analytics:** the Google Analytics measurement ID is `google_analytics_id` in `_config.yml`; the visitor badge is in `_layouts/default.html`.

[`docs/SITE_ARCHITECTURE.md`](docs/SITE_ARCHITECTURE.md) explains how the pieces fit together and how the notes are built.

## Deployment

Pushing to `main` triggers a GitHub Pages build, and the live site updates within a minute or two. The generated `_site/` folder is not committed.

## License

The site's source code is released under the [MIT License](LICENSE). It is adapted from [AcadHomepage](https://github.com/RayeRen/acad-homepage.github.io) by Yi Ren, which builds on [Minimal Mistakes](https://github.com/mmistakes/minimal-mistakes) by Michael Rose.

The MIT License covers the code only. All site content, including text, photographs, the CV, research notes, and papers and their abstracts, is © 2025–2026 Yunzhou (Bobby) Zhong and co-authors where applicable. All rights reserved; it may not be reused without permission.
