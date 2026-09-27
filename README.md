# Personal Academic Website

Source code for [bobbyzhong10.github.io](https://bobbyzhong10.github.io), built with [Jekyll](https://jekyllrb.com/) and hosted on GitHub Pages.

## Getting Started

Requires Ruby and Bundler.

```bash
bundle install
bundle exec jekyll serve
```

Open `http://127.0.0.1:4000/` to preview the site locally. Restart the server after editing `_config.yml`.

## Structure

| Path | Purpose |
| --- | --- |
| `_config.yml` | Site settings and author profile |
| `_data/navigation.yml` | Top navigation menu |
| `_pages/` | Page content (home, research, education, awards, teaching, CV, notes) |
| `_layouts/`, `_includes/` | Page layout and shared components |
| `_sass/`, `assets/` | Styles, scripts, and fonts |
| `images/` | Portrait, favicons, and brand mark |
| `notes/` | Standalone note collections |
| `docs/` | Internal documentation (not published) |

## Customization

- Edit page content in `_pages/`.
- Update name, affiliation, and profile links in `_config.yml`.
- Adjust styles in `_sass/_site.scss`.
- See `docs/SITE_ARCHITECTURE.md` for a more detailed guide.

## Deployment

Pushing to `main` triggers a GitHub Pages build. Generated output in `_site/` is not committed.

## License

The site's source code is released under the [MIT License](LICENSE). It is adapted from [AcadHomepage](https://github.com/RayeRen/acad-homepage.github.io) by Yi Ren, which builds on [Minimal Mistakes](https://github.com/mmistakes/minimal-mistakes) by Michael Rose.

The MIT License covers the code only. All site content, including text, photographs, the CV, research notes, and papers and their abstracts, is © 2025–2026 Yunzhou (Bobby) Zhong and co-authors where applicable. All rights reserved; it may not be reused without permission.
