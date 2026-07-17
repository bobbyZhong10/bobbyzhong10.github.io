#!/usr/bin/env python3
"""Post-render build step: inline every <img src="assets/fig/*.svg"> into the
rendered HTML as a self-contained, responsive <svg> element, so figures travel
with the page (no dependency on the assets/fig/ folder shipping alongside).

Run AFTER pandoc, from the notes_v2/ directory:  python3 inline_figures.py
Idempotent: once a page's <img> refs are inlined there is nothing left to match.

Handles the two gotchas of inlining many SVGs on one page:
  1. id collisions  -> every id (and #ref) is namespaced per figure (g1_, g2_, ...)
  2. fixed pt sizing -> width/height attrs stripped; viewBox kept; responsive style added
"""
import re, glob, os, base64

FIGDIR = "assets/fig"
# SVGs larger than this that have a rasterized sibling .png (e.g. dense contour
# or scatter plots) are embedded as a PNG data-URI instead of inline vector,
# to keep the self-contained page from ballooning.
RASTER_THRESHOLD = 300_000
counter = [0]

def load_svg(path):
    s = open(path, encoding="utf-8").read()
    s = re.sub(r'<\?xml[^>]*\?>', '', s)
    s = re.sub(r'<!DOCTYPE[^>]*>', '', s, flags=re.IGNORECASE)
    return s.strip()

def namespace_ids(svg, prefix):
    ids = set(re.findall(r'''\bid\s*=\s*['"]([^'"]+)['"]''', svg))
    for i in sorted(ids, key=len, reverse=True):
        svg = re.sub(r'''(\bid\s*=\s*['"])''' + re.escape(i) + r'''(['"])''',
                     r'\g<1>' + prefix + i + r'\g<2>', svg)
        svg = svg.replace(f'url(#{i})', f'url(#{prefix}{i})')
        svg = re.sub(r'''((?:xlink:)?href\s*=\s*['"])#''' + re.escape(i) + r'''(['"])''',
                     r'\g<1>#' + prefix + i + r'\g<2>', svg)
    return svg

def make_responsive(svg, alt):
    m = re.match(r'<svg\b[^>]*>', svg, re.DOTALL)
    tag = m.group(0)
    newtag = re.sub(r'''\s+width\s*=\s*['"][^'"]*['"]''', '', tag)
    newtag = re.sub(r'''\s+height\s*=\s*['"][^'"]*['"]''', '', newtag)
    add = ' role="img"'
    if alt:
        add += ' aria-label="' + alt.replace('"', '&quot;') + '"'
    add += ' style="max-width:100%;height:auto;display:block;margin:0 auto"'
    newtag = newtag[:-1] + add + '>'
    return svg.replace(tag, newtag, 1)

def inline_one(match):
    src = match.group('src')
    alt = match.group('alt') or ''
    path = os.path.join(FIGDIR, os.path.basename(src))
    if not os.path.exists(path):
        return match.group(0)
    counter[0] += 1
    png = os.path.splitext(path)[0] + ".png"
    if os.path.getsize(path) > RASTER_THRESHOLD and os.path.exists(png):
        b64 = base64.b64encode(open(png, "rb").read()).decode("ascii")
        a = ' alt="' + alt.replace('"', '&quot;') + '"' if alt else ''
        return (f'<img src="data:image/png;base64,{b64}"{a} '
                'style="max-width:100%;height:auto;display:block;margin:0 auto" />')
    svg = load_svg(path)
    svg = namespace_ids(svg, f"g{counter[0]}_")
    return make_responsive(svg, alt)

img_re = re.compile(
    r'<img\s+src="(?P<src>assets/fig/[^"]+\.svg)"(?:\s+alt="(?P<alt>[^"]*)")?[^>]*?/>',
    re.DOTALL)

def main():
    total = 0
    for html in sorted(glob.glob("*.html")):
        txt = open(html, encoding="utf-8").read()
        if 'src="assets/fig/' not in txt:
            continue
        counter[0] = 0
        new = img_re.sub(inline_one, txt)
        open(html, "w", encoding="utf-8").write(new)
        total += counter[0]
        print(f"{html}: inlined {counter[0]} figure(s)")
    print(f"done: {total} figures inlined")

if __name__ == "__main__":
    main()
