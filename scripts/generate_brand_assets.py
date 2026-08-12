#!/usr/bin/env python3
"""Generate the site's BZ monogram assets from one reproducible design."""

from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
IMAGES = ROOT / "images"
FONT_PATH = Path("/System/Library/Fonts/Supplemental/Georgia Bold.ttf")

CANVAS = 2048
SCALE = CANVAS / 512
NAVY = "#011F5B"
CRIMSON = "#990000"
WHITE = "#FFFFFF"


def scaled(value: float) -> int:
    return round(value * SCALE)


def draw_master() -> Image.Image:
    image = Image.new("RGBA", (CANVAS, CANVAS), (0, 0, 0, 0))
    draw = ImageDraw.Draw(image)

    draw.rounded_rectangle(
        (scaled(24), scaled(24), scaled(488), scaled(488)),
        radius=scaled(112),
        fill=NAVY,
    )
    draw.rounded_rectangle(
        (scaled(72), scaled(72), scaled(440), scaled(90)),
        radius=scaled(9),
        fill=CRIMSON,
    )

    font = ImageFont.truetype(str(FONT_PATH), scaled(244))
    draw.text(
        (scaled(256), scaled(350)),
        "BZ",
        font=font,
        fill=WHITE,
        anchor="ms",
        stroke_width=0,
    )
    return image


def save_png(master: Image.Image, filename: str, size: int) -> None:
    output = master.resize((size, size), Image.Resampling.LANCZOS)
    output.save(IMAGES / filename, optimize=True)


def save_svg() -> None:
    svg = """<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 512 512" role="img" aria-labelledby="brand-title brand-description">
  <title id="brand-title">BZ monogram</title>
  <desc id="brand-description">White BZ initials on a navy rounded square with a crimson rule.</desc>
  <rect x="24" y="24" width="464" height="464" rx="112" fill="#011f5b"/>
  <rect x="72" y="72" width="368" height="18" rx="9" fill="#990000"/>
  <text x="256" y="350" fill="#ffffff" font-family="Georgia, 'Times New Roman', serif" font-size="244" font-weight="700" text-anchor="middle">BZ</text>
</svg>
"""
    (IMAGES / "brand-mark.svg").write_text(svg, encoding="utf-8")


def main() -> None:
    if not FONT_PATH.exists():
        raise SystemExit(f"Required font not found: {FONT_PATH}")

    master = draw_master()
    save_svg()

    for filename, size in (
        ("favicon-16x16.png", 16),
        ("favicon-32x32.png", 32),
        ("apple-touch-icon.png", 180),
        ("android-chrome-192x192.png", 192),
        ("android-chrome-512x512.png", 512),
        ("mstile-150x150.png", 150),
    ):
        save_png(master, filename, size)

    master.save(
        IMAGES / "favicon.ico",
        format="ICO",
        sizes=[(16, 16), (32, 32), (48, 48), (64, 64)],
    )


if __name__ == "__main__":
    main()
