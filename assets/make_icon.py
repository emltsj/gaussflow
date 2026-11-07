"""Generate the application icon (assets/icon.ico) — run once before packaging."""
from __future__ import annotations

from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

OUT = Path(__file__).resolve().parent / "icon.ico"


def make(size: int = 256) -> Image.Image:
    img = Image.new("RGBA", (size, size), (0, 0, 0, 0))
    d = ImageDraw.Draw(img)
    r = int(size * 0.18)
    d.rounded_rectangle([0, 0, size - 1, size - 1], radius=r, fill=(31, 78, 120, 255))
    # molecule motif: three nodes + bonds
    cx, cy = size * 0.42, size * 0.52
    nodes = [(cx, cy), (cx + size * 0.22, cy - size * 0.18),
             (cx + size * 0.05, cy + size * 0.24)]
    for a in range(len(nodes)):
        for b in range(a + 1, len(nodes)):
            d.line([nodes[a], nodes[b]], fill=(214, 228, 240, 255),
                   width=max(2, int(size * 0.022)))
    for (x, y) in nodes:
        rr = size * 0.055
        d.ellipse([x - rr, y - rr, x + rr, y + rr], fill=(124, 227, 139, 255))
    # "G" letter
    try:
        font = ImageFont.truetype("arialbd.ttf", int(size * 0.42))
    except OSError:
        font = ImageFont.load_default()
    d.text((size * 0.60, size * 0.12), "G", font=font, fill="white")
    return img


def main() -> None:
    img = make(256)
    img.save(OUT, sizes=[(16, 16), (24, 24), (32, 32), (48, 48),
                         (64, 64), (128, 128), (256, 256)])
    print(f"icon written: {OUT} ({OUT.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
