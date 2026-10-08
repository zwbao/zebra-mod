"""Build the galloping zebra drawn above the prompt (hooks/ui-sprite.ts) from Muybridge's photographs.

Source: Eadweard Muybridge, "The Horse in Motion" (1878), public domain — eleven consecutive
positions of one gallop stride. Each silhouette is cut from its panel, the rider removed (the hump
above the back is cut along a line from the croup to the withers), the frames aligned on the ground
line and the body's centre, reduced to braille dots (2 x 4 per terminal cell) and given zebra
stripes. The sprite is one colour (the terminal's own text colour), so it reads on light and dark
themes alike.

usage: python3 tools/zebra_sprite/build.py [--rows 6] [--out frames.json]
needs: numpy, scipy, pillow (not the mod: only this build step)
"""

from __future__ import annotations

import argparse
import json
import os
import urllib.request

import numpy as np
from PIL import Image
from scipy import ndimage as ndi

URL = "https://upload.wikimedia.org/wikipedia/commons/d/d2/The_Horse_in_Motion_high_res.jpg"
CACHE = os.path.expanduser("~/.cache/zebra-mod/muybridge_horse_in_motion.jpg")
# panel boxes measured on the plate scaled to 1138 x 702 (four columns, three rows; the 12th is standing)
XS = [(27, 289), (300, 561), (571, 835), (846, 1108)]
YS = [(30, 195), (206, 371), (384, 555)]
STRIPE_PERIOD = 5
BITS = [(0, 0, 0x01), (0, 1, 0x02), (0, 2, 0x04), (1, 0, 0x08), (1, 1, 0x10), (1, 2, 0x20), (0, 3, 0x40), (1, 3, 0x80)]


def disk(r: int) -> np.ndarray:
    y, x = np.ogrid[-r:r + 1, -r:r + 1]
    return x * x + y * y <= r * r


def largest(mask: np.ndarray) -> np.ndarray:
    lab, n = ndi.label(mask)
    if n == 0:
        return mask
    return lab == (1 + int(np.argmax(ndi.sum(mask, lab, range(1, n + 1)))))


def silhouettes(image: np.ndarray) -> list:
    h, w = image.shape
    fx, fy = w / 1138, h / 702
    out = []
    for k in range(11):
        x0, x1 = XS[k % 4]
        y0, y1 = YS[k // 4]
        g = image[int(y0 * fy) + 25:int(y1 * fy) - 25, int(x0 * fx) + 25:int(x1 * fx) - 25]
        p = g < 100
        ph, pw = p.shape
        # the ground lines run the full width of the panel; a horse never does
        frac = p.mean(axis=1)
        rows = [y for y in range(int(ph * 0.60), ph) if frac[y] > 0.90]
        gy = rows[0] if rows else ph
        p[gy:, :] = False
        core = largest(ndi.binary_opening(p, structure=disk(6)))       # drops the thin grid lines and numbers
        keep = p & ndi.binary_dilation(core, structure=disk(7))
        holes = ndi.binary_fill_holes(keep) & ~keep                      # small holes only, not the space between legs
        hl, hn = ndi.label(holes)
        if hn:
            sizes = ndi.sum(holes, hl, range(1, hn + 1))
            keep |= np.isin(hl, 1 + np.where(sizes < 2500)[0])
        keep = largest(ndi.binary_opening(keep, structure=disk(2)))
        # the rider: the hump above the back, between the croup and the withers
        cols = np.where(keep.any(axis=0))[0]
        xa, xb = cols[0], cols[-1]
        span = xb - xa
        top = np.full(pw, float(ph))
        for x in cols:
            top[x] = np.argmax(keep[:, x])
        sm = ndi.uniform_filter1d(top, 31, mode="nearest")
        mid = slice(xa + int(span * 0.35), xa + int(span * 0.72))
        xr = mid.start + int(np.argmin(sm[mid]))
        right = slice(xr, xa + int(span * 0.85))
        xrr = right.start + int(np.argmax(sm[right]))
        xl = xr
        while xl > xa + int(span * 0.15):
            xl -= 1
            if sm[xl] - sm[xr] > 40 and abs(sm[xl] - sm[xl - 40]) < 6:
                break
        yl, yr = top[xl], top[xrr]
        for x in range(xl, xrr + 1):
            t = (x - xl) / max(1, xrr - xl)
            keep[:int(yl + (yr - yl) * t + 12 * np.sin(np.pi * t)), x] = False
        keep = largest(ndi.binary_opening(keep, structure=disk(3)))
        # what is left of the ground under the hooves
        ys, xs = np.where(keep)
        low = np.zeros_like(keep)
        low[max(0, ys.max() - 25):ys.max() + 1] = True
        flat = ndi.binary_opening(keep, structure=np.ones((1, 45))) & low & ~ndi.binary_opening(keep, structure=np.ones((9, 9)))
        keep = largest(keep & ~ndi.binary_dilation(flat, structure=np.ones((5, 3))))
        ys, xs = np.where(keep)
        body = keep[ys.min():ys.min() + int((gy - ys.min()) * 0.55)]
        _, bx = np.where(body)
        out.append({"mask": keep, "gy": gy, "cx": float(bx.mean())})
    return out


def dots(frames: list, rows: int) -> list:
    hd = rows * 4
    cx, gy = 2000, 1200
    boxes = []
    for f in frames:
        ys, xs = np.where(f["mask"])
        boxes.append((xs.min() + cx - f["cx"], ys.min() + gy - f["gy"], xs.max() + cx - f["cx"]))
    x0 = min(b[0] for b in boxes)
    y0 = min(b[1] for b in boxes)
    x1 = max(b[2] for b in boxes)
    wd = int(round((x1 - x0) * hd / (gy - y0)))
    wd += wd % 2
    out = []
    for f in frames:
        canvas = np.zeros((int(gy - y0) + 1, int(x1 - x0) + 1), np.float32)
        ys, xs = np.where(f["mask"])
        yy = ys + int(round(gy - f["gy"] - y0))
        xx = xs + int(round(cx - f["cx"] - x0))
        ok = (yy >= 0) & (yy < canvas.shape[0]) & (xx >= 0) & (xx < canvas.shape[1])
        canvas[yy[ok], xx[ok]] = 1
        small = Image.fromarray((canvas * 255).astype(np.uint8)).resize((wd, hd), Image.BOX)
        out.append(np.asarray(small, np.float32) / 255.0 > 0.38)
    return out


def clean(m: np.ndarray) -> np.ndarray:
    m = m.copy()
    h, w = m.shape
    for y in range(h - 3, h):                    # ground runs under the hooves
        x = 0
        while x < w:
            if m[y, x]:
                j = x
                while j < w and m[y, j]:
                    j += 1
                if j - x > 3 and m[y - 1, x:j].sum() <= (j - x) // 3:
                    m[y, x:j] = False
                x = j
            else:
                x += 1
    lab, n = ndi.label(m, structure=np.ones((3, 3)))   # specks not joined to the horse
    if n > 1:
        sizes = ndi.sum(m, lab, range(1, n + 1))
        big = 1 + int(np.argmax(sizes))
        for k in range(1, n + 1):
            if k != big and sizes[k - 1] < 6:
                m[lab == k] = False
    return m


def zebra(m: np.ndarray) -> np.ndarray:
    m = clean(m)
    h, w = m.shape
    on = m.copy()
    cols = np.where(m.any(axis=0))[0]
    xa, xb = cols[0], cols[-1]
    span = xb - xa + 1
    top = np.array([np.argmax(m[:, x]) if m[:, x].any() else h for x in range(w)])
    back = min(top[x] for x in range(xa + int(span * 0.30), xa + int(span * 0.65)))
    belly = back + int(round(0.50 * (h - back)))
    neck0, head0, tail1 = xa + int(span * 0.66), xa + int(span * 0.84), xa + int(span * 0.20)
    for y in range(h):
        for x in range(w):
            if not m[y, x] or x < tail1 or x >= head0 or y >= belly:   # tail, head and legs stay solid
                continue
            if x >= neck0:
                phase = (x + (y * 2) // 3) % STRIPE_PERIOD             # neck stripes lean forward
            elif x < xa + int(span * 0.40):
                phase = (x + y // 2) % STRIPE_PERIOD                   # rump stripes lean back
            else:
                phase = x % STRIPE_PERIOD
            if phase == 0:
                on[y, x] = False
    head = [x for x in range(head0, xb + 1) if m[:, x].any()]
    if len(head) >= 6:                                                 # the eye
        ex = head[-6]
        if top[ex] + 2 < h:
            on[top[ex] + 2, ex] = False
    for x in range(neck0 + 1, head0 - 1, 2):                           # the mane, above the neck
        if top[x] - 1 >= 0:
            on[top[x] - 1, x] = True
    return on


def braille(on: np.ndarray) -> list:
    h, w = on.shape
    rows = []
    for r in range(h // 4):
        row = ""
        for c in range(w // 2):
            b = 0
            for dx, dy, bit in BITS:
                if on[r * 4 + dy, c * 2 + dx]:
                    b |= bit
            row += chr(0x2800 + b)
        rows.append(row)
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rows", type=int, default=6)
    ap.add_argument("--out", default="zebra_frames.json")
    a = ap.parse_args()
    if not os.path.exists(CACHE):
        os.makedirs(os.path.dirname(CACHE), exist_ok=True)
        req = urllib.request.Request(URL, headers={"User-Agent": "zebra-mod sprite build"})
        with urllib.request.urlopen(req) as r, open(CACHE, "wb") as f:
            f.write(r.read())
    image = np.asarray(Image.open(CACHE).convert("L")).astype(np.float32)
    frames = [braille(zebra(m)) for m in dots(silhouettes(image), a.rows)]
    with open(a.out, "w", encoding="utf-8") as f:
        json.dump(frames, f, ensure_ascii=False, indent=0)
    print(f"{len(frames)} frames, {len(frames[0])} rows x {len(frames[0][0])} cells -> {a.out}")


if __name__ == "__main__":
    main()
