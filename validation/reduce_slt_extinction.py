"""Reduce the four SLT/SN2024ggi airmass-sweep nights to per-frame zero points.

`slt.py` reduced one night (2024-04-14) and found it could not measure absolute
extinction: the target sets all night, so airmass and time run together and a
transparency drift is indistinguishable from an airmass term. The obvious next
move is more nights — 2024-04-12/13/14/16 each sweep airmass 1.81 to ~3.8 in
five bands. This script reduces all four so `analyze_slt_extinction.py` can ask
whether a *joint* fit (one extinction shared across nights, a free transparency
term per night) breaks the degeneracy one night could not.

It does not (see `analyze_slt_extinction.py` and `slt.MULTINIGHT`). The value of
running it anyway is that the failure is now quantified and reproducible rather
than argued from a single night.

Method, following slt.py where it can:

* Reference magnitudes are SkyMapper DR4 PSF mags (`_psf`), AB-ish, for stellar
  sources (`class_star` > 0.9, `flags` == 0). The field is at dec -32.8, outside
  Pan-STARRS DR2, so SkyMapper is used with NO colour-term transform applied —
  the same caveat slt.py carries. Cached to `skymapper_sn2024ggi_refstars.csv`.
* Per frame: project the catalogue onto the frame's WCS, keep references in the
  magnitude window 14-18.5 (below saturation, above the noise), aperture
  photometry with a per-frame FWHM-scaled aperture, per-star ZP = m_ref - m_inst
  where m_inst = -2.5 log10(net / exptime), frame ZP = 2.5-sigma-clipped mean.
* Frame FWHM is a coarse second-moment estimate, deliberately cheap. It sets
  only the aperture *scale*; the extinction result is a ZP-vs-airmass slope and
  a constant aperture bias cancels into the per-frame zero point, so this
  crudeness does not touch the conclusion. It is not the FWHM slt.py quotes.

Frames are ~4.6 GB and stay out of the repository (`data/raw/*` is gitignored).
What is committed is the output, `slt_sn2024ggi_zp_per_frame.csv`, so the fit in
`analyze_slt_extinction.py` reproduces without them.

Run: `uv run --with numpy --with pandas --with astropy --with photutils \
      python validation/reduce_slt_extinction.py`
"""
import glob
import io
import os
import re
import urllib.parse
import urllib.request
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
from astropy.io import fits
from astropy.stats import sigma_clip, sigma_clipped_stats
from astropy.wcs import WCS
from photutils.aperture import (CircularAnnulus, CircularAperture,
                                aperture_photometry)

HERE = os.path.dirname(os.path.abspath(__file__))
FRAMES = os.path.join(HERE, "data", "raw", "slt_sn2024ggi")
REFSTARS = os.path.join(HERE, "data", "skymapper_sn2024ggi_refstars.csv")
OUT = os.path.join(HERE, "data", "slt_sn2024ggi_zp_per_frame.csv")

FIELD = dict(ra=169.5762, dec=-32.8470, radius=0.15)
NIGHTS = ["2024-04-12", "2024-04-13", "2024-04-14", "2024-04-16"]
BAND_OF = {"up": "u", "gp": "g", "rp": "r", "ip": "i", "zp": "z"}
MAG_WINDOW = (14.0, 18.5)   # ref mag: avoid saturation (bright) and low SNR (faint)
EDGE_MARGIN = 25            # px kept clear of the frame edge
SATURATION_ADU = 58000      # peak-pixel guard, below the 16-bit ceiling


def fetch_refstars():
    """SkyMapper DR4 cone search, cached. Returns the reference-star table."""
    if os.path.exists(REFSTARS):
        return pd.read_csv(REFSTARS)
    url = "https://skymapper.anu.edu.au/sm-cone/public/query?" + urllib.parse.urlencode(
        dict(RA=FIELD["ra"], DEC=FIELD["dec"], SR=FIELD["radius"], RESPONSEFORMAT="CSV"))
    raw = urllib.request.urlopen(url, timeout=90).read().decode()
    df = pd.read_csv(io.StringIO(raw))
    keep = (df["class_star"] > 0.90) & (df["flags"] == 0) & (df["ngood"] >= 1)
    sel = df[keep]
    out = pd.DataFrame({"ra": sel["raj2000"], "dec": sel["dej2000"],
                        "class_star": sel["class_star"]})
    for b in "ugriz":
        out[b] = sel[f"{b}_psf"]
        out[f"e_{b}"] = sel[f"e_{b}_psf"]
    out = out[out[["g", "r", "i"]].notna().all(axis=1)]
    out.to_csv(REFSTARS, index=False)
    return out


def _fwhm(data, xs, ys, bkg):
    """Coarse second-moment FWHM in px, median over the reference stars."""
    sig = []
    for x, y in zip(xs, ys):
        xi, yi = int(round(x)), int(round(y))
        if xi < 8 or yi < 8 or xi > data.shape[1] - 9 or yi > data.shape[0] - 9:
            continue
        cut = data[yi - 8:yi + 9, xi - 8:xi + 9] - bkg
        if cut.max() <= 0:
            continue
        yy, xx = np.mgrid[0:cut.shape[0], 0:cut.shape[1]]
        w = np.clip(cut, 0, None)
        tot = w.sum()
        if tot <= 0:
            continue
        cx = (w * xx).sum() / tot
        cy = (w * yy).sum() / tot
        sx = np.sqrt(max((w * (xx - cx) ** 2).sum() / tot, 0.1))
        sy = np.sqrt(max((w * (yy - cy) ** 2).sum() / tot, 0.1))
        s = 0.5 * (sx + sy)
        if 0.7 < s < 6:
            sig.append(s)
    if len(sig) < 3:
        return None
    return float(np.median(sig)) * 2.3548


def _night_of(path):
    m = re.search(r"D(\d{4}-\d{2}-\d{2})T", os.path.basename(path))
    return m.group(1) if m else "?"


def reduce_frame(path, band, refs):
    with fits.open(path) as hd:
        h = hd[0].header
        data = hd[0].data.astype(float)
    airmass = h.get("AIRMASS")
    mjd = h.get("MJD-OBS") or h.get("JD-OBS")
    if airmass is None:
        return None
    airmass = float(airmass)
    exp = float(h.get("EXPTIME", 1))
    try:
        wcs = WCS(h)
    except Exception:
        return None

    col = BAND_OF[band]
    r = refs[refs[col].notna() & refs[f"e_{col}"].notna()]
    r = r[(r[col] > MAG_WINDOW[0]) & (r[col] < MAG_WINDOW[1]) & (r[f"e_{col}"] < 0.05)]
    if len(r) < 4:
        return None
    x, y = wcs.all_world2pix(r["ra"].values, r["dec"].values, 0)
    inb = ((x > EDGE_MARGIN) & (x < data.shape[1] - EDGE_MARGIN)
           & (y > EDGE_MARGIN) & (y < data.shape[0] - EDGE_MARGIN))
    x, y, mref = x[inb], y[inb], r[col].values[inb]
    if len(x) < 4:
        return None

    _, med, _ = sigma_clipped_stats(data, sigma=3.0)
    fwhm = _fwhm(data, x, y, med)
    if fwhm is None or not np.isfinite(fwhm):
        return None

    pos = np.transpose((x, y))
    ap = CircularAperture(pos, r=2.5 * fwhm)
    an = CircularAnnulus(pos, r_in=4.0 * fwhm, r_out=6.0 * fwhm)
    flux = aperture_photometry(data, ap)["aperture_sum"].value
    bkg = np.array([
        np.median(m.multiply(data)[m.data > 0]) if m is not None else med
        for m in an.to_mask(method="center")])
    net = flux - bkg * ap.area
    peak = np.array([
        np.nanmax(m.multiply(data)) if m is not None else np.nan
        for m in ap.to_mask(method="center")])

    good = (net > 0) & (peak < SATURATION_ADU) & np.isfinite(net)
    net, mref = net[good], mref[good]
    if len(net) < 4:
        return None
    zp = mref - (-2.5 * np.log10(net / exp))
    clip = sigma_clip(zp, sigma=2.5, maxiters=3)
    zp = zp[~clip.mask]
    if len(zp) < 4:
        return None
    return dict(night=_night_of(path), band=band, mjd=float(mjd),
                airmass=airmass, fwhm_px=round(fwhm, 2), n_ref=len(zp),
                zp=round(float(np.mean(zp)), 4),
                zp_err=round(float(np.std(zp) / np.sqrt(len(zp))), 4))


def main():
    refs = fetch_refstars()
    rows = []
    for night in NIGHTS:
        for band in BAND_OF:
            for f in sorted(glob.glob(f"{FRAMES}/*{night}*_{band}_*.fits")):
                try:
                    row = reduce_frame(f, band, refs)
                except Exception:
                    row = None
                if row:
                    rows.append(row)
        print(f"{night}: {sum(r['night'] == night for r in rows)} frames", flush=True)
    df = pd.DataFrame(rows).sort_values(["band", "mjd"])
    df.to_csv(OUT, index=False)
    print(f"wrote {OUT}: {len(df)} frames", flush=True)


if __name__ == "__main__":
    main()
