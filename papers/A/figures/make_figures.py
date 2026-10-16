"""Paper A figures (journal format): TIFF (LZW, 600 dpi, RGB) plus PNG previews in papers/A/figures/.

Usage (from the repository root, after `cd models && Rscript A/run_all.R`):
    python papers/A/figures/make_figures.py

Reads only the figure-data CSVs written by models/A/run_all.R (papers/A/figures/data/Figure_<n>.csv),
plus the GADM boundaries and the unit frame for Supplementary Figure S1 (ruling P19: the uniform unit-IRR
map moved from the main text to the supplement; main figures are 1-3). Journal sizes: double column 17.5 cm wide,
single column 8.5 cm, maximum height 21 cm.
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import TwoSlopeNorm  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
DATA = HERE / "data"
CM = 1 / 2.54
DPI = 600

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
                     "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7, "axes.spines.top": False,
                     "axes.spines.right": False})
RED, BLUE, GREY, DARK = "#b2182b", "#2166ac", "#9e9e9e", "#333333"


def save(fig, n):
    tif = HERE / f"Figure_{n}.tif"
    png = HERE / f"Figure_{n}.png"
    fig.savefig(png, dpi=200, bbox_inches="tight", facecolor="white")
    fig.savefig(tif, dpi=DPI, bbox_inches="tight", facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    # force RGB (no alpha) and record dpi
    from PIL import Image
    with Image.open(tif) as im:
        im.convert("RGB").save(tif, compression="tiff_lzw", dpi=(DPI, DPI))
    plt.close(fig)
    print("wrote", tif.relative_to(ROOT), "and", png.name)


def fig1():
    d = pd.read_csv(DATA / "Figure_1.csv").sort_values("year")
    fig, ax = plt.subplots(3, 1, figsize=(17.5 * CM, 15 * CM), sharex=True,
                           gridspec_kw={"height_ratios": [1.1, 1, 1]})
    # fold periods
    folds = d.dropna(subset=["fold_label"]).groupby("fold_label")["year"].agg(["min", "max"])
    for a in ax:
        for _, r in folds.iterrows():
            a.axvspan(r["min"] - 0.5, r["max"] + 0.5, color="#f3e3c3", zorder=0, lw=0)
    a = ax[0]
    a.bar(d["year"], d["oni_jjas"], color=np.where(d["oni_jjas"] >= 0, RED, BLUE), width=0.7, label="ONI (JJAS)")
    a.plot(d["year"], d["roni_jjas"], color=DARK, marker="o", ms=2.5, lw=1, label="RONI (JJAS)")
    a.axhline(0, color=DARK, lw=0.5)
    a.set_ylabel("Index, °C")
    a.legend(loc="upper left", frameon=False, ncol=2)
    a.set_title("(a) June–September ENSO indices; shaded: leave-one-El Niño-out folds", loc="left")
    a = ax[1]
    a.plot(d["year"], d["mal_inc_per_1000"], color=DARK, marker="o", ms=3, lw=1.2)
    a.set_ylabel("Cases per 1,000")
    a.set_title("(b) National malaria incidence (2006 masked: partial-year source)", loc="left")
    a.set_ylim(bottom=0)
    a = ax[2]
    dd = d[d["year"] >= 2007]
    a.plot(dd["year"], dd["den_inc_per_100k"], color=DARK, marker="s", ms=3, lw=1.2)
    a.set_ylabel("Cases per 100,000")
    a.set_title("(c) National dengue incidence (units with missing counts excluded in that year)", loc="left")
    a.set_ylim(bottom=0)
    a.set_xlim(2001.4, 2025.6)
    a.set_xticks(range(2002, 2026, 2))
    a.set_xlabel("Year")
    fig.align_ylabels(ax)
    fig.tight_layout()
    save(fig, 1)


def unit_polygons():
    sys.path.insert(0, str(ROOT))
    from pipeline.boundaries import load
    st = load("state")
    frame = pd.read_csv(ROOT / "pipeline" / "reference" / "state_frame.csv")
    m = {gid: r.unit_id for r in frame.itertuples() for gid in r.gid1_members.split("|")}
    st["unit_id"] = st["region_id"].map(m)
    st = st.dropna(subset=["unit_id"])
    u = st.dissolve(by="unit_id", as_index=False)
    u["geometry"] = u.geometry.simplify(0.02)
    return u


def figS1():
    """Supplementary Figure S1 (ruling P19; was Figure 2): unit-specific IRR maps."""
    d = pd.read_csv(DATA / "Figure_S1.csv")
    u = unit_polygons().merge(d, on="unit_id", how="left")
    vals = np.log(pd.concat([d["irr_lag0"], d["irr_lag1"]]))
    lim = max(float(np.nanmax(np.abs(vals))) * 1.05, np.log(1.2))
    norm = TwoSlopeNorm(vmin=np.exp(-lim), vcenter=1.0, vmax=np.exp(lim))
    fig, ax = plt.subplots(1, 2, figsize=(17.5 * CM, 11 * CM))
    for a, col, title in ((ax[0], "irr_lag0", "(a) Lag 0 (same year)"), (ax[1], "irr_lag1", "(b) Lag 1 (following year)")):
        missing = u[u[col].isna()]
        u.dropna(subset=[col]).plot(column=col, cmap="RdBu_r", norm=norm, ax=a, edgecolor="#555555", linewidth=0.25)
        if len(missing):
            missing.plot(ax=a, facecolor="white", edgecolor="#555555", hatch="////", linewidth=0.25)
        a.set_title(title, loc="left")
        a.set_axis_off()
    sm = plt.cm.ScalarMappable(cmap="RdBu_r", norm=norm)
    cb = fig.colorbar(sm, ax=ax, orientation="horizontal", fraction=0.05, pad=0.04, shrink=0.6)
    cb.set_label("Unit-specific incidence rate ratio per 1 °C ONI (1 = no association)")
    save(fig, "S1")


def fig2():
    """Figure 2 (was Figure 3): Model C / Model B MALE ratio by fold and origin."""
    d = pd.read_csv(DATA / "Figure_2.csv", dtype={"label": str, "split_id": str})
    fig, ax = plt.subplots(2, 1, figsize=(8.5 * CM, 14 * CM))
    styles = {"malaria": dict(color=DARK, marker="o"), "dengue": dict(color=RED, marker="s", mfc="white")}
    for a, scheme, title in ((ax[0], "loeo", "(a) Leave-one-El Niño-out folds"),
                             (ax[1], "roll", "(b) Rolling-origin forecast years")):
        s = d[d["scheme"] == scheme].copy()
        s["key"] = pd.to_numeric(s["split_id"]) if scheme == "roll" else pd.to_numeric(s["label"].str[:4])
        labs = s.drop_duplicates("key").sort_values("key")
        pos = {k: i for i, k in enumerate(labs["key"])}
        for o, st in styles.items():
            so = s[s["outcome"] == o].sort_values("key")
            off = -0.12 if o == "malaria" else 0.12
            a.plot([pos[k] + off for k in so["key"]], so["ratio"], ls="none", ms=4, label=o.capitalize(), **st)
        a.axhline(1, color=GREY, ls="--", lw=0.8)
        a.set_xticks(range(len(labs)))
        a.set_xticklabels(labs["label"] if scheme == "loeo" else labs["split_id"], rotation=45 if scheme == "roll" else 0)
        a.set_ylabel("MALE ratio, Model C / Model B")
        a.set_title(title, loc="left")
        a.legend(frameon=False, loc="best")
    fig.tight_layout()
    save(fig, 2)


def fig3():
    """Figure 3 (was Figure 4): expected excess cases by unit, 2026 and 2027."""
    d = pd.read_csv(DATA / "Figure_3.csv")
    order = d[(d["index"] == "ONI") & (d["year"] == 2027)].sort_values("excess_mean")["unit_name"].tolist()
    ypos = {u: i for i, u in enumerate(order)}
    fig, ax = plt.subplots(1, 2, figsize=(17.5 * CM, 18 * CM), sharey=True)
    for a, yr in zip(ax, (2026, 2027)):
        for ix, off, mfc in (("ONI", 0.15, DARK), ("RONI", -0.15, "white")):
            s = d[(d["year"] == yr) & (d["index"] == ix)]
            y = np.array([ypos[u] for u in s["unit_name"]]) + off
            a.hlines(y, s["excess_lo"], s["excess_hi"], color=DARK if ix == "ONI" else GREY, lw=0.9)
            a.plot(s["excess_mean"], y, ls="none", marker="o", ms=3.2, mec=DARK, mfc=mfc,
                   label=f"{ix} model")
        a.axvline(0, color=GREY, lw=0.6)
        a.set_xscale("symlog", linthresh=10)
        a.set_title(f"({'a' if yr == 2026 else 'b'}) {yr} ({'lag 0' if yr == 2026 else 'lag 1'})", loc="left")
        a.grid(axis="x", color="#e5e5e5", lw=0.5)
    ax[0].set_yticks(range(len(order)))
    ax[0].set_yticklabels(order)
    ax[0].legend(frameon=False, loc="lower right")
    fig.supxlabel("Expected excess malaria cases relative to ONI = 0 (symmetric log scale)", fontsize=8)
    fig.tight_layout()
    save(fig, 3)


if __name__ == "__main__":
    for f in (fig1, fig2, fig3, figS1):
        f()
