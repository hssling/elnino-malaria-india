"""Paper A revision: maps, Moran's I / LISA and modified Mann-Kendall trends (reviewer analyses R2.2 and R2.3).

Run from the repository root:
    python -I papers/A/revision/spatial_trends.py [--unit-irr models/A/output/rev_unit_irr_unpooled.csv]

Writes (papers/A/revision/): unit_weights_queen.csv, moran_results.csv, lisa_clusters.csv, trends_mk.csv,
sir_by_unit.csv, spatial_trend_keys.json; (papers/A/figures/): Figure_S6..S10 (TIFF 600 dpi RGB + PNG) and
papers/A/figures/data/figure_S*.csv. Permutations: 9,999, seed 20261117.
"""
import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from scipy import stats  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
REV = ROOT / "papers" / "A" / "revision"
FIG = ROOT / "papers" / "A" / "figures"
sys.path.insert(0, str(ROOT))

SEED = 20261117
NPERM = 9999
CM = 1 / 2.54
DPI = 600
WEIGHTS_CSV = REV / "unit_weights_queen.csv"
MAL_PERIODS = [(2002, 2009), (2010, 2017), (2018, 2025)]
DEN_PERIODS = [(2007, 2017), (2018, 2025)]
# constants added before the log in the trend analysis (incidence units per 1,000 / per 100,000)
LOG_CONST = {"malaria": 0.01, "dengue": 0.1}
SCALE = {"malaria": 1000.0, "dengue": 1e5}

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
                     "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7, "axes.spines.top": False,
                     "axes.spines.right": False})


# ----------------------------------------------------------------------------------- data
def frame():
    return pd.read_csv(ROOT / "pipeline" / "reference" / "state_frame.csv")


def panel(outcome):
    d = pd.read_csv(ROOT / "data" / "derived" / "model" / f"{outcome}_model.csv",
                    usecols=["unit_id", "year", "cases", "population"])
    return d.sort_values(["unit_id", "year"]).reset_index(drop=True)


def unit_polygons(simplify=None):
    """GADM 4.1 level-1 polygons dissolved to the 34 pooled units (WGS84)."""
    from pipeline.boundaries import load
    st = load("state")
    fr = frame()
    m = {g: r.unit_id for r in fr.itertuples() for g in r.gid1_members.split("|")}
    st["unit_id"] = st["region_id"].map(m)
    st = st.dropna(subset=["unit_id"])
    u = st.dissolve(by="unit_id", as_index=False)[["unit_id", "geometry"]]
    if simplify:
        u["geometry"] = u.geometry.simplify(simplify)
    return u


# ----------------------------------------------------------------------------------- weights
def build_weights(buffer_m=1000, crs="EPSG:7755"):
    """Queen contiguity of dissolved units. Each polygon is buffered by buffer_m/2 in a metric CRS so that gaps
    (slivers) up to buffer_m close; returns the edge list (unit_i < unit_j alphabetically by frame order)."""
    u = unit_polygons().to_crs(crs)
    g = u.geometry.buffer(buffer_m / 2.0)
    ids = u["unit_id"].tolist()
    edges = []
    for a in range(len(ids)):
        for b in range(a + 1, len(ids)):
            if g.iloc[a].intersects(g.iloc[b]):
                edges.append((ids[a], ids[b]))
    return pd.DataFrame(edges, columns=["unit_i", "unit_j"]).sort_values(["unit_i", "unit_j"]).reset_index(drop=True)


def read_weights(weights_csv):
    e = pd.read_csv(weights_csv)
    return list(zip(e["unit_i"], e["unit_j"]))


# ----------------------------------------------------------------------------------- Moran
def _weight_matrix(units, edges):
    idx = {u: k for k, u in enumerate(units)}
    W = np.zeros((len(units), len(units)))
    for a, b in edges:
        if a in idx and b in idx:
            W[idx[a], idx[b]] = W[idx[b], idx[a]] = 1.0
    return W


def _prepare(values, weights_csv):
    """Drop units with missing values, then units without neighbours among the rest; row-standardise."""
    vals = {k: float(v) for k, v in values.items() if v is not None and np.isfinite(v)}
    units = sorted(vals)
    W = _weight_matrix(units, read_weights(weights_csv))
    keep = W.sum(axis=1) > 0
    units = [u for u, k in zip(units, keep) if k]
    W = W[np.ix_(keep, keep)]
    Wr = W / W.sum(axis=1, keepdims=True)
    x = np.array([vals[u] for u in units])
    return units, Wr, x


def _moran_stat(z, Wr):
    return (len(z) / Wr.sum()) * (z @ Wr @ z) / (z @ z)


def moran_i(values, weights_csv=WEIGHTS_CSV, nperm=NPERM, seed=SEED, detail=False):
    """Global Moran's I (row-standardised queen weights; units with missing values or no neighbours dropped).

    Two-sided permutation p = (1 + #{|I* - E| >= |I - E|}) / (nperm + 1), E = -1/(n-1), over nperm random
    reassignments of the values to the units. Returns (I, p) or, with detail=True, (I, p, n_used, units_used).
    """
    units, Wr, x = _prepare(values, weights_csv)
    n = len(x)
    z = x - x.mean()
    obs = _moran_stat(z, Wr)
    rng = np.random.default_rng(seed)
    perm = np.argsort(rng.random((nperm, n)), axis=1)
    zp = z[perm]
    num = np.einsum("pi,ij,pj->p", zp, Wr, zp)
    sim = (n / Wr.sum()) * num / (z @ z)
    e = -1.0 / (n - 1)
    p = (1 + np.sum(np.abs(sim - e) >= abs(obs - e) - 1e-12)) / (nperm + 1)
    return (obs, float(p), n, units) if detail else (obs, float(p))


def lisa(values, weights_csv=WEIGHTS_CSV, nperm=NPERM, seed=SEED):
    """Local Moran with conditional permutation (unit i held fixed, its neighbour values drawn from the others)."""
    units, Wr, x = _prepare(values, weights_csv)
    n = len(x)
    z = (x - x.mean()) / x.std()
    lag = Wr @ z
    Ii = z * lag
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(n):
        k = int((Wr[i] > 0).sum())
        others = np.delete(z, i)
        draw = np.argsort(rng.random((nperm, n - 1)), axis=1)[:, :k]
        sim = z[i] * others[draw].mean(axis=1)
        p = (1 + np.sum(np.abs(sim) >= abs(Ii[i]) - 1e-12)) / (nperm + 1)
        if z[i] >= 0 and lag[i] >= 0:
            q = "HH"
        elif z[i] < 0 and lag[i] < 0:
            q = "LL"
        elif z[i] >= 0:
            q = "HL"
        else:
            q = "LH"
        rows.append((units[i], x[i], z[i], lag[i], Ii[i], k, p, q, q if p < 0.05 else "NS"))
    return pd.DataFrame(rows, columns=["unit_id", "value", "z", "lag_z", "local_I", "n_neighbours", "p_perm", "quadrant",
                                       "cluster"])


# ----------------------------------------------------------------------------------- incidence / SIR
def period_incidence(d, y0, y1, scale):
    s = d[(d.year >= y0) & (d.year <= y1) & d.cases.notna()]
    g = s.groupby("unit_id").agg(cases=("cases", "sum"), pop=("population", "sum"), n_years=("year", "size"))
    g["inc"] = g["cases"] / g["pop"] * scale
    return g


def sir_table(d):
    s = d[d.cases.notna()].copy()
    tot = s.groupby("year")[["cases", "population"]].sum()
    s["expected"] = s["population"] * s["year"].map(tot["cases"] / tot["population"])
    g = s.groupby("unit_id").agg(observed=("cases", "sum"), expected=("expected", "sum"))
    g["sir"] = g["observed"] / g["expected"]
    return g


# ----------------------------------------------------------------------------------- modified Mann-Kendall
def sen_slope(t, x):
    t = np.asarray(t, float)
    x = np.asarray(x, float)
    i, j = np.triu_indices(len(x), 1)
    dt = t[j] - t[i]
    ok = dt != 0
    slope = float(np.median((x[j] - x[i])[ok] / dt[ok]))
    return slope, float(np.median(x - slope * t))


def _mk_s_var(x):
    n = len(x)
    i, j = np.triu_indices(n, 1)
    s = float(np.sum(np.sign(x[j] - x[i])))
    _, cnt = np.unique(x, return_counts=True)
    tie = np.sum(cnt * (cnt - 1) * (2 * cnt + 5))
    return s, (n * (n - 1) * (2 * n + 5) - tie) / 18.0


def _rank_acf(r, nlags):
    r = r - r.mean()
    den = np.sum(r * r)
    return np.array([np.sum(r[: len(r) - k] * r[k:]) / den for k in range(1, nlags + 1)])


def mann_kendall(t, x, modified=True):
    """Mann-Kendall test with Sen's slope; modified=True applies the Hamed and Rao (1998) variance correction.

    The series is detrended with Sen's slope; the autocorrelations rho_k of the ranks of the detrended series, at
    lags whose 95% limits {-1 +/- 1.96 sqrt(n-k-1)}/(n-k) exclude the estimate, give the factor
    n/n* = 1 + 2/(n(n-1)(n-2)) * sum_k (n-k)(n-k-1)(n-k-2) rho_k and Var*(S) = Var(S) * n/n*.
    Returns dict(n, S, tau, z, p, slope, intercept, correction, sig_lags). p is two-sided, normal reference with
    continuity correction.
    """
    t = np.asarray(t, float)
    x = np.asarray(x, float)
    ok = np.isfinite(x)
    t, x = t[ok], x[ok]
    n = len(x)
    s, var = _mk_s_var(x)
    slope, intercept = sen_slope(t, x)
    corr, sig = 1.0, []
    if modified and n >= 4:
        ranks = stats.rankdata(x - slope * t)
        rho = _rank_acf(ranks, n - 1)
        acc = 0.0
        for k, r in enumerate(rho, start=1):
            lo = (-1 - 1.96 * np.sqrt(n - k - 1)) / (n - k)
            hi = (-1 + 1.96 * np.sqrt(n - k - 1)) / (n - k)
            if r < lo or r > hi:
                sig.append(k)
                acc += (n - k) * (n - k - 1) * (n - k - 2) * r
        corr = 1 + 2.0 / (n * (n - 1) * (n - 2)) * acc
        if not corr > 0:  # degenerate correction: fall back to the uncorrected test
            corr, sig = 1.0, []
    z = 0.0 if s == 0 else (s - np.sign(s)) / np.sqrt(var * corr)
    p = 2 * (1 - stats.norm.cdf(abs(z)))
    return {"n": n, "S": s, "tau": s / (n * (n - 1) / 2), "z": float(z), "p": float(p), "slope": slope,
            "intercept": intercept, "correction": float(corr), "sig_lags": sig}


def trends(outcome, d, names):
    c = LOG_CONST[outcome]
    rows = []
    for u, g in d.groupby("unit_id"):
        g = g[g.cases.notna()]
        inc = g.cases / g.population * SCALE[outcome]
        y = np.log(inc.values + c)
        r = mann_kendall(g.year.values, y)
        rows.append({"outcome": outcome, "unit_id": u, "unit_name": names[u], "n_years": r["n"], "S": r["S"],
                     "tau": r["tau"], "z_mod": r["z"], "p_mod": r["p"],
                     "p_unmodified": mann_kendall(g.year.values, y, modified=False)["p"],
                     "variance_factor": r["correction"], "sig_rank_lags": ";".join(map(str, r["sig_lags"])),
                     "sen_slope_log": r["slope"], "sen_pct_per_year": 100 * (np.exp(r["slope"]) - 1),
                     "mean_incidence": inc.mean()})
    out = pd.DataFrame(rows)
    out["direction"] = np.where(out.p_mod < 0.05, np.where(out.sen_slope_log < 0, "decreasing", "increasing"), "none")
    return out


# ----------------------------------------------------------------------------------- figures
def save_fig(fig, name):
    from PIL import Image
    tif, png = FIG / f"{name}.tif", FIG / f"{name}.png"
    fig.savefig(png, dpi=200, bbox_inches="tight", facecolor="white")
    fig.savefig(tif, dpi=DPI, bbox_inches="tight", facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    with Image.open(tif) as im:
        im.convert("RGB").save(tif, compression="tiff_lzw", dpi=(DPI, DPI))
    print("wrote", tif.relative_to(ROOT))


def _style_map(ax, title):
    ax.set_axis_off()
    ax.set_aspect(1.0)
    ax.set_title(title, loc="left", fontsize=8)


def shade(i, k):
    return plt.get_cmap("Greys")(0.12 + 0.8 * i / max(k - 1, 1))


def choropleth_panel(ax, poly, vals, bins, title):
    """bins: increasing edges (len k+1); class i shaded from lighter to darker grey; NaN hatched."""
    k = len(bins) - 1
    cls = np.clip(np.digitize(vals, bins[1:-1]), 0, k - 1)
    ok = np.isfinite(vals)
    for i in range(k):
        sub = poly[ok & (cls == i)]
        if len(sub):
            sub.plot(ax=ax, facecolor=shade(i, k), edgecolor="black", linewidth=0.3)
    if (~ok).any():
        poly[~ok].plot(ax=ax, facecolor="white", edgecolor="black", hatch="////", linewidth=0.3)
    _style_map(ax, title)


def class_legend(fig, bins, fmt, title, anchor, ncol=5, missing=False):
    k = len(bins) - 1
    hs = [Patch(facecolor=shade(i, k), edgecolor="black", linewidth=0.4,
                label=f"{fmt.format(bins[i])} to {fmt.format(bins[i + 1])}") for i in range(k)]
    if missing:
        hs.append(Patch(facecolor="white", edgecolor="black", hatch="////", linewidth=0.4, label="No data"))
    fig.legend(handles=hs, loc="upper center", bbox_to_anchor=anchor, ncol=ncol, frameon=False, title=title,
               handlelength=1.2, columnspacing=1.0)


def quantile_bins(arrs, k, decimals):
    v = np.concatenate([a[np.isfinite(a)] for a in arrs])
    b = np.round(np.quantile(v, np.linspace(0, 1, k + 1)), decimals)
    f = 10 ** decimals
    b[0], b[-1] = np.floor(v.min() * f) / f, np.ceil(v.max() * f) / f
    return np.unique(b)


def fig_incidence(poly, tables, periods, name, label, outcome, sir=None):
    vals = [tables[p]["inc"].reindex(poly["unit_id"]).values for p in periods]
    dec = 2 if outcome == "malaria" else 1
    bins = quantile_bins(vals, 5, dec)
    n_p = len(periods) + (1 if sir is not None else 0)
    nrow = int(np.ceil(n_p / 2))
    fig, axes = plt.subplots(nrow, 2, figsize=(17.5 * CM, (7.6 * nrow + 2.4) * CM))
    axes = np.atleast_1d(axes).ravel()
    fmt = "{:.2f}" if outcome == "malaria" else "{:.1f}"
    miss = any(np.isnan(v).any() for v in vals)
    for k, (p, v) in enumerate(zip(periods, vals)):
        choropleth_panel(axes[k], poly, v, bins, f"({'abcd'[k]}) {p[0]}–{p[1]}")
    k = len(periods)
    rows = []
    for p, v in zip(periods, vals):
        rows += [{"panel": f"mean_incidence_{p[0]}_{p[1]}", "unit_id": u, "value": x} for u, x in zip(poly["unit_id"], v)]
    if sir is not None:
        sv = sir["sir"].reindex(poly["unit_id"]).values
        sb = np.array([0, 0.5, 0.8, 1.25, 2.0, float(np.ceil(np.nanmax(sv)))])
        choropleth_panel(axes[k], poly, sv, sb, f"({'abcd'[k]}) Standardised incidence ratio, 2002–2025")
        rows += [{"panel": "sir_whole_period", "unit_id": u, "value": x} for u, x in zip(poly["unit_id"], sv)]
        k += 1
    for a in axes[k:]:
        a.set_axis_off()
    fig.subplots_adjust(hspace=0.15, wspace=0.02, bottom=0.2)
    if sir is not None:
        class_legend(fig, bins, fmt, f"Mean incidence, {label} (panels a–c)", (0.27, 0.15), ncol=2, missing=miss)
        class_legend(fig, sb, "{:g}", "Standardised incidence ratio (panel d)", (0.74, 0.15), ncol=2)
    else:
        class_legend(fig, bins, fmt, f"Mean incidence, {label}", (0.5, 0.17), ncol=3, missing=miss)
    save_fig(fig, name)
    pd.DataFrame(rows).to_csv(FIG / "data" / f"figure_{name.split('_')[1]}.csv", index=False)


LISA_STYLE = {"HH": ("0.05", ""), "LL": ("white", "////"), "HL": ("0.45", ""), "LH": ("0.8", "xxxx"), "NS": ("0.95", "")}
LISA_LAB = {"HH": "High–high", "LL": "Low–low", "HL": "High–low", "LH": "Low–high", "NS": "Not significant"}


def fig_lisa(poly, lisas):
    fig, axes = plt.subplots(1, 2, figsize=(17.5 * CM, 9.2 * CM))
    for ax, (key, title) in zip(axes, [("malaria", "(a) Malaria, 2002–2025"), ("dengue", "(b) Dengue, 2007–2025")]):
        cl = lisas[key].set_index("unit_id")["cluster"]
        for c, (fc, h) in LISA_STYLE.items():
            sub = poly[poly["unit_id"].map(cl).eq(c)]
            if len(sub):
                sub.plot(ax=ax, facecolor=fc, hatch=h, edgecolor="black", linewidth=0.3)
        excl = poly[~poly["unit_id"].isin(cl.index)]
        if len(excl):
            excl.plot(ax=ax, facecolor="white", edgecolor="black", linestyle=":", linewidth=0.5)
        _style_map(ax, title)
    hs = [Patch(facecolor=LISA_STYLE[c][0], hatch=LISA_STYLE[c][1], edgecolor="black", linewidth=0.4, label=LISA_LAB[c])
          for c in ["HH", "LL", "HL", "LH", "NS"]]
    hs.append(Patch(facecolor="white", edgecolor="black", linestyle=":", label="No neighbours"))
    fig.legend(handles=hs, loc="upper center", bbox_to_anchor=(0.5, 0.1), ncol=6, frameon=False, handlelength=1.2,
               columnspacing=0.9)
    fig.subplots_adjust(wspace=0.02, bottom=0.14)
    save_fig(fig, "Figure_S8")
    rows = []
    for key in ["malaria", "dengue"]:
        rows += [{"panel": key, "unit_id": r.unit_id, "log_mean_incidence": r.value, "cluster": r.cluster}
                 for r in lisas[key].itertuples()]
    pd.DataFrame(rows).to_csv(FIG / "data" / "figure_S8.csv", index=False)


def fig_trends(tr):
    fig, axes = plt.subplots(1, 2, figsize=(17.5 * CM, 13.5 * CM))
    for ax, (o, title) in zip(axes, [("malaria", "(a) Malaria, 2002–2025 (2006 missing)"),
                                     ("dengue", "(b) Dengue, 2007–2025")]):
        d = tr[tr.outcome == o].sort_values("sen_pct_per_year").reset_index(drop=True)
        y = np.arange(len(d))
        sig = (d.p_mod < 0.05).values
        ax.axvline(0, color="0.3", lw=0.6)
        ax.scatter(d.sen_pct_per_year[~sig], y[~sig], facecolor="white", edgecolor="black", s=22, zorder=3,
                   label="p ≥ 0.05")
        ax.scatter(d.sen_pct_per_year[sig], y[sig], facecolor="black", edgecolor="black", s=22, zorder=3,
                   label="p < 0.05")
        ax.set_yticks(y)
        ax.set_yticklabels(d.unit_id, fontsize=6.5)
        ax.set_ylim(-0.7, len(d) - 0.3)
        ax.set_xlabel("Sen's slope, % change in incidence per year")
        ax.set_title(title, loc="left")
        ax.grid(axis="y", color="0.9", lw=0.4)
        ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    save_fig(fig, "Figure_S9")
    tr[["outcome", "unit_id", "sen_pct_per_year", "p_mod"]].to_csv(FIG / "data" / "figure_S9.csv", index=False)


def fig_unit_irr(poly, irr, moran):
    v = irr.set_index("unit_id")["log_irr_lag0"].reindex(poly["unit_id"]).values
    lim = float(np.nanmax(np.abs(v)))
    bins = np.linspace(-lim, lim, 6)
    fig, ax = plt.subplots(figsize=(9 * CM, 10.5 * CM))
    choropleth_panel(ax, poly, v, bins, "Unpooled unit-specific ONI IRR, lag 0")
    hs = [Patch(facecolor=shade(i, 5), edgecolor="black", linewidth=0.4,
                label=f"{np.exp(bins[i]):.2f} to {np.exp(bins[i + 1]):.2f}") for i in range(5)]
    hs.append(Patch(facecolor="white", edgecolor="black", hatch="////", linewidth=0.4, label="Not estimated"))
    fig.legend(handles=hs, loc="upper center", bbox_to_anchor=(0.5, 0.08), ncol=3, frameon=False, handlelength=1.2,
               title=f"IRR per 1 °C ONI (Moran's I = {moran[0]:.2f}, p = {moran[1]:.3f})")
    fig.subplots_adjust(bottom=0.2)
    save_fig(fig, "Figure_S10")
    pd.DataFrame({"unit_id": poly["unit_id"], "log_irr_lag0": v, "irr_lag0": np.exp(v)}).to_csv(
        FIG / "data" / "figure_S10.csv", index=False)


# ----------------------------------------------------------------------------------- main
def fmt_p(p):
    return "p < 0.001" if p < 0.001 else f"p = {p:.3f}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--unit-irr", default=None)
    args = ap.parse_args()
    (FIG / "data").mkdir(parents=True, exist_ok=True)
    fr = frame()
    names = dict(zip(fr.unit_id, fr.unit_name))

    W = build_weights()
    W.to_csv(WEIGHTS_CSV, index=False)
    nb = {u: [] for u in fr.unit_id}
    for a, b in W.itertuples(index=False):
        nb[a].append(b)
        nb[b].append(a)
    print("queen neighbours:")
    for u in fr.unit_id:
        print(f"  {u:7s} {len(nb[u]):2d}  {','.join(sorted(nb[u]))}")
    islands = [u for u in nb if not nb[u]]
    print("islands (no neighbours):", islands)

    panels = {o: panel(o) for o in ["malaria", "dengue"]}
    periods = {"malaria": MAL_PERIODS, "dengue": DEN_PERIODS}
    tables = {o: {p: period_incidence(panels[o], *p, SCALE[o]) for p in periods[o]} for o in panels}
    whole = {"malaria": period_incidence(panels["malaria"], 2002, 2025, SCALE["malaria"]),
             "dengue": period_incidence(panels["dengue"], 2007, 2025, SCALE["dengue"])}
    sir = {o: sir_table(panels[o]) for o in panels}
    pd.concat([sir[o].assign(outcome=o) for o in sir]).reset_index().to_csv(REV / "sir_by_unit.csv", index=False)

    keys = {}

    def key(name, value, display, desc):
        keys[f"A.rev.{name}"] = {"value": value, "display": display, "description": desc}

    mrows = []
    for o in panels:
        for p in periods[o]:
            v = {u: np.log(x) for u, x in tables[o][p]["inc"].items() if x > 0}
            I, pp, n, used = moran_i(v, WEIGHTS_CSV, detail=True)
            mrows.append({"outcome": o, "period": f"{p[0]}-{p[1]}", "statistic": "log mean incidence", "I": I,
                          "expected_I": -1 / (n - 1), "p_perm": pp, "n_units": n, "n_permutations": NPERM,
                          "excluded_units": ";".join(sorted(set(fr.unit_id) - set(used)))})
            tag = f"moran.{o}_{p[0]}_{p[1]}"
            key(f"{tag}.I", I, f"I = {I:.2f}", f"Global Moran's I of log mean {o} incidence {p[0]}-{p[1]} (queen, row-standardised, {n} units; islands and units without data excluded)")
            key(f"{tag}.p", pp, fmt_p(pp), f"Two-sided permutation p ({NPERM} permutations, seed {SEED}) for Moran's I, {o} {p[0]}-{p[1]}")
            key(f"{tag}.n_units", n, str(n), f"Units in Moran's I, {o} {p[0]}-{p[1]}")

    lisas, lrows = {}, []
    for o in panels:
        v = {u: np.log(x) for u, x in whole[o]["inc"].items() if x > 0}
        L = lisa(v, WEIGHTS_CSV)
        L["unit_name"] = L.unit_id.map(names)
        L.insert(0, "outcome", o)
        lisas[o] = L
        lrows.append(L)
        for c, lab in [("HH", "high-high"), ("LL", "low-low"), ("HL", "high-low"), ("LH", "low-high")]:
            mem = L[L.cluster == c].unit_name.tolist()
            key(f"lisa.{o}.{c.lower()}_units", mem, ", ".join(mem) if mem else "none",
                f"Units in the {lab} LISA cluster for {o} (whole period, conditional permutation p < 0.05)")
            key(f"lisa.{o}.{c.lower()}_n", len(mem), str(len(mem)), f"Number of {lab} LISA units, {o}")
    pd.concat(lrows).to_csv(REV / "lisa_clusters.csv", index=False)

    irr, irr_moran = None, None
    if args.unit_irr:
        f = Path(args.unit_irr)
        f = f if f.is_absolute() else ROOT / f
        if f.exists():
            irr = pd.read_csv(f)
            v = {r.unit_id: r.log_irr_lag0 for r in irr.itertuples() if np.isfinite(r.log_irr_lag0)}
            I, pp, n, used = moran_i(v, WEIGHTS_CSV, detail=True)
            irr_moran = (I, pp)
            mrows.append({"outcome": "malaria", "period": "2002-2025", "statistic": "unpooled log IRR (lag 0)", "I": I,
                          "expected_I": -1 / (n - 1), "p_perm": pp, "n_units": n, "n_permutations": NPERM,
                          "excluded_units": ";".join(sorted(set(fr.unit_id) - set(used)))})
            key("moran.unpooled_irr_lag0.I", I, f"I = {I:.2f}", f"Moran's I of unpooled unit-specific lag-0 log ONI IRRs ({n} units)")
            key("moran.unpooled_irr_lag0.p", pp, fmt_p(pp), f"Two-sided permutation p ({NPERM}) for Moran's I of unpooled lag-0 log IRRs")
            print(f"unpooled IRR Moran's I = {I:.3f}, p = {pp:.4f}, n = {n}")
        else:
            print(f"NOTE: {f} not found: unpooled IRR Moran's I and Figure_S10 skipped")
    pd.DataFrame(mrows).to_csv(REV / "moran_results.csv", index=False)

    tr = pd.concat([trends(o, panels[o], names) for o in panels])
    tr.to_csv(REV / "trends_mk.csv", index=False)
    for o in panels:
        t = tr[tr.outcome == o]
        dec, inc = int((t.direction == "decreasing").sum()), int((t.direction == "increasing").sum())
        med = float(t.sen_pct_per_year.median())
        nund = int(((t.p_unmodified < 0.05) & (t.sen_slope_log < 0)).sum())
        key(f"trend.{o}.n_decreasing", dec, str(dec), f"Units (of {len(t)}) with a significant (modified Mann-Kendall p < 0.05) decreasing trend in log {o} incidence")
        key(f"trend.{o}.n_increasing", inc, str(inc), f"Units with a significant increasing trend, {o}")
        key(f"trend.{o}.n_none", len(t) - dec - inc, str(len(t) - dec - inc), f"Units with no significant trend, {o}")
        key(f"trend.{o}.median_sen_pct", med, f"{med:.1f}%", f"Median across units of Sen's slope as % change per year, {o} (100*(exp(slope)-1) of log(incidence + {LOG_CONST[o]}))")
        key(f"trend.{o}.n_units", len(t), str(len(t)), f"Units in the trend analysis, {o}")
        dn = t[t.direction == "decreasing"].unit_name.tolist()
        up = t[t.direction == "increasing"].unit_name.tolist()
        key(f"trend.{o}.decreasing_units", dn, ", ".join(dn) or "none", f"Units with significant decreasing trend, {o}")
        key(f"trend.{o}.increasing_units", up, ", ".join(up) or "none", f"Units with significant increasing trend, {o}")
        key(f"trend.{o}.n_unmodified_decreasing", nund, str(nund), f"Units with unmodified Mann-Kendall p < 0.05 and decreasing slope, {o} (comparison)")
    for o in panels:
        key(f"sir.{o}.max", float(sir[o].sir.max()), f"{sir[o].sir.max():.1f}", f"Largest {o} standardised incidence ratio ({names[sir[o].sir.idxmax()]})")
    key("weights.n_islands", len(islands), str(len(islands)), f"Units without queen neighbours ({', '.join(islands)})")
    key("weights.n_edges", len(W), str(len(W)), "Queen-contiguity edges among the 34 units (1 km buffer)")
    (REV / "spatial_trend_keys.json").write_text(json.dumps(keys, indent=1, ensure_ascii=False), encoding="utf-8")

    poly = unit_polygons(0.02)
    fig_incidence(poly, tables["malaria"], MAL_PERIODS, "Figure_S6", "per 1,000", "malaria", sir=sir["malaria"])
    fig_incidence(poly, tables["dengue"], DEN_PERIODS, "Figure_S7", "per 100,000", "dengue")
    fig_lisa(poly, lisas)
    fig_trends(tr)
    if irr is not None:
        fig_unit_irr(poly, irr, irr_moran)
    print(pd.DataFrame(mrows)[["outcome", "period", "statistic", "I", "p_perm", "n_units"]].to_string())
    for o in panels:
        print(o, "LISA:", {c: lisas[o][lisas[o].cluster == c].unit_id.tolist() for c in ["HH", "LL", "HL", "LH"]})
        t = tr[tr.outcome == o]
        print(o, "decreasing:", t[t.direction == "decreasing"].unit_id.tolist(), "increasing:",
              t[t.direction == "increasing"].unit_id.tolist(), "median %/yr", round(t.sen_pct_per_year.median(), 2))


if __name__ == "__main__":
    main()
