"""Display-ready CSVs for the revised supplement (Supplementary Tables S3, S4a-d, S7, S8, S9).

Reads the revision outputs (models/A/output/rev_*.csv, papers/A/revision/*.csv) and writes
papers/A/revision/supp/*.csv, which papers/A/actatropica/supplement.md renders with <<CSV:...>> markers.
Formatting only: no estimate is recomputed. Run: python -I papers/A/revision/make_supp_tables.py
"""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "papers/A/revision/supp"
RO = ROOT / "models/A/output"
RV = ROOT / "papers/A/revision"
MINUS = "−"


def rd(p):
    with open(p, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def wr(name, header, rows):
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / name, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def num(x, d=2):
    s = f"{float(x):.{d}f}"
    return s.replace("-", MINUS)


def ci(e, lo, hi):
    return f"{num(e)} ({num(lo)}–{num(hi)})"


def pval(p):
    p = float(p)
    return "<0.001" if p < 0.001 else f"{p:.3f}"


def intc(x):
    s = f"{float(x):,.0f}"
    return s.replace("-", MINUS)


def s3():
    rows = []
    for r in rd(RO / "rev_inference.csv"):
        if r["variant"] == "glm_model":
            continue
        lab = {"primary": "Primary model, Wald", "year_iid": "+ independent year random effect",
               "year_ar1": "+ AR(1) year random effect (random intercepts only)",
               "cluster_year": "NB GLM, unit fixed effects, clustered by year (HC1, t(22))",
               "cluster_twoway": "NB GLM, clustered by year and unit (not part of the rule)",
               "boot": "Primary model, moving-block year bootstrap (3-year blocks, 999 resamples)"}[r["variant"]]
        rows.append(["Inference", lab, "782", r["display_lag0"].replace("95% CI ", "").replace(" (random intercepts only)", ""),
                     r["display_lag1"].replace("95% CI ", "").replace(" (random intercepts only)", "")])
    for r in rd(RO / "rev_leave_one_episode_out.csv"):
        rows.append(["Episode removed", f"{r['label']} (years {r['years_dropped'].replace(';', ', ')})", r["n"],
                     r["display_lag0"].replace("95% CI ", ""), r["display_lag1"].replace("95% CI ", "")])
    labels = {"trend_ns3": "Trend ns(year, 3)", "trend_ns6": "Trend ns(year, 6)", "trend_ns8": "Trend ns(year, 8)",
              "trend_unit_linear": "Trend ns(year, 4) + unit-specific linear trends",
              "aber_none": "No testing-effort term", "aber_lag1": "log(ABER) of the previous year",
              "aber_ns3": "ns(log ABER, 3)", "aber_positivity": "Test positivity (offset log persons tested)"}
    for r in rd(RO / "rev_trend_aber.csv"):
        grp = "Trend" if r["family"] == "trend" else "Testing effort"
        rows.append([grp, labels.get(r["fit"], r["label"]), r["n"], r["display_lag0"].replace("95% CI ", ""),
                     r["display_lag1"].replace("95% CI ", "")])
    wr("S3_temporal_specification.csv", ["Group", "Analysis", "Unit-years", "ONI lag 0, IRR (95% CI)",
                                          "ONI lag 1, IRR (95% CI)"], rows)


def s4():
    rows = [[r["outcome"].capitalize(), r["period"].replace("-", "–"), r["statistic"], num(r["I"]), pval(r["p_perm"]),
             r["n_units"]] for r in rd(RV / "moran_results.csv")]
    wr("S4a_moran.csv", ["Outcome", "Period", "Statistic", "Moran's I", "Permutation p", "Units"], rows)
    rows = []
    for r in rd(RV / "trends_mk.csv"):
        rows.append([r["outcome"].capitalize(), r["unit_name"], r["n_years"], num(r["sen_pct_per_year"], 1),
                     pval(r["p_mod"]), r["direction"]])
    wr("S4b_trends.csv", ["Outcome", "Unit", "Years", "Sen's slope, % per year", "Modified Mann–Kendall p",
                          "Trend (p < 0.05)"], rows)
    nm = {r["unit_id"]: r["unit_name"] for r in rd(RV / "trends_mk.csv")}
    rows = [[r["zone"], ", ".join(nm.get(u, u) for u in r["units"].split(";")), ci(r["irr_lag0"], r["lo_lag0"], r["hi_lag0"]),
             ci(r["irr_lag1"], r["lo_lag1"], r["hi_lag1"])] for r in rd(RO / "rev_zone_interaction.csv")]
    wr("S4c_zones.csv", ["Zone", "Units", "ONI lag 0, IRR (95% CI)", "ONI lag 1, IRR (95% CI)"], rows)
    names = {r["unit_id"]: r["unit_name"] for r in rd(RV / "trends_mk.csv")}
    rows = []
    for r in rd(RO / "rev_unit_irr_unpooled.csv"):
        if r["status"] == "ok":
            rows.append([names.get(r["unit_id"], r["unit_id"]), r["n"], ci(r["irr_lag0"], r["lo"], r["hi"]),
                         ci(r["irr_lag1"], r["lo_lag1"], r["hi_lag1"])])
        else:
            rows.append([names.get(r["unit_id"], r["unit_id"]), r["n"], "not estimated", "not estimated"])
    wr("S4d_unpooled.csv", ["Unit", "Unit-years", "ONI lag 0, IRR (95% CI)", "ONI lag 1, IRR (95% CI)"], rows)


def s7():
    rows = []
    for r in rd(RO / "rev_nonlinear_climate.csv"):
        rows.append([{"linear climate": "Linear", "spline climate (ns, 3 df)": "Natural splines, 3 df"}.get(r["model"], r["model"]), r["n"], ci(r["irr_lag0"], r["lo_lag0"], r["hi_lag0"]),
                     ci(r["irr_lag1"], r["lo_lag1"], r["hi_lag1"]), pval(r["lrt_p"]) if "spline" in r["model"] else "–"])
    wr("S7a_nonlinear.csv", ["Local-climate terms", "Unit-years", "ONI lag 0, IRR (95% CI)", "ONI lag 1, IRR (95% CI)",
                             "Likelihood-ratio p vs linear"], rows)
    rows = [[r["tercile"].capitalize(), r["lag"], ci(r["irr"], r["lo"], r["hi"]), pval(r["wald_p_equal"])]
            for r in rd(RO / "rev_tmax_tercile.csv")]
    wr("S7b_tmax_tercile.csv", ["Climatological Tmax tercile", "Lag", "Tmax anomaly, IRR per 1 °C (95% CI)",
                                "Wald p, equal slopes"], rows)


def s8():
    lab = {"mal_loeo": "Malaria, LOEO", "mal_roll": "Malaria, rolling origin",
           "den_loeo": "Dengue, LOEO", "den_roll": "Dengue, rolling origin"}
    rows = [[lab[r["job"]], r["model"], r["n"], num(r["male"], 3), num(r["rmse_log1p"], 3), num(r["log_score"], 3),
             f"{100 * float(r['coverage95']):.1f}%", f"{100 * float(r['coverage95_yearvar']):.1f}%"]
            for r in rd(RO / "rev_early_warning_summary.csv")]
    wr("S8_early_warning.csv", ["Outcome, validation", "Model", "Scored unit-years", "MALE", "RMSE, log(1 + y)",
                                "Log score", "Coverage, nominal 95%", "Coverage with year variance"], rows)


def s9():
    rows = []
    for r in rd(RO / "rev_scenarios.csv"):
        oni = float(r["oni"])
        rows.append(["Primary" if r["model"] == "primary" else "Primary + year random effect",
                     f"+{oni:.1f}" if abs(oni - 1.9525) > 1e-3 else "+1.95 (observed)",
                     r["year"], f"{intc(r['excess_mean'])} ({intc(r['excess_lo'])} to {intc(r['excess_hi'])})"])
    wr("S9_scenarios.csv", ["Model", "2026 JJAS ONI, °C", "Year", "Difference in expected cases vs ONI = 0, "
                            "mean (nominal 95% interval)"], rows)


if __name__ == "__main__":
    s3(); s4(); s7(); s8(); s9()
    print("wrote", sorted(p.name for p in OUT.glob("*.csv")))
