# Paper A, Acta Tropica ACTROP-D-26-02276: revision analysis plan (major revision)

Written 17 November 2026, after reading the two reviews and before any revision analysis was run. It fixes the
additional analyses, their specifications and how each result will be reported. Every analysis below is reported
whatever its result. The primary model, its estimates and the 191 existing keys of `papers/A/results.json` do not
change; new numbers enter `results.json` as new keys under a dated contract amendment in `papers/A/results_keys.md`.
Decision rule for the revised conclusions: the ONI association is described as robust only if its 95% interval
excludes 1 under the year random effect, the AR(1) year effect, the cluster-robust (year) and the year-block bootstrap
inference; otherwise the text states that it is not robust to shared temporal dependence.

## R2.1 Shared temporal dependence (Reviewer 2, major 1)

1. Primary model + iid year random effect `(1 | year)`.
2. Primary model + AR(1) year random effect shared by all units (glmmTMB `ar1(yearf + 0 | one)`).
3. Cluster-robust inference: NB GLM (MASS::glm.nb) with unit fixed effects and the primary fixed terms; `sandwich::vcovCL`
   clustered by year (HC1) with t(G − 1) reference (G = number of years), and two-way clustered by year and unit.
4. Year-block bootstrap of the primary glmmTMB model: moving blocks of 3 consecutive years resampled with replacement
   to the original length, 999 resamples, seed 20261117, percentile 95% intervals; failures counted and reported.
5. Residual diagnostics of the primary model: randomized quantile (Dunn–Smyth) residuals (uniform QQ, dispersion), the
   year-mean residual series with its autocorrelation function and Ljung–Box test (lags 1–3), the intraclass
   correlation of residuals by year. Supplementary figure.
6. How much of ONI's variance the trend can absorb: R² of JJAS ONI regressed on ns(year, k) for k = 3, 4, 6, 8.
7. Leave-one-ENSO-episode-out estimation: episodes from the NOAA rule on the ONI series (|ONI| ≥ 0.5 °C for ≥ 5
   consecutive overlapping seasons, El Niño and La Niña), restricted to moderate or stronger episodes (peak |ONI| ≥
   1.0 °C). For each, every calendar year the episode touches is dropped, plus the following year (which carries it
   through the lag-1 term); the primary model is refitted. Forest plot (supplementary).
8. Figure 1 already overlays ONI on national malaria and dengue incidence; a supplementary panel adds national
   detrended incidence (observed / primary-model expectation without ONI) against ONI.

## R2.2 Geography (Reviewer 2, major 2; Reviewer 1, major 2)

1. Maps of mean malaria incidence (per 1,000) and dengue incidence (per 100,000) for 2002–2009, 2010–2017 and
   2018–2025 (dengue from 2007), and of the standardised incidence ratio (observed / expected from national rates by
   year) over the whole period. GADM 4.1 state polygons pooled to the 34 units.
2. Global Moran's I (queen contiguity of the pooled units, row-standardised, islands without neighbours dropped from
   the statistic and reported) on log mean incidence per period, permutation p (9,999); local Moran (LISA) clusters
   (high–high, low–low, high–low, low–high at permutation p < 0.05) mapped for the whole period.
3. Spatial variation in the ENSO effect, kept separate from clustering of disease: (a) ONI × zone interaction in the
   primary model, zones = the five regional zonal councils plus the North Eastern Council (Northern, Central, Eastern,
   Western, Southern with the island territories, North-Eastern), Wald test of the interaction; (b) unpooled unit-specific
   ONI IRRs (separate NB GLMs per unit with the primary fixed terms; units with too few years or non-convergence listed)
   mapped, with Moran's I of the unpooled lag-0 log IRRs.
4. Power for slope heterogeneity is reported qualitatively from the profile/boundary estimate and the number of
   episodes; no new test is added.

## R2.3 Trends and diagnostic effort (Reviewer 2, major 3)

1. State-level trends: modified Mann–Kendall test allowing for serial correlation (Hamed and Rao, 1998) and Sen's slope
   of log annual incidence per unit (malaria 2002–2025, dengue 2007–2025, 2006 malaria missing); table and small multiples.
2. Alternative trend specifications for the primary model: ns(year, 3), ns(year, 6), ns(year, 8), and unit-specific linear
   trends (random slope on centred year, uncorrelated).
3. Alternative adjustment for diagnostic effort: no ABER; log ABER of the previous year; ns(log ABER, 3); test positivity
   (offset log persons tested, no ABER term). Plus the existing testing-definition sensitivities.

## R2.4 Climate nonlinearity and ecology (Reviewer 2, major 4)

1. Primary model + natural splines (3 df) of March–June Tmax anomaly and JJAS rainfall departure at lags 0 and 1;
   likelihood-ratio test against the linear-climate model; ONI IRRs with and without the splines; partial-response
   curves (supplementary).
2. If unit climatological March–June Tmax (1981–2010) is available from the existing derived climate data, the Tmax
   anomaly slope by climatology tercile (cool, middle, hot units); otherwise not done and stated.
3. Mediation is relabelled exploratory; no new mediation analysis.

## R2.5 Early-warning validation (Reviewer 2, major 5)

Information available by 1 June of year t only: JJAS ONI of t − 1, March–May ONI of year t (published early June),
local climate at lag 1 only, log ABER of t − 1, previous-year cases, indicators known in advance. Nested models
A_ew (no climate), B_ew (+ lag-1 local climate), C_ew (+ ONI terms). Same LOEO folds and rolling origins, same
metrics and identical-row scoring as the original validation. Coverage additionally reported for intervals that add
the year-level variance from the iid year-effect model.

## R2.6 Outlook as conditional scenarios (Reviewer 2, major 6; Reviewer 1, major 1)

Expected national excess for 2026 (lag 0) and 2027 (lag 1) under JJAS ONI = +0.5, +1.0, +1.5 and the observed +1.95 °C,
with the primary model and with the iid year-effect model; relative to ONI = 0; no "attributable" wording. Title and
Abstract revised so the outlook is a conditional scenario, not a headline.

## Text-only items

Lag choice and malaria seasonality; annual time step; association (not attribution) language throughout; ecological
and aggregation limits (sub-state heterogeneity, e.g. Odisha, Madhya Pradesh, Maharashtra); boundary random-slope
variances explained; vector ecology (An. culicifacies, An. stephensi, An. fluviatilis; Mordecai et al., 2013;
Martens, 1998); Pirani et al. (2024) cited; conclusions reflect the small predictive contribution and outlook
uncertainty; reference formatting checked.
