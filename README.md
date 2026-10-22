# El Niño and malaria across Indian states, 2002–2025

Data and code for the study of the El Niño–Southern Oscillation (ENSO) and annual malaria (and dengue) counts in 34 Indian state-level units, 2002–2025, with an outlook for 2026–27.

Repository: https://github.com/hssling/elnino-malaria-india

**Paper (manuscript under review):**
Siddalingaiah H S. El Niño and malaria in Indian states, 2002–2025, with a 2026–27 outlook. Manuscript under review.

Please cite the paper once it is published, and this repository through `CITATION.cff`.

## Contents

```
.
├── data/
│   ├── raw_files_manifest.csv         public URL, retrieval time, SHA-256 and size of every raw file used
│   └── derived/
│       ├── outcomes/                  state x year outcome panels and their provenance files
│       ├── model/                     model-ready inputs read by the R analysis
│       ├── climate_unit_year.csv      climate exposures pooled to the 34 units, 1981–2025
│       └── enso_{monthly,annual,episodes}.csv
├── models/                            R analysis (renv-locked)
│   ├── common/                        shared model functions (glmmTMB negative binomial) and testthat suite
│   ├── A/                             the analysis of this study (run_all.R, steps 01–06)
│   │   └── output/                    supplementary tables written by the analysis (CSV)
│   ├── renv.lock, renv/, .Rprofile, DESCRIPTION, setup.R
├── papers/A/
│   ├── results.json                   every number reported in the paper (value, display string, script, description)
│   └── figures/                       figure data (CSV) and the plotting script
├── pipeline/                          Python pipeline: raw sources -> derived tables
│   ├── outcomes/                      parsers, source precedence, gates for the outcome panels
│   └── reference/                     unit frame, source list, reviewed exceptions, QC reviews
├── tests/                             pytest suite for the pipeline (and, with -m r, the R suite)
├── tools/build_subdivision_map.py     builds pipeline/reference/district_to_imd_subdivision.csv
├── run_all.py                         pipeline orchestrator
├── pyproject.toml, requirements.txt
├── LICENSE (code, MIT), DATA_LICENSE.md (derived data, CC BY 4.0), CITATION.cff, .zenodo.json
```

No manuscript text is included. Raw source files are not redistributed (see Provenance and Reproduction).

## Data dictionary

All tables are UTF-8 CSV. Missing is an empty cell (NaN), never 0; a printed 0 is 0.

### The 34-unit frame

`pipeline/reference/state_frame.csv` defines 34 analysis units (`unit_id`, `unit_name`, member GADM state ids, name aliases). Three units pool states, because the outcome sources report them together in some years and separately in others:

| unit_id | members | reason |
|---|---|---|
| AP_TG | Andhra Pradesh + Telangana | Telangana was formed in 2014; sources before and after report one or two states |
| JK | Jammu and Kashmir + Ladakh | Ladakh became a UT in 2019; GADM and the pre-2019 sources have one unit |
| DNH_DD | Dadra and Nagar Haveli + Daman and Diu | merged in 2020; earlier sources report two rows |

Pooling rule (applied within one source document): member counts are summed with NaN skipped, so a pooled unit-year is NaN only if every member is NaN; a pooled value with a NaN member is recorded as `partial_pool`. Totals and sub-totals ("All India", "Grand Total") are recognised as `TOTAL` and never summed into units. Climate exposures for the same units (`climate_unit_year.csv`) are pooled as an area-weighted mean over non-NaN members.

### Outcome panels (`data/derived/outcomes/`)

Each panel has one row per unit-year (NaN rows included): `unit_id`, `unit_name`, `year`, the variables, `provisional_any` and `sources` (winning source per variable, `variable:source_id;...`; `masked` or `non_reporting` where a reviewed disposition blanked the value). Panels hold counts; `population` is merged for rate construction and malaria `aber` (annual blood examination rate) is tested / population x 100.

| File | Units x years | Variables | NaN cells |
|---|---|---|---|
| `malaria_state_year.csv` | 34 x 2002–2025 | `tested`, `cases`, `pf` (P. falciparum incl. mixed), `deaths`, `population`, `aber` | 34 per variable: every unit in 2006 (partial-year count, masked; see Missing-data rules) |
| `dengue_state_year.csv` | 34 x 2007–2025 | `cases`, `deaths`, `population` | `cases` 16, `deaths` 9 (14 and 7 reviewed non-reporting; West Bengal 2023 and 2024 partial-year counts, masked) |
| `population_unit_year.csv` | 34 x 2001–2027 | `population`, `source` | none |

Other files:

- `outcome_records_all.csv.gz`: every parsed malaria and dengue record from every source in long format (`outcome, unit_id, raw_label, year, variable, value, source_id, provisional`), including values later overridden, printed zeros and totals. The two panels are rebuilt exactly from this file plus `pipeline/reference/outcomes_reviewed.csv`.
- `cross_source_report.csv`: every pair of sources that disagree on a unit-year-variable (malaria 142, dengue 270), with `provisional_pair`, `lower_bound_pair`, `partial_pool`, `reviewed` and `gating`.
- `missing_sources.csv`: sources that could not be obtained (`mal_sit_0509`, `mal_sit_1317`, `den_nvbdcp14`) and why.

Source precedence. `pipeline/reference/outcomes_sources.csv` lists every source with its public URL, years, `rank` and notes. For each unit-year-variable the value comes from the highest-ranked source that has one, except that a final value beats a provisional one at any rank. Rank 0 means cross-check only.

| Outcome | Sources, low to high precedence |
|---|---|
| malaria | `mal_sit_0206` (2002–06), `mal_sit_0509` (missing), `mal_upto2011` (2007–10), `mal_sit_1317` (missing), `mal_sit_1519` (2015–19), `mal_ar_2010` to `mal_ar_2025` (NCVBDC annual reports; 2017 and 2018 give deaths only), `mal_sit_2022` (2022–25) |
| dengue | `den_nvbdcp12`, `den_nvbdcp` (2007–13), `den_ar1415` (2010–14), `den_nvbdcp14` (missing), `den_nhp2018`, `den_nhp2019`, `den_nhp2021`, `den_nhp2023` (CBHI National Health Profile, 2013–2020), `den_opendengue` (2015–24), `den_ncvbdc` (2021–25) |
| population | `census_a02` (Census 2001 and 2011), `interpolated_census` (2002–10), `rgi_proj2020` (RGI/NCP projections, 2012–27) |

Malaria `tested` is "B.S.E." to 2019, "Total (BSE+RDT)" in 2020–21 and "Tested" (slides + RDT + others) from 2022; it is treated as one series. OpenDengue (2015–24) copies the NVBDCP/MoHFW state tables, so it is not an independent check of them.

### Climate and ENSO tables (`data/derived/`)

`enso_monthly.csv`: one row per ONI season (keyed to the central month): `season`, `year`, `month`, `oni_total`, `oni_anom` (ONI, °C), `phase` (`el_nino`, `la_nina`, `neutral`; NOAA rule: ONI rounded to 0.1 °C ≥ +0.5 or ≤ −0.5 for at least 5 consecutive overlapping seasons), `episode_id`, `provisional`, `ongoing`, `nino34`, `nino34_anom` (ERSSTv5, 1991–2020 base), `dmi` (Dipole Mode Index, HadISST).

`enso_episodes.csv`: one row per El Niño / La Niña episode: `episode_id`, `phase`, `start_season`, `start_year`, `end_season`, `end_year`, `n_seasons`, `peak_anom`, `intensity` (weak < 1.0, moderate 1.0–1.4, strong 1.5–1.9, very_strong ≥ 2.0), `provisional`, `ongoing`.

`enso_annual.csv`: one row per calendar year.

| Column | Definition |
|---|---|
| `oni_jjas` | mean ONI of the four seasons centred on June–September (MJJ, JJA, JAS, ASO); NaN unless all four exist |
| `oni_djf_next` | ONI of the DJF season centred on January of year + 1 |
| `el_nino_developing` / `la_nina_developing` | the JAS season is in an El Niño / La Niña episode |
| `el_nino_decaying` | the DJF season centred on January of this year is in an El Niño episode that started in an earlier year |
| `roni_jjas`, `roni_djf_next` | as `oni_jjas` / `oni_djf_next` for the Relative ONI |
| `nino34_jjas` | June–September mean of the monthly Niño-3.4 anomaly |
| `dmi_jjas`, `dmi_son` | June–September and September–November means of the DMI |

`climate_unit_year.csv`: 34 units x 1981–2025 (1,530 rows), from India Meteorological Department (IMD) gridded data, area-weighted to states and pooled to units. Columns: `rain_annual_mm`, `rain_jjas_mm` (June–September), `rain_jjas_normal_mm` (1981–2010 mean), `rain_jjas_dep_pct` (100 x (JJAS − normal) / normal), `spi3_sep` (3-month SPI ending September), `spi12_dec`, `onset_doy`, `onset_anom_days` (monsoon onset and its anomaly vs 1981–2010), `tmax_mamj_mean`, `tmax_mamj_anom` (March–June mean Tmax and its anomaly vs 1981–2010, °C), `tmax_ge40_days`, `heatwave_days` (Tmax ≥ 45 °C, or ≥ 40 °C and ≥ 4.5 °C above the day-of-year normal), `tmin_mean`, `n_members`, `n_members_used`. SPI: gamma fit per calendar month on 1981–2010 non-zero accumulations, zero-rain probability q (exact zeros at q/2).

### Model inputs (`data/derived/model/`)

`malaria_model.csv` (34 x 2002–2025) and `dengue_model.csv` (34 x 2007–2025): the panel columns plus the climate columns `rain_jjas_dep_pct`, `spi3_sep`, `tmax_mamj_anom`, `heatwave_days`, `onset_anom_days` and their `_lag1` (previous year, merged on year − 1, never shifted); ENSO exposures `oni_jjas`, `roni_jjas`, `nino34_jjas`, `dmi_jjas` and their `_lag1`, `oni_djf_prev`, `roni_djf_prev`, `dmi_son_prev`, the phase flags; programme indicators `nfme2016` (year ≥ 2016), `web_report2019`, `tested_def2022`, `covid2020_21`, `label2014`, `covid2020`; `log_pop`; malaria `log_aber`; `log1p_cases_lag1`. `population_2001_2027.csv`: `unit_id`, `year`, `population`, `log_pop`. `folds_developing.csv` / `folds_decaying.csv`: `fold_id`, `year`, `label` (cross-validation folds anchored on El Niño developing or decaying years).

### Results

`papers/A/results.json` holds every number reported in the paper as `{"<key>": {"value", "display", "script", "description"}}`. `models/A/output/` holds the supplementary tables: `table_S1_sources.csv` (sources), `dropped_rows.csv`, `compare_per_split_*.csv` (validation per split), `mediation_runs.csv`, `model_fits.csv`, `asym_iod.csv`, `index_comparison.csv`, `outlook_by_unit*.csv`, `unit_slopes_malaria.csv`, and `sessionInfo.txt` / `session_packages.csv`. `papers/A/figures/data/` holds the data behind each figure.

## Missing-data rules

- Missing is NaN, never 0. In NVBDCP/NCVBDC tables "-", "NR" or blank means not reported (NaN). A printed 0 is kept as 0 in `outcome_records_all.csv.gz`; whether it reaches a panel as 0 or NaN follows the reviewed rules below.
- Reviewed exceptions: `pipeline/reference/outcomes_reviewed.csv` (`outcome, unit_id, year, variable, source_id, check, disposition, reason, evidence`). `disposition` is `mask` (NaN, source `masked`), `non_reporting` (NaN), `keep` or `use_source:<id>`. The completeness gate fails on any NaN cell that is not listed there. (The file also holds reviewed rows for other outcomes built by the same pipeline; see Pipeline scope.)
- Plausibility: a unit-year is flagged when its value is below 0.2 or above 5 times the centred rolling 7-year median of the unit's series (medians below 20 not tested). A flagged zero is judged non-reporting; a flagged non-zero value is kept only when a second source corroborates it or the printed cell matches the parse.
- Partial years: a count that the document says stops before the end of the year is masked, never used as an annual total.
  - Malaria 2006 (`mal_sit_0206`, "Report received till 25th November, 2006"): the only 2006 source, so all 34 units x `tested`, `cases`, `pf`, `deaths` are masked. The full-year table (`mal_sit_0509`) was not obtainable (Internet Archive rate limit).
  - Dengue West Bengal 2023 and 2024 ("Till 13.9.2023", "Till 21.2.2024"): `cases` and `deaths` masked.
- Provisional values (`provisional_any` True): malaria 2005 (all units; marked provisional in the 2002–06 table) and West Bengal 2017.
- Climate: Andaman and Nicobar rain is masked in every year (no valid IMD grid cell); Lakshadweep rain is a Kerala-coast proxy and is set to NaN in the models (`RAIN_NA_UNITS` in `models/A/common_a.R`); Tamil Nadu and Puducherry 2011–12 rain is masked (IMD grid artefact); Arunachal Pradesh 1999 is masked. Reviewed QC masks are in `pipeline/reference/qc_reviewed.csv`.
- The R analysis drops rows with a missing model variable and lists every dropped unit-year in `models/A/output/dropped_rows.csv`.

## Provenance and licences of the sources

Every raw file is listed in `data/raw_files_manifest.csv` (public URL, retrieval time in UTC, SHA-256, bytes). The files themselves are not redistributed. Source-level URLs and notes are in `pipeline/reference/outcomes_sources.csv` and `models/A/output/table_S1_sources.csv`.

| Source | Used for | Publisher, URL | Terms and attribution |
|---|---|---|---|
| Malaria situation tables and annual reports, 2002–2025 | malaria panel | National Center for Vector Borne Diseases Control (NCVBDC, formerly NVBDCP), Ministry of Health and Family Welfare, https://ncvbdc.mohfw.gov.in (2018–2025 reports: `pipeline/reference/ncvbdc_source_urls.csv`; the analysed files are byte-identical, by SHA-256, to the official PDFs at these URLs, accessed 16 Oct 2026) | Public Government of India documents; aggregated counts reproduced with attribution to NCVBDC |
| Dengue cases and deaths | dengue panel | NCVBDC (as above); Central Bureau of Health Intelligence, National Health Profile, https://cbhidghs.mohfw.gov.in | Public Government of India documents; attribution to NCVBDC and CBHI |
| OpenDengue | dengue panel (2015–24) | OpenDengue, https://opendengue.org; Clarke J et al. 2024, Sci Data 11:296, https://doi.org/10.1038/s41597-024-03120-7 | CC BY 4.0; attribution to the OpenDengue project |
| Census of India 2001/2011 (Table A-02); RGI/NCP population projections 2011–2036 | population | Office of the Registrar General and Census Commissioner, India (https://censusindia.gov.in); National Commission on Population projections (https://nhm.gov.in/New_Updates_2018/Report_Population_Projection_2019.pdf) | Public Government of India data; attribution to RGI and the National Commission on Population |
| Oceanic Niño Index and Relative ONI (ERSSTv6); monthly Niño-3.4 (ERSSTv5) | ENSO exposure | NOAA Climate Prediction Center, https://www.cpc.ncep.noaa.gov/data/indices/ | US Government (NOAA) public data; attribution to NOAA CPC |
| Dipole Mode Index (HadISST) | IOD exposure | NOAA Physical Sciences Laboratory, https://psl.noaa.gov/gcos_wgsp/Timeseries/Data/dmi.had.long.data | NOAA public data; attribution to NOAA PSL |
| IMD 0.25° daily gridded rainfall; 1° daily gridded Tmax/Tmin | local climate | India Meteorological Department, Pune, https://imdpune.gov.in/cmpg/Griddata/rainfall.php and .../maxtemp.php | IMD gridded-data terms: raw grids are not redistributed here; only area-weighted aggregates are. Cite Pai DS et al. 2014, MAUSAM 65(1):1–18 (rainfall) and Srivastava AK et al. 2009, Atmos Sci Lett 10:249–254 (temperature) |
| IMD sub-divisional rainfall 1901–2017; IMD district list by subdivision | climate validation only | data.gov.in (https://data.gov.in/resource/sub-divisional-monthly-rainfall-1901-2017); IMD (https://mausam.imd.gov.in) | Government Open Data License – India; IMD terms |
| GADM 4.1 India boundaries | area weights, maps | https://gadm.org | Obtained from GADM; not redistributed (GADM licence: free for academic and non-commercial use) |
| NASA POWER (MERRA-2) | climate cross-check only | NASA Langley Research Center, https://power.larc.nasa.gov | Obtained from NASA POWER; not redistributed |

The pipeline also reads a few user-supplied files from `data/external/` (listed with their origins in `pipeline/copy_external.py`): the OpenDengue admin-1 annual aggregation, GADM, NASA POWER points, ONI/DMI snapshots used only as a download fallback, and two national series used only as cross-checks (an NCVBDC-based national malaria series and the WHO SEARO national dengue series). They are obtained from their publishers and not redistributed.

## Reproduction

### (a) The R analysis from the derived tables

Requirements: R 4.5.3 (checked by `setup.R`), internet access for the first `renv` restore, about 16 GB RAM.

```
cd models
Rscript setup.R          # restores the package library from renv.lock (renv)
Rscript A/run_all.R      # full analysis
```

The full run takes about 10 hours (1,000 validation and projection draws, B = 1,000 mediation bootstrap resamples). It uses 6 parallel (PSOCK) workers by default; on a 16 GB machine keep `--cores` at 6 or fewer (`Rscript A/run_all.R --cores=4`). Results do not depend on the number of workers (fixed seeds). `--reuse` reuses finished steps from `models/A/_cache/`; `--step=fits|validation|mediation|outlook|results` runs one step.

The analysis reads only `data/derived/` (model inputs, ENSO tables, outcome provenance files), `pipeline/reference/` (unit frame, source list, NCVBDC URLs) and `data/raw_files_manifest.csv` (retrieval dates for Table S1). It writes `papers/A/results.json`, `models/A/output/*.csv` and `papers/A/figures/data/*.csv`; figures are then drawn with `python papers/A/figures/make_figures.py` (Supplementary Figure S1 also needs the GADM file in `data/external/`).

Shared R functions have a testthat suite: `cd models && Rscript common/tests/run.R` (or `python -m pytest -m r` from the repository root).

### (b) Rebuilding the derived tables from the raw sources

A raw → panel rebuild requires downloading the raw files from the official URLs listed in `data/raw_files_manifest.csv` and `pipeline/reference/outcomes_sources.csv`; most downloads are automatic. Compare the SHA-256 of each downloaded file with `data/raw_files_manifest.csv`: documents re-published by the agencies under the same URL may differ.

```
python -m pip install -r requirements.txt   # Python 3.10 or later (tested with 3.14); pdftotext (Poppler) must be on PATH
python run_all.py                           # ENSO download, IMD download and aggregation (1981-2025),
                                            # climate tables, outcome parsing, model inputs, all gates
```

Steps can also be run one at a time: `python -m pipeline.fetch_enso`, `python -m pipeline.outcomes.assemble` (parses the cached raw outcome files and downloads missing ones), `python -m pipeline.climate_units`, `python -m pipeline.model_inputs`. `python run_all.py --skip-download` rebuilds from files already in `data/raw/`. Notes:

- The IMD gridded download is large (45 years of daily grids); the build is resumable.
- Some sources are Internet Archive captures, which may refuse automated requests (HTTP 429); the parsers then record the source as missing.
- The two NCVBDC files with opaque URLs (2025 annual report, situation table from 2022) may be replaced when NCVBDC updates them.
- The files of `data/external/` (see Provenance) must be supplied by the user (`ELNINO_EXTERNAL_SOURCE_DIR`, see `pipeline/copy_external.py`).
- Without any downloads, the malaria and dengue panels can be checked against their records and reviewed exceptions, and the model inputs rebuilt from the distributed tables:

```
python -c "from pipeline.outcomes import assemble, gates; spec = {k: assemble.GATE_SPEC[k] for k in ('malaria', 'dengue')}; [print(r.passed, r.name, r.detail[:100]) for r in gates.run_outcome_gates(spec=spec)]"
python -m pipeline.model_inputs
```

(`run_outcome_gates` rewrites `cross_source_report.csv` with identical content; `model_inputs` rewrites `data/derived/model/` with identical content.)

Python tests: `python -m pytest` (network tests and the R suite are deselected by default).

### Pipeline scope

The outcome pipeline was written for a wider programme and also contains parsers, reference rows and tests for three outcome panels of a separate study (NCRB heat-stroke deaths and farmer suicides, RBI/ICRISAT foodgrain production). Those panels and their model inputs are not part of this study and are not distributed here; `pipeline/model_inputs.py` builds their model inputs only when a full rebuild has produced them. Code comments cite internal decision identifiers (for example "Ruling P12"); the decisions that affect the data are described in this README.

## Software

- R 4.5.3; packages pinned in `models/renv.lock` (renv; glmmTMB, data.table, jsonlite, mgcv, ggplot2, testthat). The session of the reported run is in `models/A/output/sessionInfo.txt`.
- Python 3.14 (tested; 3.10 or later required) with pandas 2.3, numpy 2.3, scipy 1.16, geopandas 1.1, shapely 2.1, pyogrio 0.12, requests 2.32, matplotlib 3.10, pdfplumber 0.11, Pillow, pytest; pdftotext (Poppler) 4.00.

## Licences

- Code: MIT (`LICENSE`).
- Derived data (`data/derived/`, `data/raw_files_manifest.csv`, `pipeline/reference/`, `papers/A/results.json`, `models/A/output/`, `papers/A/figures/data/`): CC BY 4.0 (`DATA_LICENSE.md`), subject to the third-party source terms above.

## Contact

Dr Siddalingaiah H S, Department of Community Medicine, Shridevi Institute of Medical Sciences and Research Hospital, Tumakuru, Karnataka, India. E-mail: hssling@yahoo.com. ORCID: https://orcid.org/0000-0002-4771-8285.
