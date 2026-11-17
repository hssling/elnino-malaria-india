# Paper A revision helpers (models/A/revision_a.R): episode derivation, dropped years, moving-block
# bootstrap resampler, early-warning predictor set, conditional scenarios.

local({
  root <- if (exists("MODELS_ROOT")) MODELS_ROOT else normalizePath(file.path(getwd(), "..", "..", ".."), winslash = "/")
  if (!exists("A_CACHE")) source(file.path(root, "A", "common_a.R"), local = globalenv())
  if (!exists("outlook_newdata")) source(file.path(root, "A", "04_outlook.R"), local = globalenv())
  source(file.path(root, "A", "revision_a.R"), local = globalenv())
})

# synthetic consecutive seasonal ONI starting at DJF of `y0`
syn_oni <- function(values, y0 = 2000L) {
  n <- length(values)
  m <- (seq_len(n) - 1L) %% 12L + 1L
  y <- y0 + (seq_len(n) - 1L) %/% 12L
  data.table::data.table(season = REV_SEASONS[m], year = y, month = m, oni = values)
}

test_that("episodes follow the NOAA rule on a synthetic ONI series", {
  v <- rep(0, 48)
  v[5:10] <- c(0.5, 0.8, 1.1, 1.24, 0.9, 0.6)   # El Nino, 6 seasons (MAM-AUG 2000), peak 1.2
  v[15:18] <- c(0.7, 0.9, 1.0, 0.8)             # 4 seasons only: not an episode
  v[22:26] <- -c(0.45, 0.6, 0.8, 0.7, 0.5)      # La Nina, 5 seasons; -0.45 rounds to -0.5; peak -0.8 (weak)
  v[34:41] <- -c(0.6, 1.0, 1.5, 1.66, 1.4, 1.0, 0.7, 0.5)  # La Nina, 8 seasons, peak -1.7, OCT 2002-MAY 2003
  e <- derive_episodes(syn_oni(v))
  expect_equal(nrow(e), 3L)
  expect_identical(e$phase, c("el_nino", "la_nina", "la_nina"))
  expect_identical(e$n_seasons, c(6L, 5L, 8L))
  expect_equal(e$peak, c(1.2, -0.8, -1.7))
  expect_identical(e$start_season[1], "AMJ")
  expect_identical(e$start_year, c(2000L, 2001L, 2002L))
  # moderate-or-stronger filter drops the weak La Nina
  m <- derive_episodes(syn_oni(v), min_peak = 1.0)
  expect_identical(m$episode_id, c("en2000", "ln2002"))
  # a qualifying run reaching the last season is flagged ongoing
  w <- c(rep(0, 10), rep(0.9, 6))
  expect_true(derive_episodes(syn_oni(w))$ongoing)
})

test_that("dropped years are the years an episode touches plus the following year", {
  # OND 2002 (centre Nov) to MAM 2003: months Oct 2002 - May 2003
  expect_identical(episode_years_touched(2002L, 11L, 2003L, 4L), c(2002L, 2003L))
  # a DJF start touches the previous December; an NDJ end touches the next January
  expect_identical(episode_years_touched(2005L, 1L, 2005L, 12L), c(2004L, 2005L, 2006L))
  expect_identical(episode_drop_years(c(2002L, 2003L)), c(2002L, 2003L, 2004L))
  v <- rep(0, 48); v[34:41] <- -c(0.6, 1.0, 1.5, 1.66, 1.4, 1.0, 0.7, 0.5)
  e <- derive_episodes(syn_oni(v))
  expect_identical(episode_years(e$years_touched), c(2002L, 2003L))
  expect_identical(episode_years(e$years_dropped), c(2002L, 2003L, 2004L))
})

test_that("the moving-block resampler returns contiguous 3-year blocks of the original length", {
  set.seed(7)
  yrs <- 2002:2025
  for (r in 1:50) {
    s <- block_resample_years(yrs, 3L)
    expect_length(s, length(yrs))
    expect_true(all(s %in% yrs))
    blocks <- split(s, ceiling(seq_along(s) / 3))
    expect_true(all(vapply(blocks, function(b) all(diff(b) == 1L), logical(1))))
  }
  # a length that is not a multiple of the block length is cut; the last partial block is still contiguous
  s <- block_resample_years(2002:2011, 3L)
  expect_length(s, 10L)
  expect_true(all(diff(s[1:3]) == 1L) && all(diff(s[4:6]) == 1L) && all(diff(s[7:9]) == 1L))
  expect_error(block_resample_years(c(2002L, 2004L, 2005L), 3L), "consecutive")
  # boot_panel stacks the drawn years' rows (duplicates allowed)
  d <- data.table::data.table(unit_id = rep(c("a", "b"), each = 4), year = rep(2002:2005, 2), y = 1:8)
  bp <- boot_panel(d, c(2003L, 2004L, 2003L))
  expect_equal(nrow(bp), 6L)
  expect_identical(bp$year, rep(c(2003L, 2004L, 2003L), each = 2))
})

test_that("early-warning models use no year-t information except March-May ONI and known indicators", {
  for (o in c("malaria", "dengue")) {
    sp <- ew_specs(ew_covars(o))
    for (m in names(sp)) expect_length(ew_year_t_vars(sp[[m]]), 0L)
    expect_true(all(c("oni_mam", "oni_jjas_lag1") %in% sp$C$fixed))
    expect_false(any(c("oni_jjas", "log_aber", LOCAL_CLIMATE[!grepl("_lag1$", LOCAL_CLIMATE)]) %in%
                       unlist(lapply(sp, spec_vars))))
  }
  # a year-t covariate would be caught
  bad <- ew_specs(ew_covars("malaria"))$C; bad$fixed <- c(bad$fixed, "oni_jjas")
  expect_identical(ew_year_t_vars(bad), "oni_jjas")
  # oni_mam is the MAM season of the same year; log_aber_lag1 is the unit's previous-year log ABER
  oni <- syn_oni(seq(-1, 1, length.out = 36), y0 = 2010L)
  d <- data.table::data.table(unit_id = factor(rep(c("U1", "U2"), each = 3)), year = rep(2010:2012, 2),
                              log_aber = c(1, 2, 3, 4, 5, 6))
  e <- add_ew_terms(d, oni)
  expect_equal(e$oni_mam, rep(oni[season == "MAM"]$oni, 2))
  expect_equal(e$log_aber_lag1, c(NA, 1, 2, NA, 4, 5))
})

test_that("the conditional-scenario excess is zero at ONI 0 (primary and year-effect machinery)", {
  units <- c("BR", "OD", "JH", "MP", "CG", "AS", "MZ", "TR", "GJ", "RJ", "KA", "WB")
  d <- sim_panel(n_units = length(units), seed = 3)
  d[, unit_id := factor(units[as.integer(unit_id)])]
  data.table::setorder(d, unit_id, year)
  d[, `:=`(oni_jjas = x, oni_jjas_lag1 = data.table::shift(x, fill = 0)), by = unit_id]
  ys <- stats::setNames(stats::rnorm(24, 0, 0.2), 2002:2025)
  d[, cases := stats::rnbinom(.N, size = 8, mu = mu * exp(ys[as.character(year)]))]
  spec <- nb_spec("cases", fixed = INDEX_TERMS$oni, re = "int")
  f <- fit_nb(spec, d)
  nd_fun <- function(fit, terms, v) {
    base <- data.table::CJ(unit_id = levels(fit$data$unit_id), year = 2026:2027)
    base[, `:=`(log_pop = 16, lag25 = 0.3)]
    obs <- data.table::copy(base); cf <- data.table::copy(base)
    obs[year == 2026, `:=`(oni_jjas = v, oni_jjas_lag1 = lag25)]; cf[year == 2026, `:=`(oni_jjas = 0, oni_jjas_lag1 = lag25)]
    obs[year == 2027, `:=`(oni_jjas = 0, oni_jjas_lag1 = v)]; cf[year == 2027, `:=`(oni_jjas = 0, oni_jjas_lag1 = 0)]
    list(obs = obs[, !"lag25"], cf = cf[, !"lag25"])
  }
  # primary path: project() on the scenario pair, as rev_scenarios_primary does with outlook_newdata
  for (v in c(0, 1)) {
    nd <- nd_fun(f, INDEX_TERMS$oni, v)
    p <- project(f, nd$obs, nd$cf, INDEX_TERMS$oni, n_draws = 200, seed = 1)
    if (v == 0) expect_equal(max(abs(p$total$excess_mean)), 0) else expect_gt(max(abs(p$total$excess_mean)), 0)
  }
  # year-effect path
  fy <- fit_nb_extra(spec, d, REV_YEAR_IID)
  s <- rev_scenarios_year(fy, c(0, 1), n_draws = 200, seed = 1, newdata_fun = nd_fun)
  expect_equal(max(abs(s[oni == 0, c(excess_mean, excess_lo, excess_hi)])), 0)
  expect_gt(max(abs(s[oni == 1, excess_mean])), 0)
  expect_identical(sort(unique(s$year)), 2026:2027)
})

test_that("Ljung-Box with gaps equals Box.test on a complete series", {
  set.seed(2)
  x <- stats::arima.sim(list(ar = 0.4), 30)
  for (k in 1:3) {
    a <- ljung_box_na(as.numeric(x), k); b <- stats::Box.test(x, lag = k, type = "Ljung-Box")
    expect_equal(a$statistic, unname(b$statistic)); expect_equal(a$p, b$p.value)
  }
})
