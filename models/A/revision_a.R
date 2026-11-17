# Paper A revision (Acta Tropica ACTROP-D-26-02276, major revision): helper functions for the reviewer
# analyses fixed in papers/A/revision/revision_plan.md (sections R2.1 items 1-7, R2.2 item 3, R2.3 items 2-3,
# R2.4 items 1-2, R2.5 and R2.6). Functions only: nothing runs when this file is sourced. The driver is
# models/A/07_revision.R. The shared code in models/common/R and the existing models/A scripts are reused
# unchanged (fit_nb, compare_models, project, outlook_newdata, primary_spec, ...).
# Non-ASCII characters are written as escapes.

REV_SEED <- 20261117L                    # revision seed (plan R2.1 item 4)
REV_SEASONS <- c("DJF", "JFM", "FMA", "MAM", "AMJ", "MJJ", "JJA", "JAS", "ASO", "SON", "OND", "NDJ")

# ------------------------------------------------------------------ ENSO episodes (plan R2.1 item 7)

# NOAA rounding of the ONI to 0.1 degC, half away from zero (as pipeline/enso.py round_half_away).
rev_round_half_away <- function(x) sign(x) * floor(abs(x) * 10 + 0.5 + 1e-9) / 10

# Seasonal ONI from the A pipeline's ENSO table (data/derived/enso_monthly.csv, built from
# data/external/oni_snapshot.ascii.txt and the CPC refresh): season, year, month (centre month), oni.
rev_oni_seasonal <- function(path = file.path(REPO_ROOT, "data", "derived", "enso_monthly.csv")) {
  m <- data.table::fread(path)
  out <- m[, .(season, year = as.integer(year), month = as.integer(month), oni = as.numeric(oni_anom))]
  data.table::setorder(out, year, month)
  out[]
}

# Calendar years touched by an episode running from the season centred on (y0, m0) to the season centred
# on (y1, m1): every month from m0 - 1 to m1 + 1 (a 3-month season spans its centre month +- 1).
episode_years_touched <- function(y0, m0, y1, m1) {
  first <- y0 * 12L + (m0 - 1L) - 1L                    # month index (0-based months) of the first month
  last <- y1 * 12L + (m1 - 1L) + 1L
  sort(unique(as.integer(seq(first, last) %/% 12L)))
}

# Years dropped for an episode: every calendar year it touches plus the following year (carried by lag 1).
episode_drop_years <- function(touched) sort(unique(as.integer(c(touched, max(touched) + 1L))))

# NOAA rule on the seasonal ONI: |ONI| >= threshold (after rounding to 0.1) for >= min_run consecutive
# overlapping seasons, El Nino and La Nina separately; peak = largest |ONI| (rounded) in the run.
# Episodes with rounded |peak| < min_peak are dropped. A qualifying run that reaches the last season is
# flagged ongoing. oni: data.table(season, year, month, oni), one row per consecutive season.
derive_episodes <- function(oni, threshold = 0.5, min_run = 5L, min_peak = 0) {
  oni <- data.table::as.data.table(oni)
  data.table::setorder(oni, year, month)
  idx <- oni$year * 12L + oni$month
  if (any(diff(idx) != 1L)) stop("derive_episodes: seasons must be consecutive", call. = FALSE)
  r <- rev_round_half_away(oni$oni)
  sgn <- ifelse(r >= threshold, 1L, ifelse(r <= -threshold, -1L, 0L))
  sgn[is.na(sgn)] <- 0L
  rl <- rle(sgn)
  ends <- cumsum(rl$lengths); starts <- ends - rl$lengths + 1L
  keep <- rl$values != 0L & rl$lengths >= min_run
  out <- data.table::rbindlist(lapply(which(keep), function(k) {
    i <- starts[k]; j <- ends[k]
    pk <- r[i:j][which.max(abs(r[i:j]))]
    touched <- episode_years_touched(oni$year[i], oni$month[i], oni$year[j], oni$month[j])
    data.table::data.table(
      type = if (rl$values[k] > 0) "El Ni\u00f1o" else "La Ni\u00f1a",
      phase = if (rl$values[k] > 0) "el_nino" else "la_nina",
      start_season = oni$season[i], start_year = oni$year[i], end_season = oni$season[j], end_year = oni$year[j],
      n_seasons = as.integer(j - i + 1L), peak = pk, ongoing = j == nrow(oni),
      years_touched = paste(touched, collapse = ";"),
      years_dropped = paste(episode_drop_years(touched), collapse = ";"))
  }))
  if (!nrow(out)) return(out)
  out <- out[abs(peak) >= min_peak - 1e-9]
  out[, episode_id := sprintf("%s%d", ifelse(phase == "el_nino", "en", "ln"), start_year)]
  out[, label := sprintf("%s %d\u2013%s", type, start_year, substr(as.character(end_year), 3, 4))]
  out[start_year == end_year, label := sprintf("%s %d", type, start_year)]
  data.table::setcolorder(out, c("episode_id", "label", "type"))
  out[]
}

# Years in an episode's years_dropped string, as integers.
episode_years <- function(s) as.integer(strsplit(s, ";", fixed = TRUE)[[1]])

# ------------------------------------------------------------------ moving-block bootstrap (plan R2.1 item 4)

# One moving-block resample of a run of consecutive years: blocks of `block` consecutive years, start
# drawn uniformly from the length(years) - block + 1 possible starts, with replacement, concatenated and
# cut to n_out years. Uses the caller's RNG stream.
block_resample_years <- function(years, block = 3L, n_out = length(years)) {
  years <- sort(as.integer(years))
  if (any(diff(years) != 1L)) stop("block_resample_years: years must be consecutive", call. = FALSE)
  starts <- years[seq_len(length(years) - block + 1L)]
  nb <- ceiling(n_out / block)
  s <- starts[sample.int(length(starts), nb, replace = TRUE)]
  as.integer(unlist(lapply(s, function(a) a + seq_len(block) - 1L)))[seq_len(n_out)]
}

# Panel rows for a resampled year sequence (a year drawn twice contributes its rows twice).
boot_panel <- function(d, yseq) {
  d <- data.table::as.data.table(d)
  rows <- split(seq_len(nrow(d)), d$year)
  d[unlist(rows[as.character(yseq)], use.names = FALSE)]
}

# ------------------------------------------------------------------ fitter hooks

# A glmmTMB fitter for fit_nb() that adds extra formula terms (e.g. a year random effect) and lets `prep`
# add the columns they need. The data passed by fit_nb are already complete on the spec's variables;
# extra variables must be complete too (use fit_nb_extra()).
fitter_plus <- function(extra, prep = NULL) {
  force(extra); force(prep)
  function(formula, data, family = glmmTMB::nbinom2()) {
    if (!is.null(prep)) data <- prep(data)
    f <- stats::as.formula(paste(paste(deparse(formula, width.cutoff = 500L), collapse = " "), "+",
                                 paste(extra, collapse = " + ")), env = globalenv())
    glmmTMB::glmmTMB(f, data = data, family = family)
  }
}

# yearf over the full span of years in the data (a year with no rows keeps its level, so the AR(1)
# process stays equally spaced in calendar time); one = a single-level factor for the shared AR(1) term;
# year_c = centred year per decade.
prep_year <- function(d) {
  d <- data.table::copy(d)
  yr <- seq(min(d$year), max(d$year))
  d[, yearf := factor(year, levels = yr)]
  d[, one := factor(1L)]
  d[, year_c := (year - 2013.5) / 10]
  d
}

# fit_nb with extra formula terms: rows are first restricted to complete extra_vars so fit$n is the
# number of rows the model used.
fit_nb_extra <- function(spec, data, extra, extra_vars = character(), prep = prep_year) {
  d <- data.table::as.data.table(data)
  if (length(extra_vars)) d <- complete_rows(d, extra_vars)$data
  fit_nb(spec, droplevels(d), fitter = fitter_plus(extra, prep))
}

REV_YEAR_IID <- "(1 | yearf)"
REV_YEAR_AR1 <- "ar1(yearf + 0 | one)"
REV_UNIT_TREND <- "(0 + year_c | unit_id)"

# ------------------------------------------------------------------ IRRs and displays

# IRR row for one term; adds the RE-structure suffix used by the sensitivity keys when the structure
# differs from the primary model's.
rev_irr <- function(fit, term, primary_re = "diag") {
  r <- irr_table(fit, term)
  suf <- if (!identical(fit$re_used, primary_re)) RE_SUFFIX[[fit$re_used]] else ""
  list(irr = r$irr, lo = r$lo, hi = r$hi, p = r$p, display = paste0(disp_irr(r), suf), re_used = fit$re_used)
}

# Variance-component SD of a grouping term in a glmmTMB fit (NA if absent).
rev_group_sd <- function(fit, group) {
  vc <- glmmTMB::VarCorr(fit$model)$cond[[group]]
  if (is.null(vc)) return(NA_real_)
  unname(attr(vc, "stddev")[1])
}

rev_ar1_rho <- function(fit, group = "one") {
  vc <- glmmTMB::VarCorr(fit$model)$cond[[group]]
  if (is.null(vc)) return(NA_real_)
  cr <- attr(vc, "correlation")
  if (is.null(cr) || nrow(cr) < 2) return(NA_real_)
  cr[1, 2]
}

# Wald test that a set of fixed effects is zero.
rev_wald <- function(fit, terms) {
  b <- glmmTMB::fixef(fit$model)$cond[terms]
  V <- stats::vcov(fit$model)$cond[terms, terms, drop = FALSE]
  stat <- as.numeric(t(b) %*% solve(V, b))
  list(chisq = stat, df = length(terms), p = stats::pchisq(stat, length(terms), lower.tail = FALSE))
}

# Linear combination sum(beta[terms]) as an IRR with a Wald CI.
rev_lincomb_irr <- function(fit, terms) {
  b <- glmmTMB::fixef(fit$model)$cond[terms]
  V <- stats::vcov(fit$model)$cond[terms, terms, drop = FALSE]
  est <- sum(b); se <- sqrt(sum(V))
  list(irr = exp(est), lo = exp(est - 1.959964 * se), hi = exp(est + 1.959964 * se), log_irr = est, se = se)
}

# "x to y" display of a range of point estimates.
rev_range_disp <- function(v, digits = 2) paste(fmt_num(min(v), digits), "to", fmt_num(max(v), digits))

# ------------------------------------------------------------------ cluster-robust NB GLM (plan R2.1 item 3)

# NB GLM with unit fixed effects and the primary fixed terms (ns(year, 4) trend, offset log population).
glm_nb_unitfe <- function(data, terms = c(INDEX_TERMS$oni, MAL_COVARS), trend_df = 4L) {
  d <- complete_rows(data.table::as.data.table(data), c("cases", terms, "log_pop", "year", "unit_id"))$data
  d <- droplevels(d)
  f <- stats::as.formula(paste("cases ~", paste(terms, collapse = " + "),
                               sprintf("+ splines::ns(year, df = %d) + unit_id + offset(log_pop)", trend_df)))
  m <- MASS::glm.nb(f, data = d, control = stats::glm.control(maxit = 100))
  list(model = m, data = d)
}

# IRR with a t(df) reference from a coefficient and its (robust) variance.
rev_t_irr <- function(b, v, df) {
  q <- if (is.finite(df)) stats::qt(0.975, df) else 1.959964
  se <- sqrt(v)
  list(irr = exp(b), lo = exp(b - q * se), hi = exp(b + q * se), se = se,
       p = 2 * (if (is.finite(df)) stats::pt(-abs(b / se), df) else stats::pnorm(-abs(b / se))))
}

# ------------------------------------------------------------------ residual diagnostics (plan R2.1 item 5)

# Randomized quantile (Dunn-Smyth) residuals of an NB2 fit, conditional on the random effects:
# u ~ U(F(y - 1), F(y)), r = qnorm(u). Seed applied locally.
rqr_nb <- function(y, mu, size, seed = REV_SEED) {
  .with_seed(seed, {
    a <- stats::pnbinom(y - 1, size = size, mu = mu)
    b <- stats::pnbinom(y, size = size, mu = mu)
    u <- a + stats::runif(length(y)) * (b - a)
    list(u = u, r = stats::qnorm(u))
  })
}

# Ljung-Box Q on a series that may contain NA (gaps kept in place): autocorrelations from
# acf(na.action = na.pass), n = number of non-missing values, chi-square on `lag` df. Without NA this
# equals stats::Box.test(type = "Ljung-Box").
ljung_box_na <- function(x, lag) {
  n <- sum(is.finite(x))
  rho <- stats::acf(x, lag.max = lag, plot = FALSE, na.action = stats::na.pass)$acf[-1][seq_len(lag)]
  q <- n * (n + 2) * sum(rho^2 / (n - seq_len(lag)))
  list(statistic = q, df = lag, p = stats::pchisq(q, lag, lower.tail = FALSE), acf = rho)
}

# ------------------------------------------------------------------ joint draws for a fit with a year random effect

# Mean and covariance of (identified beta, unit random effects) from the joint precision of a glmmTMB fit
# with a unit term and possibly other terms (e.g. (1 | yearf)). The marginal covariance of that block is
# the inverse of the Schur complement of the other random effects (exact, and stable when a variance is
# near zero). The unit block is checked against ranef() modes.
rev_joint_unit <- function(fit) {
  m <- fit$model
  b_hat <- glmmTMB::fixef(m)$cond; ok <- !is.na(b_hat)
  rs <- m$modelInfo$reStruc$condReStruc
  sizes <- vapply(rs, function(z) as.integer(z$blockSize * z$blockReps), integer(1))
  grp <- sub("^.*\\|\\s*", "", names(rs))
  ku <- which(grp == "unit_id")
  if (length(ku) != 1L) stop("rev_joint_unit: expected one unit_id term", call. = FALSE)
  modes <- as.matrix(glmmTMB::ranef(m)$cond$unit_id)
  U <- nrow(modes); K <- ncol(modes)
  if (sizes[ku] != U * K) stop("rev_joint_unit: unit block size mismatch", call. = FALSE)
  Q <- .joint_precision(fit)
  if (is.null(Q)) stop("rev_joint_unit: joint precision unavailable", call. = FALSE)
  Q <- as.matrix(Q)
  nm <- rownames(Q); ib <- which(nm == "beta"); ibb <- which(nm == "b")
  if (length(ib) != sum(ok) || length(ibb) != sum(sizes)) stop("rev_joint_unit: joint precision layout", call. = FALSE)
  off <- c(0L, cumsum(sizes))[ku]
  iu <- ibb[off + seq_len(sizes[ku])]
  par <- m$obj$env$last.par.best
  b_vec <- par[m$obj$env$random]
  bu <- b_vec[off + seq_len(sizes[ku])]
  if (max(abs(matrix(bu, U, K, byrow = TRUE) - modes)) > 1e-8)
    stop("rev_joint_unit: unit random effects do not map unit-major onto ranef() modes", call. = FALSE)
  # condition on the variance parameters of the other (year) terms: a year SD at the boundary has a
  # degenerate curvature that makes the joint precision indefinite; the unit-term variance parameters,
  # the dispersion and the other random effects are integrated over as in .joint_beta_b()
  it <- which(nm == "theta")
  nth <- vapply(rs, function(z) as.integer(z$blockNumTheta), integer(1))
  if (length(it) != sum(nth)) stop("rev_joint_unit: theta layout", call. = FALSE)
  drop_th <- it[-(c(0L, cumsum(nth))[ku] + seq_len(nth[ku]))]
  A <- c(ib, iu); B <- setdiff(seq_len(nrow(Q)), c(A, drop_th))
  S <- Q[A, A, drop = FALSE]
  if (length(B)) {
    RB <- chol(Q[B, B, drop = FALSE])
    W <- backsolve(RB, Q[B, A, drop = FALSE], transpose = TRUE)
    S <- S - crossprod(W)
  }
  C <- .invert_joint_precision((S + t(S)) / 2)
  list(mean = c(b_hat[ok], bu), cov = (C + t(C)) / 2, ok = ok, modes = modes, U = U, K = K)
}

# Joint parameter draws (beta, unit RE array) for such a fit; same layout as .draw_params().
rev_draw_unit <- function(fit, n_draws, seed) {
  b_hat <- glmmTMB::fixef(fit$model)$cond
  .with_seed(seed, {
    jb <- rev_joint_unit(fit)
    p <- sum(jb$ok)
    D <- matrix(MASS::mvrnorm(n_draws, jb$mean, jb$cov), nrow = n_draws)
    B <- matrix(0, n_draws, length(b_hat), dimnames = list(NULL, names(b_hat)))
    B[, jb$ok] <- D[, seq_len(p), drop = FALSE]
    bb <- D[, p + seq_len(jb$U * jb$K), drop = FALSE]
    R <- array(NA_real_, c(jb$U, jb$K, n_draws), dimnames = c(dimnames(jb$modes), list(NULL)))
    for (k in seq_len(jb$K)) R[, k, ] <- t(bb[, (seq_len(jb$U) - 1L) * jb$K + k, drop = FALSE])
    list(beta_draws = B, re_draws = R)
  })
}

# mu draws (rows x draws) for newdata from (beta, unit RE) draws plus a new-year effect per draw:
# year_draws is a years x draws matrix (rownames = years) of year effects, or NULL.
rev_mu_draws <- function(fit, newdata, pd, year_draws = NULL) {
  p <- .fixed_design(fit, newdata)
  eta <- p$X %*% t(pd$beta_draws[, colnames(p$X), drop = FALSE]) + p$off + .re_contrib(pd$re_draws, newdata)
  if (!is.null(year_draws)) eta <- eta + year_draws[as.character(newdata$year), , drop = FALSE]
  mu <- exp(eta); dimnames(mu) <- NULL
  mu
}

# ------------------------------------------------------------------ early-warning validation (plan R2.5)

# Information available by 1 June of year t: JJAS ONI of t - 1, March-May ONI of year t (published early
# June), local climate at lag 1, log ABER of t - 1, previous-year cases, the indicators known in advance
# (programme and COVID steps), the trend in calendar year and the population offset.
EW_LOCAL_CLIMATE <- c("rain_jjas_dep_pct_lag1", "spi3_sep_lag1", "tmax_mamj_anom_lag1")
EW_ONI <- c("oni_jjas_lag1", "oni_mam")
EW_KNOWN_T <- c("oni_mam", "nfme2016", "covid2020_21", "log_pop", "year")   # year-t columns allowed

ew_covars <- function(outcome = c("malaria", "dengue")) {
  outcome <- match.arg(outcome)
  if (outcome == "malaria") c("log_aber_lag1", "nfme2016", "covid2020_21") else "covid2020_21"
}

ew_specs <- function(covars) {
  A <- nb_spec("cases", fixed = c(covars, "log1p_cases_lag1"), re = "int")
  B <- nb_spec("cases", fixed = c(A$fixed, EW_LOCAL_CLIMATE), re = "int")
  C <- nb_spec("cases", fixed = c(B$fixed, EW_ONI), slopes = EW_ONI, re = "us")
  list(A = A, B = B, C = C)
}

# Variables of a spec that carry year-t information other than the allowed set (should be empty).
ew_year_t_vars <- function(spec) {
  v <- setdiff(spec_vars(spec), spec$outcome)
  setdiff(v[!grepl("_lag1$", v)], EW_KNOWN_T)
}

# Adds oni_mam (MAM season ONI of the same year, from the seasonal ONI table) and log_aber_lag1 (the
# unit's log ABER of the previous calendar year; NA when that year is absent or has no ABER).
add_ew_terms <- function(d, oni = rev_oni_seasonal()) {
  d <- data.table::copy(data.table::as.data.table(d))
  mam <- oni[season == "MAM", .(year, oni_mam = oni)]
  d[mam, oni_mam := i.oni_mam, on = "year"]
  if ("log_aber" %in% names(d)) {
    prev <- d[, .(unit_id, year = year + 1L, log_aber_lag1 = log_aber)]
    d[prev, log_aber_lag1 := i.log_aber_lag1, on = c("unit_id", "year")]
  }
  d[]
}

ew_validation_jobs <- function() {
  folds <- read_folds("developing")
  loeo <- loeo_splits(folds)
  list(
    mal_loeo = list(data = "malaria", specs = ew_specs(ew_covars("malaria")), splits = loeo, seed = SEED + 1L),
    mal_roll = list(data = "malaria", specs = ew_specs(ew_covars("malaria")),
                    splits = rolling_splits(2002:2025, 2015), seed = SEED + 2L),
    den_loeo = list(data = "dengue", specs = ew_specs(ew_covars("dengue")), splits = loeo, seed = SEED + 3L),
    den_roll = list(data = "dengue", specs = ew_specs(ew_covars("dengue")),
                    splits = rolling_splits(2007:2025, 2015), seed = SEED + 4L))
}

ew_data <- function(outcome) add_ew_terms(if (outcome == "malaria") malaria_data() else dengue_data())

# Coverage of nominal 95% prediction intervals that add the year-level variance: per split and model, the
# same spec plus an iid year random effect is fitted to the same training rows as in compare_models();
# (beta, unit RE) are drawn jointly; each test year gets a new year effect ~ N(0, sd_year^2) per draw,
# shared by all units; NB2 noise is added. Scored on exactly the rows scored in `cmp`. Seeds as cv_predict.
ew_yearvar_coverage <- function(specs, data, splits, cmp, n_draws = 1000L, seed = 1L, log_file = NULL) {
  common <- complete_rows(data.table::as.data.table(data), unique(unlist(lapply(specs, spec_vars))))$data
  lag_buffer <- isTRUE(cmp$lag_buffer)
  rows <- list(); info <- list()
  for (i in seq_along(splits)) {
    sp <- splits[[i]]; sid <- as.character(sp$id)
    train <- common[!year %in% sp$test_years]
    if (!is.null(sp$train_max)) train <- train[year <= sp$train_max]
    if (lag_buffer) train <- train[!(year - 1L) %in% sp$test_years]
    for (mname in names(specs)) {
      keys <- cmp$per_split[model == mname & scored == TRUE & split_id == sid, .(unit_id, year, y, lo, hi)]
      if (!nrow(keys)) next
      test <- common[keys[, .(unit_id, year)], on = c("unit_id", "year")]
      test[, unit_id := as.character(unit_id)]
      split_seed <- seed + i - 1L
      fit <- tryCatch(fit_nb_extra(specs[[mname]], train, REV_YEAR_IID), error = function(e) e)
      if (inherits(fit, "error")) {
        info[[length(info) + 1L]] <- data.table::data.table(model = mname, split_id = sid, status = "fit_failed",
                                                            sd_year = NA_real_, re_used = NA_character_)
        next
      }
      sdy <- rev_group_sd(fit, "yearf")
      pd <- tryCatch(rev_draw_unit(fit, n_draws, split_seed), error = function(e) e)
      if (inherits(pd, "error")) {
        info[[length(info) + 1L]] <- data.table::data.table(model = mname, split_id = sid, status = "draw_failed",
                                                            sd_year = sdy, re_used = fit$re_used)
        next
      }
      # unidentified coefficients (ruling P5) stay 0 in every draw (rev_draw_unit draws the identified block)
      yrs <- sort(unique(test$year))
      yd <- .with_seed(split_seed + 200000L,
                       matrix(stats::rnorm(length(yrs) * n_draws, 0, sdy), length(yrs), n_draws,
                              dimnames = list(as.character(yrs), NULL)))
      mu <- rev_mu_draws(fit, test, pd, yd)
      ydr <- draw_y(fit, mu, seed = split_seed + 100000L)
      q <- apply(ydr, 1, stats::quantile, probs = c(0.025, 0.975), names = FALSE)
      rows[[length(rows) + 1L]] <- data.table::data.table(model = mname, split_id = sid,
                                                          unit_id = as.character(test$unit_id), year = test$year,
                                                          y = test$cases, lo_yv = q[1, ], hi_yv = q[2, ])
      info[[length(info) + 1L]] <- data.table::data.table(model = mname, split_id = sid, status = "ok",
                                                          sd_year = sdy, re_used = fit$re_used)
      if (!is.null(log_file)) alog(sprintf("  yearvar %s split %s sd_year %.4f", mname, sid, sdy), file = log_file)
    }
  }
  pr <- data.table::rbindlist(rows); inf <- data.table::rbindlist(info)
  summ <- pr[, .(coverage_yearvar = mean(y >= lo_yv & y <= hi_yv), n = .N), by = model]
  list(pred = pr, info = inf, summary = summ)
}

# ------------------------------------------------------------------ conditional scenarios (plan R2.6)

# Expected national excess under JJAS ONI = v relative to ONI = 0, for each v, with the primary-model
# machinery of 04_outlook.R (outlook_newdata + project, same draws and seed for every v).
rev_scenarios_primary <- function(fit, values, terms = INDEX_TERMS$oni, n_draws = N_DRAWS, seed = SEED + 7L) {
  data.table::rbindlist(lapply(values, function(v) {
    nd <- outlook_newdata(fit, terms, v)
    pr <- project(fit, nd$obs, nd$cf, terms, n_draws = n_draws, seed = seed)
    cbind(data.table::data.table(model = "primary", oni = v), pr$total)
  }))
}

# The same with a fit that has an iid year random effect: (beta, unit RE) drawn jointly, plus a new-year
# effect ~ N(0, sd_year^2) per draw and year shared by all units and by both scenarios.
rev_scenarios_year <- function(fit, values, terms = INDEX_TERMS$oni, n_draws = N_DRAWS, seed = SEED + 7L,
                               newdata_fun = outlook_newdata) {
  pd <- rev_draw_unit(fit, n_draws, seed)
  sdy <- rev_group_sd(fit, "yearf")
  data.table::rbindlist(lapply(values, function(v) {
    nd <- newdata_fun(fit, terms, v)
    obs <- data.table::copy(nd$obs); cf <- data.table::copy(nd$cf)
    yrs <- sort(unique(obs$year))
    yd <- .with_seed(seed + 200000L, matrix(stats::rnorm(length(yrs) * n_draws, 0, sdy), length(yrs), n_draws,
                                            dimnames = list(as.character(yrs), NULL)))
    ex <- rev_mu_draws(fit, obs, pd, yd) - rev_mu_draws(fit, cf, pd, yd)
    data.table::rbindlist(lapply(yrs, function(y) {
      tot <- colSums(ex[obs$year == y, , drop = FALSE])
      q <- stats::quantile(tot, c(0.025, 0.975), names = FALSE)
      data.table::data.table(model = "year_iid", oni = v, year = y, n_rows = sum(obs$year == y),
                             excess_mean = mean(tot), excess_lo = q[1], excess_hi = q[2])
    }))
  }))
}

# ------------------------------------------------------------------ zones (plan R2.2 item 3a)

# Zonal councils (States Reorganisation Act 1956) plus the North Eastern Council; the island territories
# go with the Southern zone (plan R2.2 item 3).
REV_ZONES <- c(
  CH = "Northern", DL = "Northern", HR = "Northern", HP = "Northern", JK = "Northern", PB = "Northern",
  RJ = "Northern",
  CG = "Central", MP = "Central", UK = "Central", UP = "Central",
  BR = "Eastern", JH = "Eastern", OD = "Eastern", WB = "Eastern",
  DNH_DD = "Western", GA = "Western", GJ = "Western", MH = "Western",
  AN = "Southern", AP_TG = "Southern", KA = "Southern", KL = "Southern", LD = "Southern", PY = "Southern",
  TN = "Southern",
  AR = "North-Eastern", AS = "North-Eastern", MN = "North-Eastern", ML = "North-Eastern", MZ = "North-Eastern",
  NL = "North-Eastern", SK = "North-Eastern", TR = "North-Eastern")
REV_ZONE_REF <- "Eastern"
zone_key <- function(z) gsub("[^a-z]", "", tolower(z))

# Adds zone indicators (z_<zone>) and ONI x zone columns for the non-reference zones.
add_zone_terms <- function(d) {
  d <- data.table::copy(d)
  d[, zone := REV_ZONES[as.character(unit_id)]]
  if (anyNA(d$zone)) stop("add_zone_terms: unit without a zone", call. = FALSE)
  for (z in setdiff(unique(REV_ZONES), REV_ZONE_REF)) {
    k <- zone_key(z)
    data.table::set(d, j = paste0("z_", k), value = as.numeric(d$zone == z))
    data.table::set(d, j = paste0("oni_x_", k), value = d$oni_jjas * (d$zone == z))
    data.table::set(d, j = paste0("oni_lag1_x_", k), value = d$oni_jjas_lag1 * (d$zone == z))
  }
  d[]
}

# ------------------------------------------------------------------ unpooled unit fits (plan R2.2 item 3b)

REV_MIN_YEARS <- 15L   # a unit needs at least this many complete years for its own NB GLM

unpooled_unit_fit <- function(du, terms = c(INDEX_TERMS$oni, MAL_COVARS), trend_df = 4L) {
  du <- complete_rows(data.table::as.data.table(du), c("cases", terms, "log_pop", "year"))$data
  out <- data.table::data.table(unit_id = as.character(du$unit_id[1]), n = nrow(du), irr_lag0 = NA_real_,
                                lo = NA_real_, hi = NA_real_, log_irr_lag0 = NA_real_, se = NA_real_,
                                irr_lag1 = NA_real_, lo_lag1 = NA_real_, hi_lag1 = NA_real_, theta = NA_real_,
                                status = "ok")
  if (nrow(du) < REV_MIN_YEARS) { out$status <- sprintf("too few years (%d < %d)", nrow(du), REV_MIN_YEARS); return(out) }
  f <- stats::as.formula(paste("cases ~", paste(terms, collapse = " + "),
                               sprintf("+ splines::ns(year, df = %d) + offset(log_pop)", trend_df)))
  warn <- character()
  m <- tryCatch(withCallingHandlers(MASS::glm.nb(f, data = du, control = stats::glm.control(maxit = 100)),
                                    warning = function(w) { warn <<- c(warn, conditionMessage(w)); invokeRestart("muffleWarning") }),
                error = function(e) e)
  if (inherits(m, "error")) { out$status <- paste("error:", conditionMessage(m)); return(out) }
  co <- summary(m)$coefficients
  if (!all(INDEX_TERMS$oni %in% rownames(co)) || !isTRUE(m$converged) || length(warn) ||
      any(!is.finite(co[INDEX_TERMS$oni, 2]))) {
    out$status <- paste("not converged:", paste(unique(c(if (!isTRUE(m$converged)) "glm.nb did not converge", warn)),
                                                collapse = "; "))
    if (!nzchar(sub("not converged: ", "", out$status))) out$status <- "not converged: non-finite or aliased ONI terms"
    return(out)
  }
  b0 <- co["oni_jjas", 1]; s0 <- co["oni_jjas", 2]; b1 <- co["oni_jjas_lag1", 1]; s1 <- co["oni_jjas_lag1", 2]
  out[, `:=`(irr_lag0 = exp(b0), lo = exp(b0 - 1.959964 * s0), hi = exp(b0 + 1.959964 * s0), log_irr_lag0 = b0,
             se = s0, irr_lag1 = exp(b1), lo_lag1 = exp(b1 - 1.959964 * s1), hi_lag1 = exp(b1 + 1.959964 * s1),
             theta = m$theta)]
  out
}

# ------------------------------------------------------------------ nonlinear climate (plan R2.4)

REV_CLIM_VARS <- c("tmax_mamj_anom", "tmax_mamj_anom_lag1", "rain_jjas_dep_pct", "rain_jjas_dep_pct_lag1")
rev_ns_terms <- function(vars = REV_CLIM_VARS, df = 3L) sprintf("splines::ns(%s, df = %d)", vars, df)

# Partial response of one spline (or linear) term relative to the variable's value 0, with Wald bands.
partial_response <- function(fit, var, grid, base) {
  nd <- base[rep(1L, length(grid) + 1L)]
  data.table::set(nd, j = var, value = c(0, grid))
  X <- .fixed_design(fit, nd)$X
  D <- sweep(X[-1, , drop = FALSE], 2, X[1, ])
  b <- glmmTMB::fixef(fit$model)$cond[colnames(X)]
  b[is.na(b)] <- 0
  V <- stats::vcov(fit$model)$cond
  keep <- colnames(X)[colSums(abs(D)) > 0]
  est <- as.numeric(D[, keep, drop = FALSE] %*% b[keep])
  se <- sqrt(pmax(rowSums((D[, keep, drop = FALSE] %*% V[keep, keep, drop = FALSE]) * D[, keep, drop = FALSE]), 0))
  data.table::data.table(variable = var, x = grid, log_rr = est, rr = exp(est), lo = exp(est - 1.959964 * se),
                         hi = exp(est + 1.959964 * se))
}

# Unit climatological March-June Tmax 1981-2010 (from data/derived/climate_unit_year.csv) and terciles.
tmax_climatology <- function(path = file.path(REPO_ROOT, "data", "derived", "climate_unit_year.csv")) {
  cy <- data.table::fread(path)
  cl <- cy[year %in% 1981:2010, .(tmax_clim = mean(tmax_mamj_mean, na.rm = TRUE), n_years = sum(is.finite(tmax_mamj_mean))),
           by = unit_id]
  q <- stats::quantile(cl$tmax_clim, c(1 / 3, 2 / 3), names = FALSE)
  cl[, tercile := as.character(cut(tmax_clim, c(-Inf, q, Inf), labels = c("cool", "middle", "hot")))]
  list(table = cl[], cut = q)
}

add_tercile_terms <- function(d, clim) {
  d <- data.table::copy(d)
  d[clim, tercile := i.tercile, on = "unit_id"]
  for (t in c("cool", "middle", "hot")) {
    data.table::set(d, j = paste0("tmax_", t), value = d$tmax_mamj_anom * (d$tercile == t))
    data.table::set(d, j = paste0("tmax_lag1_", t), value = d$tmax_mamj_anom_lag1 * (d$tercile == t))
  }
  d[]
}
