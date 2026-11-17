# Paper A revision (Acta Tropica ACTROP-D-26-02276): reviewer analyses fixed in
# papers/A/revision/revision_plan.md before they were run (R2.1 items 1-7, R2.2 item 3, R2.3 items 2-3,
# R2.4 items 1-2, R2.5, R2.6). Run from the models/ folder after A/run_all.R has filled models/A/_cache:
#   Rscript A/07_revision.R --step=fits      # single refits (temporal, episodes, trend, ABER, climate, zones)
#   Rscript A/07_revision.R --step=boot      # 999 year-block bootstrap refits
#   Rscript A/07_revision.R --step=ew        # early-warning validation (compare_models) + year-variance coverage
#   Rscript A/07_revision.R --step=results   # tables (models/A/output/rev_*.csv), scenarios, Figures S2-S5,
#                                            # and papers/A/results.json = the 191 original keys + A.rev.* keys
#   Rscript A/07_revision.R                  # all four steps in order
#   --cores=N  PSOCK workers (default 4). Every fit is cached in models/A/_cache/rev_*.rds (or
#   _cache/rev_boot/) and reused, so an interrupted run resumes; --fresh ignores the revision caches.
# The original keys are rebuilt from the cached original fits by build_results() (05_results.R), exactly as
# run_all.R --step=results does, so they stay byte-identical. NOTE: run_all.R --step=results on its own
# writes the 191 original keys only; run this script's results step afterwards to add the A.rev.* keys.

.rev_main <- sys.nframe() == 0L     # TRUE when run with Rscript; FALSE when sourced (e.g. by PSOCK workers)
if (.rev_main) {
  .rev_args <- commandArgs(trailingOnly = TRUE)
  options(paperA.rev_fresh = "--fresh" %in% .rev_args)
  Sys.setenv(OMP_NUM_THREADS = "1")
  source("common/R/load.R")
  source("A/common_a.R")
  for (.f in c("01_fits.R", "02_validation.R", "03_mediation.R", "04_outlook.R", "05_results.R", "conclusion_text.R",
               "revision_a.R"))
    source(file.path("A", .f))
}

REV_SCRIPT <- "models/A/07_revision.R"
REV_FIG_DIR <- file.path(REPO_ROOT, "papers", "A", "figures")
REV_BOOT_DIR <- file.path(A_CACHE, "rev_boot")
REV_B <- 999L
REV_FRESH <- isTRUE(getOption("paperA.rev_fresh", FALSE))

rev_cached <- function(name, expr) {
  f <- file.path(A_CACHE, paste0("rev_", name, ".rds"))
  if (!REV_FRESH && file.exists(f)) return(readRDS(f))
  v <- expr
  saveRDS(v, f)
  v
}
rev_read <- function(name) {
  f <- file.path(A_CACHE, paste0("rev_", name, ".rds"))
  if (!file.exists(f)) stop("revision cache missing: ", f, " (run the step that makes it)", call. = FALSE)
  readRDS(f)
}

# ------------------------------------------------------------------ episodes

rev_episodes <- function() {
  ep <- derive_episodes(rev_oni_seasonal(), min_peak = 1.0)
  # moderate or stronger episodes that remove at least one panel year (2002-2025)
  ep[, in_panel := vapply(years_dropped, function(s) any(episode_years(s) %in% 2002:2025), logical(1))]
  ep[in_panel == TRUE][, in_panel := NULL][]
}

# ------------------------------------------------------------------ single refits (one job = one fit)

rev_fit_jobs <- function() {
  P <- primary_spec()
  mal <- function() malaria_data()
  J <- list()
  add <- function(name, spec, data, extra = NULL, extra_vars = character())
    J[[name]] <<- list(spec = spec, data = data, extra = extra, extra_vars = extra_vars)
  # R2.1 items 1-2: shared year effect
  add("year_iid", P, mal, REV_YEAR_IID)
  add("year_ar1", P, mal, REV_YEAR_AR1)
  # R2.1 item 7: leave one moderate-or-stronger ENSO episode out
  ep <- rev_episodes()
  for (i in seq_len(nrow(ep))) local({
    yrs <- episode_years(ep$years_dropped[i])
    add(paste0("loepo_", ep$episode_id[i]), P, function() malaria_data()[!year %in% yrs])
  })
  # R2.3 item 2: trend specifications
  for (k in c(3L, 6L, 8L)) local({ s <- P; s$trend_df <- k; add(sprintf("trend_ns%d", k), s, mal) })
  add("trend_unit_linear", P, mal, REV_UNIT_TREND)
  # R2.3 item 3: diagnostic effort
  add("aber_none", primary_spec(drop = "log_aber"), mal)
  add("aber_lag1", primary_spec(covars = c("log_aber_lag1", "nfme2016", "covid2020_21")),
      function() add_ew_terms(malaria_data()))
  add("aber_ns3", primary_spec(drop = "log_aber"), mal, "splines::ns(log_aber, df = 3)", "log_aber")
  pos <- primary_spec(drop = "log_aber"); pos$offset <- "log_tested"
  add("aber_positivity", pos, function() {
    d <- malaria_data(); d[, log_tested := ifelse(tested > 0, log(tested), NA_real_)]; d })
  # R2.4 item 1: linear vs spline local climate, same rows
  add("nonlin_linear", primary_spec(extra = REV_CLIM_VARS), mal)
  add("nonlin_spline", P, mal, paste(rev_ns_terms(), collapse = " + "), REV_CLIM_VARS)
  # R2.4 item 2: Tmax anomaly slope by climatology tercile
  add("tmax_tercile", primary_spec(extra = c(paste0("tmax_", c("cool", "middle", "hot")),
                                             paste0("tmax_lag1_", c("cool", "middle", "hot")))),
      function() add_tercile_terms(malaria_data(), tmax_climatology()$table))
  # R2.2 item 3a: ONI x zone
  zk <- zone_key(setdiff(unique(REV_ZONES), REV_ZONE_REF))
  add("zone", primary_spec(extra = c(paste0("z_", zk), paste0("oni_x_", zk), paste0("oni_lag1_x_", zk))),
      function() add_zone_terms(malaria_data()))
  J
}

.rev_run_fit_job <- function(name) {
  f <- file.path(A_CACHE, paste0("rev_fit_", name, ".rds"))
  if (!REV_FRESH && file.exists(f)) return(name)
  j <- rev_fit_jobs()[[name]]
  t0 <- Sys.time()
  alog("fit ", name, file = "rev_fits.log")
  d <- j$data()
  fit <- tryCatch(if (is.null(j$extra)) fit_nb(j$spec, d) else fit_nb_extra(j$spec, d, j$extra, j$extra_vars),
                  error = function(e) e)
  saveRDS(fit, f)
  alog(if (inherits(fit, "error")) paste("  FAILED:", conditionMessage(fit))
       else sprintf("  %s: n = %d, re_used = %s (%.1f min)", name, fit$n, fit$re_used,
                    as.numeric(difftime(Sys.time(), t0, units = "mins"))), file = "rev_fits.log")
  name
}

# ------------------------------------------------------------------ year-block bootstrap

rev_boot_sequences <- function(B = REV_B, years = 2002:2025) {
  .with_seed(REV_SEED, lapply(seq_len(B), function(i) block_resample_years(years, 3L)))
}

.rev_run_boot_rep <- function(i) {
  f <- file.path(REV_BOOT_DIR, sprintf("rep_%04d.rds", i))
  if (file.exists(f)) return(i)
  yseq <- rev_boot_sequences()[[i]]
  d <- boot_panel(malaria_data(), yseq)
  s <- primary_spec(); s$re <- "diag"          # the structure the primary model converged with
  fit <- tryCatch(fit_nb(s, d), error = function(e) e)
  out <- if (inherits(fit, "error")) list(i = i, ok = FALSE, error = conditionMessage(fit), re_used = NA_character_,
                                          beta = c(oni_jjas = NA_real_, oni_jjas_lag1 = NA_real_), n = nrow(d))
         else list(i = i, ok = TRUE, error = NA_character_, re_used = fit$re_used,
                   beta = glmmTMB::fixef(fit$model)$cond[INDEX_TERMS$oni], n = fit$n)
  saveRDS(out, f)
  i
}

# ------------------------------------------------------------------ early warning

.rev_run_ew_job <- function(name) {
  f <- file.path(A_CACHE, paste0("rev_ew_", name, ".rds"))
  if (!REV_FRESH && file.exists(f)) return(name)
  job <- ew_validation_jobs()[[name]]
  t0 <- Sys.time()
  alog("start ew ", name, file = paste0("rev_ew_", name, ".log"))
  res <- compare_models(job$specs, ew_data(job$data), job$splits, n_draws = N_DRAWS, seed = job$seed)
  res$job <- name; res$seed <- job$seed
  saveRDS(res, f)
  alog(sprintf("done ew %s in %.1f min", name, as.numeric(difftime(Sys.time(), t0, units = "mins"))),
       file = paste0("rev_ew_", name, ".log"))
  name
}

.rev_run_ew_yv_job <- function(name) {
  f <- file.path(A_CACHE, paste0("rev_ewyv_", name, ".rds"))
  if (!REV_FRESH && file.exists(f)) return(name)
  job <- ew_validation_jobs()[[name]]
  cmp <- readRDS(file.path(A_CACHE, paste0("rev_ew_", name, ".rds")))
  t0 <- Sys.time()
  alog("start ew year-variance ", name, file = paste0("rev_ewyv_", name, ".log"))
  res <- ew_yearvar_coverage(job$specs, ew_data(job$data), job$splits, cmp, n_draws = N_DRAWS, seed = job$seed,
                             log_file = paste0("rev_ewyv_", name, ".log"))
  saveRDS(res, f)
  alog(sprintf("done ew year-variance %s in %.1f min", name, as.numeric(difftime(Sys.time(), t0, units = "mins"))),
       file = paste0("rev_ewyv_", name, ".log"))
  name
}

# ------------------------------------------------------------------ parallel runner

.rev_worker_init <- function(root, libs, fresh) {
  .libPaths(libs)
  setwd(root)
  Sys.setenv(OMP_NUM_THREADS = "1")
  source("common/R/load.R")
  source("A/common_a.R")
  source("A/04_outlook.R")
  source("A/revision_a.R")
  options(paperA.rev_fresh = fresh)
  source("A/07_revision.R")
  options(glmmTMB.cores = 1L)
  invisible(TRUE)
}

rev_parallel <- function(items, fun, cores) {
  if (!length(items)) return(invisible(NULL))
  cl <- parallel::makePSOCKcluster(min(cores, length(items)))
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterCall(cl, .rev_worker_init, MODELS_ROOT, .libPaths(), REV_FRESH)
  parallel::clusterApplyLB(cl, items, fun)
}

run_rev_fits <- function(cores) {
  todo <- names(rev_fit_jobs())
  if (!REV_FRESH) todo <- todo[!file.exists(file.path(A_CACHE, paste0("rev_fit_", todo, ".rds")))]
  alog("revision fits to run: ", paste(todo, collapse = ", "), file = "rev_fits.log")
  rev_parallel(todo, .rev_run_fit_job, cores)
}

run_rev_boot <- function(cores) {
  dir.create(REV_BOOT_DIR, showWarnings = FALSE)
  if (REV_FRESH) unlink(list.files(REV_BOOT_DIR, full.names = TRUE))
  todo <- which(!file.exists(file.path(REV_BOOT_DIR, sprintf("rep_%04d.rds", seq_len(REV_B)))))
  alog(sprintf("bootstrap: %d of %d replicates to run", length(todo), REV_B), file = "rev_boot.log")
  t0 <- Sys.time()
  # chunks so progress is logged
  for (ch in split(todo, ceiling(seq_along(todo) / 100))) {
    rev_parallel(ch, .rev_run_boot_rep, cores)
    alog(sprintf("bootstrap: %d done (%.1f min)", sum(file.exists(file.path(REV_BOOT_DIR, sprintf("rep_%04d.rds", seq_len(REV_B))))),
                 as.numeric(difftime(Sys.time(), t0, units = "mins"))), file = "rev_boot.log")
  }
}

run_rev_ew <- function(cores) {
  jobs <- names(ew_validation_jobs())
  rev_parallel(jobs[REV_FRESH | !file.exists(file.path(A_CACHE, paste0("rev_ew_", jobs, ".rds")))], .rev_run_ew_job, cores)
  rev_parallel(jobs[REV_FRESH | !file.exists(file.path(A_CACHE, paste0("rev_ewyv_", jobs, ".rds")))], .rev_run_ew_yv_job, cores)
}

# ------------------------------------------------------------------ summaries (results step)

rev_fit <- function(name) rev_read(paste0("fit_", name))

rev_summaries <- function(fits, outlook) {
  fp <- fits$fits$malaria_primary
  S <- list(primary = fp)
  terms <- INDEX_TERMS$oni
  mal <- malaria_data()

  # ---- R2.1 inference variants
  inf <- list()
  inf_row <- function(variant, label, r0, r1, extra = list())
    data.table::data.table(variant = variant, label = label, irr_lag0 = r0$irr, lo_lag0 = r0$lo, hi_lag0 = r0$hi,
                           irr_lag1 = r1$irr, lo_lag1 = r1$lo, hi_lag1 = r1$hi,
                           display_lag0 = r0$display, display_lag1 = r1$display,
                           note = if (length(extra)) paste(names(extra), unlist(extra), sep = " = ", collapse = "; ") else "")
  wr <- function(r) { r$display <- disp_irr(r); r }
  inf$primary <- inf_row("primary", "Primary model (Wald)", rev_irr(fp, "oni_jjas"), rev_irr(fp, "oni_jjas_lag1"))
  for (v in c("year_iid", "year_ar1")) {
    f <- rev_fit(v)
    if (inherits(f, "error")) {
      inf[[v]] <- data.table::data.table(variant = v, label = v, irr_lag0 = NA_real_, lo_lag0 = NA_real_, hi_lag0 = NA_real_,
                                         irr_lag1 = NA_real_, lo_lag1 = NA_real_, hi_lag1 = NA_real_,
                                         display_lag0 = "not estimated (did not converge)",
                                         display_lag1 = "not estimated (did not converge)",
                                         note = conditionMessage(f))
    } else {
      ex <- list(sd_year = signif(rev_group_sd(f, if (v == "year_iid") "yearf" else "one"), 4), re_used = f$re_used, n = f$n)
      if (v == "year_ar1") ex$rho <- signif(rev_ar1_rho(f), 4)
      inf[[v]] <- inf_row(v, if (v == "year_iid") "Primary + iid year effect" else "Primary + AR(1) year effect",
                          rev_irr(f, "oni_jjas"), rev_irr(f, "oni_jjas_lag1"), ex)
    }
  }
  S$year_fits <- list(year_iid = rev_fit("year_iid"), year_ar1 = rev_fit("year_ar1"))
  # cluster-robust NB GLM with unit fixed effects
  g <- rev_cached("glm_unitfe", glm_nb_unitfe(mal))
  co <- stats::coef(g$model)
  G_year <- data.table::uniqueN(g$data$year); G_unit <- data.table::uniqueN(g$data$unit_id)
  V0 <- stats::vcov(g$model)
  Vy <- sandwich::vcovCL(g$model, cluster = g$data$year, type = "HC1")
  V2 <- sandwich::vcovCL(g$model, cluster = data.frame(year = g$data$year, unit = as.character(g$data$unit_id)),
                         type = "HC1", multi0 = FALSE)
  mk <- function(V, df) lapply(terms, function(t) wr(rev_t_irr(co[[t]], V[t, t], df)))
  m0 <- mk(V0, Inf); my <- mk(Vy, G_year - 1L); m2 <- mk(V2, min(G_year, G_unit) - 1L)
  inf$glm_model <- inf_row("glm_model", "NB GLM, unit fixed effects (model-based)", m0[[1]], m0[[2]],
                           list(theta = signif(g$model$theta, 4), n = nrow(g$data)))
  inf$cluster_year <- inf_row("cluster_year", "NB GLM, cluster-robust by year (HC1, t(G-1))", my[[1]], my[[2]],
                              list(G = G_year, df = G_year - 1L))
  inf$cluster_twoway <- inf_row("cluster_twoway", "NB GLM, two-way cluster-robust (year and unit)", m2[[1]], m2[[2]],
                                list(G_year = G_year, G_unit = G_unit, df = min(G_year, G_unit) - 1L))
  S$glm <- list(G_year = G_year, G_unit = G_unit, n = nrow(g$data), theta = g$model$theta,
                rows = list(model = m0, year = my, twoway = m2))
  # year-block bootstrap
  bf <- file.path(REV_BOOT_DIR, sprintf("rep_%04d.rds", seq_len(REV_B)))
  if (!all(file.exists(bf))) stop("bootstrap incomplete: ", sum(file.exists(bf)), " of ", REV_B, call. = FALSE)
  bt <- data.table::rbindlist(lapply(bf, function(p) {
    r <- readRDS(p)
    data.table::data.table(rep = r$i, ok = r$ok, re_used = r$re_used, b_lag0 = r$beta[["oni_jjas"]],
                           b_lag1 = r$beta[["oni_jjas_lag1"]], n = r$n, error = r$error)
  }))
  S$boot <- bt
  okb <- bt[ok == TRUE]
  bq <- function(b) exp(stats::quantile(b, c(0.025, 0.975), names = FALSE))
  pe <- irr_table(fp, terms)
  q0 <- bq(okb$b_lag0); q1 <- bq(okb$b_lag1)
  b0 <- list(irr = pe$irr[1], lo = q0[1], hi = q0[2]); b0$display <- disp_irr(b0)
  b1 <- list(irr = pe$irr[2], lo = q1[1], hi = q1[2]); b1$display <- disp_irr(b1)
  inf$boot <- inf_row("boot", "Primary model, year-block bootstrap (blocks of 3 years, percentile)", b0, b1,
                      list(B = REV_B, ok = nrow(okb), failed = sum(!bt$ok), int_fallback = sum(okb$re_used == "int")))
  S$inference <- data.table::rbindlist(inf, fill = TRUE)

  # ---- R2.1 item 5 residual diagnostics
  rq <- rqr_nb(fp$data$cases, stats::fitted(fp$model), stats::sigma(fp$model))
  rd <- data.table::data.table(unit_id = as.character(fp$data$unit_id), year = fp$data$year, y = fp$data$cases,
                               mu = stats::fitted(fp$model), u = rq$u, r = rq$r,
                               pearson = stats::residuals(fp$model, type = "pearson"))
  ym <- data.table::data.table(year = 2002:2025)[rd[, .(mean_r = mean(r), n = .N), by = year], on = "year",
                                                 `:=`(mean_r = i.mean_r, n = i.n)]
  lb <- lapply(1:3, function(k) ljung_box_na(ym$mean_r, k))
  acf_all <- stats::acf(ym$mean_r, lag.max = 8, plot = FALSE, na.action = stats::na.pass)$acf[-1]
  lm_icc <- lme4::lmer(r ~ 1 + (1 | yearf), data = rd[, yearf := factor(year)], REML = TRUE)
  vc <- as.data.frame(lme4::VarCorr(lm_icc))
  icc <- vc$vcov[vc$grp == "yearf"] / sum(vc$vcov)
  pdisp <- sum(rd$pearson^2) / stats::df.residual(fp$model)
  ks <- suppressWarnings(stats::ks.test(rd$u, "punif"))
  S$resid <- list(rows = rd, year_means = ym, ljung_box = lb, acf = acf_all, icc = icc, pearson_dispersion = pdisp,
                  ks_p = ks$p.value, sd_r = stats::sd(rd$r), n_years = sum(is.finite(ym$mean_r)))

  # ---- R2.1 item 6: R2 of JJAS ONI on ns(year, k)
  ea <- data.table::fread(file.path(REPO_ROOT, "data", "derived", "enso_annual.csv"))
  r2 <- data.table::rbindlist(lapply(c(3L, 4L, 6L, 8L), function(k) {
    a <- ea[year %in% 2002:2025]; a1 <- ea[year %in% 2001:2024]
    data.table::data.table(k = k, n_years = nrow(a),
                           r2_lag0 = summary(stats::lm(oni_jjas ~ splines::ns(year, df = k), data = a))$r.squared,
                           r2_lag1 = summary(stats::lm(a1$oni_jjas ~ splines::ns(a$year, df = k)))$r.squared)
  }))
  S$r2 <- r2

  # ---- R2.1 item 7: episodes and leave-one-episode-out
  ep <- rev_episodes()
  S$episodes <- ep
  lo <- data.table::rbindlist(lapply(seq_len(nrow(ep)), function(i) {
    f <- rev_fit(paste0("loepo_", ep$episode_id[i]))
    if (inherits(f, "error")) stop("leave-one-episode-out fit failed: ", ep$episode_id[i], ": ", conditionMessage(f))
    a <- rev_irr(f, "oni_jjas"); b <- rev_irr(f, "oni_jjas_lag1")
    data.table::data.table(episode_id = ep$episode_id[i], label = ep$label[i], years_dropped = ep$years_dropped[i],
                           n = f$n, re_used = f$re_used, irr_lag0 = a$irr, lo_lag0 = a$lo, hi_lag0 = a$hi,
                           irr_lag1 = b$irr, lo_lag1 = b$lo, hi_lag1 = b$hi, display_lag0 = a$display, display_lag1 = b$display)
  }))
  S$loepo <- lo

  # ---- R2.3 items 2-3: trend and ABER
  sens_tab <- function(names, labels) data.table::rbindlist(lapply(seq_along(names), function(i) {
    f <- rev_fit(names[i])
    if (inherits(f, "error"))
      return(data.table::data.table(fit = names[i], label = labels[i], n = NA_integer_, re_used = NA_character_,
                                    irr_lag0 = NA_real_, lo_lag0 = NA_real_, hi_lag0 = NA_real_, irr_lag1 = NA_real_,
                                    lo_lag1 = NA_real_, hi_lag1 = NA_real_,
                                    display_lag0 = "not estimated (did not converge)",
                                    display_lag1 = "not estimated (did not converge)"))
    a <- rev_irr(f, "oni_jjas"); b <- rev_irr(f, "oni_jjas_lag1")
    data.table::data.table(fit = names[i], label = labels[i], n = f$n, re_used = f$re_used, irr_lag0 = a$irr,
                           lo_lag0 = a$lo, hi_lag0 = a$hi, irr_lag1 = b$irr, lo_lag1 = b$lo, hi_lag1 = b$hi,
                           display_lag0 = a$display, display_lag1 = b$display)
  }))
  S$trend <- sens_tab(c("trend_ns3", "trend_ns6", "trend_ns8", "trend_unit_linear"),
                      c("ns(year, 3)", "ns(year, 6)", "ns(year, 8)", "ns(year, 4) + unit linear trends (random slope)"))
  S$aber <- sens_tab(c("aber_none", "aber_lag1", "aber_ns3", "aber_positivity"),
                     c("No ABER term", "log ABER of the previous year", "ns(log ABER, 3)",
                       "Test positivity (offset log persons tested, no ABER term)"))
  ut <- rev_fit("trend_unit_linear")
  S$unit_trend_sd <- if (inherits(ut, "error")) NA_real_ else {
    vc <- glmmTMB::VarCorr(ut$model)$cond; unname(attr(vc[[length(vc)]], "stddev")[1]) }

  # ---- R2.4 item 1: nonlinear climate
  fl <- rev_fit("nonlin_linear"); fs <- rev_fit("nonlin_spline")
  for (f in list(fl, fs)) if (inherits(f, "error")) stop("climate fit failed: ", conditionMessage(f))
  if (fl$n != fs$n) stop("nonlinear LRT: different rows (", fl$n, " vs ", fs$n, ")")
  l1 <- stats::logLik(fs$model); l0 <- stats::logLik(fl$model)
  st <- as.numeric(2 * (l1 - l0)); df <- as.integer(attr(l1, "df") - attr(l0, "df"))
  S$nonlin <- list(linear = fl, spline = fs, lrt = list(stat = st, df = df, p = stats::pchisq(max(st, 0), df, lower.tail = FALSE)),
                   same_re = identical(fl$re_used, fs$re_used))
  base <- data.table::data.table(unit_id = fs$data$unit_id[1], oni_jjas = 0, oni_jjas_lag1 = 0,
                                 log_aber = stats::median(fs$data$log_aber), nfme2016 = 0, covid2020_21 = 0,
                                 year = 2013L, log_pop = 0)
  for (v in REV_CLIM_VARS) data.table::set(base, j = v, value = 0)
  S$curves <- data.table::rbindlist(lapply(REV_CLIM_VARS, function(v) {
    g <- seq(stats::quantile(fs$data[[v]], 0.01), stats::quantile(fs$data[[v]], 0.99), length.out = 101)
    rbind(cbind(model = "spline", partial_response(fs, v, g, base)),
          cbind(model = "linear", partial_response(fl, v, g, base)))
  }))
  S$curve_data <- fs$data[, c("unit_id", "year", REV_CLIM_VARS), with = FALSE]

  # ---- R2.4 item 2: Tmax tercile
  ft <- rev_fit("tmax_tercile")
  cl <- tmax_climatology()
  S$tercile <- list(fit = ft, clim = cl, wald_lag0 = if (!inherits(ft, "error")) {
    b <- paste0("tmax_", c("cool", "middle", "hot"))
    L <- rbind(c(1, -1, 0), c(0, 1, -1))
    be <- glmmTMB::fixef(ft$model)$cond[b]; V <- stats::vcov(ft$model)$cond[b, b]
    w <- as.numeric(t(L %*% be) %*% solve(L %*% V %*% t(L), L %*% be)); list(chisq = w, df = 2L, p = stats::pchisq(w, 2, lower.tail = FALSE))
  }, wald_lag1 = if (!inherits(ft, "error")) {
    b <- paste0("tmax_lag1_", c("cool", "middle", "hot"))
    L <- rbind(c(1, -1, 0), c(0, 1, -1))
    be <- glmmTMB::fixef(ft$model)$cond[b]; V <- stats::vcov(ft$model)$cond[b, b]
    w <- as.numeric(t(L %*% be) %*% solve(L %*% V %*% t(L), L %*% be)); list(chisq = w, df = 2L, p = stats::pchisq(w, 2, lower.tail = FALSE))
  })

  # ---- R2.2 item 3a: zones
  fz <- rev_fit("zone")
  if (inherits(fz, "error")) stop("zone fit failed: ", conditionMessage(fz))
  zk <- zone_key(setdiff(unique(REV_ZONES), REV_ZONE_REF))
  zones <- c("Northern", "Central", "Eastern", "Western", "Southern", "North-Eastern")
  zt <- data.table::rbindlist(lapply(zones, function(z) {
    k <- zone_key(z)
    t0 <- if (z == REV_ZONE_REF) "oni_jjas" else c("oni_jjas", paste0("oni_x_", k))
    t1 <- if (z == REV_ZONE_REF) "oni_jjas_lag1" else c("oni_jjas_lag1", paste0("oni_lag1_x_", k))
    a <- rev_lincomb_irr(fz, t0); b <- rev_lincomb_irr(fz, t1)
    data.table::data.table(zone = z, n_units = sum(REV_ZONES == z), units = paste(names(REV_ZONES)[REV_ZONES == z], collapse = ";"),
                           irr_lag0 = a$irr, lo_lag0 = a$lo, hi_lag0 = a$hi, irr_lag1 = b$irr, lo_lag1 = b$lo, hi_lag1 = b$hi)
  }))
  S$zone <- list(fit = fz, table = zt, wald = rev_wald(fz, c(paste0("oni_x_", zk), paste0("oni_lag1_x_", zk))),
                 wald_lag0 = rev_wald(fz, paste0("oni_x_", zk)), wald_lag1 = rev_wald(fz, paste0("oni_lag1_x_", zk)))

  # ---- R2.2 item 3b: unpooled unit fits
  S$unpooled <- rev_cached("unpooled", data.table::rbindlist(lapply(split(mal, mal$unit_id), unpooled_unit_fit)))

  # ---- R2.6 scenarios
  e <- outlook$enso
  vals <- c(0.5, 1.0, 1.5, e$oni)
  sp <- rev_cached("scen_primary", rev_scenarios_primary(fp, vals))
  # the observed-value scenario must reproduce the published outlook
  for (y in 2026:2027) {
    a <- sp[oni == e$oni & year == y]; b <- outlook$oni$proj$total[year == y]
    if (abs(a$excess_mean - b$excess_mean) > 1e-6 || abs(a$excess_lo - b$excess_lo) > 1e-6)
      stop("scenario at the observed ONI does not reproduce the outlook for ", y)
  }
  yfit <- S$year_fits$year_iid
  sy <- if (inherits(yfit, "error")) NULL else rev_cached("scen_year_iid", rev_scenarios_year(yfit, vals))
  S$scen <- data.table::rbindlist(list(sp, sy), fill = TRUE)
  S$oni_obs <- e$oni

  # ---- R2.5 early warning
  ew <- list(); yv <- list()
  for (j in names(ew_validation_jobs())) { ew[[j]] <- rev_read(paste0("ew_", j)); yv[[j]] <- rev_read(paste0("ewyv_", j)) }
  S$ew <- ew; S$ewyv <- yv
  S
}

# ------------------------------------------------------------------ CSV outputs

rev_write_csvs <- function(S) {
  w <- function(x, f) data.table::fwrite(x, file.path(A_OUT, paste0("rev_", f, ".csv")))
  w(S$inference, "inference")
  w(S$boot[order(rep)], "bootstrap_replicates")
  w(S$episodes, "episodes")
  w(S$loepo, "leave_one_episode_out")
  w(S$r2, "oni_trend_r2")
  rs <- S$resid
  w(rs$year_means, "residual_year_means")
  w(data.table::data.table(statistic = c("ljung_box_q_lag1", "ljung_box_q_lag2", "ljung_box_q_lag3", "ljung_box_p_lag1",
                                         "ljung_box_p_lag2", "ljung_box_p_lag3", "acf_lag1", "acf_lag2", "acf_lag3",
                                         "icc_year", "pearson_dispersion", "ks_uniform_p", "sd_quantile_residuals",
                                         "n_years"),
                           value = c(vapply(rs$ljung_box, `[[`, numeric(1), "statistic"),
                                     vapply(rs$ljung_box, `[[`, numeric(1), "p"), rs$acf[1:3], rs$icc,
                                     rs$pearson_dispersion, rs$ks_p, rs$sd_r, rs$n_years)), "residual_diagnostics")
  w(rbind(cbind(family = "trend", S$trend), cbind(family = "aber", S$aber)), "trend_aber")
  fl <- S$nonlin$linear; fs <- S$nonlin$spline
  nl <- data.table::data.table(model = c("linear climate", "spline climate (ns, 3 df)"), n = c(fl$n, fs$n),
                               re_used = c(fl$re_used, fs$re_used), logLik = c(as.numeric(stats::logLik(fl$model)), as.numeric(stats::logLik(fs$model))),
                               irr_lag0 = c(rev_irr(fl, "oni_jjas")$irr, rev_irr(fs, "oni_jjas")$irr),
                               lo_lag0 = c(rev_irr(fl, "oni_jjas")$lo, rev_irr(fs, "oni_jjas")$lo),
                               hi_lag0 = c(rev_irr(fl, "oni_jjas")$hi, rev_irr(fs, "oni_jjas")$hi),
                               irr_lag1 = c(rev_irr(fl, "oni_jjas_lag1")$irr, rev_irr(fs, "oni_jjas_lag1")$irr),
                               lo_lag1 = c(rev_irr(fl, "oni_jjas_lag1")$lo, rev_irr(fs, "oni_jjas_lag1")$lo),
                               hi_lag1 = c(rev_irr(fl, "oni_jjas_lag1")$hi, rev_irr(fs, "oni_jjas_lag1")$hi),
                               lrt_chisq = S$nonlin$lrt$stat, lrt_df = S$nonlin$lrt$df, lrt_p = S$nonlin$lrt$p)
  w(nl, "nonlinear_climate")
  w(S$curves, "nonlinear_curves")
  tc <- S$tercile
  if (!inherits(tc$fit, "error")) {
    tt <- irr_table(tc$fit, c(paste0("tmax_", c("cool", "middle", "hot")), paste0("tmax_lag1_", c("cool", "middle", "hot"))))
    tt[, `:=`(tercile = rep(c("cool", "middle", "hot"), 2), lag = rep(0:1, each = 3), n = tc$fit$n, re_used = tc$fit$re_used,
              wald_p_equal = rep(c(tc$wald_lag0$p, tc$wald_lag1$p), each = 3))]
    w(tt, "tmax_tercile")
  }
  w(tc$clim$table[, .(unit_id, unit_name = UNIT_NAMES[unit_id], tmax_clim_1981_2010 = tmax_clim, n_years, tercile)], "tmax_climatology")
  z <- S$zone
  w(cbind(z$table, wald_p_joint = z$wald$p, wald_df_joint = z$wald$df, wald_p_lag0 = z$wald_lag0$p,
          wald_p_lag1 = z$wald_lag1$p, n = z$fit$n, re_used = z$fit$re_used), "zone_interaction")
  up <- data.table::copy(S$unpooled)
  data.table::setcolorder(up, c("unit_id", "n", "irr_lag0", "lo", "hi", "log_irr_lag0", "se"))
  w(up, "unit_irr_unpooled")
  w(S$scen, "scenarios")
  ew <- data.table::rbindlist(lapply(names(S$ew), function(j) {
    s <- data.table::copy(S$ew[[j]]$summary)
    yv <- S$ewyv[[j]]$summary
    s[yv, `:=`(coverage95_yearvar = i.coverage_yearvar, n_yearvar = i.n), on = "model"]
    cbind(job = j, s, n_folds_scored = data.table::uniqueN(S$ew[[j]]$per_split[scored == TRUE]$split_id))
  }), fill = TRUE)
  w(ew, "early_warning_summary")
  w(data.table::rbindlist(lapply(names(S$ew), function(j) cbind(job = j, S$ew[[j]]$split_metrics))), "early_warning_per_split")
  w(data.table::rbindlist(lapply(names(S$ewyv), function(j) cbind(job = j, S$ewyv[[j]]$info))), "early_warning_yearvar_fits")
}

# ------------------------------------------------------------------ results keys

build_revision_results <- function(S) {
  R <- list()
  add <- function(key, value, display, description) {
    stopifnot(grepl("^A\\.rev\\.", key))
    if (key %in% names(R)) stop("duplicate key ", key)
    if (length(value) != 1L) value <- paste(value, collapse = "; ")
    R[[key]] <<- res_entry(value, display, REV_SCRIPT, description)
  }
  fp <- S$primary
  rev_small <- function(x) if (is.na(x)) "not estimated" else if (x < 0.0005) "<0.001" else fmt_num(x, 3)
  irr_keys <- function(prefix, row, what) {
    add(paste0(prefix, ".irr_lag0"), row$irr_lag0, row$display_lag0, paste0("ONI lag-0 IRR per 1 C, ", what))
    add(paste0(prefix, ".irr_lag1"), row$irr_lag1, row$display_lag1, paste0("ONI lag-1 IRR per 1 C, ", what))
  }
  rng_keys <- function(prefix, tab, what) {
    for (lg in c("lag0", "lag1")) {
      v <- tab[[paste0("irr_", lg)]]; v <- v[is.finite(v)]
      add(sprintf("%s.irr_%s_range", prefix, lg), max(v) - min(v), rev_range_disp(v),
          sprintf("Range of ONI %s IRR point estimates across %s", sub("lag", "lag-", lg), what))
    }
  }

  # ---- R2.1 inference
  I <- S$inference
  desc <- c(year_iid = "primary model + iid year random effect (1 | year), Wald 95% CI",
            year_ar1 = "primary model + AR(1) year random effect shared by all units, Wald 95% CI",
            glm_model = "NB GLM with unit fixed effects and the primary fixed terms, model-based Wald 95% CI",
            cluster_year = "NB GLM with unit fixed effects, cluster-robust by year (vcovCL HC1), t(G - 1) reference",
            cluster_twoway = "NB GLM with unit fixed effects, two-way cluster-robust by year and unit (vcovCL HC1), t(min(G) - 1) reference",
            boot = sprintf("primary model point estimate with percentile 95%% CI from %d moving-block (3-year) year bootstrap resamples, seed %d",
                           REV_B, REV_SEED))
  for (v in names(desc)) irr_keys(paste0("A.rev.temporal.", v), I[variant == v], desc[[v]])
  fy <- S$year_fits$year_iid; fa <- S$year_fits$year_ar1
  if (!inherits(fy, "error")) {
    s <- rev_group_sd(fy, "yearf")
    add("A.rev.temporal.year_iid.sd_year", s, rev_small(s), "SD of the iid year random effect (log-rate scale)")
  }
  if (!inherits(fa, "error")) {
    s <- rev_group_sd(fa, "one"); rho <- rev_ar1_rho(fa)
    add("A.rev.temporal.year_ar1.sd_year", s, rev_small(s), "Marginal SD of the AR(1) year random effect (log-rate scale)")
    add("A.rev.temporal.year_ar1.rho", rho, if (is.na(rho)) "not estimated" else fmt_num(rho, 2), "AR(1) correlation of the year effect")
  }
  add("A.rev.temporal.cluster.n_years", S$glm$G_year, as.character(S$glm$G_year), "Number of year clusters G (years with malaria data)")
  bt <- S$boot
  add("A.rev.temporal.boot.n_failed", sum(!bt$ok), as.character(sum(!bt$ok)),
      sprintf("Bootstrap refits that failed to converge (of %d)", REV_B))
  add("A.rev.temporal.boot.n_ok", sum(bt$ok), fmt_int(sum(bt$ok)), sprintf("Bootstrap refits used (of %d)", REV_B))
  add("A.rev.temporal.boot.n_int_fallback", sum(bt$re_used == "int", na.rm = TRUE),
      as.character(sum(bt$re_used == "int", na.rm = TRUE)), "Bootstrap refits that converged only with random intercepts")
  Iv <- I[variant %in% c("year_iid", "year_ar1", "cluster_year", "boot")]
  rob <- all(is.finite(Iv$lo_lag0)) && all(Iv$lo_lag0 > 1 | Iv$hi_lag0 < 1)
  verdict <- if (rob) "robust to shared temporal dependence" else "not robust to shared temporal dependence"
  add("A.rev.temporal.lag0_verdict", verdict, verdict,
      paste0("Plan decision rule: robust only if the lag-0 95% interval excludes 1 under the iid year effect, AR(1) year effect, ",
             "cluster-robust (year) and year-block bootstrap inference; excluding 1 in: ",
             paste(Iv$variant[Iv$lo_lag0 > 1 | Iv$hi_lag0 < 1], collapse = ", ")))
  n_ex <- sum(Iv$lo_lag0 > 1 | Iv$hi_lag0 < 1, na.rm = TRUE)
  add("A.rev.temporal.lag0_n_excluding_1", n_ex, sprintf("%d of %d", n_ex, nrow(Iv)),
      "Inference variants of the decision rule whose lag-0 95% interval excludes 1")
  rng_keys("A.rev.temporal", I[variant != "primary"], "the R2.1 inference variants (year effects, NB GLM, bootstrap)")
  for (lg in c("lag0", "lag1")) {
    v <- I[variant != "primary"][[paste0("lo_", lg)]]; v <- v[is.finite(v)]
    add(sprintf("A.rev.temporal.lower_%s_range", lg), max(v) - min(v), rev_range_disp(v),
        sprintf("Range of the lower 95%% limits of the ONI %s IRR across the R2.1 inference variants", sub("lag", "lag-", lg)))
  }
  # residuals
  rs <- S$resid
  for (k in 1:3) {
    p <- rs$ljung_box[[k]]$p
    add(sprintf("A.rev.resid.ljung_box_p_lag%d", k), p, disp_p(p),
        sprintf("Ljung-Box test up to lag %d of the year-mean randomized quantile residual series, primary model (Q = %.2f, %d years)",
                k, rs$ljung_box[[k]]$statistic, rs$n_years))
  }
  add("A.rev.resid.acf_lag1", rs$acf[1], fmt_num(rs$acf[1], 2), "Lag-1 autocorrelation of the year-mean residual series")
  add("A.rev.resid.icc_year", rs$icc, fmt_num(rs$icc, 3),
      "Intraclass correlation of the randomized quantile residuals by year (lmer, REML)")
  add("A.rev.resid.pearson_dispersion", rs$pearson_dispersion, fmt_num(rs$pearson_dispersion, 2),
      "Pearson chi-square / residual df, primary model")
  add("A.rev.resid.ks_uniform_p", rs$ks_p, disp_p(rs$ks_p), "Kolmogorov-Smirnov p, PIT values vs uniform (randomized, seed 20261117)")
  for (i in seq_len(nrow(S$r2))) {
    r <- S$r2[i]
    add(sprintf("A.rev.trend_r2.k%d", r$k), r$r2_lag0, fmt_num(r$r2_lag0, 2),
        sprintf("R-squared of JJAS ONI (2002-2025, %d years) regressed on ns(year, %d)", r$n_years, r$k))
  }
  # episodes
  ep <- S$episodes; lo <- S$loepo
  add("A.rev.episodes.n", nrow(ep), as.character(nrow(ep)),
      "Moderate or stronger ENSO episodes (NOAA rule, peak |ONI| >= 1.0 C) removing at least one 2002-2025 year")
  add("A.rev.episodes.list", paste(ep$episode_id, collapse = ";"),
      join_and(sprintf("%s (peak %s)", ep$label, signed2(ep$peak))), "Episodes used for leave-one-episode-out (rev_episodes.csv)")
  for (i in seq_len(nrow(lo))) {
    r <- lo[i]; k <- paste0("A.rev.loepo.", r$episode_id)
    irr_keys(k, r, sprintf("primary model refitted without %s (years %s dropped)", r$label, gsub(";", ", ", r$years_dropped)))
    yrs <- episode_years(r$years_dropped)
    add(paste0(k, ".years_dropped"), r$years_dropped, sprintf("%d%s%d", min(yrs), EN_DASH, max(yrs)), paste0("Years dropped, ", r$label))
    add(paste0(k, ".n"), r$n, fmt_int(r$n), paste0("Unit-years, refit without ", r$label))
  }
  rng_keys("A.rev.loepo", lo, "leave-one-episode-out refits")
  m0 <- lo[which.min(irr_lag0)]; m1 <- lo[which.max(irr_lag0)]
  add("A.rev.loepo.irr_lag0_range_text", max(lo$irr_lag0) - min(lo$irr_lag0),
      sprintf("%s (without %s) to %s (without %s)", fmt_num(m0$irr_lag0, 2), m0$label, fmt_num(m1$irr_lag0, 2), m1$label),
      "Range of lag-0 IRRs across leave-one-episode-out refits, with the episode left out")
  nex <- sum(lo$lo_lag0 > 1)
  add("A.rev.loepo.n_lag0_excluding_1", nex, sprintf("%d of %d", nex, nrow(lo)),
      "Leave-one-episode-out refits whose lag-0 95% CI lies above 1")

  # ---- R2.2 item 3
  z <- S$zone
  for (i in seq_len(nrow(z$table))) {
    r <- z$table[i]; k <- paste0("A.rev.zone.", zone_key(r$zone))
    add(paste0(k, ".irr_lag0"), r$irr_lag0, disp_irr(list(irr = r$irr_lag0, lo = r$lo_lag0, hi = r$hi_lag0)),
        sprintf("Zone-specific ONI lag-0 IRR, %s zone (%d units), ONI x zone model", r$zone, r$n_units))
    add(paste0(k, ".irr_lag1"), r$irr_lag1, disp_irr(list(irr = r$irr_lag1, lo = r$lo_lag1, hi = r$hi_lag1)),
        sprintf("Zone-specific ONI lag-1 IRR, %s zone, ONI x zone model", r$zone))
  }
  add("A.rev.zone.wald_p", z$wald$p, disp_p(z$wald$p),
      sprintf("Wald test of the ONI x zone interactions at both lags (chi-square %.2f on %d df); structure %s", z$wald$chisq, z$wald$df, z$fit$re_used))
  add("A.rev.zone.wald_p_lag0", z$wald_lag0$p, disp_p(z$wald_lag0$p),
      sprintf("Wald test of the lag-0 ONI x zone interactions (chi-square %.2f on %d df)", z$wald_lag0$chisq, z$wald_lag0$df))
  add("A.rev.zone.wald_p_lag1", z$wald_lag1$p, disp_p(z$wald_lag1$p),
      sprintf("Wald test of the lag-1 ONI x zone interactions (chi-square %.2f on %d df)", z$wald_lag1$chisq, z$wald_lag1$df))
  add("A.rev.zone.wald_chisq", z$wald$chisq, sprintf("%s on %d df", fmt_num(z$wald$chisq, 2), z$wald$df), "Joint Wald chi-square, ONI x zone")
  rng_keys("A.rev.zone", z$table, "the six zones (ONI x zone model)")
  up <- S$unpooled
  okk <- up[status == "ok"]
  add("A.rev.unpooled.n_units_fitted", nrow(okk), as.character(nrow(okk)), "Units with a converged unit-specific NB GLM")
  add("A.rev.unpooled.n_units_dropped", nrow(up) - nrow(okk), as.character(nrow(up) - nrow(okk)),
      sprintf("Units without a unit-specific estimate (fewer than %d complete years or no convergence)", REV_MIN_YEARS))
  dr <- up[status != "ok"]
  add("A.rev.unpooled.units_dropped", if (nrow(dr)) paste(dr$unit_id, collapse = ";") else "none",
      if (nrow(dr)) join_and(UNIT_NAMES[dr$unit_id]) else "none", paste("Units dropped:", paste(dr$unit_id, dr$status, sep = ": ", collapse = " | ")))
  add("A.rev.unpooled.irr_lag0_range", max(okk$irr_lag0) - min(okk$irr_lag0), rev_range_disp(okk$irr_lag0),
      "Range of unpooled unit-specific lag-0 IRRs")
  add("A.rev.unpooled.irr_lag0_median", stats::median(okk$irr_lag0), fmt_num(stats::median(okk$irr_lag0), 2),
      "Median of unpooled unit-specific lag-0 IRRs")

  # ---- R2.3 items 2-3
  for (i in seq_len(nrow(S$trend))) {
    r <- S$trend[i]; k <- paste0("A.rev.trend.", sub("^trend_", "", r$fit))
    irr_keys(k, r, paste0("primary model with trend ", r$label, "; structure ", r$re_used))
    add(paste0(k, ".n"), r$n, fmt_int(r$n), paste0("Unit-years, trend ", r$label))
  }
  rng_keys("A.rev.trend", S$trend, "the alternative trend specifications")
  if (is.finite(S$unit_trend_sd))
    add("A.rev.trend.unit_linear.sd_slope", S$unit_trend_sd, fmt_num(S$unit_trend_sd, 3),
        "SD of the unit random slope on centred year (log-rate per decade)")
  for (i in seq_len(nrow(S$aber))) {
    r <- S$aber[i]; k <- paste0("A.rev.aber.", sub("^aber_", "", r$fit))
    irr_keys(k, r, paste0("diagnostic-effort variant: ", r$label, "; structure ", r$re_used))
    add(paste0(k, ".n"), r$n, fmt_int(r$n), paste0("Unit-years, ", r$label))
  }
  rng_keys("A.rev.aber", S$aber, "the alternative diagnostic-effort adjustments")

  # ---- R2.4
  nl <- S$nonlin
  for (m in c("linear", "spline")) {
    f <- nl[[m]]; a <- rev_irr(f, "oni_jjas"); b <- rev_irr(f, "oni_jjas_lag1")
    w <- if (m == "linear") "primary + linear March-June Tmax anomaly and JJAS rainfall departure (lags 0 and 1)"
         else "primary + natural splines (3 df) of March-June Tmax anomaly and JJAS rainfall departure (lags 0 and 1)"
    add(sprintf("A.rev.nonlin.%s.irr_lag0", m), a$irr, a$display, paste0("ONI lag-0 IRR, ", w, "; structure ", f$re_used))
    add(sprintf("A.rev.nonlin.%s.irr_lag1", m), b$irr, b$display, paste0("ONI lag-1 IRR, ", w, "; structure ", f$re_used))
  }
  add("A.rev.nonlin.n", nl$spline$n, fmt_int(nl$spline$n), "Unit-years in both climate models (AN and LD rainfall excluded)")
  add("A.rev.nonlin.lrt_p", nl$lrt$p, disp_p(nl$lrt$p),
      sprintf("LRT spline vs linear climate (chi-square %.2f on %d df; same rows; same RE structure: %s)", nl$lrt$stat, nl$lrt$df, nl$same_re))
  add("A.rev.nonlin.lrt_chisq", nl$lrt$stat, sprintf("%s on %d df", fmt_num(nl$lrt$stat, 2), nl$lrt$df), "LRT statistic, spline vs linear climate")
  tc <- S$tercile
  if (!inherits(tc$fit, "error")) {
    tt <- irr_table(tc$fit, c(paste0("tmax_", c("cool", "middle", "hot")), paste0("tmax_lag1_", c("cool", "middle", "hot"))))
    for (i in seq_len(nrow(tt))) {
      t <- sub("^tmax_(lag1_)?", "", tt$term[i]); lg <- if (grepl("lag1", tt$term[i])) "lag1" else "lag0"
      add(sprintf("A.rev.tmax_tercile.%s.irr_%s", t, lg), tt$irr[i], disp_irr(tt[i]),
          sprintf("IRR per 1 C March-June Tmax anomaly (%s), units in the %s tercile of 1981-2010 March-June Tmax; structure %s",
                  sub("lag", "lag ", lg), t, tc$fit$re_used))
    }
    add("A.rev.tmax_tercile.wald_p_lag0", tc$wald_lag0$p, disp_p(tc$wald_lag0$p),
        sprintf("Wald test of equal lag-0 Tmax slopes across terciles (chi-square %.2f on 2 df)", tc$wald_lag0$chisq))
    add("A.rev.tmax_tercile.wald_p_lag1", tc$wald_lag1$p, disp_p(tc$wald_lag1$p),
        sprintf("Wald test of equal lag-1 Tmax slopes across terciles (chi-square %.2f on 2 df)", tc$wald_lag1$chisq))
    add("A.rev.tmax_tercile.cutpoints", paste(round(tc$clim$cut, 3), collapse = ";"),
        sprintf("%s and %s \u00b0C", fmt_num(tc$clim$cut[1], 1), fmt_num(tc$clim$cut[2], 1)),
        "Tercile cut points of unit 1981-2010 March-June mean Tmax")
  }

  # ---- R2.5 early warning
  for (j in names(S$ew)) {
    o <- if (startsWith(j, "mal")) "malaria" else "dengue"; sch <- sub("^[a-z]+_", "", j)
    pre <- sprintf("A.rev.ew.%s.%s", o, sch)
    cm <- S$ew[[j]]; s <- cm$summary; yv <- S$ewyv[[j]]$summary
    for (m in c("A", "B", "C")) {
      r <- s[model == m]; fl <- !is.na(r$flag)
      add(sprintf("%s.%s.male", pre, m), r$male, if (fl) FLAG_TEXT else fmt_num(r$male, 3), paste0("Early-warning MALE, Model ", m, "_ew"))
      add(sprintf("%s.%s.logscore", pre, m), r$log_score, if (fl) FLAG_TEXT else fmt_num(r$log_score, 3),
          paste0("Early-warning NB log score, Model ", m, "_ew"))
      add(sprintf("%s.%s.coverage", pre, m), r$coverage95, if (fl) FLAG_TEXT else fmt_pct(r$coverage95, 1),
          sprintf("Coverage of nominal 95%% prediction intervals, Model %s_ew (n_draws %d, seed %d)", m, N_DRAWS, cm$seed))
      cy <- yv[model == m]
      add(sprintf("%s.%s.coverage_yearvar", pre, m), if (nrow(cy)) cy$coverage_yearvar else NA_real_,
          if (nrow(cy) && !fl) fmt_pct(cy$coverage_yearvar, 1) else FLAG_TEXT,
          sprintf("Coverage of nominal 95%% intervals adding the year-level variance (iid year-effect refit per split), Model %s_ew; %d rows",
                  m, if (nrow(cy)) cy$n else 0L))
    }
    ratio <- s[model == "C"]$male / s[model == "B"]$male
    add(paste0(pre, ".c_vs_b.male_ratio"), ratio, if (is.na(ratio)) FLAG_TEXT else fmt_num(ratio, 3), "Early-warning MALE ratio C_ew / B_ew, same rows")
    ratio_ba <- s[model == "B"]$male / s[model == "A"]$male
    add(paste0(pre, ".b_vs_a.male_ratio"), ratio_ba, if (is.na(ratio_ba)) FLAG_TEXT else fmt_num(ratio_ba, 3), "Early-warning MALE ratio B_ew / A_ew, same rows")
    dl <- s[model == "C"]$log_score - s[model == "B"]$log_score
    add(paste0(pre, ".c_vs_b.logscore_diff"), dl, if (is.na(dl)) FLAG_TEXT else fmt_num(dl, 3), "Early-warning mean log score, C_ew minus B_ew")
    nf <- data.table::uniqueN(cm$per_split[scored == TRUE]$split_id)
    add(paste0(pre, ".n_splits"), nf, as.character(nf), "Folds (LOEO) or origins (rolling) with scored rows")
    add(paste0(pre, ".common_n"), cm$common_n, fmt_int(cm$common_n), "Scored unit-years common to A_ew, B_ew, C_ew")
    if (sch == "loeo") {
      wv <- data.table::dcast(cm$split_metrics[model %in% c("B", "C")], split_id ~ model, value.var = "male")
      kb <- sum(wv$C < wv$B)
      add(paste0(pre, ".c_vs_b.n_folds_better"), kb, sprintf("%d of %d", kb, nrow(wv)), "LOEO folds with C_ew MALE < B_ew MALE")
    }
  }

  # ---- R2.6 scenarios
  sc <- S$scen
  add("A.rev.scen.oni_observed", S$oni_obs, signed2(S$oni_obs), "Observed JJAS ONI 2026 used as the fourth scenario")
  lab <- function(v) if (isTRUE(all.equal(v, S$oni_obs))) "oni_obs" else sprintf("oni%03d", round(100 * v))
  for (i in seq_len(nrow(sc))) {
    r <- sc[i]
    add(sprintf("A.rev.scen.%s.%s.excess_%d", r$model, lab(r$oni), r$year), r$excess_mean,
        fmt_ci(r$excess_mean, r$excess_lo, r$excess_hi, digits = 0, label = "nominal 95% interval", big.mark = ","),
        sprintf("All-unit expected excess malaria cases %d, JJAS ONI %s %s vs ONI = 0 (%s), %s model; mean and 2.5-97.5 percentiles of %d joint draws (seed %d)%s",
                r$year, if (r$year == 2026) "2026 (lag 0)" else "2026 (lag 1 in 2027)", formatC(r$oni, format = "f", digits = 4),
                if (r$year == 2026) "lag-0 term" else "lag-1 term", if (r$model == "primary") "primary" else "iid year-effect",
                N_DRAWS, SEED + 7L, if (r$model == "year_iid") "; new-year effect drawn from N(0, sd_year^2), shared by units and scenarios" else ""))
  }
  R
}

# ------------------------------------------------------------------ figures

rev_save_fig <- function(p, name, width_mm = 174, height_mm = 120) {
  tif <- file.path(REV_FIG_DIR, paste0(name, ".tif")); png <- file.path(REV_FIG_DIR, paste0(name, ".png"))
  grDevices::tiff(tif, width = width_mm, height = height_mm, units = "mm", res = 600, compression = "lzw",
                  type = "cairo", bg = "white")
  print(p); grDevices::dev.off()
  grDevices::png(png, width = width_mm, height = height_mm, units = "mm", res = 200, type = "cairo", bg = "white")
  print(p); grDevices::dev.off()
  invisible(c(tif, png))
}

rev_theme <- function() ggplot2::theme_bw(base_size = 8) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), strip.background = ggplot2::element_rect(fill = "grey92"),
                 plot.title = ggplot2::element_text(size = 8, face = "bold"))

rev_figures <- function(S) {
  library(ggplot2)
  # ---- S2 residual diagnostics
  rd <- S$resid$rows; ym <- S$resid$year_means
  n <- nrow(rd)
  qq <- data.table::data.table(expected = (seq_len(n) - 0.5) / n, observed = sort(rd$u))
  pa <- ggplot(qq, aes(expected, observed)) + geom_abline(linetype = 2, colour = "grey50") + geom_point(size = 0.5) +
    labs(x = "Expected uniform quantile", y = "Observed PIT value", title = "a  Uniform Q-Q plot") + coord_equal() + rev_theme()
  pb <- ggplot(rd, aes(log(mu), r)) + geom_hline(yintercept = 0, colour = "grey50") + geom_point(size = 0.5, alpha = 0.5) +
    labs(x = "log fitted cases", y = "Quantile residual", title = "b  Residuals vs fitted") + rev_theme()
  pc <- ggplot(ym[is.finite(mean_r)], aes(year, mean_r)) + geom_hline(yintercept = 0, colour = "grey50") +
    geom_line() + geom_point(size = 1) +
    labs(x = "Year", y = "Mean quantile residual", title = "c  Year-mean residual series") + rev_theme()
  ac <- data.table::data.table(lag = seq_along(S$resid$acf), acf = S$resid$acf)
  ci <- 1.96 / sqrt(S$resid$n_years)
  lb <- S$resid$ljung_box
  pd <- ggplot(ac, aes(lag, acf)) + geom_hline(yintercept = 0) + geom_segment(aes(xend = lag, yend = 0)) +
    geom_hline(yintercept = c(-ci, ci), linetype = 2, colour = "grey40") +
    scale_x_continuous(breaks = ac$lag) +
    labs(x = "Lag (years)", y = "Autocorrelation", title = "d  ACF of the year-mean residuals",
         subtitle = sprintf("Ljung\u2013Box p = %s, %s, %s (lags 1\u20133); ICC = %s",
                         fmt_p(lb[[1]]$p), fmt_p(lb[[2]]$p), fmt_p(lb[[3]]$p), fmt_num(S$resid$icc, 3))) + rev_theme()
  rev_save_fig(cowplot::plot_grid(pa, pb, pc, pd, ncol = 2), "Figure_S2", 174, 150)

  # ---- S3 leave-one-episode-out forest
  lo <- S$loepo
  I <- S$inference[variant == "primary"]
  fr <- rbind(data.table::data.table(label = "All years (primary model)", lag = "Lag 0", irr = I$irr_lag0, lo = I$lo_lag0, hi = I$hi_lag0),
              data.table::data.table(label = "All years (primary model)", lag = "Lag 1", irr = I$irr_lag1, lo = I$lo_lag1, hi = I$hi_lag1),
              lo[, .(label = paste0("Without ", label, " (", gsub(";", ", ", years_dropped), ")"), lag = "Lag 0", irr = irr_lag0, lo = lo_lag0, hi = hi_lag0)],
              lo[, .(label = paste0("Without ", label, " (", gsub(";", ", ", years_dropped), ")"), lag = "Lag 1", irr = irr_lag1, lo = lo_lag1, hi = hi_lag1)])
  fr[, label := factor(label, levels = rev(unique(label)))]
  fr[, primary := grepl("^All years", label)]
  p3 <- ggplot(fr, aes(irr, label)) + geom_vline(xintercept = 1, linetype = 2, colour = "grey40") +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.25, orientation = "y") + geom_point(aes(shape = primary), size = 1.8) +
    scale_shape_manual(values = c(`TRUE` = 15, `FALSE` = 16), guide = "none") +
    facet_wrap(~ lag) + scale_x_log10(breaks = c(0.9, 1, 1.1, 1.2, 1.3)) +
    labs(x = "IRR per 1 \u00b0C JJAS ONI (log scale, 95% CI)", y = NULL) + rev_theme()
  p3 <- p3 + theme(panel.spacing = grid::unit(1.2, "lines"))
  rev_save_fig(p3, "Figure_S3", 174, 110)

  # ---- S4 nonlinear partial responses
  cv <- data.table::copy(S$curves)
  labs_v <- c(tmax_mamj_anom = "March\u2013June Tmax anomaly (\u00b0C), lag 0", tmax_mamj_anom_lag1 = "March\u2013June Tmax anomaly (\u00b0C), lag 1",
              rain_jjas_dep_pct = "JJAS rainfall departure (%), lag 0", rain_jjas_dep_pct_lag1 = "JJAS rainfall departure (%), lag 1")
  cv[, panel := factor(labs_v[variable], levels = labs_v)]
  cv[, model := factor(model, levels = c("spline", "linear"), labels = c("Natural spline (3 df)", "Linear"))]
  rug <- data.table::rbindlist(lapply(names(labs_v), function(v) {
    x <- S$curve_data[[v]]; q <- stats::quantile(x, c(0.01, 0.99))
    data.table::data.table(panel = factor(labs_v[[v]], levels = labs_v), x = x[x >= q[1] & x <= q[2]]) }))
  p4 <- ggplot(cv, aes(x, rr)) + geom_hline(yintercept = 1, colour = "grey50") +
    geom_ribbon(data = cv[model == "Natural spline (3 df)"], aes(ymin = lo, ymax = hi), fill = "grey80") +
    geom_line(aes(linetype = model)) + geom_rug(data = rug, aes(x = x), inherit.aes = FALSE, alpha = 0.2, length = grid::unit(0.02, "npc")) +
    facet_wrap(~ panel, scales = "free_x") + scale_y_log10() +
    scale_linetype_manual(values = c("Natural spline (3 df)" = 1, "Linear" = 2), name = NULL) +
    labs(x = NULL, y = "Rate ratio relative to 0 (log scale)",
         title = sprintf("Likelihood-ratio test, splines vs linear: %s", disp_p(S$nonlin$lrt$p))) +
    rev_theme() + theme(legend.position = "bottom")
  rev_save_fig(p4, "Figure_S4", 174, 140)

  # ---- S5 conditional outlook across ONI
  sc <- data.table::copy(S$scen)
  sc[, model := factor(model, levels = c("primary", "year_iid"), labels = c("Primary model", "With iid year effect"))]
  sc[, yr := factor(year, labels = c("2026 (ONI at lag 0)", "2027 (2026 ONI at lag 1)"))]
  p5 <- ggplot(sc, aes(oni, excess_mean / 1000, shape = model)) + geom_hline(yintercept = 0, colour = "grey50") +
    geom_errorbar(aes(ymin = excess_lo / 1000, ymax = excess_hi / 1000), width = 0.06, position = position_dodge(width = 0.12)) +
    geom_point(size = 1.8, position = position_dodge(width = 0.12)) + facet_wrap(~ yr) +
    scale_shape_manual(values = c(16, 1), name = NULL) +
    scale_x_continuous(breaks = sort(unique(sc$oni)), labels = function(x) formatC(x, format = "f", digits = 2)) +
    labs(x = "2026 JJAS ONI scenario (\u00b0C)", y = "Expected excess malaria cases vs ONI = 0 (thousands)") +
    rev_theme() + theme(legend.position = "bottom")
  rev_save_fig(p5, "Figure_S5", 174, 100)
}

# ------------------------------------------------------------------ results step

run_rev_results <- function() {
  fits <- readRDS(file.path(A_CACHE, "fits.rds"))
  cmp <- run_validation(reuse = TRUE)
  med <- run_mediation(reuse = TRUE)
  outlook <- readRDS(file.path(A_CACHE, "outlook.rds"))
  S <- rev_summaries(fits, outlook)
  rev_write_csvs(S)
  alog("wrote models/A/output/rev_*.csv", file = "rev_results.log")
  rev_figures(S)
  alog("wrote Figures S2-S5", file = "rev_results.log")
  R0 <- build_results(fits, cmp, med, outlook)
  R1 <- build_revision_results(S)
  if (length(intersect(names(R0), names(R1)))) stop("revision keys collide with original keys")
  write_results(c(R0, R1), file.path(PAPER_A, "results.json"))
  alog(sprintf("wrote papers/A/results.json (%d original + %d revision keys)", length(R0), length(R1)), file = "rev_results.log")
  saveRDS(S, file.path(A_CACHE, "rev_summaries.rds"))
  invisible(list(S = S, R0 = R0, R1 = R1))
}

# ------------------------------------------------------------------ main (only when run as a script)

if (.rev_main) {
  step <- sub("^--step=", "", grep("^--step=", .rev_args, value = TRUE))
  cores <- sub("^--cores=", "", grep("^--cores=", .rev_args, value = TRUE))
  cores <- if (length(cores)) as.integer(cores) else 4L
  want <- function(s) !length(step) || s %in% step
  t0 <- Sys.time()
  alog("Paper A revision: steps ", if (length(step)) paste(step, collapse = ",") else "all", "; cores ", cores,
       "; fresh ", REV_FRESH, file = "rev_run.log")
  if (want("fits")) run_rev_fits(cores)
  if (want("boot")) run_rev_boot(cores)
  if (want("ew")) run_rev_ew(cores)
  if (want("results")) run_rev_results()
  alog(sprintf("Paper A revision finished (%.1f min)", as.numeric(difftime(Sys.time(), t0, units = "mins"))),
       file = "rev_run.log")
}
