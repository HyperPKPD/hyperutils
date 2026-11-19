par1 <- function(res, name, row = 1L) res$parameters[[name]][row]

# ------------------------------------------------------------ registry
test_that("the parameter registry is complete and consistent", {
  p <- hp_nca_parameters()
  expect_named(p, c("NAME", "KIND", "LABEL"))
  expect_false(anyDuplicated(p$NAME) > 0)
  expect_true(all(c("CMAX", "TMAX", "AUCLAST", "AUCINF_OBS", "LAMZ", "THALF",
                    "CL_F_OBS", "CL_OBS", "VZ_OBS", "VSS_OBS", "AUCTAU", "CAVG") %in% p$NAME))
  expect_true(all(nzchar(p$LABEL)))
})

test_that("units are built from the unit labels", {
  o <- hp_nca_options(dose_unit = "mg", time_unit = "h", conc_unit = "ng/mL")
  expect_equal(hp_nca_unit("CMAX", o), "ng/mL")
  expect_equal(hp_nca_unit("TMAX", o), "h")
  expect_equal(hp_nca_unit("AUCLAST", o), "ng/mL*h")
  expect_equal(hp_nca_unit("LAMZ", o), "1/h")
  expect_equal(hp_nca_unit("CL_F_OBS", o), "mg/(ng/mL*h)")
  expect_equal(hp_nca_unit("VZ_F_OBS", o), "mg/ng/mL")
  expect_equal(hp_nca_unit("AUC_0_4", o), "ng/mL*h")
  expect_equal(hp_nca_unit("FLUCT_PCT", o), "%")
  expect_equal(hp_nca_unit("PTR", o), "")
})

# ------------------------------------------------------------- options
test_that("options have documented defaults and normalise names", {
  o <- hp_nca_options()
  expect_s3_class(o, "hp_nca_options")
  expect_equal(o$auc_method, "linear-up-log-down")
  expect_equal(o$blq_rule, "tmax-zero-missing")
  expect_equal(o$lz$method, "best-fit")
  expect_equal(o$lz$min_points, 3)
  expect_equal(hp_nca_options(auc_method = "Linear_Up_Log_Down")$auc_method, "linear-up-log-down")
  expect_equal(hp_nca_options(blq_rule = "winnonlin")$blq_rule, "tmax-zero-missing")
  expect_equal(hp_nca_options(blq_rule = "monolix")$blq_rule, "first-zero-mid-half-last-missing")
  expect_equal(hp_nca_options(lz_method = "last_n")$lz$method, "last-n")
})

test_that("options reject unknown choices with a helpful message", {
  expect_error(hp_nca_options(auc_method = "simpson"), "AUC method")
  expect_error(hp_nca_options(blq_rule = "guess"), "BLQ rule")
  expect_error(hp_nca_options(steady_state = "maybe"), "steady state")
  expect_error(hp_nca_options(lz_method = "magic"), "lambda-z")
  expect_error(hp_nca_options(duplicates = "x"), "duplicates")
  expect_error(hp_nca_options(c0 = "x"), "C0")
})

# --------------------------------------------------- analytic checks, IV
test_that("an intravenous bolus without noise reproduces the model exactly", {
  D <- 100; V <- 40; k <- 0.15
  r <- hp_nca(iv_bolus(dose = D, v = V, k = k), route = "bolus")
  expect_equal(par1(r, "AUCINF_OBS"), D / (V * k), tolerance = 1e-9)
  expect_equal(par1(r, "AUCINF_PRED"), D / (V * k), tolerance = 1e-9)
  expect_equal(par1(r, "CL_OBS"), V * k, tolerance = 1e-9)
  expect_equal(par1(r, "VZ_OBS"), V, tolerance = 1e-9)
  expect_equal(par1(r, "LAMZ"), k, tolerance = 1e-9)
  expect_equal(par1(r, "THALF"), log(2) / k, tolerance = 1e-9)
  expect_equal(par1(r, "C0"), D / V, tolerance = 1e-9)
  expect_equal(par1(r, "MRTINF_OBS"), 1 / k, tolerance = 1e-6)
  expect_equal(par1(r, "VSS_OBS"), V, tolerance = 1e-6)
  expect_equal(par1(r, "CMAX"), D / V * exp(-k * 0.25))
  expect_match(r$parameters$FLAGS[1], "C0_BACKEXTRAP")
  expect_equal(r$parameters$ROUTE[1], "iv bolus")
  # the oral-only names are empty for an intravenous dose
  expect_true(is.na(par1(r, "CL_F_OBS")))
})

test_that("the linear rule overestimates a convex decline, the default does not", {
  D <- 100; V <- 40; k <- 0.15
  d <- iv_bolus(dose = D, v = V, k = k)
  exact <- D / (V * k)
  lud <- hp_nca(d, route = "bolus")$parameters
  lin <- hp_nca(d, options = hp_nca_options(auc_method = "linear"), route = "bolus")$parameters
  expect_equal(lud$AUCLAST + lud$CLAST / lud$LAMZ, exact, tolerance = 1e-9)
  expect_gt(lin$AUCLAST + lin$CLAST / lin$LAMZ, exact * 1.005)
})

test_that("c0 settings change how the bolus start is found", {
  d <- iv_bolus()
  obs <- hp_nca(d, options = hp_nca_options(c0 = "observed"), route = "bolus")
  expect_true(is.na(par1(obs, "C0")))
  expect_match(obs$parameters$FLAGS[1], "AUC_FROM_FIRST_SAMPLE")
  first <- hp_nca(d, options = hp_nca_options(c0 = "first"), route = "bolus")
  expect_equal(par1(first, "C0"), 100 / 40 * exp(-0.15 * 0.25))
  # an observed C0 is used when it is there
  d0 <- iv_bolus(times = c(0, 0.5, 1, 2, 4, 8, 12, 24))
  r0 <- hp_nca(d0, route = "bolus")
  expect_equal(par1(r0, "C0"), 100 / 40)
  expect_false(grepl("C0_BACKEXTRAP", r0$parameters$FLAGS[1]))
})

# -------------------------------------------------- analytic checks, oral
test_that("an oral profile without noise gives the model values", {
  D <- 100; ka <- 2; ke <- 0.1; V <- 50
  r <- hp_nca(oral_1cmt(dose = D, ka = ka, ke = ke, v = V), route = "oral")
  # the best fit may reach back to where absorption still matters a little
  expect_equal(par1(r, "LAMZ"), ke, tolerance = 0.01)
  expect_equal(par1(r, "THALF"), log(2) / ke, tolerance = 0.01)
  # the last points are purely mono-exponential
  tail4 <- hp_nca(oral_1cmt(dose = D, ka = ka, ke = ke, v = V),
                  options = hp_nca_options(lz_method = "last-n", lz_points = 4), route = "oral")
  expect_equal(par1(tail4, "LAMZ"), ke, tolerance = 1e-8)
  expect_equal(par1(r, "AUCINF_OBS"), D / (V * ke), tolerance = 0.01)
  expect_equal(par1(r, "CL_F_OBS"), V * ke, tolerance = 0.01)
  expect_equal(par1(r, "VZ_F_OBS"), V, tolerance = 0.01)
  expect_equal(r$parameters$ROUTE[1], "extravascular")
  d <- oral_1cmt(dose = D, ka = ka, ke = ke, v = V)
  obs <- d[d$EVID == 0, ]
  expect_equal(par1(r, "CMAX"), max(obs$DV))
  expect_equal(par1(r, "TMAX"), obs$TIME[which.max(obs$DV)])
  expect_equal(par1(r, "CMAX_D"), max(obs$DV) / D)
  expect_equal(par1(r, "TLAST"), 48)
})

test_that("areas equal plain R trapezoids on the profile with a zero at the dose", {
  d <- oral_1cmt()
  obs <- d[d$EVID == 0, ]
  for (m in c("linear", "linear-log", "linear-up-log-down")) {
    r <- hp_nca(d, options = hp_nca_options(auc_method = m), route = "oral")
    expect_equal(par1(r, "AUCLAST"), ref_auc(c(0, obs$TIME), c(0, obs$DV), m),
                 tolerance = 1e-10, info = m)
  }
})

test_that("a hand-worked profile gives the expected areas", {
  d <- make_profile(1, c(1, 2, 4, 8), c(10, 8, 4, 1))
  lin <- hp_nca(d, options = hp_nca_options(auc_method = "linear", lz_min_points = 3),
                route = "oral")
  # 0-1: 5, 1-2: 9, 2-4: 12, 4-8: 10
  expect_equal(par1(lin, "AUCLAST"), 36)
  expect_equal(par1(lin, "CMAX"), 10)
  expect_equal(par1(lin, "TMAX"), 1)
  expect_equal(par1(lin, "TLAG"), 0)
  # AUMC with linear trapezoids: t * C = 0, 10, 16, 16, 8
  expect_equal(par1(lin, "AUMCLAST"), (0 + 10) / 2 * 1 + (10 + 16) / 2 * 1 +
                 (16 + 16) / 2 * 2 + (16 + 8) / 2 * 4)
  lud <- hp_nca(d, route = "oral")
  expect_equal(par1(lud, "AUCLAST"), 5 + (10 - 8) / log(10 / 8) + (8 - 4) / log(2) * 2 + (4 - 1) / log(4) * 4)
})

test_that("TLAG is the last time before the first measurable concentration", {
  d <- make_profile(1, c(0.5, 1, 2, 4, 8, 12), c(0, 0, 3, 5, 2, 1))
  r <- hp_nca(d, route = "oral")
  expect_equal(par1(r, "TLAG"), 1)
})

test_that("CLAST and TLAST are the last measurable concentration", {
  d <- make_profile(1, c(1, 2, 4, 8, 12), c(5, 8, 4, 1, 0))
  r <- hp_nca(d, route = "oral")
  expect_equal(par1(r, "TLAST"), 8)
  expect_equal(par1(r, "CLAST"), 1)
})

# ---------------------------------------------------------- terminal slope
test_that("best-fit picks the window with the highest adjusted R-squared", {
  # a biphasic tail: the last three points are mono-exponential
  t <- c(1, 2, 4, 8, 12, 24, 36)
  c <- c(8, 10, 7, 4, 3, 1, 0.5)
  d <- make_profile(1, t, c)
  r <- hp_nca(d, route = "oral")
  cand <- r$lambda_z[r$lambda_z$CHOSEN == 0 | r$lambda_z$CHOSEN == 1, ]
  expect_equal(sum(cand$CHOSEN), 1)
  best <- cand[which.max(cand$R2ADJ), ]
  chosen <- cand[cand$CHOSEN == 1, ]
  expect_equal(chosen$R2ADJ, best$R2ADJ, tolerance = 1e-4)
  expect_equal(par1(r, "LAMZ"), chosen$LAMBDA_Z)
  expect_equal(par1(r, "THALF"), log(2) / chosen$LAMBDA_Z)
})

test_that("last-n and time-range select the stated points", {
  d <- oral_1cmt(ke = 0.1)
  ln <- hp_nca(d, options = hp_nca_options(lz_method = "last-n", lz_points = 4), route = "oral")
  expect_equal(par1(ln, "LAMZ_N"), 4)
  tr <- hp_nca(d, options = hp_nca_options(lz_method = "time-range", lz_from = 12, lz_to = 36),
               route = "oral")
  expect_equal(par1(tr, "LAMZ_START"), 12)
  expect_equal(par1(tr, "LAMZ_END"), 36)
  none <- hp_nca(d, options = hp_nca_options(lz_method = "none"), route = "oral")
  expect_true(is.na(par1(none, "LAMZ")))
  expect_true(is.na(par1(none, "AUCINF_OBS")))
  expect_match(none$parameters$FLAGS[1], "NO_LZ")
})

test_that("per-subject overrides and limits on the terminal fit work", {
  d <- rbind(oral_1cmt(id = 1), oral_1cmt(id = 2))
  o <- hp_nca_options(lz_overrides = list(list(id = "2", from = 8, to = 48)))
  r <- hp_nca(d, options = o, route = "oral")
  expect_equal(r$parameters$LAMZ_START[2], 8)
  expect_true(r$parameters$LAMZ_START[1] != 8)
  skip <- hp_nca(d, options = hp_nca_options(lz_overrides = list(list(id = "1", skip = TRUE))),
                 route = "oral")
  expect_true(is.na(skip$parameters$LAMZ[1]))
  expect_false(is.na(skip$parameters$LAMZ[2]))
  excl <- hp_nca(d, options = hp_nca_options(lz_overrides = list(list(id = "1", points = 3,
                                                                       exclude = 48))), route = "oral")
  expect_equal(excl$parameters$LAMZ_END[1], 36)
  # a minimum number of points
  few <- hp_nca(oral_1cmt(), options = hp_nca_options(lz_min_points = 5), route = "oral")
  expect_gte(par1(few, "LAMZ_N"), 5)
})

test_that("a poor terminal fit is refused and flagged", {
  d <- make_profile(1, c(1, 2, 4, 8, 12), c(5, 9, 3, 8, 2))
  r <- hp_nca(d, route = "oral")
  expect_true(is.na(par1(r, "LAMZ")))
  expect_match(r$parameters$FLAGS[1], "NO_LZ")
  expect_true(is.na(par1(r, "AUCINF_OBS")))
  expect_false(is.na(par1(r, "AUCLAST")))
})

test_that("a short terminal phase and a large extrapolated area are flagged", {
  d <- make_profile(1, c(0.5, 1, 2, 4, 6), c(2, 5, 4, 3, 2.5))
  r <- hp_nca(d, route = "oral")
  fl <- r$parameters$FLAGS[1]
  expect_match(fl, "LZ_SPAN_SHORT")
  expect_match(fl, "AUC_EXTRAP_GT_20")
  om <- hp_nca(d, options = hp_nca_options(lz_extrap_action = "omit"), route = "oral")
  expect_true(is.na(par1(om, "AUCINF_OBS")))
  expect_match(om$parameters$FLAGS[1], "INF_OMITTED")
  req <- hp_nca(d, options = hp_nca_options(lz_require_span = TRUE), route = "oral")
  expect_true(is.na(par1(req, "LAMZ")))
})

# --------------------------------------------------------------- BLQ rules
blq_data <- function() {
  make_profile(1, c(0.5, 1, 2, 4, 8, 12, 24), c(0.1, 3, 2, 1, 0.5, 0.2, 0.05))
}
blq_status <- function(rule, lloq = 0.3, d = blq_data()) {
  r <- hp_nca(d, options = hp_nca_options(lloq = lloq, blq_rule = rule), route = "oral")
  p <- r$profiles
  p[p$STATUS != "added zero at dose", ]
}

test_that("every BLQ rule treats the samples as defined", {
  z <- blq_status("zero")
  expect_equal(z$STATUS[z$TIME %in% c(0.5, 12, 24)], rep("blq set to zero", 3))
  expect_equal(z$CONC_USED[z$TIME %in% c(0.5, 12, 24)], c(0, 0, 0))

  m <- blq_status("missing")
  expect_equal(m$STATUS[m$TIME %in% c(0.5, 12, 24)], rep("blq left out", 3))
  expect_true(all(is.na(m$CONC_USED[m$TIME %in% c(0.5, 12, 24)])))

  l <- blq_status("lloq")
  expect_equal(l$CONC_USED[l$TIME %in% c(0.5, 12, 24)], rep(0.3, 3))

  h <- blq_status("half-lloq")
  expect_equal(h$CONC_USED[h$TIME %in% c(0.5, 12, 24)], rep(0.15, 3))

  w <- blq_status("tmax-zero-missing")
  expect_equal(w$STATUS[w$TIME == 0.5], "blq set to zero")
  expect_equal(w$STATUS[w$TIME %in% c(12, 24)], rep("blq left out", 2))

  mo <- blq_status("first-zero-mid-half-last-missing")
  expect_equal(mo$STATUS[mo$TIME == 0.5], "blq set to zero")
  expect_equal(mo$STATUS[mo$TIME %in% c(12, 24)], rep("blq left out", 2))
})

test_that("the Monolix rule uses half the limit between quantifiable samples", {
  d <- make_profile(1, c(1, 2, 4, 8, 12, 24), c(3, 2, 0.1, 1, 0.5, 0.2))
  s <- blq_status("first-zero-mid-half-last-missing", lloq = 0.3, d = d)
  expect_equal(s$STATUS[s$TIME == 4], "blq set to lloq/2")
  expect_equal(s$CONC_USED[s$TIME == 4], 0.15)
  expect_equal(s$STATUS[s$TIME == 24], "blq left out")
})

test_that("BLQ rules change the areas as expected", {
  d <- blq_data()
  area <- function(rule) {
    hp_nca(d, options = hp_nca_options(lloq = 0.3, blq_rule = rule), route = "oral")$parameters$AUCLAST
  }
  # the sample of 0.1 at 0.5 h: set to zero it lowers the first segment by
  # 0.75 relative to leaving it out (0 -> 3 over 0.5 h against 0 -> 3 over 1 h)
  expect_equal(area("missing") - area("zero"), 0.75, tolerance = 1e-10)
  # BLQ samples after the last measurable concentration add no area
  expect_equal(area("zero"), area("tmax-zero-missing"), tolerance = 1e-12)
  expect_gt(area("lloq"), area("zero"))
  expect_gt(area("half-lloq"), area("zero"))
})

test_that("a rule that needs a limit stops without one", {
  d <- make_profile(1, c(0.5, 1, 2, 4), c(0, 3, 2, 1))
  expect_error(hp_nca(d, options = hp_nca_options(blq_rule = "half-lloq", blq_flagged = TRUE),
                      route = "oral", blq_flag = "MDV"), "limit of quantification")
})

test_that("a BLQ flag column marks samples without a limit", {
  d <- blq_data()
  d$BLQ <- c(0, 1, 0, 0, 0, 0, 1, 1)
  r <- hp_nca(d, blq_flag = "BLQ", route = "oral")
  p <- r$profiles
  expect_true(all(p$STATUS[p$TIME %in% c(0.5, 24)] %in% c("blq set to zero", "blq left out")))
  # without lloq and without a flag nothing is BLQ
  r0 <- hp_nca(d, route = "oral")
  expect_false(any(grepl("blq", r0$profiles$STATUS)))
})

test_that("log-transformed DV is back-transformed", {
  d <- oral_1cmt()
  a <- hp_nca(d, route = "oral")
  d2 <- d
  d2$DV <- log(d2$DV)
  b <- hp_nca(d2, route = "oral", dv_log = TRUE)
  expect_equal(a$parameters$AUCLAST, b$parameters$AUCLAST, tolerance = 1e-10)
})

# ------------------------------------------------------------ duplicates etc.
test_that("duplicate samples are combined or refused", {
  d <- make_profile(1, c(1, 2, 2, 4, 8, 12), c(5, 8, 6, 4, 2, 1))
  m <- hp_nca(d, route = "oral")
  expect_match(paste(m$notes, collapse = " "), "combined")
  mx <- hp_nca(d, options = hp_nca_options(duplicates = "max"), route = "oral")
  expect_equal(mx$profiles$CONC_USED[mx$profiles$TIME == 2], 8)
  fi <- hp_nca(d, options = hp_nca_options(duplicates = "first"), route = "oral")
  expect_equal(fi$profiles$CONC_USED[fi$profiles$TIME == 2], 8)
  la <- hp_nca(d, options = hp_nca_options(duplicates = "last"), route = "oral")
  expect_equal(la$profiles$CONC_USED[la$profiles$TIME == 2], 6)
  expect_equal(m$profiles$CONC_USED[m$profiles$TIME == 2], 7)
  expect_error(hp_nca(d, options = hp_nca_options(duplicates = "error"), route = "oral"),
               "distinct times")
})

test_that("profiles with too few samples are skipped with a reason", {
  d <- rbind(make_profile(1, c(1, 2), c(5, 3)), oral_1cmt(id = 2))
  r <- hp_nca(d, route = "oral")
  expect_match(r$parameters$SKIPPED[r$parameters$ID == "1"], "usable sample")
  expect_equal(r$parameters$SKIPPED[r$parameters$ID == "2"], "")
  expect_true(is.na(r$parameters$CMAX[r$parameters$ID == "1"]))
  expect_equal(nrow(r$parameters), 2)
  expect_match(capture.output(print(r))[1], "1 profile")
})

test_that("a dataset where nothing can be analysed is an error", {
  d <- make_profile(1, c(1, 2), c(5, 3))
  expect_error(hp_nca(d, route = "oral"), "none of the")
  expect_error(hp_nca(data.frame(ID = 1, TIME = 0, DV = 1), route = "oral"), "amt")
})

test_that("samples before the first dose are not analysed", {
  d <- rbind(
    data.frame(ID = 1, TIME = -1, AMT = 0, DV = 0.2, EVID = 0, MDV = 0, CMT = 2),
    oral_1cmt())
  r <- hp_nca(d, route = "oral")
  expect_match(paste(r$notes, collapse = " "), "before a subject's first dose")
  expect_false(-1 %in% r$profiles$TIME)
})

test_that("max_time and min_observations limit the profile", {
  d <- oral_1cmt()
  r <- hp_nca(d, options = hp_nca_options(max_time = 12), route = "oral")
  expect_equal(par1(r, "TLAST"), 12)
  expect_error(hp_nca(d, options = hp_nca_options(min_observations = 100), route = "oral"),
               "usable sample")
})

# -------------------------------------------------------------- partial AUC
test_that("partial areas are cut with the interpolation of the method", {
  d <- oral_1cmt()
  o <- hp_nca_options(partial = list(c(0, 4), c(2, 12)))
  r <- hp_nca(d, options = o, route = "oral")
  obs <- d[d$EVID == 0, ]
  tt <- c(0, obs$TIME); cc <- c(0, obs$DV)
  expect_equal(par1(r, "AUC_0_4"), ref_auc(tt[tt <= 4], cc[tt <= 4]), tolerance = 1e-10)
  keep <- tt >= 2 & tt <= 12
  expect_equal(par1(r, "AUC_2_12"), ref_auc(tt[keep], cc[keep]), tolerance = 1e-10)
  # beyond the last sample the tail is extrapolated
  far <- hp_nca(d, options = hp_nca_options(partial = list(c(0, 200))), route = "oral")
  expect_equal(par1(far, "AUC_0_200"), par1(far, "AUCINF_OBS"), tolerance = 1e-6)
  noex <- hp_nca(d, options = hp_nca_options(partial = list(c(0, 200)), partial_extrapolate = FALSE),
                 route = "oral")
  expect_true(is.na(par1(noex, "AUC_0_200")))
})

# ------------------------------------------------------------ infusion
test_that("an intravenous infusion reproduces the model", {
  D <- 100; V <- 40; k <- 0.15; Tinf <- 1; CL <- V * k
  tt <- c(seq(0.1, 1, by = 0.1), 1.25, 1.5, 2, 3, 4, 6, 8, 12, 16, 24)
  cp <- ifelse(tt <= Tinf, D / Tinf / CL * (1 - exp(-k * tt)),
               D / Tinf / CL * (1 - exp(-k * Tinf)) * exp(-k * (tt - Tinf)))
  d <- make_profile(1, tt, cp, D, obs_cmt = 1)
  d$RATE <- c(D / Tinf, rep(0, length(tt)))
  r <- hp_nca(d, route = "infusion")
  expect_equal(r$parameters$ROUTE[1], "iv infusion")
  expect_equal(par1(r, "INF_DURATION"), Tinf)
  expect_equal(par1(r, "LAMZ"), k, tolerance = 1e-8)
  expect_equal(par1(r, "AUCINF_OBS"), D / CL, tolerance = 0.005)
  expect_equal(par1(r, "CL_OBS"), CL, tolerance = 0.005)
  expect_equal(par1(r, "MRTINF_OBS"), 1 / k, tolerance = 0.02)
})

# ---------------------------------------------------------- steady state
test_that("steady state gives the dosing-interval measures", {
  D <- 100; ka <- 1.5; ke <- 0.1; V <- 50; tau <- 12
  tt <- c(seq(0.25, 4, by = 0.25), 5, 6, 8, 10, 12)
  cc <- D * ka / (V * (ka - ke)) *
    (exp(-ke * tt) / (1 - exp(-ke * tau)) - exp(-ka * tt) / (1 - exp(-ka * tau)))
  d <- make_profile(1, tt, cc, D, extra = list(II = c(tau, rep(0, length(tt))),
                                                SS = c(1, rep(0, length(tt)))))
  r <- hp_nca(d, route = "oral")
  expect_equal(par1(r, "STEADY_STATE"), 1)
  expect_equal(par1(r, "TAU"), tau)
  CL <- V * ke
  expect_equal(par1(r, "AUCTAU"), D / CL, tolerance = 0.005)
  expect_equal(par1(r, "CL_F_SS"), CL, tolerance = 0.005)
  expect_equal(par1(r, "CAVG"), D / CL / tau, tolerance = 0.005)
  expect_equal(par1(r, "CMIN"), min(cc))
  expect_equal(par1(r, "PTR"), max(cc) / min(cc))
  expect_equal(par1(r, "FLUCT_PCT"), 100 * (max(cc) - min(cc)) / par1(r, "CAVG"))
  expect_equal(par1(r, "SWING"), (max(cc) - min(cc)) / min(cc))
  expect_equal(par1(r, "RACC_PRED"), 1 / (1 - exp(-ke * tau)), tolerance = 0.02)
  # infinity-based measures are not reported at steady state
  expect_true(is.na(par1(r, "AUCINF_OBS")))
})

test_that("steady_state = none treats the interval as a single dose", {
  tau <- 12
  tt <- c(0.5, 1, 2, 4, 8, 12)
  d <- make_profile(1, tt, c(3, 4, 3.5, 2.5, 1.5, 1), extra = list(II = c(tau, rep(0, 6)),
                                                                    SS = c(1, rep(0, 6))))
  r <- hp_nca(d, options = hp_nca_options(steady_state = "none"), route = "oral")
  expect_equal(par1(r, "STEADY_STATE"), 0)
  a <- hp_nca(d, route = "oral")
  expect_equal(par1(a, "STEADY_STATE"), 1)
})

test_that("a steady-state flag without an interval is reported", {
  d <- make_profile(1, c(0.5, 1, 2, 4, 8), c(3, 4, 3.5, 2.5, 1.5), extra = list(SS = c(1, rep(0, 5))))
  r <- hp_nca(d, options = hp_nca_options(steady_state = "last"), route = "oral")
  expect_match(paste(r$notes, collapse = " "), "dosing interval")
  expect_match(r$parameters$FLAGS[1], "SS_NO_TAU")
})

# ----------------------------------------------------- repeated doses
test_that("ADDL and II are expanded into dosing intervals", {
  conc_at <- function(t, doses = c(0, 24, 48), ka = 1.2, ke = 0.1, v = 50) {
    sum(vapply(doses, function(dd) {
      if (t < dd) 0 else 100 * ka / (v * (ka - ke)) * (exp(-ke * (t - dd)) - exp(-ka * (t - dd)))
    }, 1))
  }
  tm <- c(0.5, 1, 2, 4, 8, 12, 24, 24.5, 25, 26, 28, 32, 36, 48, 48.5, 49, 50, 52, 56, 60, 72)
  d <- rbind(
    data.frame(ID = 1, TIME = 0, AMT = 100, DV = NA, EVID = 1, MDV = 1, CMT = 1, ADDL = 2, II = 24),
    data.frame(ID = 1, TIME = tm, AMT = 0, DV = vapply(tm, conc_at, 1), EVID = 0, MDV = 0,
               CMT = 2, ADDL = 0, II = 0))
  r <- hp_nca(d, route = "oral")
  expect_equal(r$parameters$OCCASION, 1:3)
  expect_equal(r$parameters$TAU[1:2], c(24, 24))
  expect_equal(r$parameters$DOSE, rep(100, 3))
  expect_gt(r$parameters$CMAX[3], r$parameters$CMAX[1])
  # an interval that follows an earlier dose includes its tail
  expect_match(r$parameters$FLAGS[2], "INF_AFTER_EARLIER_DOSES")
  # intervals = "first" and "last"
  f <- hp_nca(d, options = hp_nca_options(intervals = "first"), route = "oral")
  expect_equal(f$parameters$OCCASION, 1)
  l <- hp_nca(d, options = hp_nca_options(intervals = "last"), route = "oral")
  expect_equal(l$parameters$OCCASION, 3)
})

# ----------------------------------------------------------- route
test_that("the route is decided from the data when auto", {
  oral <- hp_nca(oral_1cmt(), route = "auto")
  expect_equal(oral$parameters$ROUTE, "extravascular")
  iv <- hp_nca(iv_bolus(), route = "auto")
  expect_equal(iv$parameters$ROUTE, "iv bolus")
  expect_equal(hp_nca(oral_1cmt(), route = "po")$parameters$ROUTE, "extravascular")
  expect_error(hp_nca(oral_1cmt(), route = "rectal"), "route")
})

test_that("mixed routes use the compartment of the dose", {
  both <- rbind(iv_bolus(id = 1), oral_1cmt(id = 2))
  both$CMT[both$ID == 1] <- 2   # the bolus goes into the observation compartment
  r <- hp_nca(both, route = "mixed")
  expect_equal(r$parameters$ROUTE, c("iv bolus", "extravascular"))
  expect_error(hp_nca(both[, setdiff(names(both), "CMT")], route = "mixed"), "CMT")
})

# ------------------------------------------------------------ groups etc.
test_that("groups, filters and column names are honoured", {
  d <- example_data()
  r <- hp_nca(d, route = "oral", lloq = 0.5, groups = "DOSEGRP")
  expect_equal(sort(unique(r$parameters$GROUP)), c("DOSEGRP=100", "DOSEGRP=200"))
  f <- hp_nca(d, route = "oral", lloq = 0.5, filter = "DOSEGRP == 100")
  expect_equal(nrow(f$parameters), 12)
  of <- hp_nca(d, route = "oral", lloq = 0.5, obs_filter = "TIME <= 12")
  expect_equal(max(of$profiles$TIME), 12)
  expect_match(paste(of$notes, collapse = " "), "observation filter")
  ren <- d
  names(ren)[names(ren) == "ID"] <- "SUBJ"
  names(ren)[names(ren) == "TIME"] <- "HOURS"
  rr <- hp_nca(ren, route = "oral", lloq = 0.5, id = "SUBJ", time = "HOURS")
  expect_equal(rr$parameters$CMAX, hp_nca(d, route = "oral", lloq = 0.5)$parameters$CMAX)
  expect_error(hp_nca(d, route = "oral", id = "NOPE"), "NOPE|id")
})

# ----------------------------------------------------- the example study
test_that("the example study gives stable results", {
  r <- hp_nca(example_data(), route = "oral", lloq = 0.5, groups = "DOSEGRP")
  p <- r$parameters
  expect_equal(nrow(p), 24)
  expect_equal(sum(p$SKIPPED != ""), 0)
  expect_equal(p$CMAX[1], 1.7568)
  expect_equal(p$TMAX[1], 4)
  expect_equal(p$AUCLAST[1], 13.01264, tolerance = 1e-6)
  expect_equal(p$THALF[1], 6.315145, tolerance = 1e-6)
  expect_true(all(p$CMAX[p$GROUP == "DOSEGRP=200"] > 2))
  expect_true(all(p$AUCINF_OBS >= p$AUCLAST, na.rm = TRUE))
  expect_true(all(abs(p$THALF * p$LAMZ - log(2)) < 1e-10, na.rm = TRUE))
})

# ---------------------------------------------------------------- summaries
test_that("the summary reproduces base R statistics", {
  r <- hp_nca(example_data(), route = "oral", lloq = 0.5, groups = "DOSEGRP")
  s <- r$summary
  expect_true(all(c("STRATUM", "PARAMETER", "UNIT", "N", "MEAN", "SD", "CV_PCT", "GEOMEAN",
                    "GEOCV_PCT", "MEDIAN", "MIN", "MAX") %in% names(s)))
  p <- r$parameters
  g <- p$CMAX[p$GROUP == "DOSEGRP=100"]
  row <- s[s$PARAMETER == "CMAX" & grepl("100", s$STRATUM), ]
  expect_equal(row$N, 12)
  expect_equal(row$MEAN, mean(g))
  expect_equal(row$SD, sd(g))
  expect_equal(row$CV_PCT, 100 * sd(g) / mean(g))
  expect_equal(row$GEOMEAN, exp(mean(log(g))))
  expect_equal(row$MEDIAN, median(g))
  expect_equal(row$UNIT, "ng/mL")
  # bookkeeping columns are not summarised unless asked for
  expect_false("N_SAMPLES" %in% s$PARAMETER)
  s2 <- hp_nca_summary(p, by = "dose", statistics = c("n", "median", "q1", "q3", "harmonic", "sem"),
                       parameters = c("CMAX", "N_SAMPLES"), options = r$options)
  expect_setequal(unique(s2$PARAMETER), c("CMAX", "N_SAMPLES"))
  expect_true(all(c("Q1", "Q3", "HARMONIC_MEAN", "SEM") %in% names(s2)))
  expect_equal(s2$MEDIAN[s2$PARAMETER == "CMAX"][1], median(p$CMAX[p$DOSE == 100]))
})

test_that("a single-group study gets one stratum", {
  s <- hp_nca(example_data(), route = "oral", lloq = 0.5, filter = "DOSEGRP == 100")$summary
  expect_equal(unique(s$STRATUM), "all")
})

# ----------------------------------------------------------- documentation tables
test_that("the lambda-z and profile tables document the analysis", {
  r <- hp_nca(example_data(), route = "oral", lloq = 0.5)
  lz <- r$lambda_z
  expect_true(all(c("ID", "CHOSEN", "POINTS", "START", "END", "LAMBDA_Z", "HALF_LIFE", "R2",
                    "R2ADJ", "SPAN_RATIO", "HOW") %in% names(lz)))
  chosen <- lz[lz$CHOSEN == 1, ]
  expect_equal(nrow(chosen), sum(!is.na(r$parameters$LAMZ)))
  expect_equal(chosen$LAMBDA_Z, r$parameters$LAMZ[!is.na(r$parameters$LAMZ)])
  pr <- r$profiles
  expect_true(all(c("ID", "TIME", "CONC", "CONC_USED", "STATUS", "IN_LAMBDA_Z", "LINE") %in% names(pr)))
  expect_equal(sum(pr$IN_LAMBDA_Z[pr$ID == "1"]), r$parameters$LAMZ_N[1])
  # LINE points at the row of the dataset (the header is row 1)
  d <- example_data()
  obs <- pr[pr$STATUS == "observed", ]
  expect_equal(d$DV[obs$LINE - 1L], obs$CONC)
  expect_equal(hp_nca_table(r$results), r$parameters)
  expect_equal(hp_nca_lambda_z_table(r$results), r$lambda_z)
  expect_equal(hp_nca_profile_table(r$results), r$profiles)
})

test_that("the step functions can be run one by one", {
  d <- example_data()
  o <- hp_nca_options(lloq = 0.5)
  x <- hp_nca_read(d, route = "oral", lloq = 0.5)
  expect_length(x$subjects, 24)
  expect_equal(nrow(x$subjects[[1]]$doses), 1)
  p <- hp_nca_profiles(x, o)
  expect_length(p, 24)
  expect_equal(p[[1]]$route, "extravascular")
  a <- hp_nca_analyse(p[[1]], o)
  expect_equal(a$params$CMAX, hp_nca(d, route = "oral", lloq = 0.5)$parameters$CMAX[1])
  expect_named(a, c("id", "group", "occasion", "route", "dose", "tau", "steady_state", "params",
                    "lz", "flags", "points", "skip"), ignore.order = TRUE)
})

test_that("dose_factor scales clearance and volume only", {
  d <- oral_1cmt()
  a <- hp_nca(d, route = "oral")
  b <- hp_nca(d, options = hp_nca_options(dose_factor = 1000), route = "oral")
  expect_equal(b$parameters$CL_F_OBS, 1000 * a$parameters$CL_F_OBS)
  expect_equal(b$parameters$VZ_F_OBS, 1000 * a$parameters$VZ_F_OBS)
  expect_equal(b$parameters$AUCINF_OBS, a$parameters$AUCINF_OBS)
})

test_that("results do not depend on the order of the subjects", {
  d <- example_data()
  a <- hp_nca(d, route = "oral", lloq = 0.5)$parameters
  ids <- rev(unique(d$ID))
  d2 <- do.call(rbind, lapply(ids, function(i) d[d$ID == i, ]))
  b <- hp_nca(d2, route = "oral", lloq = 0.5)$parameters
  b <- b[order(as.numeric(b$ID)), ]
  rownames(b) <- NULL
  expect_equal(a, b)
})
