# Small datasets with known answers, shared by the tests. Nothing here writes
# to the user's file space.

# one subject: a dose record and observation records
make_profile <- function(id, times, conc, dose = 100, obs_cmt = 2, extra = list()) {
  d <- rbind(
    data.frame(ID = id, TIME = 0, AMT = dose, DV = NA_real_, EVID = 1, MDV = 1, CMT = 1),
    data.frame(ID = id, TIME = times, AMT = 0, DV = conc, EVID = 0, MDV = 0, CMT = obs_cmt)
  )
  for (nm in names(extra)) d[[nm]] <- extra[[nm]]
  d
}

# one-compartment model, intravenous bolus, no noise
iv_bolus <- function(id = 1, dose = 100, v = 40, k = 0.15,
                     times = c(0.25, 0.5, 1, 2, 4, 6, 8, 12, 24)) {
  make_profile(id, times, dose / v * exp(-k * times), dose, obs_cmt = 1)
}

# one-compartment model, first-order absorption, no noise (F = 1)
oral_1cmt <- function(id = 1, dose = 100, ka = 2, ke = 0.1, v = 50,
                      times = c(seq(0.25, 4, by = 0.25), 5, 6, 8, 10, 12, 16, 24, 36, 48)) {
  cp <- dose * ka / (v * (ka - ke)) * (exp(-ke * times) - exp(-ka * times))
  make_profile(id, times, cp, dose)
}

example_data <- function() hp_read_data(hp_example_file())

# trapezoids written in plain R: the reference the C++ code is checked against
ref_auc <- function(t, c, method = "linear-up-log-down") {
  tot <- 0
  tmax <- t[which.max(c)]
  for (i in seq_len(length(t) - 1L)) {
    dt <- t[i + 1L] - t[i]; c1 <- c[i]; c2 <- c[i + 1L]
    use_log <- switch(method,
      linear = FALSE,
      "linear-log" = t[i] >= tmax && c1 > 0 && c2 > 0 && c1 != c2,
      "linear-up-log-down" = c2 < c1 && c1 > 0 && c2 > 0)
    tot <- tot + if (use_log) (c1 - c2) / log(c1 / c2) * dt else (c1 + c2) / 2 * dt
  }
  tot
}
