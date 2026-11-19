test_that("hp_seg_auc matches the linear and logarithmic trapezoid", {
  expect_equal(hp_seg_auc(0, 4, 1, 5, FALSE), 4.5)
  expect_equal(hp_seg_auc(2, 5, 4, 3, FALSE), 8)
  expect_equal(hp_seg_auc(2, 5, 4, 3, TRUE), (5 - 3) / log(5 / 3) * 2)
  # equal concentrations: the logarithmic rule reduces to a rectangle
  expect_equal(hp_seg_auc(0, 3, 2, 3, TRUE), 6)
  # a zero concentration is always integrated linearly
  expect_equal(hp_seg_auc(0, 4, 2, 0, TRUE), 4)
  # nearly equal concentrations lose no accuracy
  expect_equal(hp_seg_auc(0, 1, 1, 1 + 1e-12, TRUE), 1, tolerance = 1e-9)
})

test_that("hp_seg_aumc integrates t * C(t)", {
  expect_equal(hp_seg_aumc(0, 2, 2, 4, FALSE), (0 * 2 + 2 * 4) / 2 * 2)
  k <- log(5 / 3) / 2
  f <- function(t) t * 5 * exp(-k * (t - 2))
  expect_equal(hp_seg_aumc(2, 5, 4, 3, TRUE), integrate(f, 2, 4)$value, tolerance = 1e-6)
})

test_that("hp_curve_area agrees with plain R trapezoids for the three methods", {
  t <- c(0, 1, 2, 4, 8, 12)
  c <- c(0, 4, 5, 3, 1, 0.4)
  for (m in 0:2) {
    nm <- c("linear", "linear-log", "linear-up-log-down")[m + 1L]
    expect_equal(hp_curve_area(t, c, m, 2, 0, 12), ref_auc(t, c, nm), tolerance = 1e-12)
  }
})

test_that("hp_curve_area cuts partial areas with the interpolation of the method", {
  t <- c(0, 1, 2, 4, 8, 12)
  c <- c(0, 4, 5, 3, 1, 0.4)
  c6 <- hp_curve_conc_at(t, c, 2L, 2, 6)
  expect_equal(c6, 3 * exp(-log(3) / 4 * 2))
  expect_equal(hp_curve_area(t, c, 2L, 2, 0, 6),
               ref_auc(c(t[1:4], 6), c(c[1:4], c6)))
  # additivity at an arbitrary inner time
  for (m in 0:2) {
    whole <- hp_curve_area(t, c, m, 2, 0, 12)
    split <- hp_curve_area(t, c, m, 2, 0, 5.5) + hp_curve_area(t, c, m, 2, 5.5, 12)
    expect_equal(whole, split, tolerance = 1e-12)
  }
  # limits outside the curve add nothing
  expect_equal(hp_curve_area(t, c, 2L, 2, 0, 100), hp_curve_area(t, c, 2L, 2, 0, 12))
})

test_that("the first-moment area gives the mean residence time of a mono-exponential", {
  k <- 0.2
  t <- seq(0, 60, by = 0.05)
  c <- exp(-k * t)
  aumc <- hp_curve_area(t, c, 2L, 0, 0, 60, TRUE)
  auc <- hp_curve_area(t, c, 2L, 0, 0, 60)
  expect_equal(auc, (1 - exp(-k * 60)) / k, tolerance = 1e-8)
  expect_equal(aumc / auc, 1 / k, tolerance = 1e-3)
})

test_that("hp_curve_conc_at interpolates", {
  t <- c(0, 2, 4); c <- c(0, 6, 2)
  expect_equal(hp_curve_conc_at(t, c, 0L, 2, 1), 3)
  expect_equal(hp_curve_conc_at(t, c, 0L, 2, 3), 4)
  expect_equal(hp_curve_conc_at(t, c, 2L, 2, 3), sqrt(6 * 2))
  expect_equal(hp_curve_conc_at(t, c, 2L, 2, 2), 6)
})

test_that("hp_fit_line is the least-squares line", {
  x <- c(1, 2, 3, 4, 5)
  y <- c(2.1, 3.9, 6.2, 7.8, 10.3)
  f <- hp_fit_line(x, y)
  m <- lm(y ~ x)
  s <- summary(m)
  expect_equal(f[1:2], unname(coef(m)))
  expect_equal(f[3], s$r.squared)
  expect_equal(f[4], s$adj.r.squared)
  expect_equal(f[5], s$coefficients["x", "Std. Error"])
  expect_equal(f[6], 3)
})

test_that("hp_fit_line honours weights and degenerate input", {
  x <- c(1, 2, 3, 4, 5); y <- c(2.1, 3.9, 6.2, 7.8, 10.3); w <- c(1, 2, 3, 2, 1)
  expect_equal(hp_fit_line(x, y, w)[1:2], unname(coef(lm(y ~ x, weights = w))))
  expect_true(all(is.na(hp_fit_line(1, 2))))
  expect_true(all(is.na(hp_fit_line(c(2, 2, 2), c(1, 2, 3)))))
  two <- hp_fit_line(c(1, 2), c(1, 3))
  expect_equal(two[1:2], c(-1, 2))
  expect_true(is.na(two[4]))
})

test_that("hp_lz_windows regresses log concentration on the last n points", {
  t <- c(2, 4, 6, 8, 12, 24)
  c <- c(5, 3, 2, 1.4, 0.7, 0.1)
  w <- hp_lz_windows(t, c, 3:5)
  expect_equal(nrow(w), 3)
  expect_equal(unname(w[, "n"]), 3:5)
  for (i in 1:3) {
    n <- 2 + i
    idx <- (length(t) - n + 1):length(t)
    m <- lm(log(c[idx]) ~ t[idx])
    expect_equal(unname(w[i, "lambda"]), -unname(coef(m)[2]), tolerance = 1e-10)
    expect_equal(unname(w[i, "intercept"]), unname(coef(m)[1]), tolerance = 1e-10)
    expect_equal(unname(w[i, "adj_r2"]), summary(m)$adj.r.squared, tolerance = 1e-10)
    expect_equal(unname(w[i, "start"]), t[idx[1]])
    expect_equal(unname(w[i, "end"]), 24)
  }
})

test_that("hp_lz_windows weights are 1/y and 1/y^2", {
  t <- c(2, 4, 6, 8, 12, 24)
  c <- c(5, 3, 2, 1.4, 0.7, 0.1)
  idx <- 3:6
  for (code in 1:2) {
    w <- hp_lz_windows(t, c, 4L, code)
    ref <- lm(log(c[idx]) ~ t[idx], weights = 1 / c[idx]^code)
    expect_equal(unname(w[1, "lambda"]), -unname(coef(ref)[2]), tolerance = 1e-8)
  }
})

test_that("hp_lz_windows skips sizes that cannot be fitted", {
  w <- hp_lz_windows(c(1, 2, 3), c(3, 2, 1), c(2L, 3L, 6L))
  expect_equal(unname(w[, "n"]), c(2, 3))
})

test_that("hp_exposure gives the observed exposure measures", {
  t <- c(0.5, 1, 2, 4, 8, 12, 24)
  c <- c(1.2, 2.8, 3.5, 2.6, 1.3, 0.7, 0.1)
  e <- hp_exposure(t, c)
  expect_equal(unname(e["cmax"]), 3.5)
  expect_equal(unname(e["tmax"]), 2)
  expect_equal(unname(e["clast"]), 0.1)
  expect_equal(unname(e["tlast"]), 24)
  expect_equal(unname(e["nobs"]), 7)
  # a zero at time zero is assumed when there is no sample there
  expect_equal(unname(e["auclast"]), ref_auc(c(0, t), c(0, c)), tolerance = 1e-12)
})
