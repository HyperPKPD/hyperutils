nca_fixture <- function() hp_nca(example_data(), route = "oral", lloq = 0.5, groups = "DOSEGRP")

test_that("figure options have defaults and keep what is given", {
  o <- hp_nca_figure_options()
  expect_s3_class(o, "hp_nca_figures")
  expect_equal(o$error_bars, "sd")
  expect_true(o$profiles && o$mean_profile)
  o2 <- hp_nca_figure_options(boxplots = c("CMAX", "AUCLAST"), error_bars = "ci95", individual = FALSE,
                              units = list(time = "h", conc = "ng/mL", dose = "mg"))
  expect_equal(o2$boxplots, c("CMAX", "AUCLAST"))
  expect_equal(o2$error_bars, "ci95")
  expect_false(o2$individual)
})

test_that("hp_nca_plots draws the documented figures", {
  r <- nca_fixture()
  f <- hp_nca_plots(r, hp_nca_figure_options(boxplots = c("CMAX", "AUCLAST"), individual_nrow = 2,
                                             individual_ncol = 3))
  expect_true(all(c("profiles_linear", "profiles_semilog", "mean_profile_linear", "mean_profile_semilog",
                    "box_cmax", "box_auclast", "individual_01") %in% names(f)))
  # 24 profiles in panels of 2 x 3
  expect_equal(sum(grepl("^individual_", names(f))), 4)
  for (p in f) {
    expect_s3_class(p, "ggplot")
    expect_no_error(ggplot2::ggplot_build(p))
  }
})

test_that("figure options switch figures on and off", {
  r <- nca_fixture()
  f <- hp_nca_plots(r, hp_nca_figure_options(profiles = FALSE, mean_profile = TRUE, semilog = FALSE,
                                             individual = FALSE))
  expect_equal(names(f), "mean_profile_linear")
  g <- hp_nca_plots(r, hp_nca_figure_options(mean_profile = FALSE, linear = FALSE, individual = FALSE))
  expect_equal(names(g), "profiles_semilog")
})

test_that("mean profiles can use every kind of error bar and the geometric mean", {
  r <- nca_fixture()
  for (e in c("sd", "sem", "ci95", "off")) {
    f <- hp_nca_plots(r, hp_nca_figure_options(error_bars = e, geometric_mean = e == "sd", individual = FALSE,
                                               profiles = FALSE, linear = FALSE))
    expect_no_error(ggplot2::ggplot_build(f$mean_profile_semilog))
  }
  f <- hp_nca_plots(r, hp_nca_figure_options(time_bin = 1, individual = FALSE, profiles = FALSE))
  expect_no_error(ggplot2::ggplot_build(f$mean_profile_linear))
})

test_that("the individual pages show the terminal fit", {
  r <- nca_fixture()
  f <- hp_nca_plots(r, hp_nca_figure_options(individual_nrow = 2, individual_ncol = 2, max_individual = 4))
  expect_equal(sum(grepl("^individual_", names(f))), 1)
  lines <- vapply(f$individual_01$layers, function(l) class(l$geom)[1], "")
  expect_true(any(lines %in% c("GeomLine", "GeomPath", "GeomSegment", "GeomAbline")))
  none <- hp_nca_plots(r, hp_nca_figure_options(terminal_points = FALSE, individual_lines = FALSE,
                                                individual_yscale = "linear", max_individual = 2))
  expect_no_error(ggplot2::ggplot_build(none$individual_01))
})

test_that("figures of an intravenous analysis can be drawn", {
  d <- do.call(rbind, lapply(1:3, function(i) iv_bolus(id = i, k = 0.1 + 0.02 * i)))
  r <- hp_nca(d, route = "bolus")
  f <- hp_nca_plots(r, hp_nca_figure_options(units = list(time = "h", conc = "mg/L", dose = "mg")))
  expect_true("profiles_linear" %in% names(f))
  expect_no_error(ggplot2::ggplot_build(f$profiles_semilog))
})

test_that("the NCA figures can be saved", {
  r <- nca_fixture()
  f <- hp_nca_plots(r, hp_nca_figure_options(individual = FALSE))
  dir <- file.path(tempdir(), "t-nca-figs")
  files <- vapply(names(f), function(nm) hp_save_plot(f[[nm]], nm, dir, 5, 3, dpi = 50), "")
  expect_true(all(file.exists(files)))
  unlink(dir, recursive = TRUE)
})
