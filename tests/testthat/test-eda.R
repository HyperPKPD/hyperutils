two_dose <- function() {
  data.frame(
    ID = c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2),
    TIME = c(0, 1, 2, 12, 13, 14, 0, 1, 2, 4),
    AMT = c(100, 0, 0, 150, 0, 0, 200, 0, 0, 0),
    DV = c(NA, 2, 1.5, NA, 3, 2.5, NA, 0.2, 5, 3),
    EVID = c(1, 0, 0, 1, 0, 0, 1, 0, 0, 0),
    MDV = c(1, 0, 0, 1, 0, 0, 1, 0, 0, 0),
    WT = c(70, 70, 70, 70, 70, 70, 80, 80, 80, 80),
    SEX = c("M", "M", "M", "M", "M", "M", "F", "F", "F", "F"),
    stringsAsFactors = FALSE)
}

test_that("hp_eda_prepare derives the time and dose columns", {
  e <- hp_eda_prepare(two_dose(), lloq = 0.5)
  expect_s3_class(e, "hp_eda_data")
  r <- e$rec
  expect_equal(r$TAD, c(0, 1, 2, 0, 1, 2, 0, 1, 2, 4))
  expect_equal(r$TAFD, c(0, 1, 2, 12, 13, 14, 0, 1, 2, 4))
  expect_equal(r$DOSE, c(rep(100, 3), rep(150, 3), rep(200, 4)))
  expect_equal(r$EVENT, c("dose", "observation", "observation", "dose", "observation",
                          "observation", "dose", "observation", "observation", "observation"))
  expect_equal(r$PKOBS, c(0, 1, 1, 0, 1, 1, 0, 1, 1, 1))
  expect_equal(r$BLQ[r$EVENT == "observation"], c(0, 0, 0, 0, 1, 0, 0))
  expect_equal(r$DVNORM[2], 2 / 100)
  expect_equal(r$DVNORM[5], 3 / 150)
  expect_true(is.na(r$DVNORM[1]))
})

test_that("hp_eda_prepare builds the subject frame", {
  e <- hp_eda_prepare(two_dose(), lloq = 0.5)
  s <- e$subj
  expect_equal(nrow(s), 2)
  expect_equal(s$NOBS, c(4, 3))
  expect_equal(s$NDOSE, c(2, 1))
  expect_equal(s$NBLQ, c(0, 1))
  expect_equal(s$DOSE1, c(100, 200))
  expect_equal(s$TOTDOSE, c(250, 200))
  # exposure of the first dosing interval, BLQ samples left out
  expect_equal(s$CMAX, c(2, 5))
  expect_equal(s$TMAX, c(1, 2))
  expect_equal(s$CLAST, c(1.5, 3))
  expect_equal(s$TLAST, c(2, 4))
  expect_equal(s$AUCLAST[1], 1 + (2 - 1.5) / log(2 / 1.5), tolerance = 1e-10)
  expect_equal(s$AUCLAST[2], 5 + (5 - 3) / log(5 / 3) * 2, tolerance = 1e-10)
  expect_equal(s$CMAXD, s$CMAX / s$DOSE1)
  expect_equal(s$AUCLASTD, s$AUCLAST / s$DOSE1)
  expect_equal(hp_eda_subjects(e), s[, setdiff(names(s), c("EVENT", "PKOBS", "DVNORM", "TAD_CALC", "TAD", "TAFD", "BLQ"))])
})

test_that("the observed exposure agrees with the C++ kernel and the NCA", {
  d <- example_data()
  e <- hp_eda_prepare(d, dvid = "DVID", lloq = 0.5)
  s <- e$subj
  one <- d[d$ID == 3 & d$EVID == 0, ]
  one <- one[one$DV >= 0.5, ]
  x <- hp_exposure(one$TIME, one$DV)
  expect_equal(s$CMAX[s$ID == 3], unname(x["cmax"]))
  expect_equal(s$AUCLAST[s$ID == 3], unname(x["auclast"]), tolerance = 1e-10)
  # the exploratory exposure leaves BLQ samples out, like the "missing" rule
  r <- hp_nca(d, options = hp_nca_options(lloq = 0.5, blq_rule = "missing"), route = "oral")$parameters
  expect_equal(s$CMAX, r$CMAX)
  expect_equal(s$AUCLAST, r$AUCLAST, tolerance = 1e-10)
})

test_that("covariates are found and can be given", {
  e <- hp_eda_prepare(two_dose(), lloq = 0.5)
  expect_equal(e$info$categorical, c("WT", "SEX"))
  g <- hp_eda_prepare(example_data(), dvid = "DVID", lloq = 0.5)
  expect_setequal(g$info$continuous, c("WT", "AGE"))
  expect_equal(g$info$categorical, "SEX")
  h <- hp_eda_prepare(example_data(), dvid = "DVID", continuous = "WT", categorical = c("SEX", "DOSEGRP"))
  expect_equal(h$info$continuous, "WT")
  expect_equal(h$info$categorical, c("SEX", "DOSEGRP"))
})

test_that("filters keep doses and count what they remove", {
  d <- two_dose()
  e <- hp_eda_prepare(d, lloq = 0.5, filter = "ID == 1")
  expect_equal(e$info$individuals, 1)
  o <- hp_eda_prepare(d, obs_filter = "TIME < 2")
  expect_equal(sum(o$rec$EVENT == "dose"), 3)
  expect_equal(sum(o$rec$EVENT == "observation"), 2)
  expect_error(hp_eda_prepare(d[, -1]), "ID")
})

test_that("a BLQ flag column is used without a limit", {
  d <- two_dose()
  d$BLQ <- c(0, 0, 0, 0, 0, 1, 0, 1, 0, 0)
  e <- hp_eda_prepare(d)
  expect_equal(sum(e$rec$BLQ == 1, na.rm = TRUE), 2)
})

test_that("DVID separates PK from PD observations", {
  d <- example_data()
  pd <- d[d$EVID == 0, ]
  pd$DVID <- 2
  pd$DV <- 10 * pd$DV / (1 + pd$DV)
  e <- hp_eda_prepare(rbind(d, pd), dvid = "DVID", lloq = 0.5)
  expect_equal(e$info$pk_dvid, "1")
  expect_equal(e$info$pd_dvid, "2")
  expect_equal(sum(e$rec$PKOBS), 288)
  expect_equal(sum(e$rec$EVENT == "observation"), 576)
})

# ----------------------------------------------------------------- tables
test_that("hp_eda_tspec stores its fields", {
  s <- hp_eda_tspec("T1", "cov", "Covariates", "summary", "subjects", vars = c("wt", "age"), by = "dosegrp")
  expect_s3_class(s, "hp_eda_tspec")
  expect_equal(s$vars, c("WT", "AGE"))
  expect_equal(s$by, "DOSEGRP")
  expect_equal(s$kind, "summary")
})

test_that("summary tables agree with base R", {
  e <- hp_eda_prepare(example_data(), dvid = "DVID", lloq = 0.5)
  t <- hp_eda_table(e, hp_eda_tspec(kind = "summary", scope = "subjects", vars = "WT", by = "DOSEGRP",
                                    stats = c("n", "mean", "sd", "median", "min", "max"), digits = 6))
  expect_equal(names(t), c("Variable", "DOSEGRP", "n", "mean", "sd", "median", "min", "max"))
  expect_equal(nrow(t), 2)
  s <- e$subj
  w <- s$WT[s$DOSEGRP == 100]
  row <- t[grepl("100", t$DOSEGRP), ]
  expect_equal(as.numeric(row$n), 12)
  expect_equal(as.numeric(row$mean), signif(mean(w), 6))
  expect_equal(as.numeric(row$sd), signif(sd(w), 6))
  expect_equal(as.numeric(row$median), median(w))
  expect_equal(as.numeric(row$min), min(w))
  all_ <- hp_eda_table(e, hp_eda_tspec(kind = "summary", scope = "subjects", vars = "WT", by = "DOSEGRP",
                                       overall = TRUE, stats = "n"))
  expect_equal(all_$DOSEGRP[3], "All")
  expect_equal(all_$n[3], "24")
  expect_equal(attr(t, "note"), "24 individuals")
})

test_that("count tables give counts and percentages", {
  e <- hp_eda_prepare(example_data(), dvid = "DVID", lloq = 0.5)
  t <- hp_eda_table(e, hp_eda_tspec(kind = "count", scope = "subjects", vars = "SEX", by = "DOSEGRP",
                                    overall = TRUE))
  expect_equal(t$Level, c("0", "1"))
  expect_equal(t$All, c("12 (50.0%)", "12 (50.0%)"))
  expect_equal(t[["DOSEGRP=100"]], c("6 (50.0%)", "6 (50.0%)"))
})

test_that("listing tables respect filters and max_rows", {
  e <- hp_eda_prepare(example_data(), dvid = "DVID", lloq = 0.5)
  t <- hp_eda_table(e, hp_eda_tspec(kind = "listing", scope = "obs", vars = c("ID", "TIME", "DV", "BLQ"),
                                    filter = "BLQ == 1", max_rows = 5))
  expect_equal(nrow(t), 5)
  expect_true(all(t$BLQ == "1"))
  expect_match(attr(t, "note"), "first 5 of 78")
})

test_that("table errors are explicit", {
  e <- hp_eda_prepare(example_data(), dvid = "DVID", lloq = 0.5)
  expect_error(hp_eda_table(e, hp_eda_tspec(kind = "summary", vars = "NOPE")), "NOPE")
  expect_error(hp_eda_table(e, hp_eda_tspec(kind = "summary")), "no variables")
  expect_error(hp_eda_table(e, hp_eda_tspec(kind = "weird", vars = "WT", scope = "subjects")), "weird")
})

test_that("the data summary counts records, individuals and BLQ", {
  e <- hp_eda_prepare(example_data(), dvid = "DVID", lloq = 0.5)
  s <- hp_eda_data_summary(e)
  expect_named(s, c("Item", "Value"))
  val <- function(item) s$Value[s$Item == item]
  expect_equal(val("Records"), "312")
  expect_equal(val("Individuals"), "24")
  expect_equal(val("Dose records"), "24")
  expect_equal(val("Observations below the LLOQ"), "78")
  expect_equal(val("First-dose levels"), "100, 200")
})

test_that("the missing-value table counts missing values", {
  d <- two_dose()
  d$WT[c(2, 3)] <- NA
  m <- hp_eda_missing(hp_eda_prepare(d))
  row <- m[m$Column == "WT", ]
  expect_equal(row[["Missing records"]], "2")
  expect_equal(row[["% records"]], "20.0")
  expect_equal(row[["Missing observations"]], "2")
})

test_that("the correlation table is the Pearson matrix of the continuous covariates", {
  e <- hp_eda_prepare(example_data(), dvid = "DVID", lloq = 0.5)
  cm <- hp_eda_correlation(e)
  expect_equal(names(cm), c("", "WT", "AGE"))
  expect_equal(cm$WT[1], "1.000")
  r <- cor(e$subj$WT, e$subj$AGE)
  expect_equal(cm$AGE[1], sprintf("%.3f", r))
  expect_equal(cm$WT[2], cm$AGE[1])
})

test_that("hp_cell formats numbers compactly", {
  expect_equal(hp_cell(c(12, 3.14159265, NA), 3), c("12", "3.14", ""))
  expect_equal(hp_cell(0.000012345, 3), "1.23e-5")
  expect_equal(hp_cell(1234.5678, 6), "1234.57")
  expect_equal(hp_cell(0.5), "0.5")
  expect_equal(hp_cell(100), "100")
})

# ------------------------------------------------------------------ plots
test_that("hp_eda_spec stores its fields and upper-cases variables", {
  s <- hp_eda_spec(kind = "scatter", x = "wt", y = "cmax", scope = "subjects", smooth = "loess")
  expect_s3_class(s, "hp_eda_spec")
  expect_equal(s$x, "WT")
  expect_equal(s$y, "CMAX")
  expect_equal(s$smooth, "loess")
  expect_null(s$color)
})

test_that("axis labels use the units", {
  u <- list(time = "h", conc = "ng/mL", dose = "mg")
  expect_equal(hp_eda_label("TIME", u), "Time (h)")
  expect_equal(hp_eda_label("TAD", u), "Time after dose (h)")
  expect_equal(hp_eda_label("CMAX", u), "CMAX (ng/mL)")
  expect_equal(hp_eda_label("DVNORM", u), "Dose-normalised DV")
  expect_equal(hp_eda_label("WT", u), "WT")
  expect_equal(hp_eda_label("TIME", list(time = "", conc = "", dose = "")), "Time")
})

eda_fixture <- function() hp_eda_prepare(example_data(), dvid = "DVID", lloq = 0.5, color = "DOSEGRP")

draws <- function(e, spec, units = list(time = "h", conc = "ng/mL", dose = "mg")) {
  pages <- hp_eda_plot(e, spec, units)
  expect_type(pages, "list")
  expect_gte(length(pages), 1)
  for (p in pages) {
    expect_s3_class(p, "ggplot")
    expect_no_error(ggplot2::ggplot_build(p))
  }
  pages
}

test_that("every kind of exploratory plot can be drawn", {
  e <- eda_fixture()
  draws(e, hp_eda_spec(kind = "line", x = "TAFD", y = "DV", color = "DOSEGRP", pk_only = TRUE, log_y = TRUE))
  draws(e, hp_eda_spec(kind = "scatter", x = "TAFD", y = "DV", pk_only = TRUE, smooth = "loess"))
  draws(e, hp_eda_spec(kind = "scatter", x = "WT", y = "CMAX", scope = "subjects", smooth = "linear",
                       color = "DOSEGRP"))
  draws(e, hp_eda_spec(kind = "scatter", x = "DOSE1", y = "AUCLAST", scope = "subjects", smooth = "power",
                       log_x = TRUE, log_y = TRUE, identity = TRUE))
  draws(e, hp_eda_spec(kind = "summary", x = "TAFD", y = "DV", color = "DOSEGRP", centre = "median",
                       spread = "iqr", pk_only = TRUE, log_y = TRUE))
  draws(e, hp_eda_spec(kind = "summary", x = "TAFD", y = "DVNORM", color = "DOSEGRP", centre = "geomean",
                       spread = "sd", bins = 6, pk_only = TRUE))
  draws(e, hp_eda_spec(kind = "box", x = "DOSEGRP", y = "CMAX", scope = "subjects"))
  draws(e, hp_eda_spec(kind = "box", x = "DOSEGRP", y = c("CMAX", "AUCLAST"), scope = "subjects", free = "y"))
  draws(e, hp_eda_spec(kind = "histogram", x = "WT", scope = "subjects", hist_bins = 6))
  draws(e, hp_eda_spec(kind = "histogram", x = "CMAX", scope = "subjects", log_x = TRUE, color = "DOSEGRP"))
  draws(e, hp_eda_spec(kind = "bar", x = "SEX", scope = "subjects", color = "DOSEGRP"))
  draws(e, hp_eda_spec(kind = "bar", x = "DOSEGRP", y = "CMAX", scope = "subjects"))
  draws(e, hp_eda_spec(kind = "corr", x = c("WT", "AGE", "CMAX"), scope = "subjects"))
  draws(e, hp_eda_spec(kind = "pairs", x = c("WT", "AGE", "CMAX", "AUCLAST"), scope = "subjects"))
})

test_that("plot options are applied", {
  e <- eda_fixture()
  p <- draws(e, hp_eda_spec(title = "My title", kind = "line", x = "TAFD", y = "DV", pk_only = TRUE,
                            hlines = 0.5, xlab = "hours", ylab = "conc", xlim = c(0, 24),
                            ylim = c(NA, 10), legend = FALSE))[[1]]
  expect_equal(p$labels$title, "My title")
  expect_equal(p$labels$x, "hours")
  expect_equal(p$labels$y, "conc")
  expect_equal(p$theme$legend.position, "none")
  expect_match(p$labels$caption, "hyperutils eda")
})

test_that("facets are split over pages", {
  e <- eda_fixture()
  pages <- draws(e, hp_eda_spec(kind = "line", x = "TAFD", y = "DV", facet = "ID", pk_only = TRUE,
                                ncol = 4, nrow = 3))
  expect_length(pages, 2)
  one <- draws(e, hp_eda_spec(kind = "line", x = "TAFD", y = "DV", facet = "ID", pk_only = TRUE,
                              ncol = 4, nrow = 3, max_pages = 1))
  expect_length(one, 1)
})

test_that("a filter in the specification selects the rows", {
  e <- eda_fixture()
  p <- draws(e, hp_eda_spec(kind = "line", x = "TAFD", y = "DV", pk_only = TRUE, filter = "DOSEGRP == 200"))[[1]]
  expect_equal(length(unique(p$data$.BY)), 12)
})

test_that("a PK-PD plot pairs the response with the concentration", {
  d <- example_data()
  pd <- d[d$EVID == 0, ]
  pd$DVID <- 2
  pd$DV <- 10 * pd$DV / (1 + pd$DV)
  e <- hp_eda_prepare(rbind(d, pd), dvid = "DVID", lloq = 0.5)
  p <- draws(e, hp_eda_spec(kind = "pkpd", x = "DV", y = "DV", color = "DOSEGRP"))
  expect_s3_class(p[[1]], "ggplot")
  expect_error(hp_eda_plot(eda_fixture(), hp_eda_spec(kind = "pkpd", x = "DV", y = "DV")), "DVID")
})

test_that("plot errors are explicit", {
  e <- eda_fixture()
  expect_error(hp_eda_plot(e, hp_eda_spec(kind = "scatter", x = "NOPE", y = "DV")), "NOPE")
  expect_error(hp_eda_plot(e, hp_eda_spec(kind = "weird", x = "TIME", y = "DV")), "weird")
  expect_error(hp_eda_plot(e, hp_eda_spec(kind = "corr", x = "WT", scope = "subjects")), "two")
  expect_error(hp_eda_plot(e, hp_eda_spec(kind = "scatter", x = "TIME", y = "DV", filter = "DV > 1e9")),
               "no rows")
})

test_that("pages can be saved", {
  e <- eda_fixture()
  pages <- hp_eda_plot(e, hp_eda_spec(kind = "line", x = "TAFD", y = "DV", facet = "ID", pk_only = TRUE,
                                      ncol = 4, nrow = 3))
  dir <- file.path(tempdir(), "t-eda-pages")
  f <- hp_save_plots(pages, "ind", dir, 6, 4, dpi = 50)
  expect_equal(basename(f), c("ind_01.png", "ind_02.png"))
  expect_true(all(file.exists(f)))
  unlink(dir, recursive = TRUE)
})
