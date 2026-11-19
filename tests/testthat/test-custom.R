cust_obs <- function() hp_custom_prep(example_data(), filter = "EVID == 0")
cust_subj <- function() hp_custom_prep(example_data(), distinct = "ID")

# ------------------------------------------------------------------ prep
test_that("hp_custom_prep runs the steps in the documented order", {
  d <- example_data()
  p <- hp_custom_prep(d, mutate = "X = DV * 2", filter = "EVID == 0", distinct = "ID",
                      sort = "-WT", head = 3, select = c("ID", "WT", "X"))
  expect_equal(names(p), c("ID", "WT", "X"))
  expect_equal(nrow(p), 3)
  expect_true(all(diff(p$WT) <= 0))
  # the filter runs after the mutate, so it can use the new column
  q <- hp_custom_prep(d, mutate = "HIGH = AMT >= 200", filter = "HIGH")
  expect_true(all(q$AMT >= 200))
  expect_equal(rownames(q), as.character(seq_len(nrow(q))))
})

test_that("hp_custom_prep distinct keeps the first row of each combination", {
  d <- example_data()
  p <- hp_custom_prep(d, distinct = c("ID"))
  expect_equal(nrow(p), 24)
  expect_equal(p$TIME, rep(0, 24))
  p2 <- hp_custom_prep(d, distinct = c("DOSEGRP", "SEX"))
  expect_equal(nrow(p2), 4)
})

test_that("hp_custom_prep sorts up and down on several columns", {
  d <- data.frame(A = c(2, 1, 2, 1), B = c("x", "y", "z", "w"))
  expect_equal(hp_custom_prep(d, sort = c("A", "B"))$B, c("w", "y", "x", "z"))
  expect_equal(hp_custom_prep(d, sort = c("-A", "B"))$B, c("x", "z", "w", "y"))
  expect_equal(hp_custom_prep(d, sort = "-B")$B, c("z", "y", "x", "w"))
})

test_that("hp_custom_prep reports unknown columns", {
  d <- data.frame(A = 1:3)
  expect_error(hp_custom_prep(d, select = "B"), "select `B`")
  expect_error(hp_custom_prep(d, sort = "B"), "sort `B`")
  expect_error(hp_custom_prep(d, distinct = "B"), "distinct `B`")
})

test_that("hp_custom_prep without steps returns the data", {
  d <- data.frame(A = 1:3, B = c("a", "b", "c"))
  expect_equal(hp_custom_prep(d), d)
})

# ------------------------------------------------------------------ plots
built <- function(p) {
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
  p
}

test_that("every layer can be drawn", {
  obs <- cust_obs(); subj <- cust_subj()
  geoms_xy <- c("point", "jitter", "line", "path", "step")
  for (g in geoms_xy)
    built(hp_custom_plot(obs[obs$ID <= 4, ], list(geoms = g, x = "TIME", y = "DV", group = "ID")))
  built(hp_custom_plot(subj, list(geoms = c("point", "smooth"), x = "WT", y = "AGE")))
  built(hp_custom_plot(subj, list(geoms = c("point", "smooth"), x = "WT", y = "AGE",
                                  smooth_loess = TRUE, smooth_se = TRUE, span = 0.9)))
  built(hp_custom_plot(subj, list(geoms = "boxplot", x = "factor(SEX)", y = "WT", outliers = FALSE)))
  built(hp_custom_plot(subj, list(geoms = "violin", x = "factor(SEX)", y = "WT")))
  built(hp_custom_plot(subj, list(geoms = "bar", x = "factor(SEX)")))
  built(hp_custom_plot(subj, list(geoms = "bar", x = "factor(SEX)", y = "WT", summary_fun = "median")))
  built(hp_custom_plot(data.frame(g = c("a", "b"), v = c(3, 5)), list(geoms = "col", x = "g", y = "v")))
  built(hp_custom_plot(subj, list(geoms = "histogram", x = "WT", nbins = 5)))
  built(hp_custom_plot(subj, list(geoms = "histogram", x = "WT", binwidth = 5)))
  built(hp_custom_plot(subj, list(geoms = "density", x = "WT", fill = "factor(SEX)")))
  built(hp_custom_plot(subj, list(geoms = "area", x = "AGE", y = "WT")))
  built(hp_custom_plot(subj, list(geoms = "ecdf", x = "WT")))
  built(hp_custom_plot(subj, list(geoms = "qq", x = "WT")))
  for (e in c("sd", "sem", "ci95"))
    built(hp_custom_plot(obs, list(geoms = c("errorbar", "pointrange"), x = "TIME", y = "DV", summary_error = e)))
  cells <- aggregate(DV ~ TIME + DOSEGRP, data = obs, FUN = mean)
  built(hp_custom_plot(cells, list(geoms = "tile", x = "factor(TIME)", y = "factor(DOSEGRP)", fill = "DV")))
  built(hp_custom_plot(subj[1:5, ], list(geoms = c("point", "text"), x = "WT", y = "AGE", label = "ID")))
})

test_that("an unknown layer is an error", {
  expect_error(hp_custom_plot(cust_subj(), list(geoms = "sculpture", x = "WT", y = "AGE")), "sculpture")
})

test_that("aesthetics accept expressions and make discrete legends", {
  subj <- cust_subj()
  p <- built(hp_custom_plot(subj, list(geoms = "point", x = "WT", y = "log10(AGE)",
                                       color = "factor(DOSEGRP)", shape = "SEX", size = "AGE")))
  expect_equal(p$data$.y, log10(subj$AGE))
  expect_s3_class(p$data$.color, "factor")
  expect_s3_class(p$data$.shape, "factor")
  expect_error(hp_custom_plot(subj, list(geoms = "point", x = "WT", y = "NOPE")), "NOPE")
})

test_that("scales, limits and reference lines are applied", {
  obs <- cust_obs()
  p <- built(hp_custom_plot(obs, list(geoms = "point", x = "TIME", y = "DV", y_log = TRUE, x_log = FALSE,
                                      xlim = c(0, 24), hline = 0.5, vline = c(4, 8), identity = TRUE,
                                      abline = list(c(0, 1)), ref_color = "red", x_breaks = c(0, 12, 24))))
  kinds <- vapply(p$layers, function(l) class(l$geom)[1], "")
  expect_equal(sum(kinds == "GeomHline"), 1)
  expect_equal(sum(kinds == "GeomVline"), 2)
  expect_equal(sum(kinds == "GeomAbline"), 2)
  expect_true(inherits(p$scales$get_scales("y")$trans, "trans") ||
                inherits(p$scales$get_scales("y"), "ScaleContinuousPosition"))
})

test_that("titles and legends are set from the list", {
  obs <- cust_obs()
  p <- built(hp_custom_plot(obs, list(geoms = "point", x = "TIME", y = "DV", color = "factor(DOSEGRP)",
                                      title = "T", subtitle = "S", caption = "C", x_title = "X",
                                      y_title = "Y", legend_title = "L", legend = "bottom")))
  expect_equal(p$labels$title, "T")
  expect_equal(p$labels$subtitle, "S")
  expect_equal(p$labels$caption, "C")
  expect_equal(p$labels$x, "X")
  expect_equal(p$labels$y, "Y")
  expect_equal(p$theme$legend.position, "bottom")
})

test_that("facets split the data", {
  obs <- cust_obs()
  w <- built(hp_custom_plot(obs, list(geoms = "point", x = "TIME", y = "DV", facet_wrap = "DOSEGRP",
                                      facet_ncol = 2, labels_both = TRUE, scales = "free_y")))
  expect_s3_class(w$facet, "FacetWrap")
  g <- built(hp_custom_plot(obs, list(geoms = "point", x = "TIME", y = "DV", facet_rows = "SEX",
                                      facet_cols = "DOSEGRP")))
  expect_s3_class(g$facet, "FacetGrid")
})

test_that("themes and palettes are accepted", {
  obs <- cust_obs()
  for (th in c("classic", "minimal", "bw", "light", "grey", "void"))
    built(hp_custom_plot(obs, list(geoms = "point", x = "TIME", y = "DV", theme = th)))
  for (pa in c("okabe-ito", "tableau", "set1", "dark2", "pastel", "grey", "viridis"))
    built(hp_custom_plot(obs, list(geoms = "point", x = "TIME", y = "DV", color = "factor(DOSEGRP)", palette = pa)))
})

test_that("a plot is repeatable", {
  subj <- cust_subj()
  a <- hp_custom_plot(subj, list(geoms = "jitter", x = "factor(SEX)", y = "WT"), seed = 7)
  b <- hp_custom_plot(subj, list(geoms = "jitter", x = "factor(SEX)", y = "WT"), seed = 7)
  expect_equal(ggplot2::layer_data(a), ggplot2::layer_data(b))
})

test_that("a custom plot can be extended and saved", {
  obs <- cust_obs()
  p <- hp_custom_plot(obs, list(geoms = "point", x = "TIME", y = "DV")) + ggplot2::labs(caption = "extra")
  expect_equal(p$labels$caption, "extra")
  dir <- file.path(tempdir(), "t-custom-save")
  f <- hp_save_plot(p, "c", dir, 4, 3, dpi = 50)
  expect_true(file.exists(f))
  unlink(dir, recursive = TRUE)
})

# ----------------------------------------------------------------- tables
test_that("listing keeps the rows and formats numbers", {
  d <- data.frame(ID = 1:2, V = c(1 / 3, 1234.5678), S = c("a", NA))
  t <- hp_custom_table(d, list(kind = "listing"), digits = 3)
  expect_equal(t$ID, c("1", "2"))
  expect_equal(t$V, c("0.333", "1235"))
  expect_equal(t$S, c("a", ""))
})

test_that("summary tables agree with base R", {
  s <- cust_subj()
  t <- hp_custom_table(s, list(kind = "summary", vars = "WT", by = "SEX", total = TRUE,
                               stats = c("n", "mean", "sd", "median", "min", "max"), digits = 6))
  expect_equal(names(t)[1:2], c("Variable", "SEX"))
  expect_equal(nrow(t), 3)
  w0 <- s$WT[s$SEX == 0]
  r0 <- t[t$SEX == "0", ]
  expect_equal(as.numeric(r0$N), length(w0))
  expect_equal(as.numeric(r0$Mean), signif(mean(w0), 6))
  expect_equal(as.numeric(r0$SD), signif(sd(w0), 6))
  expect_equal(as.numeric(r0$Median), median(w0))
  expect_equal(t$SEX[3], "All")
  expect_equal(as.numeric(t$N[3]), 24)
  # without vars, all numeric columns are described
  all_ <- hp_custom_table(s[, c("WT", "AGE")], list(kind = "summary", stats = "n"))
  expect_equal(all_$Variable, c("WT", "AGE"))
  expect_error(hp_custom_table(s, list(kind = "summary", vars = "NOPE")), "NOPE")
  expect_error(hp_custom_table(data.frame(a = "x"), list(kind = "summary")), "numeric")
})

test_that("demographics tables have a column per group", {
  s <- cust_subj()
  t <- hp_custom_table(s, list(kind = "demographics", vars = c("WT", "SEX"), by = "DOSEGRP", total = TRUE))
  expect_equal(names(t), c("Characteristic", "Overall (N=24)", "100 (N=12)", "200 (N=12)"))
  expect_equal(t$Characteristic[1:3], c("WT", "Mean (SD)", "Median [Min, Max]"))
  expect_match(t$`Overall (N=24)`[2], sprintf("^%s \\(", format(signif(mean(s$WT), 3), nsmall = 1)))
  sex <- t[t$Characteristic == "0", ]
  expect_equal(sex$`Overall (N=24)`, "12 (50.0%)")
  no_total <- hp_custom_table(s, list(kind = "demographics", vars = "WT", by = "DOSEGRP", total = FALSE))
  expect_equal(names(no_total), c("Characteristic", "100 (N=12)", "200 (N=12)"))
  expect_error(hp_custom_table(s, list(kind = "demographics")), "needs variables")
})

test_that("demographics tables show missing values", {
  s <- cust_subj()
  s$WT[1:2] <- NA
  t <- hp_custom_table(s, list(kind = "demographics", vars = "WT", total = TRUE))
  expect_true("Missing" %in% t$Characteristic)
  expect_equal(t$`Overall (N=24)`[t$Characteristic == "Missing"], "2")
})

test_that("frequency tables count combinations", {
  s <- cust_subj()
  one <- hp_custom_table(s, list(kind = "frequency", rows = "SEX"))
  expect_equal(one$N, c("12", "12", "24"))
  expect_equal(one$`%`, c("50.0", "50.0", "100.0"))
  two <- hp_custom_table(s, list(kind = "frequency", rows = "SEX", columns = "DOSEGRP", percent = "column"))
  expect_equal(two[["DOSEGRP=100"]], c("6 (50.0%)", "6 (50.0%)", "12 (100.0%)"))
  row <- hp_custom_table(s, list(kind = "frequency", rows = "SEX", columns = "DOSEGRP", percent = "row"))
  expect_equal(row[["DOSEGRP=100"]][1:2], c("6 (50.0%)", "6 (50.0%)"))
  none <- hp_custom_table(s, list(kind = "frequency", rows = "SEX", columns = "DOSEGRP", total = FALSE))
  expect_equal(none[["DOSEGRP=200"]], c("6", "6"))
  expect_error(hp_custom_table(s, list(kind = "frequency")), "rows")
  expect_error(hp_custom_table(s, list(kind = "frequency", rows = "NOPE")), "NOPE")
})

test_that("pivot tables summarise a value by rows and columns", {
  s <- cust_subj()
  t <- hp_custom_table(s, list(kind = "pivot", rows = "SEX", columns = "DOSEGRP", value = "WT",
                               fun = "mean", total = TRUE, digits = 6))
  expect_equal(names(t), c("SEX", "100", "200", "All"))
  expect_equal(as.numeric(t[["100"]][t$SEX == "0"]), signif(mean(s$WT[s$SEX == 0 & s$DOSEGRP == 100]), 6))
  expect_equal(t$SEX[3], "All")
  counts <- hp_custom_table(s, list(kind = "pivot", rows = "SEX", fun = "n"))
  expect_equal(counts$N, c("12", "12"))
  expect_error(hp_custom_table(s, list(kind = "pivot", rows = "SEX")), "value")
  s$ARM <- ifelse(s$DOSEGRP == 100, "low", "high")
  expect_error(hp_custom_table(s, list(kind = "pivot", rows = "SEX", value = "ARM")), "not numeric")
})

test_that("correlation tables give the matrix as text", {
  s <- cust_subj()
  pe <- hp_custom_table(s, list(kind = "correlation", vars = c("WT", "AGE", "AMT")))
  expect_equal(names(pe), c("", "WT", "AGE", "AMT"))
  expect_equal(pe$AGE[1], sprintf("%.3f", cor(s$WT, s$AGE)))
  sp <- hp_custom_table(s, list(kind = "correlation", vars = c("WT", "AGE"), spearman = TRUE))
  expect_equal(sp$AGE[1], sprintf("%.3f", cor(s$WT, s$AGE, method = "spearman")))
  expect_error(hp_custom_table(s, list(kind = "correlation", vars = "WT")), "at least two")
})

test_that("finishing options rename, hide, replace and transpose", {
  s <- cust_subj()
  base <- list(kind = "summary", vars = c("WT", "AGE"), by = "SEX", stats = c("n", "mean"), digits = 3)
  r <- hp_custom_table(s, c(base, list(rename = list(c("N", "Subjects"), c("SEX", "Sex")))))
  expect_true(all(c("Subjects", "Sex") %in% names(r)))
  h <- hp_custom_table(s, c(base, list(hide = "Variable")))
  expect_false("Variable" %in% names(h))
  d <- hp_custom_table(s, c(base, list(decimals = 1)))
  expect_match(d$Mean[1], "^[0-9]+\\.[0-9]$")
  na <- hp_custom_table(data.frame(a = c("x", NA), b = c("1", "2")), list(kind = "listing", na_text = "-"))
  expect_equal(na$a, c("x", "-"))
  tr <- hp_custom_table(s, list(kind = "summary", vars = "WT", stats = c("n", "mean"), transpose = TRUE))
  expect_equal(tr$Variable, c("N", "Mean"))
  expect_equal(names(tr)[2], "WT")
})

test_that("an unknown table kind is an error", {
  expect_error(hp_custom_table(cust_subj(), list(kind = "mosaic")), "mosaic")
})

test_that("a table can be saved and read back", {
  s <- cust_subj()
  t <- hp_custom_table(s, list(kind = "demographics", vars = c("WT", "SEX"), by = "DOSEGRP", total = TRUE))
  dir <- file.path(tempdir(), "t-custom-table")
  f <- hp_save_table(t, "demo.csv", dir)
  back <- utils::read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  expect_equal(back, t, ignore_attr = TRUE)
  unlink(dir, recursive = TRUE)
})
