test_that("the version is a string and matches DESCRIPTION", {
  v <- hp_version()
  expect_type(v, "character")
  expect_equal(v, as.character(utils::packageVersion("hyperutils")))
})

test_that("the example file exists and is a NONMEM-style dataset", {
  f <- hp_example_file()
  expect_true(file.exists(f))
  d <- hp_read_data(f)
  expect_equal(names(d)[1:5], c("ID", "TIME", "AMT", "DV", "EVID"))
  expect_equal(length(unique(d$ID)), 24)
  expect_equal(sum(d$EVID == 1), 24)
  expect_true(all(is.na(d$DV[d$EVID == 1])))
})

# ---------------------------------------------------------------- reading
write_tmp <- function(lines, ext = ".txt") {
  f <- tempfile(fileext = ext)
  writeLines(lines, f)
  f
}

test_that("hp_read_data finds the delimiter", {
  for (sep in c(",", "\t", ";", " ")) {
    f <- write_tmp(c(paste("id", "time", "dv", sep = sep), paste(1, 0, 2.5, sep = sep),
                     paste(1, 1, 1.5, sep = sep)))
    d <- hp_read_data(f)
    expect_equal(names(d), c("ID", "TIME", "DV"))
    expect_equal(d$DV, c(2.5, 1.5))
    unlink(f)
  }
})

test_that("hp_read_data handles missing codes, comments and Fortran exponents", {
  f <- write_tmp(c("#ID,TIME,DV", "1,0,.", "# a comment row", "1,1,1.5D-3", "1,2,NA"))
  d <- hp_read_data(f)
  expect_equal(names(d), c("ID", "TIME", "DV"))
  expect_equal(nrow(d), 3)
  expect_true(is.numeric(d$DV))
  expect_true(is.na(d$DV[1]))
  expect_equal(d$DV[2], 1.5e-3)
  expect_true(is.na(d$DV[3]))
  unlink(f)
})

test_that("hp_read_data options: skip, upper, header, delimiter, na.strings", {
  f <- write_tmp(c("Title line", "id;wt", "1;70.5", "2;-99"))
  d <- hp_read_data(f, skip = 1)
  expect_equal(names(d), c("ID", "WT"))
  d2 <- hp_read_data(f, skip = 1, upper = FALSE)
  expect_equal(names(d2), c("id", "wt"))
  d3 <- hp_read_data(f, skip = 1, na.strings = c("-99", ".", ""))
  expect_true(is.na(d3$WT[2]))
  d4 <- hp_read_data(f, skip = 2, header = FALSE, delimiter = ";")
  expect_equal(names(d4), c("A", "B"))
  unlink(f)
})

test_that("hp_read_data reads NONMEM table files", {
  f <- write_tmp(c("TABLE NO.     1", " ID TIME DV", " 1 0 0", " 1 1 2.5",
                   "TABLE NO.     2", " ID TIME DV", " 2 0 0", " 2 1 3.1"))
  d <- hp_read_data(f)
  expect_equal(nrow(d), 4)
  expect_true("TABLE_NO" %in% names(d))
  expect_equal(d$TABLE_NO, c(1, 1, 2, 2))
  unlink(f)
  f1 <- write_tmp(c("TABLE NO.     1", " ID TIME DV", " 1 0 0", " 1 1 2.5"))
  expect_false("TABLE_NO" %in% names(hp_read_data(f1)))
  unlink(f1)
})

test_that("hp_read_data fails on a missing file", {
  expect_error(hp_read_data(file.path(tempdir(), "no-such-file.csv")))
})

# ------------------------------------------------------------- data steps
test_that("hp_mutate runs statements in order and upper-cases names", {
  d <- data.frame(WT = c(50, 70, 95), DOSE = c(100, 100, 200))
  m <- hp_mutate(d, "lndv = log(WT); X = DOSE / WT; Y = X * 2")
  expect_equal(m$LNDV, log(d$WT))
  expect_equal(m$Y, 2 * d$DOSE / d$WT)
  expect_equal(hp_mutate(d, c("A = 1", "B = A + WT"))$B, 1 + d$WT)
  expect_equal(hp_mutate(d, NULL), d)
  expect_equal(hp_mutate(d, character()), d)
  expect_equal(hp_mutate(d, "K = 5")$K, rep(5, 3))
})

test_that("hp_mutate does not split on semicolons inside quotes", {
  d <- data.frame(A = 1:2)
  m <- hp_mutate(d, "B = paste('x;y', A)")
  expect_equal(m$B, c("x;y 1", "x;y 2"))
})

test_that("hp_mutate reports a bad statement by name", {
  d <- data.frame(A = 1:3)
  expect_error(hp_mutate(d, "B = log(NOPE)"), "NOPE")
  expect_error(hp_mutate(d, "B = 1:2"), "length|rows")
})

test_that("hp_case_when takes the first true test and NA when none", {
  w <- c(50, 70, 95)
  expect_equal(hp_case_when(w < 60 ~ "low", w < 85 ~ "mid", TRUE ~ "high"),
               c("low", "mid", "high"))
  expect_equal(hp_case_when(w < 60 ~ 1, w < 85 ~ 2), c(1, 2, NA))
  expect_equal(hp_case_when(c(NA, 1) > 0 ~ "pos", TRUE ~ "other"), c("other", "pos"))
  expect_error(hp_case_when(w < 60, TRUE), "test ~ value")
  d <- data.frame(WT = w)
  expect_equal(hp_mutate(d, "G = case_when(WT < 60 ~ 1, WT < 85 ~ 2, TRUE ~ 3)")$G, c(1, 2, 3))
})

test_that("hp_filter drops NA rows and accepts empty input", {
  d <- data.frame(ID = 1:4, WT = c(50, 70, NA, 95), SEX = c(1, 0, 1, 1))
  expect_equal(hp_filter(d, "WT > 60")$ID, c(2, 4))
  expect_equal(hp_filter(d, "WT > 40 & SEX == 1")$ID, c(1, 4))
  expect_equal(hp_filter(d, "ID %in% c(2, 4)")$ID, c(2, 4))
  expect_equal(hp_filter(d, NULL), d)
  expect_equal(hp_filter(d, ""), d)
  expect_equal(rownames(hp_filter(d, "ID > 2")), c("1", "2"))
  expect_error(hp_filter(d, "NOPE > 1"), "NOPE")
})

# ------------------------------------------------------------ statistics
test_that("hp_describe agrees with base R", {
  x <- c(1.2, 2.5, 3.1, 4.8, 5.5, NA)
  s <- hp_describe(x)
  y <- x[!is.na(x)]
  expect_equal(unname(s["n"]), 5)
  expect_equal(unname(s["nmiss"]), 1)
  expect_equal(unname(s["mean"]), mean(y))
  expect_equal(unname(s["sd"]), sd(y))
  expect_equal(unname(s["cv"]), 100 * sd(y) / mean(y))
  expect_equal(unname(s["sem"]), sd(y) / sqrt(5))
  expect_equal(unname(s["median"]), median(y))
  expect_equal(unname(s[c("min", "max")]), range(y))
  expect_equal(unname(s[c("q1", "q3")]), unname(quantile(y, c(0.25, 0.75), type = 7)))
  expect_equal(unname(s[c("p5", "p95")]), unname(quantile(y, c(0.05, 0.95), type = 7)))
  expect_equal(unname(s["geomean"]), exp(mean(log(y))))
  expect_equal(unname(s["gcv"]), 100 * sqrt(exp(sd(log(y))^2) - 1))
  ci <- mean(y) + c(-1, 1) * qt(0.975, 4) * sd(y) / sqrt(5)
  expect_equal(unname(s[c("ci_lo", "ci_hi")]), ci)
})

test_that("hp_describe handles degenerate input", {
  s <- hp_describe(c(2, 2, 2))
  expect_equal(unname(s["sd"]), 0)
  expect_true(is.na(hp_describe(c(-1, 2, 3))["geomean"]))
  one <- hp_describe(5)
  expect_equal(unname(one["n"]), 1)
  expect_true(is.na(one["sd"]))
  none <- hp_describe(numeric())
  expect_equal(unname(none["n"]), 0)
  expect_true(is.na(none["mean"]))
})

test_that("hp_centre_spread gives the centre and error bars", {
  x <- c(1.2, 2.5, 3.1, 4.8, 5.5)
  m <- hp_centre_spread(x)
  expect_equal(unname(m), c(mean(x), mean(x) - sd(x), mean(x) + sd(x)))
  expect_equal(names(m), c("centre", "lower", "upper"))
  q <- hp_centre_spread(x, "median", "iqr")
  expect_equal(unname(q), c(median(x), unname(quantile(x, c(0.25, 0.75), type = 7))))
  g <- hp_centre_spread(x, "geomean", "ci95")
  expect_equal(unname(g[1]), exp(mean(log(x))))
  expect_true(g[2] > 0 && g[2] < g[1] && g[3] > g[1])
  r <- hp_centre_spread(x, "mean", "range")
  expect_equal(unname(r[2:3]), range(x))
  expect_true(all(is.na(hp_centre_spread(NA_real_)))) 
  expect_true(all(is.na(hp_centre_spread(x, "mean", "none")[2:3])))
})

# -------------------------------------------------------------- output
test_that("hp_theme is a ggplot2 theme", {
  expect_s3_class(hp_theme(), "theme")
  expect_s3_class(hp_theme(14, "none"), "theme")
})

test_that("hp_save_table writes a faithful CSV in the folder given", {
  dir <- file.path(tempdir(), "t-save-table")
  x <- data.frame(ID = c("a", "b,c"), V = c(1 / 3, NA), N = 1:2, stringsAsFactors = FALSE)
  f <- hp_save_table(x, "x.csv", dir)
  expect_true(file.exists(f))
  expect_equal(dirname(normalizePath(f)), normalizePath(dir))
  back <- utils::read.csv(f, stringsAsFactors = FALSE)
  expect_equal(back$ID, x$ID)
  expect_identical(back$V, x$V)
  expect_equal(back$N, x$N)
  expect_match(readLines(f)[3], "^\"b,c\",,2$")
  f2 <- hp_save_table(x, "y.csv", dir, digits = 3)
  expect_match(readLines(f2)[2], "^a,0.333,1$")
  unlink(dir, recursive = TRUE)
})

test_that("saving functions need a folder and create it", {
  p <- ggplot2::ggplot(data.frame(x = 1:3, y = c(2, 4, 3)), ggplot2::aes(x, y)) +
    ggplot2::geom_line()
  expect_error(hp_save_plot(p, "p"))
  expect_error(hp_save_table(data.frame(a = 1), "a.csv"))
  dir <- file.path(tempdir(), "t-save-plot", "nested")
  f <- hp_save_plot(p, "p", dir, width = 3, height = 2, dpi = 50)
  expect_true(file.exists(f))
  expect_match(basename(f), "^p\\.png$")
  unlink(file.path(tempdir(), "t-save-plot"), recursive = TRUE)
})

test_that("hp_save_plots numbers the pages", {
  p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) + ggplot2::geom_point()
  dir <- file.path(tempdir(), "t-save-plots")
  one <- hp_save_plots(list(p), "one", dir, 3, 2, 50)
  expect_equal(basename(one), "one.png")
  two <- hp_save_plots(list(p, p), "two", dir, 3, 2, 50)
  expect_equal(basename(two), c("two_01.png", "two_02.png"))
  expect_true(all(file.exists(c(one, two))))
  unlink(dir, recursive = TRUE)
})

test_that("hp_save_plot writes SVG when svglite is available", {
  skip_if_not_installed("svglite")
  p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) + ggplot2::geom_point()
  dir <- file.path(tempdir(), "t-svg")
  f <- hp_save_plot(p, "s", dir, 3, 2, 50, formats = c("png", "svg"))
  expect_equal(sort(basename(f)), c("s.png", "s.svg"))
  expect_true(all(file.exists(f)))
  unlink(dir, recursive = TRUE)
})
