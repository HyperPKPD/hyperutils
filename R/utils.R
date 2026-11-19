#' hyperutils: Fast Exploratory Analysis, Non-Compartmental Analysis and Custom Outputs for Clinical Pharmacology
#'
#' A toolbox for the day-to-day data work of clinical pharmacology and
#' pharmacometrics: read NONMEM-style datasets, explore them with
#' publication-quality plots and tables, run non-compartmental analysis (NCA) of
#' every profile, and build custom ggplot2 figures and summary tables, with the
#' numerical core written in C++ ('Rcpp') for speed.
#'
#' @section What the package does:
#' - **Read and prepare data.** [hp_read_data()] reads comma, tab,
#'   semicolon or white-space delimited files and NONMEM table files;
#'   [hp_mutate()], [hp_filter()] and [hp_case_when()] run data steps written as
#'   R expressions; [hp_custom_prep()] chains mutate, filter, distinct, sort,
#'   head and select.
#' - **Explore.** [hp_eda_prepare()] derives time after dose, dose-normalised
#'   concentrations, BLQ flags and observed exposure per subject.
#'   [hp_eda_plot()] draws concentration-time, covariate, distribution, summary
#'   and correlation figures from a specification built with [hp_eda_spec()];
#'   [hp_eda_table()] and friends build descriptive tables.
#' - **Non-compartmental analysis.** [hp_nca()] analyses every dosing
#'   interval of every subject (oral, intravenous bolus and infusion, single
#'   dose and steady state) with WinNonlin-style rules for areas, BLQ data and
#'   the terminal slope, and returns parameter, summary, terminal-slope and
#'   profile tables. [hp_nca_plots()] draws the standard NCA figures.
#' - **Custom outputs.** [hp_custom_plot()] and [hp_custom_table()] turn a
#'   plain list describing a plot or a table into a ggplot2 figure or a
#'   summary table (demographics, frequency, pivot, correlation).
#' - **Fast.** Areas, partial areas, moments, terminal-phase regressions and
#'   exposure measures are computed in C++ ([hp_curve_area()],
#'   [hp_lz_windows()], [hp_exposure()] and others), which you can also call
#'   directly.
#'
#' @section Typical workflow:
#' 1. `d <- hp_read_data("study.csv")`
#' 2. `e <- hp_eda_prepare(d, lloq = 0.5)` then plots and tables with
#'    [hp_eda_plot()] and [hp_eda_table()].
#' 3. `res <- hp_nca(d, route = "oral", lloq = 0.5)` then
#'    `res$parameters`, `res$summary` and [hp_nca_plots()].
#' 4. Save with [hp_save_table()] and [hp_save_plots()].
#'
#' The vignettes walk through each step: `browseVignettes("hyperutils")`.
#'
#' @section Conventions:
#' Column names are upper-cased on reading, so `WT`, `wt` and `Wt` are the same
#' column. Missing values in a dataset are `.`, an empty field or `NA`. Times
#' are in the unit of the data, concentrations in the unit of the data; units
#' shown with results are labels you supply, nothing is converted.
#' Functions never write files unless you give them a folder.
#'
#' @author Obinna Obianom \email{idonshayo@gmail.com}
#' @seealso Useful links: <https://github.com/hyperpkpd/hyperutils> and the
#'   issue tracker <https://github.com/hyperpkpd/hyperutils/issues>.
#' @keywords internal
#' @useDynLib hyperutils, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @import ggplot2
#' @importFrom stats aggregate approx complete.cases cor lm loess median predict qt quantile sd setNames
#' @importFrom utils read.csv read.table write.csv head tail
"_PACKAGE"

# ---------------------------------------------------------------------------
# small internal helpers
# ---------------------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L || (length(a) == 1L && !is.environment(a) && !is.function(a) && is.na(a))) b else a

hp_num <- function(x) suppressWarnings(as.numeric(x))

# a number as text, shortest form
hp_fmt <- function(x, digits = 6) {
  ifelse(is.na(x), "", formatC(signif(x, digits), digits = digits, format = "g"))
}

hp_stop <- function(...) stop(paste0(...), call. = FALSE)

#' Version of the installed hyperutils package
#'
#' Returns the version of the package as text. Useful for the header of a
#' reproducible analysis script and for checking that a script runs against
#' the version it was written for.
#'
#' @return A character string such as `"0.1.0"`.
#' @examples
#' hp_version()
#'
#' # a script can refuse to run with an older package
#' stopifnot(utils::packageVersion("hyperutils") >= "0.1.0")
#' @export
hp_version <- function() as.character(utils::packageVersion("hyperutils"))

# ---------------------------------------------------------------------------
# reading a dataset
# ---------------------------------------------------------------------------

#' Read a NONMEM-style dataset or any delimited table
#'
#' Reads a delimited text file with a header row into a data frame, finding the
#' delimiter by itself and applying the conventions of pharmacometric datasets:
#' `.` is a missing value, column names are upper-cased, and comment rows are
#' dropped. It also reads NONMEM table files (the `TABLE NO. 1` output of
#' `$TABLE`), several tables in one file included.
#'
#' @details
#' **Delimiter.** With `delimiter = NULL` the first lines decide: a tab, then a
#' comma, then a semicolon, otherwise white space. Give `"comma"`, `"tab"`,
#' `"space"` or any single character to override.
#'
#' **Missing values.** A field that is `.`, empty or `NA` becomes `NA`
#' (change with `na.strings`). A column that is entirely numeric apart from
#' such fields is returned as numeric, including numbers written with a
#' Fortran exponent such as `1.5D-3`.
#'
#' **Column names.** Names are trimmed, a leading `#` is removed, and the names
#' are upper-cased, so later code can rely on `ID`, `TIME`, `DV`, `AMT`, and so
#' on. Set `upper = FALSE` to keep the names as written (useful when the file
#' mixes `Wt` and `WT`).
#'
#' **Comments.** A row whose first field starts with `#` or `@` is a comment
#' and is dropped, like NONMEM's `IGNORE=#` and `IGNORE=@`.
#'
#' **NONMEM tables.** When the file contains `TABLE NO.` lines they and the
#' repeated header lines are removed. If the file holds several tables a
#' `TABLE_NO` column is added.
#'
#' @param file Path of the dataset (CSV, TSV, whitespace-delimited text, or a
#'   NONMEM table file).
#' @param delimiter `NULL` (find it), `"comma"`, `"tab"`, `"space"` or a single
#'   character.
#' @param skip Number of lines to skip before the header, for files with title
#'   lines.
#' @param header `TRUE` when the first line holds the column names. With
#'   `FALSE` the columns are named `A`, `B`, `C`, and so on.
#' @param na.strings Character vector of texts that mean a missing value.
#' @param upper `TRUE` (default) upper-cases the column names.
#' @return A data frame with one column per field. Numeric columns are numeric,
#'   the others are character.
#' @seealso [hp_mutate()] and [hp_filter()] for data steps, [hp_nca()] and
#'   [hp_eda_prepare()], which take the data frame this function returns.
#' @family data preparation
#' @examples
#' f <- hp_example_file()
#' d <- hp_read_data(f)
#' head(d)
#' str(d)
#'
#' # a semicolon file with a title line and lower-case names
#' tmp <- tempfile(fileext = ".csv")
#' writeLines(c("Study XYZ-001", "id;time;dv", "1;0;.", "1;1;4.2", "1;2;3.1"), tmp)
#' hp_read_data(tmp, skip = 1)
#' hp_read_data(tmp, skip = 1, upper = FALSE)
#' unlink(tmp)
#' @export
hp_read_data <- function(file, delimiter = NULL, skip = 0, header = TRUE,
                         na.strings = c(".", "", "NA"), upper = TRUE) {
  if (!file.exists(file)) hp_stop("the dataset `", file, "` was not found")
  first <- readLines(file, n = skip + 20L, warn = FALSE)
  first <- first[nzchar(trimws(first))]
  if (skip > 0L) first <- utils::tail(first, -min(skip, length(first)))
  sep <- if (is.null(delimiter)) {
    if (any(grepl("\t", first, fixed = TRUE))) "\t"
    else if (any(grepl(",", first, fixed = TRUE))) ","
    else if (any(grepl(";", first, fixed = TRUE))) ";"
    else ""
  } else switch(delimiter, comma = ",", tab = "\t", space = "", whitespace = "", delimiter)
  src <- file; skip_ <- skip
  all_lines <- readLines(file, warn = FALSE)
  if (any(grepl("^\\s*TABLE NO\\.", all_lines))) {
    # NONMEM table file: drop the TABLE NO. lines and the header repeated at the start of every
    # table; several tables get a TABLE_NO column
    l <- all_lines; if (skip > 0L) l <- utils::tail(l, -min(skip, length(l)))
    is_t <- grepl("^\\s*TABLE NO\\.", l); tno <- cumsum(is_t)
    body <- l[!is_t & nzchar(trimws(l))]; tno <- tno[!is_t & nzchar(trimws(l))]
    hdr <- body[1L]; keep <- body != hdr | seq_along(body) == 1L
    body <- body[keep]; tno <- tno[keep]
    if (length(unique(tno[-1L])) > 1L || max(tno) > 1L) {
      body <- c(paste(hdr, "TABLE_NO"), paste(body[-1L], tno[-1L]))
    }
    src <- textConnection(body); skip_ <- 0L
    if (is.null(delimiter)) sep <- ""
    on.exit(try(close(src), silent = TRUE), add = TRUE)
  }
  d <- if (identical(sep, "")) {
    utils::read.table(src, header = header, skip = skip_, na.strings = na.strings,
                      stringsAsFactors = FALSE, comment.char = "", check.names = FALSE,
                      strip.white = TRUE, blank.lines.skip = TRUE, fill = TRUE)
  } else {
    utils::read.csv(src, header = header, sep = sep, skip = skip_, na.strings = na.strings,
                    stringsAsFactors = FALSE, check.names = FALSE, strip.white = TRUE,
                    comment.char = "", blank.lines.skip = TRUE, fill = TRUE)
  }
  if (header) {
    names(d) <- trimws(sub("^#", "", names(d)))
    if (upper) names(d) <- toupper(names(d))
  } else {
    names(d) <- LETTERS[seq_len(ncol(d))]
  }
  if (nrow(d) > 0L) {
    f <- as.character(d[[1L]])
    comment <- !is.na(f) & grepl("^[#@]", f)
    d <- d[!comment, , drop = FALSE]
    d[] <- lapply(d, function(col) {
      if (is.character(col)) {
        num <- suppressWarnings(as.numeric(sub("[dD]", "e", col)))
        if (all(is.na(col) | !is.na(num))) num else col
      } else col
    })
  }
  rownames(d) <- NULL
  d
}

#' Path of the example dataset shipped with the package
#'
#' The package includes a small simulated study so that every example and
#' vignette runs without any data of your own. It is an oral
#' one-compartment study of 24 subjects: 12 receive 100 and 12 receive 200 (in
#' any dose unit), each with one dose at time 0 and 12 observations from 0.25 to
#' 48 hours. The columns are `ID`, `TIME`, `AMT`, `DV`, `EVID`, `MDV`, `CMT`,
#' `DVID`, `WT` (kg), `SEX` (0 or 1), `AGE` and `DOSEGRP`. Some concentrations
#' are small enough to be treated as below the limit of quantification when you
#' set `lloq = 0.5`.
#'
#' @return The full path of a CSV file as a character string.
#' @seealso [hp_read_data()]
#' @family data preparation
#' @examples
#' f <- hp_example_file()
#' basename(f)
#' d <- hp_read_data(f)
#' table(d$DOSEGRP[d$EVID == 1])
#' @export
hp_example_file <- function() {
  f <- system.file("extdata", "example_pk.csv", package = "hyperutils")
  if (!nzchar(f)) hp_stop("the example file is not installed")
  f
}

# ---------------------------------------------------------------------------
# data steps in R syntax: mutate and filter
# ---------------------------------------------------------------------------

#' Vectorised if-else chain (a dplyr-style case_when)
#'
#' Evaluates a series of `test ~ value` pairs and returns, for every element,
#' the value of the first test that is true. It is the clearest way to recode a
#' continuous variable into groups inside [hp_mutate()], [hp_filter()] and
#' [hp_custom_prep()], where it is available under the name `case_when`.
#'
#' @details
#' - `TRUE ~ value` as the last pair is the default for everything not matched
#'   before.
#' - When no test is true for an element the result is `NA`.
#' - A missing test counts as false.
#' - Values are recycled to the length of the tests; mixing numbers and
#'   text gives a character result.
#'
#' @param ... Two-sided formulas `test ~ value`.
#' @return A vector with one element per row of the tests.
#' @seealso [hp_mutate()]
#' @family data preparation
#' @examples
#' w <- c(50, 70, 95, NA)
#' hp_case_when(w < 60 ~ "low", w < 85 ~ "mid", TRUE ~ "high")
#'
#' # no default: unmatched elements are NA
#' hp_case_when(w < 60 ~ 1, w < 85 ~ 2)
#'
#' # inside hp_mutate, as text
#' d <- data.frame(WT = c(50, 70, 95))
#' hp_mutate(d, "WTGRP = case_when(WT < 60 ~ 1, WT < 85 ~ 2, TRUE ~ 3)")
#' @export
hp_case_when <- function(...) {
  args <- list(...)
  env <- parent.frame()
  n <- NULL
  tests <- vector("list", length(args)); vals <- vector("list", length(args))
  for (i in seq_along(args)) {
    f <- args[[i]]
    if (!inherits(f, "formula") || length(f) != 3L) hp_stop("case_when needs `test ~ value` pairs")
    tests[[i]] <- eval(f[[2L]], environment(f) %||% env)
    vals[[i]] <- eval(f[[3L]], environment(f) %||% env)
    n <- max(n, length(tests[[i]]), length(vals[[i]]))
  }
  out <- rep(NA, n)
  done <- rep(FALSE, n)
  for (i in seq_along(args)) {
    t <- rep_len(tests[[i]], n); t[is.na(t)] <- FALSE
    take <- t & !done
    v <- rep_len(vals[[i]], n)
    if (is.character(v) && !is.character(out)) out <- as.character(out)
    out[take] <- v[take]
    done <- done | take
  }
  out
}

# split "A = f(x, 1); B = 2" at the semicolons outside quotes and brackets
hp_split_statements <- function(text) {
  chars <- strsplit(paste(text, collapse = ";"), "")[[1L]]
  out <- character(); cur <- character(); depth <- 0L; q <- ""
  for (ch in chars) {
    if (nzchar(q)) {
      cur <- c(cur, ch); if (ch == q) q <- ""
    } else if (ch == "\"" || ch == "'") {
      q <- ch; cur <- c(cur, ch)
    } else if (ch %in% c("(", "[", "{")) {
      depth <- depth + 1L; cur <- c(cur, ch)
    } else if (ch %in% c(")", "]", "}")) {
      depth <- depth - 1L; cur <- c(cur, ch)
    } else if ((ch == ";" || ch == "\n") && depth == 0L) {
      s <- trimws(paste(cur, collapse = "")); if (nzchar(s)) out <- c(out, s); cur <- character()
    } else cur <- c(cur, ch)
  }
  s <- trimws(paste(cur, collapse = "")); if (nzchar(s)) out <- c(out, s)
  out
}

hp_eval_env <- function(data, parent = baseenv()) {
  env <- new.env(parent = parent)
  env$case_when <- hp_case_when
  env$ifelse <- base::ifelse
  list2env(as.list(data), envir = env)
  env
}

#' Create or change columns with statements written in R syntax
#'
#' Runs one or more `NAME = expression` statements on a data frame, in order, so
#' that a later statement can use a column made by an earlier one. This is the
#' quick way to derive covariates, transform concentrations and recode groups
#' without writing a pipeline of your own, and the same text can be reused in
#' [hp_custom_prep()].
#'
#' @details
#' Statements are separated by semicolons (or line breaks) outside quotes and
#' brackets. Column names are upper-cased, so `LNDV = log(DV)` creates `LNDV`.
#' Expressions see the columns of the data frame as variables and may use base R
#' functions, `case_when()` and `ifelse()`. A single value is recycled to every
#' row; any other length must equal the number of rows. A statement that cannot
#' be evaluated stops with a message that names it.
#'
#' @param data A data frame, usually from [hp_read_data()].
#' @param statements A character vector of statements, or one string with
#'   semicolons. `NULL` or empty returns the data unchanged.
#' @return The data frame with the new and changed columns.
#' @seealso [hp_filter()], [hp_case_when()], [hp_custom_prep()]
#' @family data preparation
#' @examples
#' d <- data.frame(WT = c(50, 70, 95), DV = c(1, 2, 3), DOSE = c(100, 100, 200))
#' hp_mutate(d, "LNDV = log(DV); DVN = DV / DOSE")
#'
#' # a later statement uses an earlier one
#' hp_mutate(d, c("BSA = 0.007184 * WT^0.425 * 170^0.725", "DOSEBSA = DOSE / BSA"))
#'
#' # recode into groups
#' hp_mutate(d, "WTGRP = case_when(WT < 60 ~ 'light', WT < 85 ~ 'medium', TRUE ~ 'heavy')")
#' @export
hp_mutate <- function(data, statements) {
  if (is.null(statements) || !length(statements)) return(data)
  for (s in hp_split_statements(statements)) {
    m <- regmatches(s, regexec("^\\s*([A-Za-z_.][A-Za-z0-9_.]*)\\s*(<-|=)\\s*(.*)$", s))[[1L]]
    if (length(m) != 4L) hp_stop("mutate: `", s, "` is not `NAME = expression`")
    nm <- toupper(m[2L])
    env <- hp_eval_env(data)
    val <- tryCatch(eval(parse(text = m[4L])[[1L]], envir = env),
                    error = function(e) hp_stop("mutate: `", s, "`: ", conditionMessage(e)))
    if (length(val) == 1L) val <- rep(val, nrow(data))
    if (length(val) != nrow(data)) hp_stop("mutate: `", s, "` gives ", length(val), " values for ", nrow(data), " rows")
    data[[nm]] <- val
  }
  data
}

#' Keep the rows that satisfy an R expression
#'
#' Subsets a data frame with a condition written as text. It behaves like
#' `subset()` with one difference that matters for clinical data: a row whose
#' test is `NA` is dropped rather than kept. The same expression text can be
#' used with [hp_custom_prep()], [hp_nca_read()] and [hp_eda_prepare()].
#'
#' @param data A data frame.
#' @param expr A condition as text, for example `"EVID == 0 & WT > 40"`. It may
#'   use `case_when()` and `ifelse()`. `NULL` or an empty string returns the
#'   data unchanged.
#' @return The rows that pass, with row names reset.
#' @seealso [hp_mutate()], [hp_custom_prep()]
#' @family data preparation
#' @examples
#' d <- data.frame(ID = 1:4, WT = c(50, 70, NA, 95), SEX = c(1, 0, 1, 1))
#' hp_filter(d, "WT > 60")            # the NA row is dropped
#' hp_filter(d, "WT > 40 & SEX == 1")
#' hp_filter(d, "ID %in% c(2, 4)")
#' @export
hp_filter <- function(data, expr) {
  if (is.null(expr) || !nzchar(trimws(expr))) return(data)
  keep <- tryCatch(eval(parse(text = expr)[[1L]], envir = hp_eval_env(data)),
                   error = function(e) hp_stop("filter `", expr, "`: ", conditionMessage(e)))
  if (length(keep) == 1L) keep <- rep(keep, nrow(data))
  keep[is.na(keep)] <- FALSE
  out <- data[as.logical(keep), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# ---------------------------------------------------------------------------
# descriptive statistics
# ---------------------------------------------------------------------------

hp_quantile7 <- function(x, p) {
  x <- sort(x[is.finite(x)])
  if (!length(x)) return(rep(NA_real_, length(p)))
  as.numeric(stats::quantile(x, p, type = 7, names = FALSE))
}

#' Descriptive statistics of a numeric vector
#'
#' One call returns the statistics usually requested in a clinical
#' pharmacology summary table: counts, location, spread, percentiles,
#' geometric statistics and the confidence interval of the mean.
#'
#' @details
#' Missing and non-finite values are left out of every statistic and are
#' counted in `nmiss`. Percentiles use R's type 7 definition. The geometric
#' mean and geometric CV are returned only when every value is positive
#' (`NA` otherwise); the geometric CV is `100 * sqrt(exp(sd(log(x))^2) - 1)`.
#' The CV is the standard deviation as a percentage of the absolute mean.
#' The confidence interval is the 95 percent interval of the mean from the t
#' distribution.
#'
#' @param x A numeric vector.
#' @return A named numeric vector: `n`, `nmiss`, `mean`, `sd`, `cv`, `sem`,
#'   `median`, `min`, `max`, `q1`, `q3`, `p5`, `p95`, `geomean`, `gcv`,
#'   `ci_lo` and `ci_hi`.
#' @seealso [hp_centre_spread()], [hp_custom_table()], [hp_eda_table()]
#' @family descriptive statistics
#' @examples
#' x <- c(1.2, 2.5, 3.1, 4.8, NA, 5.5)
#' round(hp_describe(x), 3)
#'
#' # one statistic
#' hp_describe(x)[["geomean"]]
#'
#' # by group
#' d <- hp_read_data(hp_example_file())
#' s <- d[!duplicated(d$ID), ]
#' sapply(split(s$WT, s$DOSEGRP), function(v) hp_describe(v)[c("n", "mean", "sd")])
#' @export
hp_describe <- function(x) {
  x0 <- x
  x <- sort(x[is.finite(x)])
  n <- length(x)
  m <- if (n) mean(x) else NA_real_
  s <- if (n > 1L) stats::sd(x) else NA_real_
  pos <- all(x > 0) && n > 0L
  gm <- gcv <- NA_real_
  if (pos) {
    l <- log(x); ls <- if (n > 1L) stats::sd(l) else NA_real_
    gm <- exp(mean(l)); gcv <- sqrt(exp(ls^2) - 1) * 100
  }
  sem <- s / sqrt(n)
  tq <- if (n >= 2L) stats::qt(0.975, n - 1) else NA_real_
  c(n = n, nmiss = length(x0) - n, mean = m, sd = s,
    cv = if (is.finite(m) && m != 0) s / abs(m) * 100 else NA_real_,
    sem = sem, median = hp_quantile7(x, 0.5), min = if (n) x[1L] else NA_real_,
    max = if (n) x[n] else NA_real_, q1 = hp_quantile7(x, 0.25), q3 = hp_quantile7(x, 0.75),
    p5 = hp_quantile7(x, 0.05), p95 = hp_quantile7(x, 0.95), geomean = gm, gcv = gcv,
    ci_lo = m - tq * sem, ci_hi = m + tq * sem)
}

# the statistics of the NCA summary: geometric statistics over the positive values
hp_describe_nca <- function(x) {
  x <- sort(x[is.finite(x)])
  n <- length(x)
  m <- if (n) mean(x) else NA_real_
  s <- if (n > 1L) stats::sd(x) else NA_real_
  pos <- x[x > 0]
  gm <- gcv <- NA_real_
  if (length(pos)) {
    l <- log(pos); sl <- if (length(l) > 1L) stats::sd(l) else NA_real_
    gm <- exp(mean(l)); gcv <- if (is.finite(sl)) 100 * sqrt(exp(sl^2) - 1) else NA_real_
  }
  c(n = n, mean = m, sd = s, cv = if (is.finite(m) && m != 0 && is.finite(s)) 100 * s / abs(m) else NA_real_,
    sem = if (n >= 2L) s / sqrt(n) else NA_real_, geomean = gm, geocv = gcv,
    median = hp_quantile7(x, 0.5), min = if (n) x[1L] else NA_real_, q1 = hp_quantile7(x, 0.25),
    q3 = hp_quantile7(x, 0.75), max = if (n) x[n] else NA_real_,
    harmonic = if (length(pos)) length(pos) / sum(1 / pos) else NA_real_)
}

#' Centre and error bar of a set of values
#'
#' Returns a measure of centre and the lower and upper end of an error bar, as
#' used for mean concentration-time profiles and summary plots.
#'
#' @details
#' With `centre = "geomean"` and `spread` one of `"sd"`, `"sem"` or `"ci95"`,
#' the bars are computed on the log scale and back-transformed, so they are
#' asymmetric and always positive. Missing values are ignored, and only
#' positive values are used for the geometric mean.
#'
#' @param x A numeric vector.
#' @param centre `"mean"`, `"median"` or `"geomean"`.
#' @param spread `"sd"`, `"sem"`, `"ci95"` (95 percent confidence interval of
#'   the mean), `"iqr"` (first to third quartile), `"p90"` (5th to 95th
#'   percentile), `"range"` or `"none"`.
#' @return A named numeric vector `centre`, `lower`, `upper`; all `NA` when
#'   there are no finite values, `lower` and `upper` `NA` for `"none"`.
#' @seealso [hp_describe()], [hp_eda_plot()] (summary plots use it)
#' @family descriptive statistics
#' @examples
#' x <- c(1.2, 2.5, 3.1, 4.8, 5.5)
#' hp_centre_spread(x)                         # mean and SD
#' hp_centre_spread(x, "median", "iqr")
#' hp_centre_spread(x, "geomean", "ci95")
#' hp_centre_spread(x, "mean", "range")
#' @export
hp_centre_spread <- function(x, centre = "mean", spread = "sd") {
  d <- hp_describe(x)
  if (d[["n"]] == 0) return(c(centre = NA_real_, lower = NA_real_, upper = NA_real_))
  cen <- switch(centre, mean = d[["mean"]], median = d[["median"]], geomean = d[["geomean"]], d[["mean"]])
  bar <- c(NA_real_, NA_real_)
  if (spread != "none") {
    if (centre == "geomean" && spread %in% c("sd", "sem", "ci95")) {
      l <- log(x[is.finite(x) & x > 0]); m <- mean(l); s <- stats::sd(l); k <- length(l)
      w <- switch(spread, sd = s, sem = s / sqrt(k), stats::qt(0.975, k - 1) * s / sqrt(k))
      bar <- exp(c(m - w, m + w))
    } else bar <- switch(spread, sd = cen + c(-1, 1) * d[["sd"]], sem = cen + c(-1, 1) * d[["sem"]],
                         ci95 = c(d[["ci_lo"]], d[["ci_hi"]]), iqr = c(d[["q1"]], d[["q3"]]),
                         p90 = c(d[["p5"]], d[["p95"]]), range = c(d[["min"]], d[["max"]]), c(NA_real_, NA_real_))
  }
  c(centre = cen, lower = bar[1L], upper = bar[2L])
}

#' The plot theme used by the package
#'
#' A clean, publication-style ggplot2 theme: `theme_classic()` with the legend
#' at the bottom, light grey facet strips and a small grey caption. Every figure
#' drawn by the package uses it, and you can add it to your own plots to match.
#'
#' @param base_size Base font size in points.
#' @param legend_position Position of the legend: `"bottom"`, `"top"`,
#'   `"left"`, `"right"` or `"none"`.
#' @return A ggplot2 theme object.
#' @seealso [hp_custom_plot()], [hp_eda_plot()]
#' @family output
#' @examples
#' library(ggplot2)
#' d <- hp_read_data(hp_example_file())
#' d <- d[d$EVID == 0, ]
#' ggplot(d, aes(TIME, DV, group = ID, colour = factor(DOSEGRP))) +
#'   geom_line() +
#'   labs(colour = "Dose group") +
#'   hp_theme(base_size = 12)
#' @export
hp_theme <- function(base_size = 11, legend_position = "bottom") {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(legend.position = legend_position,
                   strip.background = ggplot2::element_rect(fill = "grey95", colour = NA),
                   plot.caption = ggplot2::element_text(colour = "grey40", size = base_size * 0.75, hjust = 0))
}

#' Save one plot as a PNG and/or SVG file
#'
#' Writes a ggplot to the folder you name, creating it if necessary, in one or
#' more formats. The package never writes files by default; you always choose
#' the folder.
#'
#' @details
#' SVG output needs the 'svglite' package; when it is not installed the SVG file
#' is skipped with a warning and the PNG is still written. The size is in
#' inches; for a figure of 17 by 12 cm use `width = 17 / 2.54` and
#' `height = 12 / 2.54`.
#'
#' @param plot A ggplot object.
#' @param file File name without extension.
#' @param dir Folder to write to (created when missing).
#' @param width,height Size in inches.
#' @param dpi Resolution of raster formats, in dots per inch.
#' @param formats Character vector with `"png"`, `"svg"` or both.
#' @return The paths written, invisibly.
#' @seealso [hp_save_plots()] for multi-page figures, [hp_save_table()]
#' @family output
#' @examples
#' library(ggplot2)
#' p <- ggplot(data.frame(x = 1:3, y = c(2, 4, 3)), aes(x, y)) + geom_line()
#' out <- hp_save_plot(p, "demo", tempdir(), width = 5, height = 3, dpi = 100)
#' file.exists(out)
#' unlink(out)
#' @export
hp_save_plot <- function(plot, file, dir, width = 7, height = 5, dpi = 300,
                         formats = "png") {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  out <- character()
  for (f in formats) {
    p <- file.path(dir, paste0(file, ".", f))
    if (f == "svg" && !requireNamespace("svglite", quietly = TRUE)) {
      warning("package 'svglite' is needed for SVG files; skipped", call. = FALSE); next
    }
    ggplot2::ggsave(p, plot, width = width, height = height, dpi = dpi)
    out <- c(out, p)
  }
  invisible(out)
}

#' Save every page of a multi-page figure
#'
#' Individual-profile figures and large facet grids are split into pages
#' by the package. This function writes all the pages with a numbered suffix
#' (`name_01.png`, `name_02.png`, and so on), or a single file when the figure
#' has one page.
#'
#' @param pages A ggplot or a list of ggplots (the pages), as returned by
#'   [hp_eda_plot()].
#' @param name File name without extension.
#' @param dir,width,height,dpi,formats As in [hp_save_plot()].
#' @return The paths written, invisibly.
#' @seealso [hp_save_plot()], [hp_eda_plot()], [hp_nca_plots()]
#' @family output
#' @examples
#' d <- hp_read_data(hp_example_file())
#' e <- hp_eda_prepare(d, dvid = "DVID")
#' s <- hp_eda_spec("P01", "conc_time", "Concentration over time", "line",
#'                  "obs", x = "TAFD", y = "DV", pk_only = TRUE)
#' pages <- hp_eda_plot(e, s)
#' files <- hp_save_plots(pages, "conc_time", tempdir(), 6, 4, dpi = 72)
#' basename(files)
#' unlink(files)
#' @export
hp_save_plots <- function(pages, name, dir, width = 7, height = 5, dpi = 300, formats = "png") {
  if (inherits(pages, "ggplot")) pages <- list(pages)
  n <- length(pages); w <- max(nchar(as.character(n)), 2L)
  out <- character()
  for (i in seq_len(n)) {
    nm <- if (n <= 1L) name else sprintf(paste0("%s_%0", w, "d"), name, i)
    out <- c(out, hp_save_plot(pages[[i]], nm, dir, width, height, dpi, formats))
  }
  invisible(out)
}

# the shortest decimal text that reads back as the same number, without an
# exponent; "" for a missing value
hp_shortest <- function(x) {
  out <- character(length(x)); todo <- is.finite(x)
  for (d in 1:17) {
    if (!any(todo)) break
    s <- sprintf(paste0("%.", d, "g"), x[todo])
    ok <- as.numeric(s) == x[todo]
    idx <- which(todo)[ok]
    if (length(idx)) out[idx] <- vapply(x[idx], function(z) format(z, digits = d, scientific = FALSE, trim = TRUE), "")
    todo[idx] <- FALSE
  }
  out
}

#' Write a table as a CSV file
#'
#' Saves a data frame as a UTF-8 CSV file in the folder you name. The layout is
#' plain and predictable: a header row, no row names, missing values as empty
#' fields, and text quoted only when it contains a comma, a quote or a line
#' break.
#'
#' @details
#' Numeric columns are written in the shortest form that reads back as exactly
#' the same number, so no precision is lost, unless you set `digits`, which
#' rounds to that many significant digits and removes trailing zeros.
#'
#' @param x A data frame.
#' @param file File name, for example `"nca_parameters.csv"`.
#' @param dir Folder to write to (created when missing).
#' @param digits `NULL` for full precision, or the number of significant digits
#'   for numeric columns.
#' @return The path written, invisibly.
#' @seealso [hp_save_plot()], [hp_nca_table()]
#' @family output
#' @examples
#' r <- hp_nca(hp_read_data(hp_example_file()), route = "oral", lloq = 0.5)
#' f <- hp_save_table(r$summary, "summary.csv", tempdir(), digits = 4)
#' head(readLines(f), 3)
#' unlink(f)
#' @export
hp_save_table <- function(x, file, dir, digits = NULL) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  p <- file.path(dir, file)
  q <- function(z) {
    z <- ifelse(is.na(z), "", as.character(z))
    u <- Encoding(z) == "unknown" & validUTF8(z); Encoding(z[u]) <- "UTF-8"   # text is written as UTF-8
    need <- grepl("[,\"\n]", z)
    z[need] <- paste0("\"", gsub("\"", "\"\"", z[need], fixed = TRUE), "\"")
    z
  }
  cols <- lapply(x, function(v) {
    if (is.numeric(v) && !is.null(digits)) { r <- hp_cell(v, digits); r[!is.finite(v)] <- ""; q(r) }
    else if (is.numeric(v)) q(hp_shortest(v))
    else q(v)
  })
  lines <- c(paste(q(names(x)), collapse = ","),
             if (nrow(x)) do.call(paste, c(unname(cols), sep = ",")))
  con <- file(p, "wb"); on.exit(close(con))
  writeBin(charToRaw(enc2utf8(paste0(paste(lines, collapse = "\n"), "\n"))), con)
  invisible(p)
}
