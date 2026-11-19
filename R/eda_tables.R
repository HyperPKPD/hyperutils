# Exploratory tables: descriptive statistics, category counts, listings and the
# data summary / correlation / missing-value tables.

#' Format a number for a table cell
#'
#' Rounds numbers to a number of significant digits and writes them in the
#' compact style of the tables: whole numbers without decimals, very small or
#' very large numbers in scientific notation, trailing zeros removed and
#' non-finite values as an empty string.
#'
#' @param v A numeric vector.
#' @param digits Significant digits (default 4).
#' @return A character vector of the same length as `v`.
#' @seealso [hp_eda_table()], [hp_custom_table()]
#' @family exploratory data analysis
#' @examples
#' hp_cell(c(12, 3.14159265, 0.000012345, 123456789, NA, 0.5), digits = 3)
#' hp_cell(1234.5678, digits = 6)
#' @export
hp_cell <- function(v, digits = 4) {
  vapply(v, function(x) {
    if (!is.finite(x)) return("")
    if (x == round(x) && abs(x) < 1e12) return(format(x, scientific = FALSE, trim = TRUE))
    a <- abs(x)
    if (a != 0 && (a < 1e-4 || a >= 1e7)) {
      s <- sprintf(paste0("%.", max(digits - 1L, 0L), "e"), x)
      s <- sub("e([+-])0*([0-9])", "e\\1\\2", s)
      return(sub("e\\+", "e", s))
    }
    mag <- if (a == 0) 0 else floor(log10(a))
    dec <- min(max(digits - 1 - mag, 0), 10)
    s <- sprintf(paste0("%.", dec, "f"), x)
    if (grepl(".", s, fixed = TRUE)) s <- sub("\\.$", "", sub("0+$", "", s))
    s
  }, "")
}

#' Specification of an exploratory table
#'
#' A table is described by a specification, just as a plot is: the kind of
#' table, the variables, the grouping and the statistics. Pass it to
#' [hp_eda_table()].
#'
#' @details
#' **Kinds.**
#'
#' - `"summary"`: descriptive statistics of numeric variables, by group.
#' - `"count"`: counts and percentages of the levels of categorical variables,
#'   by group.
#' - `"listing"`: the rows themselves, for example a listing of the
#'   observations below the limit of quantification.
#'
#' **Statistics** of a summary table: `n`, `nmiss`, `mean`, `sd`, `cv`, `sem`,
#' `median`, `min`, `max`, `q1`, `q3`, `p5`, `p95`, `geomean`, `gcv`, `ci_lo`
#' and `ci_hi`. Percentiles use the type 7 definition of [stats::quantile()].
#'
#' @param id Catalogue identifier of the table, or `"custom"`.
#' @param name File name (without extension) used when the table is saved.
#' @param title Title.
#' @param kind `"summary"`, `"count"` or `"listing"`.
#' @param scope Rows used: `"subjects"`, `"obs"`, `"dose"` or `"all"`.
#' @param vars Variables to summarise, count or list.
#' @param by Grouping columns.
#' @param stats Statistics of a summary table, see Details.
#' @param filter R expression that selects rows.
#' @param pk_only Use observations of the PK `DVID` only.
#' @param digits Significant digits.
#' @param max_rows Most rows of a listing.
#' @param overall Add an overall group ("All") next to the groups.
#' @return A list of class `hp_eda_tspec`.
#' @seealso [hp_eda_table()], [hp_eda_spec()], [hp_custom_table()]
#' @family exploratory data analysis
#' @examples
#' # continuous covariates by dose group, with an overall column
#' hp_eda_tspec("custom", "cov", "Continuous covariates", "summary", "subjects",
#'              vars = c("WT", "AGE"), by = "DOSEGRP", overall = TRUE)
#'
#' # counts of a categorical covariate
#' hp_eda_tspec(kind = "count", scope = "subjects", vars = "SEX", by = "DOSEGRP")
#' @export
hp_eda_tspec <- function(id = "custom", name = "table", title = "", kind = "summary", scope = "subjects",
                         vars = character(), by = character(),
                         stats = c("n", "mean", "sd", "cv", "median", "min", "max"), filter = NULL,
                         pk_only = FALSE, digits = 4, max_rows = 1000, overall = FALSE) {
  s <- list(id = id, name = name, title = title, kind = kind, scope = scope, vars = toupper(vars),
            by = toupper(by), stats = stats, filter = filter, pk_only = pk_only, digits = digits,
            max_rows = max_rows, overall = overall)
  class(s) <- c("hp_eda_tspec", "list")
  s
}

# natural order of labels: numbers by value, then text
hp_sort_labels <- function(x) {
  n <- suppressWarnings(as.numeric(x))
  if (all(!is.na(n) | x == "")) x[order(is.na(n), n)] else sort(x)
}

#' Build an exploratory table
#'
#' Computes the table described by a specification (see [hp_eda_tspec()]) from
#' a prepared dataset (see [hp_eda_prepare()]): descriptive statistics, counts
#' and percentages, or a listing of rows.
#'
#' @param d An `hp_eda_data` object from [hp_eda_prepare()].
#' @param spec A specification from [hp_eda_tspec()].
#' @param max_categories Most distinct values of a numeric grouping column
#'   before it is cut into groups.
#' @return A data frame of text (formatted values), with the attributes
#'   `"title"` and `"note"` (how many rows were used). Write it with
#'   [hp_save_table()].
#' @seealso [hp_eda_tspec()], [hp_save_table()], [hp_eda_data_summary()]
#' @family exploratory data analysis
#' @examples
#' d <- hp_read_data(hp_example_file())
#' e <- hp_eda_prepare(d, dvid = "DVID", lloq = 0.5)
#'
#' # body weight and age by dose group
#' hp_eda_table(e, hp_eda_tspec(kind = "summary", scope = "subjects",
#'                              vars = c("WT", "AGE"), by = "DOSEGRP",
#'                              stats = c("n", "mean", "sd", "median", "min", "max")))
#'
#' # distribution of a categorical covariate
#' hp_eda_table(e, hp_eda_tspec(kind = "count", scope = "subjects",
#'                              vars = "SEX", by = "DOSEGRP", overall = TRUE))
#'
#' # listing of the first rows of concentrations
#' hp_eda_table(e, hp_eda_tspec(kind = "listing", scope = "obs",
#'                              vars = c("ID", "TIME", "DV"), max_rows = 5))
#' @export
hp_eda_table <- function(d, spec, max_categories = 20) {
  subj <- spec$scope == "subjects"
  f <- if (subj) d$subj else d$rec
  keep <- rep(TRUE, nrow(f))
  if (!is.null(spec$filter) && nzchar(spec$filter)) {
    k <- eval(parse(text = spec$filter)[[1L]], envir = hp_eval_env(f)); k[is.na(k)] <- FALSE; keep <- keep & k
  }
  if (!subj) {
    ok <- switch(spec$scope, obs = f$EVENT == "observation", dose = f$EVENT == "dose", rep(TRUE, nrow(f)))
    keep <- keep & ok & (!spec$pk_only | f$PKOBS == 1)
  }
  miss <- setdiff(c(spec$vars, spec$by), names(f))
  if (length(miss)) hp_stop("the column `", miss[1L], "` is not here")
  if (!length(spec$vars)) hp_stop("no variables")
  used <- sum(keep)
  what <- if (subj) "individuals" else "records"
  # group labels in the product order of the levels
  if (length(spec$by)) {
    lv <- lapply(spec$by, function(b) hp_eda_levels(f[[b]][keep], max_categories, 6))
    codes <- do.call(cbind, lapply(lv, function(z) as.integer(z)))
    gi <- rep(NA_integer_, nrow(f))
    kk <- which(keep)
    ok <- stats::complete.cases(codes)
    keys <- unique(codes[ok, , drop = FALSE]); keys <- keys[do.call(order, as.data.frame(keys)), , drop = FALSE]
    labs <- vapply(seq_len(nrow(keys)), function(i) paste(vapply(seq_along(spec$by), function(j)
      paste0(spec$by[j], "=", levels(lv[[j]])[keys[i, j]]), ""), collapse = ", "), "")
    m <- match(apply(codes, 1L, paste, collapse = "\r"), apply(keys, 1L, paste, collapse = "\r"))
    m[!ok] <- NA
    gi[kk] <- m
    by_name <- paste(spec$by, collapse = ", ")
  } else {
    labs <- "All"; gi <- ifelse(keep, 1L, NA_integer_); by_name <- "Group"
  }
  note <- paste(used, what)
  if (spec$kind == "summary") {
    rows <- list()
    for (v in spec$vars) {
      x <- hp_num(f[[v]])
      groups <- lapply(seq_along(labs), function(g) list(labs[g], x[which(gi == g)]))
      if (spec$overall && length(spec$by)) groups[[length(groups) + 1L]] <- list("All", x[keep])
      for (g in groups) {
        ds <- hp_describe(g[[2L]])
        rows[[length(rows) + 1L]] <- c(v, g[[1L]], vapply(spec$stats, function(s) hp_cell(ds[[s]], spec$digits), ""))
      }
    }
    out <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
    names(out) <- c("Variable", by_name, spec$stats)
  } else if (spec$kind == "count") {
    show_g <- length(spec$by) > 0L
    head <- c("Variable", "Level", if (show_g) labs, if (spec$overall || !show_g) "All")
    pct <- function(n, den) if (den == 0) "0" else sprintf("%d (%.1f%%)", n, 100 * n / den)
    rows <- list()
    gtot <- vapply(seq_along(labs), function(g) sum(gi == g, na.rm = TRUE), 1L)
    for (v in spec$vars) {
      txt <- ifelse(is.na(f[[v]]), "", if (is.numeric(f[[v]])) hp_num_text(f[[v]]) else as.character(f[[v]]))
      levs <- hp_sort_labels(unique(txt[keep]))
      for (lv in levs) {
        row <- c(v, if (lv == "") "(missing)" else lv)
        if (show_g) row <- c(row, vapply(seq_along(labs), function(g) pct(sum(gi == g & txt == lv, na.rm = TRUE), gtot[g]), ""))
        if (spec$overall || !show_g) row <- c(row, pct(sum(keep & txt == lv), used))
        rows[[length(rows) + 1L]] <- row
      }
    }
    out <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE); names(out) <- head
  } else if (spec$kind == "listing") {
    idx <- which(keep)
    if (length(idx) > spec$max_rows) { idx <- idx[seq_len(spec$max_rows)]; note <- paste0("the first ", spec$max_rows, " of ", used, " rows") }
    out <- as.data.frame(lapply(spec$vars, function(v) {
      x <- f[[v]][idx]; n <- hp_num(x)
      ifelse(is.finite(n), hp_cell(n, spec$digits), ifelse(is.na(x), "", as.character(x)))
    }), stringsAsFactors = FALSE)
    names(out) <- spec$vars
  } else hp_stop("table kind `", spec$kind, "` is not known (summary | count | listing)")
  rownames(out) <- NULL
  attr(out, "note") <- note
  attr(out, "title") <- spec$title
  out
}

#' Overview of a dataset
#'
#' A short table of the facts a reader wants first: number of records,
#' individuals, observations, doses, observations below the limit of
#' quantification, the dose levels and the number of observations per
#' individual.
#'
#' @param d An `hp_eda_data` object from [hp_eda_prepare()].
#' @return A data frame with the columns `Item` and `Value`.
#' @seealso [hp_eda_prepare()], [hp_eda_missing()]
#' @family exploratory data analysis
#' @examples
#' e <- hp_eda_prepare(hp_read_data(hp_example_file()), dvid = "DVID", lloq = 0.5)
#' hp_eda_data_summary(e)
#' @export
hp_eda_data_summary <- function(d) {
  i <- d$info
  rows <- list(
    c("Records", i$records), c("Individuals", i$individuals), c("Observation records", i$observations),
    c("PK observations (concentration plots)", i$pk_observations), c("Dose records", i$doses),
    c("Observations below the LLOQ", i$blq))
  if (i$filtered_out > 0) rows[[length(rows) + 1L]] <- c("Records removed by the filters", i$filtered_out)
  if (length(i$dvid_values)) {
    rows[[length(rows) + 1L]] <- c("DVID values", paste(i$dvid_values, collapse = ", "))
    rows[[length(rows) + 1L]] <- c("PK DVID", i$pk_dvid %||% "")
    rows[[length(rows) + 1L]] <- c("PD DVID", i$pd_dvid %||% "\u2013")
  }
  nobs <- hp_num(d$subj$NOBS); ds <- hp_describe(nobs)
  rows <- c(rows, list(
    c("First-dose levels", paste(hp_cell(i$dose_levels, 6), collapse = ", ")),
    c("Continuous covariates", paste(i$continuous, collapse = ", ")),
    c("Categorical covariates", paste(i$categorical, collapse = ", ")),
    c("Time-varying covariates", paste(i$time_varying, collapse = ", ")),
    c("Colour group", i$color %||% "\u2013"), c("Strata (facets)", i$strata %||% "\u2013"),
    c("Observations per individual (median, range)",
      sprintf("%s (%s to %s)", hp_cell(ds[["median"]], 4), hp_cell(ds[["min"]], 4), hp_cell(ds[["max"]], 4)))))
  out <- as.data.frame(do.call(rbind, lapply(rows, as.character)), stringsAsFactors = FALSE)
  names(out) <- c("Item", "Value")
  out
}

#' Correlation matrix of the continuous covariates
#'
#' Pearson correlations between the continuous covariates of the individuals,
#' computed with the pairwise complete observations. Strong correlations (for
#' example weight with body surface area) warn that covariate effects cannot be
#' told apart.
#'
#' @param d An `hp_eda_data` object from [hp_eda_prepare()].
#' @return A data frame of text: the first column names the covariate and the
#'   others hold the correlations with three decimals.
#' @seealso [hp_eda_prepare()], [hp_custom_table()] with `kind = "correlation"`
#' @family exploratory data analysis
#' @examples
#' e <- hp_eda_prepare(hp_read_data(hp_example_file()), dvid = "DVID")
#' hp_eda_correlation(e)
#' @export
hp_eda_correlation <- function(d) {
  nm <- intersect(d$info$continuous, names(d$subj))
  m <- stats::cor(as.data.frame(lapply(d$subj[nm], hp_num)), use = "pairwise.complete.obs")
  out <- data.frame(nm, stringsAsFactors = FALSE)
  names(out) <- ""
  for (j in seq_along(nm)) out[[nm[j]]] <- ifelse(is.finite(m[, j]), sprintf("%.3f", m[, j]), "")
  out
}

#' Missing values of every column
#'
#' For every column of the dataset it counts the missing values among all
#' records and among the observations, and the individuals whose first record
#' is missing the value. Use it before modelling to see which covariates can be
#' used and which need imputation.
#'
#' @param d An `hp_eda_data` object from [hp_eda_prepare()].
#' @return A data frame with the columns `Column`, `Missing records`,
#'   `% records`, `Missing observations`, `% observations`,
#'   `Individuals missing at first record` and `Text values` (values that are
#'   not numbers in a numeric-looking column).
#' @seealso [hp_eda_data_summary()]
#' @family exploratory data analysis
#' @examples
#' e <- hp_eda_prepare(hp_read_data(hp_example_file()), dvid = "DVID")
#' m <- hp_eda_missing(e)
#' head(m, 6)
#' @export
hp_eda_missing <- function(d) {
  rec <- d$rec; nobs <- sum(rec$EVENT == "observation")
  skip <- c("TAFD", "DVNORM", "EVENT", "PKOBS", "GROUP_ALL", "TAD_CALC")
  rows <- list()
  for (cn in setdiff(names(rec), skip)) {
    v <- rec[[cn]]
    txt <- ifelse(is.na(v), "", if (is.numeric(v)) hp_num_text(v) else as.character(v))
    mis <- txt == ""
    mo <- sum(mis & rec$EVENT == "observation")
    nt <- sum(!mis & is.na(hp_num(v)))
    ids <- if (cn %in% names(d$subj)) { sv <- d$subj[[cn]]; sum(is.na(sv) | as.character(sv) == "") } else 0L
    rows[[length(rows) + 1L]] <- c(cn, sum(mis), sprintf("%.1f", 100 * sum(mis) / max(nrow(rec), 1)), mo,
                                   sprintf("%.1f", 100 * mo / max(nobs, 1)), ids, nt)
  }
  out <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(out) <- c("Column", "Missing records", "% records", "Missing observations", "% observations",
                  "Individuals missing at first record", "Text values")
  out
}
