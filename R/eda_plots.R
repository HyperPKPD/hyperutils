# Exploratory plots: a plot is described by a specification (the PlotSpec of
# package) and drawn with ggplot2.

#' Specification of an exploratory plot
#'
#' An exploratory plot is described by a specification: a plain list that says
#' what to draw (the kind of plot, the variables, the records) and how (axes,
#' smoothers, facets, pages). The specification is separate from the drawing so
#' that a figure can be stored, edited, repeated for other variables and
#' reproduced. Pass it to [hp_eda_plot()].
#'
#' @details
#' **Kinds of plots (`kind`).**
#'
#' - `"scatter"`: points of `y` against `x`, with an optional smoother.
#' - `"line"`: one line per individual (or per `by` column), the usual
#'   spaghetti plot of concentration against time.
#' - `"summary"`: the centre and spread of `y` in bins of `x`, by colour group.
#' - `"box"`: box plots of `y` by the categories of `x`.
#' - `"histogram"`: distribution of `x`.
#' - `"bar"`: counts of the categories of `x`, or the centre of `y` by `x`.
#' - `"corr"`: a correlation heat map of the variables in `x` and `y`.
#' - `"pairs"`: scatter plots of all pairs of the variables in `x` and `y`.
#' - `"pkpd"`: the response against the time-matched concentration (needs
#'   `DVID`).
#'
#' **Scope (`scope`).** Which rows are drawn: `"obs"` observations, `"dose"`
#' dose records, `"all"` every record, or `"subjects"` the one-row-per-subject
#' frame (use it for covariates and exposure).
#'
#' **Several variables.** Giving several variables in `x` and/or `y` for the
#' kinds `"box"`, `"histogram"`, `"bar"` and `"scatter"` draws one panel per
#' combination, with free axes.
#'
#' @param id Catalogue identifier of the plot, or `"custom"`.
#' @param name File name (without extension) used when the figure is saved.
#' @param title Title.
#' @param kind Kind of plot, see Details.
#' @param scope Records used: `"obs"`, `"dose"`, `"all"` or `"subjects"`.
#' @param x,y Variables: columns of the record or subject frame. Names are
#'   put in upper case.
#' @param by Columns that define the lines (default: the individual).
#' @param color,facet Columns that colour the data and split it into panels.
#' @param filter R expression that selects rows, for example
#'   `"DOSEGRP == 2"`.
#' @param pk_only Use observations of the PK `DVID` only.
#' @param log_x,log_y Logarithmic axes.
#' @param smooth `"none"`, `"linear"`, `"loess"` or `"power"` (a straight line
#'   on log-log axes).
#' @param smooth_groups One smoother per colour group instead of one overall.
#' @param span Span of the loess smoother.
#' @param centre,spread Centre (`"mean"`, `"median"` or `"geomean"`) and error
#'   bars (`"sd"`, `"sem"`, `"ci95"`, `"iqr"`, `"p90"`, `"range"` or `"none"`)
#'   of a summary plot; see [hp_centre_spread()].
#' @param bins Number of x bins of a summary plot (0 = automatic). Use bins
#'   when sampling times differ between individuals.
#' @param points,lines Draw the points and/or the lines.
#' @param hist_bins Number of histogram bins (0 = automatic).
#' @param hlines,vlines Positions of horizontal and vertical reference lines.
#' @param identity Draw the line of identity (y = x).
#' @param ncol,nrow,max_pages Panels per page and the most pages drawn.
#' @param free Free axes of facet panels: `"none"`, `"x"`, `"y"` or `"both"`.
#' @param xlab,ylab Axis titles (default: derived from the variable and the
#'   units).
#' @param xlim,ylim Axis limits as `c(lo, hi)`; `NA` leaves a side open.
#' @param jitter Horizontal jitter of the points.
#' @param legend Show the legend.
#' @param max_groups Draw at most this many line groups (0 = all).
#' @return A list of class `hp_eda_spec`.
#' @seealso [hp_eda_plot()], [hp_eda_tspec()], [hp_custom_plot()]
#' @family exploratory data analysis
#' @examples
#' # concentration against time, one line per individual, log axis
#' s <- hp_eda_spec(kind = "line", x = "TIME", y = "DV", color = "DOSEGRP",
#'                  log_y = TRUE, pk_only = TRUE)
#' s$kind
#'
#' # covariate against exposure with a loess smoother
#' hp_eda_spec(kind = "scatter", x = "WT", y = "CMAX", scope = "subjects",
#'             smooth = "loess")
#'
#' # median and quartile range of concentrations in time bins
#' hp_eda_spec(kind = "summary", x = "TAD", y = "DV", color = "DOSEGRP",
#'             centre = "median", spread = "iqr", bins = 8, log_y = TRUE)
#' @export
hp_eda_spec <- function(id = "custom", name = "plot", title = "", kind = "scatter", scope = "obs",
                        x = character(), y = character(), by = NULL, color = NULL, facet = NULL,
                        filter = NULL, pk_only = FALSE, log_x = FALSE, log_y = FALSE, smooth = "none",
                        smooth_groups = TRUE, span = 0.75, centre = "mean", spread = "sd", bins = 0,
                        points = TRUE, lines = TRUE, hist_bins = 0, hlines = numeric(), vlines = numeric(),
                        identity = FALSE, ncol = 3, nrow = 3, max_pages = 10, free = "none", xlab = NULL,
                        ylab = NULL, xlim = c(NA, NA), ylim = c(NA, NA), jitter = 0, legend = TRUE,
                        max_groups = 0) {
  s <- list(id = id, name = name, title = title, kind = kind, scope = scope, x = toupper(x), y = toupper(y),
            by = if (is.null(by)) NULL else toupper(by), color = if (is.null(color)) NULL else toupper(color),
            facet = if (is.null(facet)) NULL else toupper(facet), filter = filter, pk_only = pk_only,
            log_x = log_x, log_y = log_y, smooth = smooth, smooth_groups = smooth_groups, span = span,
            centre = centre, spread = spread, bins = bins, points = points, lines = lines, hist_bins = hist_bins,
            hlines = hlines, vlines = vlines, identity = identity, ncol = ncol, nrow = nrow, max_pages = max_pages,
            free = free, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, jitter = jitter, legend = legend,
            max_groups = max_groups)
  class(s) <- c("hp_eda_spec", "list")
  s
}

#' Axis title of a variable
#'
#' Produces the axis title used for a standard variable: "Time (h)" for time,
#' "Time after dose (h)" for `TAD`, "Dose-normalised DV" for `DVNORM`, a
#' concentration with its unit for `CMAX`, `CMIN` and `CLAST`, and so on. Other
#' names are returned unchanged.
#'
#' @param var Variable name in upper case.
#' @param units A list with the elements `time`, `conc` and `dose`, each a
#'   unit text (empty for none).
#' @return A character string.
#' @seealso [hp_eda_plot()]
#' @family exploratory data analysis
#' @examples
#' hp_eda_label("TIME", list(time = "h", conc = "ng/mL", dose = "mg"))
#' hp_eda_label("CMAX", list(time = "h", conc = "ng/mL", dose = "mg"))
#' hp_eda_label("WT")
#' @export
hp_eda_label <- function(var, units = list(time = "h", conc = "", dose = "")) {
  w <- function(n, u) if (is.null(u) || !nzchar(u)) n else paste0(n, " (", u, ")")
  switch(var, TIME = w("Time", units$time), TAD = , TAD_CALC = w("Time after dose", units$time),
         TAFD = w("Time after first dose", units$time), TMAX = , TLAST = w(var, units$time),
         CMAX = , CMIN = , CLAST = w(var, units$conc),
         AUCLAST = if (is.null(units$conc) || !nzchar(units$conc)) "AUClast" else w("AUClast", paste0(units$conc, "\u00b7", units$time)),
         AMT = , DOSE = , DOSE1 = , TOTDOSE = w(var, units$dose), DVNORM = "Dose-normalised DV",
         CMAXD = "Dose-normalised Cmax", AUCLASTD = "Dose-normalised AUClast", var)
}

# a fixed colour for a layer, or nothing when the colour is mapped (a NULL
# parameter would be taken as a zero-length colour by some ggplot2 versions)
hp_fixed_colour <- function(fixed, colour) if (isTRUE(fixed)) list(colour = colour) else list()

# levels of a grouping column over the kept rows: a numeric column with more
# than `max` values is cut into `bins` quantile groups
hp_eda_levels <- function(v, max = 20, bins = 4) {
  txt <- ifelse(is.na(v), "", if (is.numeric(v)) hp_num_text(v) else as.character(v))
  uniq <- unique(txt)
  numeric <- is.numeric(v) || all(is.na(suppressWarnings(as.numeric(txt))) == (txt == ""))
  vn <- hp_num(v)
  if (max > 0 && length(uniq) > max && numeric && bins >= 2) {
    s <- sort(vn[is.finite(vn)])
    edges <- hp_quantile7(s, (0:bins) / bins)
    keep <- c(TRUE, diff(edges) >= 1e-12); edges <- edges[keep]
    k <- max(length(edges) - 1L, 1L)
    labels <- vapply(seq_len(k), function(i) sprintf("%s%s, %s]", if (i == 1L) "[" else "(", hp_short(edges[i]),
                                                     hp_short(edges[min(i + 1L, length(edges))])), "")
    idx <- vapply(vn, function(x) { if (!is.finite(x)) return(NA_integer_); g <- 1L; while (g < k && x > edges[g + 1L]) g <- g + 1L; g }, 1L)
    return(factor(labels[idx], levels = labels))
  }
  nums <- suppressWarnings(as.numeric(uniq))
  lv <- if (all(!is.na(nums) | uniq == "")) uniq[order(is.na(nums), nums)] else sort(uniq)
  factor(txt, levels = lv)
}

# rows of the record or subject frame a spec uses
hp_eda_rows <- function(d, spec) {
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
  list(frame = f, keep = keep, subjects = subj)
}

# x bins of a summary plot: few distinct values are kept, otherwise quantile bins
hp_eda_xbins <- function(x, nbins = 0) {
  u <- sort(x[is.finite(x)]); d <- u[c(TRUE, diff(u) > 1e-9 * (1 + abs(u[-1L])))]
  if (!length(d)) return(list(idx = rep(NA_integer_, length(x)), pos = numeric()))
  few <- if (nbins > 0) length(d) <= nbins else length(d) <= 40
  if (few) {
    idx <- vapply(x, function(v) if (!is.finite(v)) NA_integer_ else which(abs(d - v) <= 1e-9 * (1 + abs(d)))[1L], 1L)
    return(list(idx = idx, pos = d))
  }
  k <- if (nbins > 0) nbins else min(max(length(d) %/% 3L, 4L), 10L)
  edges <- hp_quantile7(u, (0:k) / k); edges <- edges[c(TRUE, diff(edges) >= 1e-12)]
  nb <- max(length(edges) - 1L, 1L)
  idx <- vapply(x, function(v) { if (!is.finite(v)) return(NA_integer_); g <- 1L; while (g < nb && v > edges[g + 1L]) g <- g + 1L; g }, 1L)
  pos <- vapply(seq_len(nb), function(g) { v <- x[which(idx == g)]; if (length(v)) hp_quantile7(v, 0.5) else NA_real_ }, 1)
  list(idx = idx, pos = pos)
}

hp_eda_hist_edges <- function(v, n = 0, log = FALSE) {
  t <- v[is.finite(v) & (!log | v > 0)]; if (log) t <- log10(t)
  if (!length(t)) return(numeric())
  s <- sort(t); lo <- s[1L]; hi <- s[length(s)]
  k <- if (n > 0) n else {
    w <- 2 * (hp_quantile7(s, 0.75) - hp_quantile7(s, 0.25)) / length(s)^(1 / 3)
    if (w > 0 && hi > lo) min(max(ceiling((hi - lo) / w), 5), 40) else min(max(ceiling(log2(length(s))) + 1, 1), 30)
  }
  if (!(hi > lo)) { lo <- lo - 0.5; hi <- hi + 0.5 }
  e <- lo + (hi - lo) * (0:k) / k
  if (log) 10^e else e
}

hp_eda_pages <- function(p, panels, spec) {
  per <- max(spec$ncol * spec$nrow, 1L)
  lv <- levels(panels)
  if (length(lv) <= per) return(list(p))
  chunks <- split(lv, ceiling(seq_along(lv) / per))
  chunks <- chunks[seq_len(min(length(chunks), spec$max_pages))]
  lapply(chunks, function(ch) p %+% p$data[p$data$.PANEL %in% ch, , drop = FALSE])
}

#' Draw an exploratory plot
#'
#' Draws the figure described by a specification (see [hp_eda_spec()]) from a
#' prepared dataset (see [hp_eda_prepare()]). Plots of many individuals or many
#' panels are split over pages.
#'
#' @details
#' The plot uses the package theme ([hp_theme()]) and a caption that records the
#' catalogue identifier and the date. Lines and points are coloured by the
#' `color` column; a numeric colour or facet column with more than
#' `max_categories` distinct values is cut into four groups at its quartiles.
#' Errors are explicit: asking for a column the data do not have stops with a
#' message that lists the available columns.
#'
#' @param d An `hp_eda_data` object from [hp_eda_prepare()].
#' @param spec A specification from [hp_eda_spec()].
#' @param units A list with the elements `time`, `conc` and `dose`, used in the
#'   axis titles.
#' @param max_categories Most distinct values of a colour, facet or axis
#'   column shown as they are before a numeric column is grouped.
#' @param base_size Base font size of the theme.
#' @return A list of ggplot objects, one per page. Print an element to see it
#'   or write the list with [hp_save_plots()].
#' @seealso [hp_eda_spec()], [hp_eda_prepare()], [hp_save_plots()]
#' @family exploratory data analysis
#' @examples
#' d <- hp_read_data(hp_example_file())
#' e <- hp_eda_prepare(d, dvid = "DVID", lloq = 0.5, color = "DOSEGRP")
#'
#' # spaghetti plot on a log axis
#' p <- hp_eda_plot(e, hp_eda_spec(name = "spaghetti", kind = "line", x = "TIME",
#'                                 y = "DV", color = "DOSEGRP", log_y = TRUE,
#'                                 pk_only = TRUE),
#'                  units = list(time = "h", conc = "ng/mL", dose = "mg"))
#' length(p)
#' p[[1]]
#'
#' # exposure against a covariate
#' q <- hp_eda_plot(e, hp_eda_spec(kind = "scatter", x = "WT", y = "AUCLAST",
#'                                 scope = "subjects", smooth = "linear"))
#'
#' # one individual per panel, three by three
#' r <- hp_eda_plot(e, hp_eda_spec(kind = "line", x = "TIME", y = "DV",
#'                                 facet = "ID", pk_only = TRUE, ncol = 3, nrow = 3))
#' length(r)
#' @export
hp_eda_plot <- function(d, spec, units = list(time = "h", conc = "", dose = ""), max_categories = 20,
                        base_size = 11) {
  today <- format(Sys.Date())
  cap <- paste("hyperutils eda \u00b7", spec$id, "\u00b7", today)
  th <- hp_theme(base_size)
  tt <- function(p) p + th + ggplot2::labs(title = if (nzchar(spec$title)) spec$title else NULL, caption = cap) +
    (if (!isTRUE(spec$legend)) ggplot2::theme(legend.position = "none") else NULL)
  rw <- hp_eda_rows(d, spec); f <- rw$frame; keep <- rw$keep
  have <- function(v) v %in% names(f)
  chk <- function(v, what) if (!have(v)) hp_stop("[plot] ", spec$name, ": ", what, " `", v, "` is not a column here (the columns are: ", paste(names(f), collapse = ", "), ")")
  for (v in spec$x) chk(v, "x variable"); for (v in spec$y) chk(v, "y variable")
  if (spec$kind == "corr") {
    vars <- unique(c(spec$x, spec$y)); if (length(vars) < 2L) hp_stop("a corr plot needs at least two numeric variables")
    m <- f[keep, vars, drop = FALSE]
    cm <- stats::cor(m, use = "pairwise.complete.obs")
    long <- data.frame(A = factor(rep(vars, times = length(vars)), levels = vars),
                       B = factor(rep(vars, each = length(vars)), levels = rev(vars)), r = as.vector(cm))
    p <- ggplot2::ggplot(long, ggplot2::aes(A, B, fill = r)) + ggplot2::geom_tile(colour = "white") +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", r)), size = 3) +
      ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-1, 1)) +
      ggplot2::labs(x = NULL, y = NULL, fill = "Pearson r")
    return(list(tt(p) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))))
  }
  df <- f[keep, , drop = FALSE]
  if (!nrow(df)) hp_stop("no rows are left for `", spec$name, "`")
  grp <- if (!is.null(spec$color)) { chk(spec$color, "colour column"); hp_eda_levels(df[[spec$color]], max_categories, 4) } else NULL
  fac <- NULL
  if (!is.null(spec$facet)) {
    chk(spec$facet, "facet column")
    is_id <- identical(spec$facet, d$names$id)
    fac <- hp_eda_levels(df[[spec$facet]], if (is_id) 0 else max_categories, 4)
    levels(fac) <- paste(spec$facet, levels(fac))
  }
  if (spec$kind == "pkpd") return(hp_eda_pkpd(d, spec, df, grp, fac, units, tt))
  xl <- function(v) spec$xlab %||% hp_eda_label(v, units)
  yl <- function(v) spec$ylab %||% hp_eda_label(v, units)
  # the variable combinations of the panels
  combos <- NULL
  if (spec$kind == "pairs") {
    vars <- c(spec$x, spec$y); if (length(vars) < 2L) hp_stop("a pairs plot needs at least two variables")
    combos <- do.call(rbind, lapply(seq_along(vars)[-1L], function(i) data.frame(x = vars[seq_len(i - 1L)], y = vars[i], stringsAsFactors = FALSE)))
  } else if (spec$kind %in% c("histogram", "bar", "box") && (length(spec$x) > 1L || length(spec$y) > 1L)) {
    xs <- if (length(spec$x)) spec$x else NA_character_; ys <- if (length(spec$y)) spec$y else NA_character_
    combos <- expand.grid(x = xs, y = ys, stringsAsFactors = FALSE)
  } else {
    xs <- if (length(spec$x)) spec$x else NA_character_; ys <- if (length(spec$y)) spec$y else NA_character_
    combos <- expand.grid(x = xs, y = ys, stringsAsFactors = FALSE)
  }
  multi <- nrow(combos) > 1L
  parts <- lapply(seq_len(nrow(combos)), function(i) {
    xv <- combos$x[i]; yv <- combos$y[i]
    lab <- if (!is.na(xv) && !is.na(yv)) paste(yv, "vs", xv) else if (!is.na(xv)) xv else yv
    data.frame(.X = if (!is.na(xv)) df[[xv]] else NA_real_, .Y = if (!is.na(yv)) df[[yv]] else NA_real_,
               .G = if (is.null(grp)) factor("") else grp, .F = if (is.null(fac)) factor("") else fac,
               .BY = if (!is.null(spec$by)) as.character(df[[spec$by]]) else if (d$names$id %in% names(df)) as.character(df[[d$names$id]]) else "1",
               .VAR = lab, .XV = xv, .YV = yv, stringsAsFactors = FALSE)
  })
  long <- do.call(rbind, parts)
  long$.PANEL <- factor(trimws(paste(long$.VAR, if (!is.null(fac) && !multi) "" else if (!is.null(fac)) as.character(long$.F) else "")),
                        levels = unique(trimws(paste(long$.VAR, if (!is.null(fac) && !multi) "" else if (!is.null(fac)) as.character(long$.F) else ""))))
  if (!is.null(fac) && !multi) long$.PANEL <- long$.F
  if (is.null(fac) && !multi) long$.PANEL <- factor("")
  has_g <- !is.null(grp)
  scales <- switch(spec$free, x = "free_x", y = "free_y", both = "free", "fixed")
  if (multi || spec$kind == "pairs") scales <- "free"
  facet_layer <- if (nlevels(long$.PANEL) > 1L) ggplot2::facet_wrap(ggplot2::vars(.PANEL), scales = scales, ncol = spec$ncol) else NULL
  fillname <- if (has_g) toupper(spec$color) else NULL
  base <- ggplot2::ggplot(long)
  lims <- function(p) {
    if (!all(is.na(spec$xlim))) p <- p + ggplot2::coord_cartesian(xlim = spec$xlim, ylim = if (all(is.na(spec$ylim))) NULL else spec$ylim)
    else if (!all(is.na(spec$ylim))) p <- p + ggplot2::coord_cartesian(ylim = spec$ylim)
    p
  }
  rl <- function(p) {
    for (h in spec$hlines) p <- p + ggplot2::geom_hline(yintercept = h, linetype = "dashed", colour = "grey40")
    for (v in spec$vlines) p <- p + ggplot2::geom_vline(xintercept = v, linetype = "dashed", colour = "grey40")
    if (isTRUE(spec$identity)) p <- p + ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40")
    p
  }
  logs <- function(p) {
    if (spec$log_x) p <- p + ggplot2::scale_x_log10()
    if (spec$log_y) p <- p + ggplot2::scale_y_log10()
    p
  }
  k <- spec$kind
  if (k %in% c("scatter", "pairs")) {
    pd <- long[is.finite(hp_num(long$.X)) & is.finite(hp_num(long$.Y)), , drop = FALSE]
    pd$.X <- hp_num(pd$.X); pd$.Y <- hp_num(pd$.Y)
    if (spec$log_x) pd <- pd[pd$.X > 0, , drop = FALSE]
    if (spec$log_y) pd <- pd[pd$.Y > 0, , drop = FALSE]
    aes_ <- if (has_g) ggplot2::aes(.X, .Y, colour = .G) else ggplot2::aes(.X, .Y)
    p <- ggplot2::ggplot(pd, aes_)
    if (spec$points) p <- p + ggplot2::geom_point(alpha = 0.65, size = 1.2,
                                                    position = if (spec$jitter > 0) ggplot2::position_jitter(width = spec$jitter, height = 0, seed = 1) else "identity")
    sm <- spec$smooth
    if (sm %in% c("linear", "loess")) {
      meth <- if (sm == "linear") "lm" else "loess"
      sa <- if (has_g && isTRUE(spec$smooth_groups)) ggplot2::aes(.X, .Y, colour = .G) else ggplot2::aes(.X, .Y)
      p <- p + suppressWarnings(do.call(ggplot2::geom_smooth, c(
        list(sa, method = meth, formula = y ~ x, se = FALSE, span = if (sm == "loess") spec$span else NULL, linewidth = 0.8),
        hp_fixed_colour(!(has_g && isTRUE(spec$smooth_groups)), "black"))))
    } else if (sm == "power") {
      pw <- do.call(rbind, lapply(split(pd, list(pd$.PANEL, if (has_g && isTRUE(spec$smooth_groups)) pd$.G else 1), drop = TRUE), function(z) {
        z <- z[z$.X > 0 & z$.Y > 0, ]; if (nrow(z) < 3L || length(unique(z$.X)) < 2L) return(NULL)
        m <- stats::lm(log(.Y) ~ log(.X), z); xs <- exp(seq(log(min(z$.X)), log(max(z$.X)), length.out = 50))
        data.frame(.X = xs, .Y = exp(stats::predict(m, data.frame(.X = xs))), .G = z$.G[1L], .PANEL = z$.PANEL[1L])
      }))
      if (!is.null(pw)) p <- p + do.call(ggplot2::geom_line, c(list(data = pw, mapping = if (has_g) ggplot2::aes(.X, .Y, colour = .G) else ggplot2::aes(.X, .Y), linewidth = 0.8), hp_fixed_colour(!has_g, "black")))
    }
    p <- p + ggplot2::labs(x = if (spec$kind == "pairs") NULL else if (multi) NULL else xl(spec$x[1L]),
                           y = if (multi) NULL else yl(spec$y[1L]), colour = fillname)
    if (spec$kind == "pairs") {
      pd$.PANEL <- factor(paste(pd$.YV, "vs", pd$.XV)); p <- p %+% pd
    }
    p <- p + facet_layer
    p <- tt(lims(rl(logs(p))))
    return(hp_eda_pages(p, long$.PANEL, spec))
  }
  if (k == "line") {
    pd <- long[is.finite(hp_num(long$.X)) & is.finite(hp_num(long$.Y)), , drop = FALSE]
    if (spec$log_y) pd <- pd[pd$.Y > 0, , drop = FALSE]
    pd <- pd[order(pd$.BY, pd$.X), , drop = FALSE]
    if (spec$max_groups > 0) { keepb <- utils::head(unique(pd$.BY), spec$max_groups); pd <- pd[pd$.BY %in% keepb, , drop = FALSE] }
    pd$.LG <- interaction(pd$.BY, pd$.PANEL, drop = TRUE)
    p <- ggplot2::ggplot(pd, if (has_g) ggplot2::aes(.X, .Y, group = .LG, colour = .G) else ggplot2::aes(.X, .Y, group = .LG))
    if (spec$lines) p <- p + do.call(ggplot2::geom_line, c(list(alpha = 0.55, linewidth = 0.45), hp_fixed_colour(!has_g, "grey25")))
    if (spec$points) p <- p + do.call(ggplot2::geom_point, c(list(alpha = 0.6, size = 0.9), hp_fixed_colour(!has_g, "grey25")))
    p <- p + ggplot2::labs(x = if (multi) NULL else xl(spec$x[1L]), y = if (multi) NULL else yl(spec$y[1L]), colour = fillname) + facet_layer
    p <- tt(lims(rl(logs(p))))
    return(hp_eda_pages(p, long$.PANEL, spec))
  }
  if (k == "summary") {
    pd <- long[is.finite(hp_num(long$.X)) & is.finite(hp_num(long$.Y)), , drop = FALSE]
    if (spec$log_y) pd <- pd[pd$.Y > 0, , drop = FALSE]
    pd$.X <- hp_num(pd$.X)
    out <- list()
    for (pn in levels(pd$.PANEL)) {
      z <- pd[pd$.PANEL == pn, , drop = FALSE]; if (!nrow(z)) next
      bn <- hp_eda_xbins(z$.X, spec$bins); z$.B <- bn$idx
      for (g in levels(droplevels(z$.G))) {
        zz <- z[z$.G == g & !is.na(z$.B), , drop = FALSE]
        for (b in sort(unique(zz$.B))) {
          v <- zz$.Y[zz$.B == b]; cs <- hp_centre_spread(v, spec$centre, spec$spread)
          out[[length(out) + 1L]] <- data.frame(.X = bn$pos[b], .C = cs[["centre"]], .L = cs[["lower"]], .U = cs[["upper"]], .G = g, .PANEL = pn)
        }
      }
    }
    sd_ <- do.call(rbind, out); if (is.null(sd_)) hp_stop("nothing to summarise in `", spec$name, "`")
    sd_$.PANEL <- factor(sd_$.PANEL, levels = levels(pd$.PANEL)); sd_$.G <- factor(sd_$.G, levels = levels(pd$.G))
    p <- ggplot2::ggplot(sd_, if (has_g) ggplot2::aes(.X, .C, colour = .G, group = .G) else ggplot2::aes(.X, .C, group = .G))
    if (spec$lines) p <- p + ggplot2::geom_line()
    if (spec$spread != "none") p <- p + ggplot2::geom_errorbar(ggplot2::aes(ymin = .L, ymax = .U), width = 0, na.rm = TRUE)
    if (spec$points) p <- p + ggplot2::geom_point(size = 1.6)
    p <- p + ggplot2::labs(x = if (multi) NULL else xl(spec$x[1L]),
                           y = if (multi) NULL else yl(spec$y[1L]), colour = fillname) + facet_layer
    return(hp_eda_pages(tt(lims(rl(logs(p)))), long$.PANEL, spec))
  }
  if (k == "box") {
    pd <- long[is.finite(hp_num(long$.Y)), , drop = FALSE]; pd$.Y <- hp_num(pd$.Y)
    if (spec$log_y) pd <- pd[pd$.Y > 0, , drop = FALSE]
    xl_ <- vapply(seq_len(nrow(pd)), function(i) NA_character_, "")
    pd$.XC <- hp_eda_levels(pd$.X, max_categories, 4)
    p <- ggplot2::ggplot(pd, if (has_g) ggplot2::aes(.XC, .Y, fill = .G) else ggplot2::aes(.XC, .Y)) +
      ggplot2::geom_boxplot(outlier.size = 0.8, alpha = 0.75) +
      ggplot2::labs(x = if (multi) NULL else xl(spec$x[1L]), y = if (multi) NULL else yl(spec$y[1L]), fill = fillname) + facet_layer
    return(hp_eda_pages(tt(lims(rl(logs(p)))), long$.PANEL, spec))
  }
  if (k == "histogram") {
    pd <- long; pd$.X <- hp_num(pd$.X); pd <- pd[is.finite(pd$.X) & (!spec$log_x | pd$.X > 0), , drop = FALSE]
    if (!nrow(pd)) hp_stop("nothing to draw")
    br <- hp_eda_hist_edges(pd$.X, spec$hist_bins, spec$log_x)
    p <- ggplot2::ggplot(pd, if (has_g) ggplot2::aes(.X, fill = .G) else ggplot2::aes(.X)) +
      ggplot2::geom_histogram(breaks = br, colour = "white", alpha = 0.8, position = "identity") +
      ggplot2::labs(x = if (multi) NULL else xl(spec$x[1L]), y = "Count", fill = fillname) + facet_layer
    if (spec$log_x) p <- p + ggplot2::scale_x_log10()
    return(hp_eda_pages(tt(lims(rl(p))), long$.PANEL, spec))
  }
  if (k == "bar") {
    pd <- long; pd$.XC <- hp_eda_levels(pd$.X, max_categories, 4)
    if (!all(is.na(pd$.Y))) {
      pd <- pd[is.finite(hp_num(pd$.Y)), , drop = FALSE]; pd$.Y <- hp_num(pd$.Y)
      ag <- do.call(rbind, lapply(split(pd, list(pd$.XC, pd$.G, pd$.PANEL), drop = TRUE), function(z)
        data.frame(.XC = z$.XC[1L], .G = z$.G[1L], .PANEL = z$.PANEL[1L], .C = hp_centre_spread(z$.Y, spec$centre, "none")[["centre"]])))
      p <- ggplot2::ggplot(ag, if (has_g) ggplot2::aes(.XC, .C, fill = .G) else ggplot2::aes(.XC, .C)) +
        ggplot2::geom_col(position = "dodge") + ggplot2::labs(x = xl(spec$x[1L]), y = yl(spec$y[1L]), fill = fillname) + facet_layer
    } else {
      p <- ggplot2::ggplot(pd, if (has_g) ggplot2::aes(.XC, fill = .G) else ggplot2::aes(.XC)) +
        ggplot2::geom_bar(position = "dodge") + ggplot2::labs(x = if (multi) NULL else xl(spec$x[1L]), y = "Count", fill = fillname) + facet_layer
    }
    return(hp_eda_pages(tt(lims(rl(p))), long$.PANEL, spec))
  }
  hp_stop("plot kind `", k, "` is not known (scatter | line | summary | box | histogram | bar | corr | pairs | pkpd)")
}

# PD response against the time-matched PK concentration
hp_eda_pkpd <- function(d, spec, df, grp, fac, units, tt) {
  nm <- d$names; f <- d$rec
  if (is.null(nm$dvid)) hp_stop("a pkpd plot needs a DVID column to tell PK from PD")
  pd <- d$info$pd_dvid; if (is.null(pd)) hp_stop("a pkpd plot needs a second DVID value (the PD response)")
  dvt <- ifelse(is.finite(hp_num(f[[nm$dvid]])), hp_num_text(hp_num(f[[nm$dvid]])), as.character(f[[nm$dvid]]))
  pkr <- f[f$PKOBS == 1 & is.finite(hp_num(f[[nm$dv]])), , drop = FALSE]
  conc_at <- function(id, t) {
    z <- pkr[as.character(pkr[[nm$id]]) == id, , drop = FALSE]; z <- z[order(z[[nm$time]]), ]
    tm <- z[[nm$time]]; if (!length(tm)) return(NA_real_)
    h <- which(abs(tm - t) <= 1e-6 * (1 + abs(t))); if (length(h)) return(z[[nm$dv]][h[1L]])
    i <- which(tm > t)[1L]; if (is.na(i) || i == 1L) return(NA_real_)
    z[[nm$dv]][i - 1L] + (z[[nm$dv]][i] - z[[nm$dv]][i - 1L]) * (t - tm[i - 1L]) / (tm[i] - tm[i - 1L])
  }
  sel <- which(df$EVENT == "observation" & dvt[match(rownames(df), rownames(f))] == pd)
  if (!length(sel)) hp_stop("no PD observation (DVID ", pd, ") could be matched with a PK concentration")
  rows <- lapply(sel, function(i) {
    id <- as.character(df[[nm$id]][i]); t <- df[[nm$time]][i]
    data.frame(.ID = id, .T = t, .C = conc_at(id, t), .Y = df[[nm$dv]][i], .G = if (is.null(grp)) factor("") else grp[i],
               .F = if (is.null(fac)) factor("") else fac[i])
  })
  z <- do.call(rbind, rows); z <- z[is.finite(z$.C) & is.finite(z$.Y), , drop = FALSE]
  if (spec$log_x) z <- z[z$.C > 0, , drop = FALSE]; if (spec$log_y) z <- z[z$.Y > 0, , drop = FALSE]
  if (!nrow(z)) hp_stop("no PD observation (DVID ", pd, ") could be matched with a PK concentration")
  z <- z[order(z$.ID, z$.T), , drop = FALSE]
  has_g <- !is.null(grp)
  p <- ggplot2::ggplot(z, if (has_g) ggplot2::aes(.C, .Y, colour = .G, group = .ID) else ggplot2::aes(.C, .Y, group = .ID))
  if (spec$lines) p <- p + do.call(ggplot2::geom_path, c(list(alpha = 0.45), hp_fixed_colour(!has_g, "grey25")))
  if (spec$points) p <- p + do.call(ggplot2::geom_point, c(list(alpha = 0.65, size = 1.2), hp_fixed_colour(!has_g, "grey25")))
  p <- p + ggplot2::labs(x = spec$xlab %||% hp_eda_label(spec$x[1L] %||% "Concentration", units), y = spec$ylab %||% "PD response",
                         colour = if (has_g) toupper(spec$color) else NULL)
  if (!is.null(fac)) p <- p + ggplot2::facet_wrap(ggplot2::vars(.F), scales = "free")
  if (spec$log_x) p <- p + ggplot2::scale_x_log10(); if (spec$log_y) p <- p + ggplot2::scale_y_log10()
  list(tt(p))
}
