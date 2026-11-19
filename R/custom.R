# Custom plots and tables of any data frame: the data steps, the ggplot2-style
# plots and the summary tables of a custom analysis.

#' Data steps of a custom analysis
#'
#' A small, predictable data-preparation pipeline for any data frame. The steps
#' always run in the same order, whatever the order of the arguments:
#' `mutate`, `filter`, `distinct`, `sort`, `head`, `select`. It is meant for
#' the step between reading data and drawing a figure or building a table, and
#' it keeps that step in one documented call.
#'
#' @param data A data frame.
#' @param mutate Statements in R syntax that create or change columns, one
#'   per element, for example `"LNDV = log(DV)"`. Each can use the columns made
#'   by the previous ones. `case_when()` (see [hp_case_when()]) and `ifelse()` are available.
#' @param filter Conditions in R syntax; only rows for which every condition
#'   holds are kept. Missing results count as FALSE.
#' @param distinct Columns; the first row of each combination is kept (for
#'   example `"ID"` keeps one row per individual).
#' @param sort Columns to sort by; a leading `-` sorts downwards.
#' @param head Number of rows kept from the top (`NULL` keeps all).
#' @param select Columns to keep, in this order.
#' @return A data frame.
#' @seealso [hp_mutate()], [hp_filter()], [hp_case_when()],
#'   [hp_custom_plot()], [hp_custom_table()]
#' @family custom plots and tables
#' @examples
#' d <- hp_read_data(hp_example_file())
#'
#' # one row per individual with a derived column
#' subj <- hp_custom_prep(d, filter = "EVID == 1", distinct = "ID",
#'                        mutate = "DOSEFLAG = case_when(AMT >= 200 ~ 'high', TRUE ~ 'low')",
#'                        select = c("ID", "AMT", "DOSEFLAG"))
#' head(subj)
#'
#' # the five highest concentrations
#' hp_custom_prep(d, filter = "EVID == 0", sort = "-DV", head = 5,
#'                select = c("ID", "TIME", "DV"))
#' @export
hp_custom_prep <- function(data, mutate = character(), filter = character(), distinct = character(),
                           sort = character(), head = NULL, select = character()) {
  if (length(mutate)) data <- hp_mutate(data, paste(mutate, collapse = "; "))
  for (f in filter) data <- hp_filter(data, f)
  chk <- function(cols, what) {
    miss <- setdiff(cols, names(data))
    if (length(miss)) hp_stop(what, " `", miss[1L], "` is not a column (the columns are: ", paste(names(data), collapse = ", "), ")")
  }
  if (length(distinct)) {
    chk(distinct, "distinct")
    key <- do.call(paste, c(lapply(data[distinct], as.character), sep = "\r"))
    data <- data[!duplicated(key), , drop = FALSE]
  }
  if (length(sort)) {
    desc <- startsWith(sort, "-"); cols <- sub("^-", "", sort); chk(cols, "sort")
    ord <- do.call(order, lapply(seq_along(cols), function(i) {
      v <- data[[cols[i]]]; if (desc[i]) { if (is.numeric(v)) -v else -xtfrm(v) } else v
    }))
    data <- data[ord, , drop = FALSE]
  }
  if (!is.null(head)) data <- utils::head(data, head)
  if (length(select)) { chk(select, "select"); data <- data[select] }
  rownames(data) <- NULL
  data
}

hp_color_values <- function(name) {
  switch(name,
    "okabe-ito" = c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#F0E442", "#999999"),
    tableau = c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC"),
    set1 = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999"),
    dark2 = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666"),
    pastel = c("#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2"),
    grey = grDevices::gray.colors(8, 0.2, 0.8),
    NULL)
}

# a column or an expression of the data, evaluated like a mutate
hp_custom_eval <- function(data, expr) {
  if (is.null(expr) || !nzchar(expr)) return(NULL)
  v <- tryCatch(eval(parse(text = expr)[[1L]], envir = hp_eval_env(data, parent = baseenv())),
                error = function(e) hp_stop("`", expr, "`: ", conditionMessage(e)))
  if (length(v) == 1L) v <- rep(v, nrow(data))
  v
}

#' Draw a custom plot from a list of settings
#'
#' Builds a ggplot2 figure from any data frame with a list of plain settings,
#' so that a figure is described in a few readable lines: which columns go on
#' the axes, the colour, the shape and the facets, which layers to draw, the
#' scales, the reference lines and the look. It covers the figures that
#' clinical pharmacology reports use most (spaghetti plots, box plots,
#' distributions, smoothers, summary plots with error bars, heat maps) and
#' returns an ordinary ggplot object that you can extend with `+`.
#'
#' @details
#' **Layers (`geoms`).** One or several of `point`, `jitter`, `line`, `path`,
#' `step`, `smooth`, `boxplot`, `violin`, `bar`, `col`, `histogram`,
#' `density`, `area`, `ecdf`, `qq`, `errorbar`, `pointrange`, `tile` and
#' `text`, drawn in the order given.
#'
#' **Aesthetics.** `x`, `y`, `color`, `fill`, `shape`, `size`, `linetype`,
#' `group` and `label` accept a column name or an R expression of the columns
#' (for example `"factor(DOSE)"` or `"log10(DV)"`). Text and few-valued
#' columns become discrete, so colours and groups are legends.
#'
#' **Fields of the list `def`.**
#'
#' - Layers and mapping: `geoms`, `x`, `y`, `color`, `fill`, `shape`, `size`,
#'   `linetype`, `group`, `label`.
#' - Layer settings: `alpha`, `line_width`, `bar_width`, `binwidth`, `nbins`,
#'   `jitter`, `position` (`"identity"`, `"stack"`, `"dodge"` or `"fill"`),
#'   `outliers`, `smooth_loess`, `smooth_se`, `smooth_level`, `span`,
#'   `summary_fun` (`"mean"`, `"median"`, `"sum"`, `"min"`, `"max"`, `"n"`,
#'   `"sd"`, `"q1"`, `"q3"` or `"geomean"`) and `summary_error`.
#' - Scales: `x_log`, `y_log`, `xlim`, `ylim`, `x_breaks`, `y_breaks`,
#'   `x_angle`, `flip`, `equal_axes`, `identity` (line of identity),
#'   `palette` (`"okabe-ito"`, `"tableau"`, `"set1"`, `"dark2"`,
#'   `"pastel"`, `"grey"` or `"viridis"`).
#' - Panels: `facet_wrap`, `facet_rows`, `facet_cols`, `facet_nrow`,
#'   `facet_ncol`, `scales` (`"fixed"`, `"free"`, `"free_x"` or `"free_y"`),
#'   `labels_both`.
#' - Reference lines: `hline`, `vline`, `abline`, `ref_color`,
#'   `ref_linetype`, `ref_width`.
#' - Text: `title`, `subtitle`, `caption`, `x_title`, `y_title`,
#'   `legend_title`, `color_title`, `fill_title`, `shape_title`,
#'   `size_title`, `linetype_title`, `title_align`, `title_size`,
#'   `axis_title_size`, `axis_text_size`, `legend_text_size`,
#'   `strip_text_size`.
#' - Look: `theme` (`"classic"`, `"minimal"`, `"bw"`, `"light"`, `"grey"` or
#'   `"void"`), `legend` (`"right"`, `"bottom"`, `"top"`, `"left"` or
#'   `"none"`), `grid_color`, `panel_background`, `plot_background`.
#'
#' Unknown fields are ignored, so one list can be shared by several plots.
#'
#' @param data A data frame, usually the result of [hp_custom_prep()].
#' @param def A named list of settings, see Details.
#' @param seed Seed of the random jitter, for reproducible figures.
#' @param base_size Base font size of the theme.
#' @return A ggplot object.
#' @seealso [hp_custom_prep()], [hp_custom_table()], [hp_save_plot()],
#'   [hp_eda_plot()]
#' @family custom plots and tables
#' @examples
#' d <- hp_read_data(hp_example_file())
#' obs <- hp_custom_prep(d, filter = "EVID == 0")
#'
#' # individual profiles on a log axis, coloured by dose group
#' hp_custom_plot(obs, list(geoms = c("point", "line"), x = "TIME", y = "DV",
#'                          group = "ID", color = "DOSEGRP", y_log = TRUE,
#'                          x_title = "Time (h)", y_title = "Concentration (ng/mL)",
#'                          legend_title = "Dose group"))
#'
#' # box plot of an exposure measure by group, with the points on top
#' subj <- hp_custom_prep(d, filter = "EVID == 1", distinct = "ID")
#' hp_custom_plot(subj, list(geoms = c("boxplot", "jitter"), x = "DOSEGRP",
#'                           y = "AMT", fill = "DOSEGRP", legend = "none"))
#'
#' # distribution of a covariate, one panel per sex
#' hp_custom_plot(subj, list(geoms = "histogram", x = "WT", fill = "SEX",
#'                           facet_wrap = "SEX", nbins = 10, theme = "minimal"))
#'
#' # the result is a ggplot, so it can be extended
#' p <- hp_custom_plot(obs, list(geoms = "point", x = "TIME", y = "DV"))
#' p + ggplot2::labs(caption = "Source: example data")
#' @export
hp_custom_plot <- function(data, def, seed = 1, base_size = 11) {
  g <- function(n, default = NULL) { v <- def[[n]]; if (is.null(v) || (length(v) == 1L && is.na(v))) default else v }
  isTRUE_ <- function(n) isTRUE(g(n, FALSE))
  geoms <- tolower(g("geoms", "point"))
  df <- data
  aesl <- list()
  for (a in c("x", "y", "color", "fill", "shape", "size", "linetype", "group", "label")) {
    v <- hp_custom_eval(data, g(a))
    if (is.null(v)) next
    if (a %in% c("color", "fill", "shape", "linetype", "group") && !is.numeric(v)) v <- factor(v, levels = unique(sort(v)))
    if (a %in% c("shape", "linetype", "group") && is.numeric(v)) v <- factor(v)
    df[[paste0(".", a)]] <- v
    aesl[[a]] <- as.name(paste0(".", a))
  }
  if (!is.null(aesl$color) && is.numeric(df$.color) && length(unique(df$.color)) <= 8L && g("palette", "okabe-ito") != "viridis")
    df$.color <- factor(df$.color)
  if (!is.null(aesl$x) && !is.numeric(df$.x)) df$.x <- factor(df$.x, levels = unique(sort(df$.x)))
  names(aesl)[names(aesl) == "color"] <- "colour"
  amap <- do.call(ggplot2::aes, aesl)
  p <- ggplot2::ggplot(df, amap)
  alpha <- g("alpha"); lw <- g("line_width", 0.6)
  opt <- function(...) { l <- list(...); l[!vapply(l, is.null, TRUE)] }
  pos <- switch(g("position", "identity"), stack = "stack", dodge = ggplot2::position_dodge(), fill = "fill", "identity")
  has <- function(n) !is.null(aesl[[n]])
  sumfun <- function(f) switch(f, mean = mean, median = stats::median, sum = sum, min = min, max = max,
                               n = length, sd = stats::sd, q1 = function(x) stats::quantile(x, 0.25, type = 7),
                               q3 = function(x) stats::quantile(x, 0.75, type = 7),
                               geomean = function(x) exp(mean(log(x[x > 0]))), mean)
  for (gm in geoms) {
    p <- p + switch(gm,
      point = do.call(ggplot2::geom_point, opt(alpha = alpha %||% 0.7, size = if (!has("size")) 1.6 else NULL)),
      jitter = ggplot2::geom_point(alpha = alpha %||% 0.7, size = 1.6,
                                   position = ggplot2::position_jitter(width = g("jitter", 0.2), height = 0, seed = seed)),
      line = do.call(ggplot2::geom_line, opt(alpha = alpha, linewidth = lw)),
      path = do.call(ggplot2::geom_path, opt(alpha = alpha, linewidth = lw)),
      step = do.call(ggplot2::geom_step, opt(alpha = alpha, linewidth = lw)),
      smooth = ggplot2::geom_smooth(method = if (isTRUE_("smooth_loess")) "loess" else "lm", formula = y ~ x,
                                    se = isTRUE_("smooth_se"), span = g("span", 0.75), level = g("smooth_level", 0.95),
                                    linewidth = lw),
      boxplot = do.call(ggplot2::geom_boxplot, opt(alpha = alpha %||% 0.75, outlier.shape = if (!isTRUE(g("outliers", TRUE))) NA else 19,
                                                  width = g("bar_width"))),
      violin = ggplot2::geom_violin(alpha = alpha %||% 0.6, scale = "width"),
      bar = if (has("y")) ggplot2::stat_summary(fun = sumfun(g("summary_fun", "mean")), geom = "col", position = pos,
                                                width = g("bar_width", 0.8), alpha = alpha %||% 0.85)
            else do.call(ggplot2::geom_bar, opt(position = pos, width = g("bar_width", 0.8), alpha = alpha %||% 0.85)),
      col = do.call(ggplot2::geom_col, opt(position = pos, width = g("bar_width", 0.8), alpha = alpha %||% 0.85)),
      histogram = do.call(ggplot2::geom_histogram, opt(bins = if (is.null(g("binwidth"))) g("nbins", 30L) else NULL,
                                                      binwidth = g("binwidth"), position = if (pos == "identity") "stack" else pos,
                                                      alpha = alpha %||% 0.85, colour = "white")),
      density = ggplot2::geom_density(alpha = alpha %||% 0.4, linewidth = lw),
      area = do.call(ggplot2::geom_area, opt(alpha = alpha %||% 0.5, position = pos)),
      ecdf = ggplot2::stat_ecdf(linewidth = lw),
      qq = list(ggplot2::geom_qq(ggplot2::aes(sample = if (has("x")) .x else .y), inherit.aes = FALSE),
                ggplot2::geom_qq_line(ggplot2::aes(sample = if (has("x")) .x else .y), colour = "grey40", inherit.aes = FALSE)),
      errorbar = ggplot2::stat_summary(fun.data = function(x) {
        m <- mean(x); s <- if (length(x) > 1L) stats::sd(x) else 0; n <- length(x)
        e <- switch(g("summary_error", "sd"), se = s / sqrt(n), sem = s / sqrt(n),
                    ci95 = if (n > 1L) stats::qt(0.975, n - 1L) * s / sqrt(n) else 0, s)
        data.frame(y = m, ymin = m - e, ymax = m + e)
      }, geom = "errorbar", width = 0.2),
      pointrange = ggplot2::stat_summary(fun.data = function(x) {
        m <- mean(x); s <- if (length(x) > 1L) stats::sd(x) else 0
        data.frame(y = m, ymin = m - s, ymax = m + s)
      }, geom = "pointrange"),
      tile = ggplot2::geom_tile(),
      text = ggplot2::geom_text(size = 3),
      hp_stop("geom `", gm, "` is not known"))
  }
  # reference lines
  rc <- g("ref_color", "grey40"); rt <- g("ref_linetype", "dashed"); rw <- g("ref_width", 0.5)
  for (h in g("hline", numeric())) p <- p + ggplot2::geom_hline(yintercept = h, colour = rc, linetype = rt, linewidth = rw)
  for (v in g("vline", numeric())) p <- p + ggplot2::geom_vline(xintercept = v, colour = rc, linetype = rt, linewidth = rw)
  ab <- g("abline", list())
  for (a in ab) p <- p + ggplot2::geom_abline(intercept = a[[1L]], slope = a[[2L]], colour = rc, linetype = rt, linewidth = rw)
  if (isTRUE_("identity")) p <- p + ggplot2::geom_abline(intercept = 0, slope = 1, colour = rc, linetype = rt, linewidth = rw)
  # facets
  lab_both <- isTRUE_("labels_both")
  lbl <- if (lab_both) ggplot2::label_both else ggplot2::label_value
  fcol <- function(cols, tag) {
    if (!length(cols)) return(NULL)
    v <- do.call(paste, c(lapply(cols, function(cn) as.character(hp_custom_eval(data, cn))), sep = ", "))
    v
  }
  sc <- switch(g("scales", "fixed"), free = "free", freex = "free_x", free_x = "free_x", freey = "free_y", free_y = "free_y", "fixed")
  fr <- g("facet_rows", character()); fc <- g("facet_cols", character()); fw <- g("facet_wrap", character())
  if (length(fw)) {
    p$data$.fw <- fcol(fw); p <- p + ggplot2::facet_wrap(ggplot2::vars(.fw), scales = sc, ncol = g("facet_ncol"), nrow = g("facet_nrow"),
                                                        labeller = ggplot2::as_labeller(function(x) if (lab_both) paste0(paste(fw, collapse = ", "), ": ", x) else x))
  } else if (length(fr) || length(fc)) {
    if (length(fr)) p$data$.fr <- fcol(fr)
    if (length(fc)) p$data$.fc <- fcol(fc)
    p <- p + ggplot2::facet_grid(rows = if (length(fr)) ggplot2::vars(.fr) else NULL, cols = if (length(fc)) ggplot2::vars(.fc) else NULL, scales = sc)
  }
  # scales
  if (isTRUE_("x_log")) p <- p + ggplot2::scale_x_log10()
  if (isTRUE_("y_log")) p <- p + ggplot2::scale_y_log10()
  xb <- g("x_breaks", numeric()); yb <- g("y_breaks", numeric())
  if (length(xb) && !isTRUE_("x_log")) p <- p + ggplot2::scale_x_continuous(breaks = xb)
  if (length(yb) && !isTRUE_("y_log")) p <- p + ggplot2::scale_y_continuous(breaks = yb)
  xl <- g("xlim", c(NA, NA)); yl <- g("ylim", c(NA, NA))
  if (isTRUE_("flip")) p <- p + ggplot2::coord_flip(xlim = if (all(is.na(xl))) NULL else xl, ylim = if (all(is.na(yl))) NULL else yl)
  else if (isTRUE_("equal_axes")) p <- p + ggplot2::coord_fixed(xlim = if (all(is.na(xl))) NULL else xl, ylim = if (all(is.na(yl))) NULL else yl)
  else if (!all(is.na(xl)) || !all(is.na(yl))) p <- p + ggplot2::coord_cartesian(xlim = if (all(is.na(xl))) NULL else xl, ylim = if (all(is.na(yl))) NULL else yl)
  pal <- hp_color_values(g("palette", "okabe-ito"))
  discrete_col <- !is.null(df$.color) && !is.numeric(df$.color)
  discrete_fil <- !is.null(df$.fill) && !is.numeric(df$.fill)
  if (!is.null(pal)) {
    if (discrete_col) p <- p + ggplot2::scale_colour_manual(values = rep_len(pal, nlevels(df$.color)))
    if (discrete_fil) p <- p + ggplot2::scale_fill_manual(values = rep_len(pal, nlevels(df$.fill)))
  } else if (g("palette", "") == "viridis") {
    if (discrete_col) p <- p + ggplot2::scale_colour_viridis_d() else if (!is.null(df$.color)) p <- p + ggplot2::scale_colour_viridis_c()
    if (discrete_fil) p <- p + ggplot2::scale_fill_viridis_d() else if (!is.null(df$.fill)) p <- p + ggplot2::scale_fill_viridis_c()
  }
  # labels
  p <- p + ggplot2::labs(title = if (nzchar(g("title", ""))) g("title") else NULL,
                         subtitle = if (nzchar(g("subtitle", ""))) g("subtitle") else NULL,
                         caption = if (nzchar(g("caption", ""))) g("caption") else NULL,
                         x = g("x_title", g("x")), y = g("y_title", g("y")),
                         colour = g("color_title", g("legend_title", g("color"))), fill = g("fill_title", g("legend_title", g("fill"))),
                         shape = g("shape_title", g("shape")), size = g("size_title", g("size")),
                         linetype = g("linetype_title", g("linetype")))
  th <- switch(g("theme", "classic"), bw = ggplot2::theme_bw(base_size), minimal = ggplot2::theme_minimal(base_size),
               gray = , grey = ggplot2::theme_gray(base_size), ggplot2::theme_classic(base_size))
  p <- p + th + ggplot2::theme(legend.position = g("legend", "right"))
  ang <- g("x_angle"); if (!is.null(ang)) p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = ang, hjust = if (ang > 0) 1 else 0.5))
  if (!is.null(g("title_size"))) p <- p + ggplot2::theme(plot.title = ggplot2::element_text(size = g("title_size")))
  if (!is.null(g("axis_title_size"))) p <- p + ggplot2::theme(axis.title = ggplot2::element_text(size = g("axis_title_size")))
  if (!is.null(g("axis_text_size"))) p <- p + ggplot2::theme(axis.text = ggplot2::element_text(size = g("axis_text_size")))
  if (!is.null(g("legend_text_size"))) p <- p + ggplot2::theme(legend.text = ggplot2::element_text(size = g("legend_text_size")))
  if (!is.null(g("strip_text_size"))) p <- p + ggplot2::theme(strip.text = ggplot2::element_text(size = g("strip_text_size")))
  ta <- g("title_align", "left")
  p <- p + ggplot2::theme(plot.title = ggplot2::element_text(hjust = switch(ta, center = 0.5, right = 1, 0)))
  if (!is.null(g("plot_background"))) p <- p + ggplot2::theme(plot.background = ggplot2::element_rect(fill = g("plot_background")))
  if (!is.null(g("panel_background"))) p <- p + ggplot2::theme(panel.background = ggplot2::element_rect(fill = g("panel_background")))
  if (identical(g("grid_color"), "none")) p <- p + ggplot2::theme(panel.grid = ggplot2::element_blank())
  p
}

# --------------------------------------------------------------------------
# tables
# --------------------------------------------------------------------------

hp_stat_label <- function(tok) {
  t <- tolower(trimws(tok))
  l <- c(n = "N", nmiss = "N missing", nmissing = "N missing", mean = "Mean", sd = "SD", se = "SE", sem = "SE", cv = "CV%",
         geomean = "Geo. mean", geocv = "Geo. CV%", geosd = "Geo. CV%", median = "Median", min = "Min", max = "Max", sum = "Sum",
         q1 = "Q1", q3 = "Q3", iqr = "IQR", lcl95 = "95% CI lower", ucl95 = "95% CI upper", mean_sd = "Mean (SD)",
         mean_se = "Mean (SE)", mean_ci = "Mean (95% CI)", median_range = "Median [Min, Max]", median_iqr = "Median [Q1, Q3]",
         geomean_cv = "Geo. mean (CV%)", range = "Range")
  if (t %in% names(l)) return(l[[t]])
  if (grepl("^p[0-9.]+$", t)) return(toupper(t))
  tok
}

# a number with `digits` significant digits, trailing zeros kept (digits >= 100:
# 100 + the number of decimals), "NA" when missing
hp_sig <- function(v, digits = 4) {
  vapply(v, function(x) {
    if (!is.finite(x)) return("NA")
    if (digits >= 100) return(sprintf(paste0("%.", digits - 100, "f"), x))
    if (x == 0) return("0")
    d <- min(max(digits, 1), 15); a <- abs(x)
    if (a < 1e-4 || a >= 1e7) {
      s <- sprintf(paste0("%.", d - 1, "e"), x)
      return(sub("e([+-])0*([0-9])", "e\\1\\2", sub("e\\+", "e", s)))
    }
    mag <- floor(log10(a)); dec <- max(d - 1 - mag, 0)
    sprintf(paste0("%.", dec, "f"), x)
  }, "")
}

# the text of one statistic of a vector of values
hp_stat_text <- function(x, tok, digits = 4, total = length(x)) {
  t <- tolower(trimws(tok)); v <- x[is.finite(x)]; n <- length(v)
  f <- function(z) hp_sig(z, digits)
  q <- function(p) if (n) hp_quantile7(v, p) else NA_real_
  m <- if (n) mean(v) else NA_real_; s <- if (n > 1L) stats::sd(v) else NA_real_
  se <- s / sqrt(n)
  tq <- if (n > 1L) stats::qt(0.975, n - 1L) else NA_real_
  gm <- if (n && all(v > 0)) exp(mean(log(v))) else NA_real_
  gcv <- if (n > 1L && all(v > 0)) sqrt(exp(stats::sd(log(v))^2) - 1) * 100 else NA_real_
  switch(t,
    n = as.character(n), nmiss = , nmissing = as.character(total - n),
    mean_sd = sprintf("%s (%s)", f(m), f(s)), mean_se = sprintf("%s (%s)", f(m), f(se)),
    mean_ci = sprintf("%s (%s, %s)", f(m), f(m - tq * se), f(m + tq * se)),
    median_range = sprintf("%s [%s, %s]", f(q(0.5)), f(q(0)), f(q(1))),
    median_iqr = sprintf("%s [%s, %s]", f(q(0.5)), f(q(0.25)), f(q(0.75))),
    geomean_cv = sprintf("%s (%s%%)", f(gm), f(gcv)), range = sprintf("%s \u2013 %s", f(q(0)), f(q(1))),
    mean = f(m), sd = f(s), se = , sem = f(se), cv = f(if (is.finite(m) && m != 0) s / abs(m) * 100 else NA_real_),
    geomean = f(gm), geocv = , geosd = f(gcv), median = f(q(0.5)), min = f(q(0)), max = f(q(1)),
    sum = f(if (n) sum(v) else NA_real_), q1 = f(q(0.25)), q3 = f(q(0.75)), iqr = f(q(0.75) - q(0.25)),
    lcl95 = f(m - tq * se), ucl95 = f(m + tq * se),
    if (grepl("^p[0-9.]+$", t)) f(q(as.numeric(sub("^p", "", t)) / 100)) else hp_stop("statistic `", tok, "` is not known"))
}

hp_group_index <- function(data, by) {
  if (!length(by)) return(list(label = "All", idx = rep(1L, nrow(data)), tuples = list(character()), names = character()))
  miss <- setdiff(by, names(data)); if (length(miss)) hp_stop("table `by`: `", miss[1L], "` is not a column")
  txt <- lapply(by, function(b) { v <- data[[b]]; ifelse(is.na(v), "NA", if (is.numeric(v)) hp_num_text(v) else as.character(v)) })
  key <- do.call(paste, c(txt, sep = "\r"))
  rk <- lapply(txt, function(tx) { u <- unique(tx); u1 <- u[u != "NA"]; n <- suppressWarnings(as.numeric(u1))
    o <- if (all(!is.na(n))) u1[order(n)] else u1[order(is.na(n), n, u1)]
    match(tx, c(o, if ("NA" %in% u) "NA")) })
  ord <- do.call(order, rk)
  lev <- unique(key[ord])
  tuples <- lapply(lev, function(k) strsplit(k, "\r", fixed = TRUE)[[1L]])
  list(idx = match(key, lev), tuples = tuples, names = by,
       label = vapply(tuples, function(tp) paste(tp, collapse = ", "), ""))
}

#' Build a custom table from a list of settings
#'
#' Builds a table of any data frame from a list of plain settings. Six kinds
#' cover the tables of a clinical pharmacology report: a listing of rows, a
#' descriptive summary, a demographics table with the study groups side by
#' side, frequency tables (cross-tabulations), a pivot table of a statistic
#' and a correlation matrix.
#'
#' @details
#' **Kinds (`kind`).**
#'
#' - `"listing"`: the rows as they are, numbers formatted.
#' - `"summary"`: statistics (`stats`) of numeric `vars`, optionally by the
#'   columns in `by`, with `total = TRUE` adding an "All" group.
#' - `"demographics"`: the table of baseline characteristics. `vars`,
#'   `continuous` and `categorical` name the characteristics (a numeric column
#'   with more than `max_levels` values is continuous), `by` the groups,
#'   `stats` the statistics of continuous variables (default
#'   `c("mean_sd", "median_range")`), `total` an overall column. Categorical
#'   variables show counts and percentages (`percent_digits`).
#' - `"frequency"`: counts of the combinations of `rows` (or `vars`), across
#'   the levels of `columns`, with `percent` of `"none"`, `"row"`, `"column"`
#'   or `"total"`.
#' - `"pivot"`: a statistic `fun` of the numeric column `value`, by `rows` and
#'   `columns`.
#' - `"correlation"`: Pearson (or Spearman with `spearman = TRUE`) correlations
#'   of `vars`.
#'
#' **Statistics** (`stats`, `fun`): `n`, `nmiss`, `mean`, `sd`, `cv`, `sem`,
#' `median`, `min`, `max`, `q1`, `q3`, `p5`, `p95`, `geomean`, `gcv`,
#' `mean_sd`, `median_range` and `sum`.
#'
#' **Finishing.** `digits` (significant digits) or `decimals` (fixed number of
#' decimals), `rename` (a list of `c(old, new)` pairs), `hide` (columns to
#' drop), `na_text` (text for empty cells) and `transpose` (swap rows and
#' columns).
#'
#' @param data A data frame, usually the result of [hp_custom_prep()].
#' @param def A named list of settings, see Details.
#' @param digits Significant digits used when `def` has none.
#' @return A data frame of text, ready to print, save with [hp_save_table()]
#'   or pass to a report.
#' @seealso [hp_custom_prep()], [hp_custom_plot()], [hp_save_table()],
#'   [hp_describe()]
#' @family custom plots and tables
#' @examples
#' d <- hp_read_data(hp_example_file())
#' subj <- hp_custom_prep(d, distinct = "ID")
#'
#' # baseline characteristics by dose group
#' hp_custom_table(subj, list(kind = "demographics", vars = c("WT", "AGE", "SEX"),
#'                            by = "DOSEGRP", total = TRUE))
#'
#' # statistics of a variable, by group
#' hp_custom_table(subj, list(kind = "summary", vars = "WT", by = "SEX",
#'                            stats = c("n", "mean", "sd", "min", "max")))
#'
#' # subjects by sex and dose group, with row percentages
#' hp_custom_table(subj, list(kind = "frequency", rows = "SEX", columns = "DOSEGRP",
#'                            percent = "row"))
#'
#' # mean weight by sex (rows) and dose group (columns)
#' hp_custom_table(subj, list(kind = "pivot", rows = "SEX", columns = "DOSEGRP",
#'                            value = "WT", fun = "mean"))
#'
#' # correlations of covariates
#' hp_custom_table(subj, list(kind = "correlation", vars = c("WT", "AGE")))
#' @export
hp_custom_table <- function(data, def, digits = 3) {
  g <- function(n, default = NULL) { v <- def[[n]]; if (is.null(v) || (length(v) == 1L && is.na(v))) default else v }
  kind <- tolower(g("kind", "listing")); dg <- g("digits", digits); by <- g("by", character())
  if (!is.null(g("decimals"))) dg <- 100 + min(g("decimals"), 12)
  out <- switch(kind,
    listing = {
      o <- as.data.frame(lapply(data, function(v) if (is.numeric(v)) ifelse(is.finite(v), hp_cell(v, dg), "") else ifelse(is.na(v), "", as.character(v))),
                         stringsAsFactors = FALSE, check.names = FALSE)
      o
    },
    summary = {
      vars <- g("vars", character()); if (!length(vars)) vars <- setdiff(names(data)[vapply(data, is.numeric, TRUE)], by)
      if (!length(vars)) hp_stop("the summary table has no numeric column to describe")
      stats <- g("stats", c("n", "mean", "sd", "median", "min", "max"))
      gi <- hp_group_index(data, by)
      rows <- list()
      for (v in vars) {
        if (!v %in% names(data)) hp_stop("table `vars`: `", v, "` is not a column")
        x <- data[[v]]; if (!is.numeric(x)) hp_stop("table `vars`: `", v, "` is not numeric")
        sets <- if (length(by)) lapply(seq_along(gi$tuples), function(i) list(gi$tuples[[i]], which(gi$idx == i))) else list(list(character(), seq_along(x)))
        if (length(by) && isTRUE(g("total"))) sets[[length(sets) + 1L]] <- list(rep("All", length(by)), seq_along(x))
        for (st in sets) {
          if (!length(st[[2L]])) next
          rows[[length(rows) + 1L]] <- c(v, st[[1L]], vapply(stats, function(s) hp_stat_text(x[st[[2L]]], s, dg, length(st[[2L]])), ""))
        }
      }
      o <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
      names(o) <- c("Variable", by, vapply(stats, hp_stat_label, ""))
      o
    },
    demographics = hp_custom_demographics(data, def, dg),
    frequency = hp_custom_frequency(data, def),
    pivot = hp_custom_pivot(data, def, dg),
    correlation = {
      vars <- g("vars", character()); if (!length(vars)) vars <- names(data)[vapply(data, is.numeric, TRUE)]
      if (length(vars) < 2L) hp_stop("a correlation table needs at least two numeric columns")
      m <- stats::cor(data[vars], use = "pairwise.complete.obs", method = if (isTRUE(g("spearman"))) "spearman" else "pearson")
      o <- data.frame(vars, stringsAsFactors = FALSE); names(o) <- ""
      for (j in seq_along(vars)) o[[vars[j]]] <- ifelse(is.finite(m[, j]), sprintf("%.3f", m[, j]), "")
      o
    },
    hp_stop("table kind `", kind, "` is not known (listing | summary | demographics | frequency | pivot | correlation)"))
  ren <- g("rename", list())
  for (r in ren) names(out)[names(out) == r[[1L]]] <- r[[2L]]
  hide <- g("hide", character()); if (length(hide) && length(hide) < ncol(out)) out <- out[setdiff(names(out), hide)]
  nat <- g("na_text"); if (!is.null(nat)) out[] <- lapply(out, function(z) ifelse(z == "", nat, z))
  if (isTRUE(g("transpose")) && ncol(out) > 1L) {
    h <- out[[1L]]; t <- as.data.frame(t(out[-1L]), stringsAsFactors = FALSE)
    out <- cbind(data.frame(names(out)[-1L], stringsAsFactors = FALSE), t); names(out) <- c(names(out)[1L], h)
    names(out)[1L] <- "Variable"; rownames(out) <- NULL
  }
  rownames(out) <- NULL
  out
}

hp_custom_demographics <- function(data, def, dg) {
  g <- function(n, default = NULL) { v <- def[[n]]; if (is.null(v) || (length(v) == 1L && is.na(v))) default else v }
  by <- g("by", character()); gi <- hp_group_index(data, by)
  max_lev <- g("max_levels", 6L)
  vars <- g("vars", character()); cont <- g("continuous", character()); cat_ <- g("categorical", character())
  order_ <- list()
  for (v in vars) {
    if (!v %in% names(data)) hp_stop("table `vars`: `", v, "` is not a column")
    is_c <- is.numeric(data[[v]]) && length(unique(data[[v]][!is.na(data[[v]])])) > max_lev
    order_[[length(order_) + 1L]] <- list(v, (is_c || v %in% cont) && !v %in% cat_)
  }
  for (v in cont) if (!v %in% vars) order_[[length(order_) + 1L]] <- list(v, TRUE)
  for (v in cat_) if (!v %in% vars) order_[[length(order_) + 1L]] <- list(v, FALSE)
  if (!length(order_)) hp_stop("a demographics table needs variables (vars, continuous or categorical)")
  stats <- g("stats", c("mean_sd", "median_range"))
  sets <- list(); overall <- g("total", TRUE)
  if (isTRUE(overall) || !length(by)) sets[[length(sets) + 1L]] <- list("Overall", seq_len(nrow(data)))
  if (length(by)) for (i in seq_along(gi$tuples)) sets[[length(sets) + 1L]] <- list(gi$label[i], which(gi$idx == i))
  head <- c("Characteristic", vapply(sets, function(s) sprintf("%s (N=%d)", s[[1L]], length(s[[2L]])), ""))
  pd <- g("percent_digits", 1L)
  rows <- list()
  for (o in order_) {
    v <- o[[1L]]; x <- data[[v]]
    rows[[length(rows) + 1L]] <- c(v, rep("", length(sets)))
    if (o[[2L]]) {
      if (!is.numeric(x)) hp_stop("the continuous variable `", v, "` is not numeric")
      for (s in stats) rows[[length(rows) + 1L]] <- c(hp_stat_label(s), vapply(sets, function(st) {
        xs <- x[st[[2L]]]; if (!any(is.finite(xs))) "NA" else hp_stat_text(xs, s, dg, length(xs)) }, ""))
      miss <- vapply(sets, function(st) sum(!is.finite(x[st[[2L]]])), 1L)
      if (any(miss > 0L)) rows[[length(rows) + 1L]] <- c("Missing", as.character(miss))
    } else {
      txt <- ifelse(is.na(x), "", if (is.numeric(x)) hp_num_text(x) else as.character(x))
      levs <- hp_sort_labels(unique(txt[txt != ""]))
      present <- vapply(sets, function(st) sum(txt[st[[2L]]] != ""), 1L)
      for (lv in levs) rows[[length(rows) + 1L]] <- c(lv, vapply(seq_along(sets), function(k) {
        n <- sum(txt[sets[[k]][[2L]]] == lv); sprintf("%d (%s%%)", n, if (present[k]) sprintf(paste0("%.", pd, "f"), 100 * n / present[k]) else "0") }, ""))
      miss <- vapply(seq_along(sets), function(k) length(sets[[k]][[2L]]) - present[k], 1L)
      if (any(miss > 0L)) rows[[length(rows) + 1L]] <- c("Missing", as.character(miss))
    }
  }
  o <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE); names(o) <- head
  o
}

hp_custom_frequency <- function(data, def) {
  g <- function(n, default = NULL) { v <- def[[n]]; if (is.null(v) || (length(v) == 1L && is.na(v))) default else v }
  rows_ <- g("rows", character()); if (!length(rows_)) rows_ <- g("vars", character()); cols_ <- g("columns", character())
  if (!length(rows_)) hp_stop("a frequency table needs `rows` (or vars)")
  chk <- setdiff(c(rows_, cols_), names(data)); if (length(chk)) hp_stop("table: `", chk[1L], "` is not a column")
  pc <- tolower(as.character(g("percent", "none"))); pd <- g("percent_digits", 1L); totals <- g("total", TRUE)
  rg <- hp_group_index(data, rows_); cg <- hp_group_index(data, cols_)
  nr <- length(rg$tuples); nc <- length(cg$tuples)
  m <- matrix(0L, nr, nc); for (r in seq_len(nrow(data))) m[rg$idx[r], cg$idx[r]] <- m[rg$idx[r], cg$idx[r]] + 1L
  rt <- rowSums(m); ct <- colSums(m); grand <- sum(m)
  pct <- function(n, den) if (den == 0) "NA" else sprintf(paste0("%.", pd, "f"), 100 * n / den)
  tup <- function(i) rg$tuples[[i]]
  if (!length(cols_)) {
    out <- as.data.frame(do.call(rbind, lapply(seq_len(nr), function(i) c(tup(i), as.character(rt[i]), pct(rt[i], grand)))), stringsAsFactors = FALSE)
    names(out) <- c(rows_, "N", "%")
    if (totals) out[nrow(out) + 1L, ] <- c("Total", rep("", length(rows_) - 1L), as.character(grand), sprintf(paste0("%.", pd, "f"), 100))
    return(out)
  }
  cell <- function(n, den) if (pc == "none") as.character(n) else sprintf("%d (%s%%)", n, pct(n, den))
  rows <- lapply(seq_len(nr), function(i) {
    cells <- vapply(seq_len(nc), function(j) cell(m[i, j], switch(pc, row = rt[i], column = ct[j], grand)), "")
    c(tup(i), cells, if (totals) cell(rt[i], if (pc == "total") grand else rt[i]))
  })
  if (totals) rows[[nr + 1L]] <- c("Total", rep("", length(rows_) - 1L),
                                   vapply(seq_len(nc), function(j) cell(ct[j], if (pc == "column") ct[j] else grand), ""), cell(grand, grand))
  out <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(out) <- c(rows_, paste0(paste(cols_, collapse = "/"), "=", cg$label), if (totals) "Total")
  out
}

hp_custom_pivot <- function(data, def, dg) {
  g <- function(n, default = NULL) { v <- def[[n]]; if (is.null(v) || (length(v) == 1L && is.na(v))) default else v }
  rows_ <- g("rows", character()); cols_ <- g("columns", character()); val <- g("value"); fun <- g("fun", "mean")
  if (!length(rows_)) hp_stop("a pivot table needs `rows`")
  if (is.null(val) && fun != "n") hp_stop("a pivot table needs `value` (the column to summarise)")
  chk <- setdiff(c(rows_, cols_, val), names(data)); if (length(chk)) hp_stop("table: `", chk[1L], "` is not a column")
  rg <- hp_group_index(data, rows_); cg <- hp_group_index(data, cols_)
  x <- if (is.null(val)) rep(0, nrow(data)) else data[[val]]
  if (!is.null(val) && !is.numeric(x)) hp_stop("the value column `", val, "` is not numeric")
  totals <- g("total", FALSE); two <- length(cols_) > 0L
  cell <- function(rows) {
    if (!length(rows)) return("")
    z <- x[rows]; if (!is.null(val) && !any(is.finite(z))) return("NA")
    hp_stat_text(z, fun, dg, length(rows))
  }
  nr <- length(rg$tuples); nc <- length(cg$tuples)
  rows <- lapply(seq_len(nr), function(i) {
    ri <- which(rg$idx == i)
    if (two) c(rg$tuples[[i]], vapply(seq_len(nc), function(j) cell(ri[cg$idx[ri] == j]), ""), if (totals) cell(ri))
    else c(rg$tuples[[i]], cell(ri))
  })
  if (totals && two) rows[[nr + 1L]] <- c("All", rep("", length(rows_) - 1L), vapply(seq_len(nc), function(j) cell(which(cg$idx == j)), ""), cell(seq_len(nrow(data))))
  out <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(out) <- c(rows_, if (two) c(cg$label, if (totals) "All") else hp_stat_label(fun))
  out
}
