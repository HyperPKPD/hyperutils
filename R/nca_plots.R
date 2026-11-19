# Figures of a non-compartmental analysis, drawn from the tables of hp_nca().

#' Options of the NCA figures
#'
#' Chooses which figures [hp_nca_plots()] draws and how: profiles of all
#' individuals, mean profiles by group with error bars, box plots of selected
#' parameters and pages of individual profiles with the terminal regression
#' highlighted.
#'
#' @param profiles Draw the concentration-time profiles of all individuals.
#' @param mean_profile Draw the mean profile by group.
#' @param semilog,linear Draw the semilogarithmic and/or linear versions of
#'   the profile and mean-profile figures.
#' @param error_bars Error bars of the mean profile: `"sd"`, `"sem"`,
#'   `"ci95"` or `"off"`.
#' @param geometric_mean Use the geometric mean and geometric bars.
#' @param time_bin Width of the time bins for the mean profile; `0` groups
#'   identical times only. Use it when sampling times differ slightly between
#'   subjects.
#' @param individual Draw pages of individual profiles.
#' @param max_individual Most profiles in the individual pages (`0` = all).
#' @param individual_nrow,individual_ncol Panels per page.
#' @param individual_yscale `"linear"` or `"log"` axis of the individual
#'   panels.
#' @param terminal_points Ring the points used for the terminal regression
#'   and draw its fitted line.
#' @param individual_lines Join the observations of the individual panels by
#'   lines.
#' @param boxplots Names of parameters to draw as box plots by group, for
#'   example `c("CMAX", "AUCINF_OBS")`.
#' @param units A list with the elements `time`, `conc` and `dose`, used for
#'   the axis titles (the default is `list(time = "h", conc = "", dose =
#'   "mg")`).
#' @return A list of class `hp_nca_figures`.
#' @seealso [hp_nca_plots()]
#' @family non-compartmental analysis
#' @examples
#' fo <- hp_nca_figure_options(boxplots = c("CMAX", "AUCLAST"), error_bars = "ci95",
#'                             geometric_mean = TRUE, individual = FALSE,
#'                             units = list(time = "h", conc = "ng/mL", dose = "mg"))
#' fo$error_bars
#' @export
hp_nca_figure_options <- function(profiles = TRUE, mean_profile = TRUE, semilog = TRUE, linear = TRUE,
                                  error_bars = "sd", geometric_mean = FALSE, time_bin = 0, individual = TRUE,
                                  max_individual = 0, individual_nrow = 3, individual_ncol = 3,
                                  individual_yscale = "log", terminal_points = TRUE, individual_lines = TRUE,
                                  boxplots = character(), units = list(time = "h", conc = "", dose = "mg")) {
  o <- list(profiles = profiles, mean_profile = mean_profile, semilog = semilog, linear = linear,
            error_bars = error_bars, geometric_mean = geometric_mean, time_bin = time_bin, individual = individual,
            max_individual = max_individual, individual_nrow = individual_nrow, individual_ncol = individual_ncol,
            individual_yscale = individual_yscale, terminal_points = terminal_points,
            individual_lines = individual_lines, boxplots = toupper(boxplots), units = units)
  class(o) <- c("hp_nca_figures", "list")
  o
}

hp_nca_samples <- function(res) {
  p <- res$profiles
  p <- p[is.finite(p$CONC_USED) & !grepl("^added", p$STATUS), , drop = FALSE]
  p$PROFILE <- paste(p$ID, p$OCCASION, sep = "\u00b7")
  p
}

#' Draw the standard figures of a non-compartmental analysis
#'
#' Creates the figures that usually accompany an NCA report from the result of
#' [hp_nca()]: all individual profiles, mean profiles by group, box plots of
#' parameters by group, and pages of individual profiles that show the points
#' used for the terminal slope and the fitted line.
#'
#' @details
#' BLQ samples are shown at the value they were assigned by the BLQ rule;
#' zero concentrations cannot be shown on a logarithmic axis and are left out
#' of the semilogarithmic versions. Groups come from the `groups` argument of
#' [hp_nca()]. Figures that cannot be drawn (for example when no profile could be
#' analysed) are omitted from the result.
#'
#' @param res The result of [hp_nca()].
#' @param options Figure options from [hp_nca_figure_options()].
#' @return A named list of ggplot objects: `profiles_linear`,
#'   `profiles_semilog`, `mean_profile_linear`, `mean_profile_semilog`,
#'   `box_<parameter>` (for example `box_cmax`) and `individual_01`,
#'   `individual_02`, and so on (one per page). Print an element to see it, or
#'   write them with [hp_save_plot()].
#' @seealso [hp_nca_figure_options()], [hp_save_plot()], [hp_nca()]
#' @family non-compartmental analysis
#' @examples
#' d <- hp_read_data(hp_example_file())
#' r <- hp_nca(d, route = "oral", groups = "DOSEGRP", lloq = 0.5)
#' f <- hp_nca_plots(r, hp_nca_figure_options(
#'   boxplots = "CMAX", individual_nrow = 2, individual_ncol = 3,
#'   units = list(time = "h", conc = "ng/mL", dose = "mg")))
#' names(f)
#' f$mean_profile_semilog
#'
#' # write every figure to a folder
#' out <- file.path(tempdir(), "nca_figs")
#' for (nm in names(f)) hp_save_plot(f[[nm]], nm, out, width = 7, height = 5, dpi = 72)
#' length(list.files(out))
#' unlink(out, recursive = TRUE)
#' @export
hp_nca_plots <- function(res, options = hp_nca_figure_options()) {
  o <- options; u <- o$units
  th <- hp_theme(11)
  today <- format(Sys.Date())
  cap <- function(n) paste0("hyperutils nca \u00b7 ", n, " profile(s) \u00b7 ", today)
  xl <- paste0("time after dose (", u$time, ")")
  yl <- paste0("concentration (", u$conc, ")")
  s <- hp_nca_samples(res)
  figs <- list()
  ok_prof <- length(unique(res$parameters$ID[!is.na(res$parameters$CMAX)])) > 0L
  if (!nrow(s) || !ok_prof) return(figs)
  nprof <- length(unique(s$PROFILE))
  labs <- sort(unique(as.character(s$GROUP)))
  pal <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#8C6D1F")
  cols <- stats::setNames(rep_len(c(pal, rep("#9A9A9A", max(0, length(labs) - 7L))), length(labs)), labs)
  one <- function(logy, mean_) {
    z <- s
    if (logy) z <- z[z$CONC_USED > 0, , drop = FALSE]
    z$GROUP <- factor(as.character(z$GROUP), levels = labs)
    z
  }
  if (o$profiles) for (lg in c(FALSE, TRUE)) {
    if ((lg && !o$semilog) || (!lg && !o$linear)) next
    z <- one(lg)
    p <- ggplot2::ggplot(z, ggplot2::aes(TIME, CONC_USED, group = PROFILE, colour = GROUP)) +
      ggplot2::geom_line(alpha = 0.6, linewidth = 0.45) + ggplot2::geom_point(alpha = 0.7, size = 1) +
      ggplot2::scale_colour_manual(values = cols, drop = FALSE) +
      ggplot2::labs(x = xl, y = yl, colour = if (length(labs) > 1L) "group" else NULL,
                    title = if (lg) "Concentration\u2013time profiles (semi-logarithmic)" else "Concentration\u2013time profiles",
                    caption = cap(nprof)) + th
    if (length(labs) < 2L) p <- p + ggplot2::theme(legend.position = "none")
    if (lg) p <- p + ggplot2::scale_y_log10()
    figs[[if (lg) "profiles_semilog" else "profiles_linear"]] <- p
  }
  if (o$mean_profile) for (lg in c(FALSE, TRUE)) {
    if ((lg && !o$semilog) || (!lg && !o$linear)) next
    z <- one(lg)
    key <- if (o$time_bin > 0) round(z$TIME / o$time_bin) * o$time_bin else round(z$TIME * 1e6) / 1e6
    rows <- list()
    for (g in labs) for (tk in sort(unique(key[z$GROUP == g]))) {
      v <- z$CONC_USED[z$GROUP == g & key == tk]
      cs <- hp_nca_mean_point(v, o$geometric_mean, o$error_bars)
      rows[[length(rows) + 1L]] <- data.frame(GROUP = g, TIME = tk, C = cs[1L], L = cs[2L], U = cs[3L])
    }
    md <- do.call(rbind, rows); if (is.null(md)) next
    md$GROUP <- factor(md$GROUP, levels = labs)
    what <- if (o$geometric_mean) "Geometric mean" else "Mean"
    bt <- switch(o$error_bars, sd = if (o$geometric_mean) " \u00b1 geometric SD" else " \u00b1 SD", sem = " \u00b1 SEM",
                 ci95 = " with 95 % CI", "")
    p <- ggplot2::ggplot(md, ggplot2::aes(TIME, C, colour = GROUP, group = GROUP)) + ggplot2::geom_line() +
      ggplot2::geom_point(size = 1.6) +
      ggplot2::scale_colour_manual(values = cols, drop = FALSE) +
      ggplot2::labs(x = xl, y = yl, colour = if (length(labs) > 1L) "group" else NULL,
                    title = paste0(what, " concentration", bt), caption = cap(nprof)) + th
    if (o$error_bars != "off") p <- p + ggplot2::geom_errorbar(ggplot2::aes(ymin = L, ymax = U), width = 0, na.rm = TRUE)
    if (length(labs) < 2L) p <- p + ggplot2::theme(legend.position = "none")
    if (lg) p <- p + ggplot2::scale_y_log10()
    figs[[if (lg) "mean_profile_semilog" else "mean_profile_linear"]] <- p
  }
  for (pn in o$boxplots) {
    if (!pn %in% names(res$parameters)) next
    z <- res$parameters; z$.Y <- z[[pn]]; z <- z[is.finite(z$.Y), , drop = FALSE]
    if (!nrow(z)) next
    z$GROUP <- factor(as.character(z$GROUP), levels = labs)
    p <- ggplot2::ggplot(z, ggplot2::aes(GROUP, .Y, fill = GROUP)) + ggplot2::geom_boxplot(outlier.size = 0.8, alpha = 0.75) +
      ggplot2::geom_jitter(width = 0.12, height = 0, size = 0.9, alpha = 0.6) +
      ggplot2::scale_fill_manual(values = cols, drop = FALSE) +
      ggplot2::labs(x = NULL, y = pn, title = paste(pn, "by group"), caption = cap(nprof)) + th +
      ggplot2::theme(legend.position = "none")
    figs[[paste0("box_", tolower(pn))]] <- p
  }
  if (o$individual) {
    ids <- unique(s$PROFILE)
    if (o$max_individual > 0) ids <- utils::head(ids, o$max_individual)
    lz <- res$lambda_z; lz$PROFILE <- paste(lz$ID, lz$OCCASION, sep = "\u00b7")
    pr <- res$parameters; pr$PROFILE <- paste(pr$ID, pr$OCCASION, sep = "\u00b7")
    per <- max(o$individual_nrow * o$individual_ncol, 1L)
    pages <- split(ids, ceiling(seq_along(ids) / per))
    for (k in seq_along(pages)) {
      z <- s[s$PROFILE %in% pages[[k]], , drop = FALSE]
      z$.PNL <- factor(paste0("ID ", z$ID, if (length(unique(s$OCCASION)) > 1L) paste0(" \u00b7 occasion ", z$OCCASION) else ""),
                        levels = unique(paste0("ID ", z$ID, if (length(unique(s$OCCASION)) > 1L) paste0(" \u00b7 occasion ", z$OCCASION) else "")))
      lg <- o$individual_yscale == "log"
      if (lg) z <- z[z$CONC_USED > 0, , drop = FALSE]
      p <- ggplot2::ggplot(z, ggplot2::aes(TIME, CONC_USED))
      if (o$individual_lines) p <- p + ggplot2::geom_line(colour = "grey40", linewidth = 0.4)
      p <- p + ggplot2::geom_point(size = 1.1, colour = "#0072B2")
      if (o$terminal_points) {
        zt <- z[z$IN_LAMBDA_Z == 1, , drop = FALSE]
        if (nrow(zt)) p <- p + ggplot2::geom_point(data = zt, shape = 1, size = 2.8, colour = "#D55E00", stroke = 0.7)
      }
      fit <- lapply(pages[[k]], function(pid) {
        l <- lz[lz$PROFILE == pid & lz$CHOSEN == 1, , drop = FALSE]
        zz <- z[z$PROFILE == pid & z$IN_LAMBDA_Z == 1, , drop = FALSE]
        if (!nrow(l) || !nrow(zz) || !is.finite(l$LAMBDA_Z[1L])) return(NULL)
        tt <- seq(min(zz$TIME), max(max(zz$TIME), max(z$TIME[z$PROFILE == pid])), length.out = 40)
        data.frame(.PNL = zz$.PNL[1L], TIME = tt, CONC_USED = exp(l$INTERCEPT[1L] - l$LAMBDA_Z[1L] * tt))
      })
      fit <- do.call(rbind, fit)
      if (!is.null(fit)) p <- p + ggplot2::geom_line(data = fit, colour = "#D55E00", linewidth = 0.6)
      p <- p + ggplot2::facet_wrap(ggplot2::vars(.PNL), scales = "free", ncol = o$individual_ncol) +
        ggplot2::labs(x = xl, y = yl, title = "Individual profiles with the terminal fit", caption = cap(length(ids))) + th
      if (lg) p <- p + ggplot2::scale_y_log10()
      figs[[sprintf("individual_%02d", k)]] <- p
    }
  }
  figs
}

# centre and bars of one time point of a mean profile
hp_nca_mean_point <- function(v, geometric, bars) {
  v <- v[is.finite(v)]; n <- length(v)
  if (!n) return(c(NA_real_, NA_real_, NA_real_))
  if (geometric) {
    l <- log(v); m <- mean(l); sdl <- if (n > 1L) stats::sd(l) else 0
    c0 <- exp(m); gsd <- exp(sdl)
    return(switch(bars, sd = c(c0, c0 / gsd, c0 * gsd), sem = c(c0, exp(m - sdl / sqrt(n)), exp(m + sdl / sqrt(n))),
                  ci95 = { tq <- if (n > 1L) stats::qt(0.975, n - 1L) else 0; c(c0, exp(m - tq * sdl / sqrt(n)), exp(m + tq * sdl / sqrt(n))) },
                  c(c0, c0, c0)))
  }
  m <- mean(v); s <- if (n > 1L) stats::sd(v) else 0
  switch(bars, sd = c(m, m - s, m + s), sem = c(m, m - s / sqrt(n), m + s / sqrt(n)),
         ci95 = { tq <- if (n > 1L) stats::qt(0.975, n - 1L) else 0; c(m, m - tq * s / sqrt(n), m + tq * s / sqrt(n)) },
         c(m, m, m))
}
