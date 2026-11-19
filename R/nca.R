# Non-compartmental analysis (WinNonlin- and
# Monolix-style), written in R with the areas and the terminal regression in C++.
#
# Pipeline: hp_nca_read() -> hp_nca_profiles() -> hp_nca_analyse() per profile
# -> hp_nca_summary(); hp_nca() runs all of it.

# ---------------------------------------------------------------------------
# options
# ---------------------------------------------------------------------------

hp_norm <- function(x) gsub("[_ ]", "-", tolower(trimws(x)))

hp_choice <- function(x, table, what) {
  k <- hp_norm(x)
  for (nm in names(table)) if (k %in% table[[nm]]) return(nm)
  hp_stop(what, " `", x, "`: expected ", paste(names(table), collapse = " | "))
}

#' Settings of a non-compartmental analysis
#'
#' Collects every setting of an NCA into one object that you pass to
#' [hp_nca()]. The defaults follow common WinNonlin practice: linear-up/log-down
#' areas, the terminal slope chosen by best fit of the adjusted R-squared, and
#' BLQ concentrations before the maximum set to zero. Change only what your
#' analysis plan specifies and keep the object with the analysis, so the
#' settings are documented in the script.
#'
#' @details
#' **Area methods (`auc_method`).**
#'
#' - `"linear"`: linear trapezoids everywhere.
#' - `"linear-log"`: linear trapezoids up to the maximum and logarithmic
#'   trapezoids afterwards.
#' - `"linear-up-log-down"`: linear trapezoids while the concentration rises
#'   and logarithmic trapezoids while it falls (the default).
#'
#' A logarithmic trapezoid is used only when both concentrations are positive.
#'
#' **BLQ rules (`blq_rule`).** The rules apply when `lloq` is given or a BLQ
#' flag column is used (see [hp_nca()]). A sample is then below the limit of
#' quantification (BLQ) when it is flagged, below `lloq` or not positive.
#' Without a limit and without flags no sample is treated as BLQ. The rule
#' decides what number a BLQ sample stands for in the analysis:
#'
#' - `"zero"`: set to 0.
#' - `"missing"`: left out.
#' - `"lloq"`: set to the limit of quantification.
#' - `"half-lloq"`: set to half of it.
#' - `"tmax-zero-missing"`: BLQ before the maximum is set to 0, after it left
#'   out (the WinNonlin rule, the default).
#' - `"first-zero-mid-half-last-missing"`: BLQ before the first quantifiable
#'   sample is 0, between quantifiable samples it is half the limit, after the
#'   last quantifiable sample it is left out (the Monolix rule).
#'
#' Rules that need a limit (`"lloq"`, `"half-lloq"` and the Monolix rule) stop
#' with a message when neither `lloq` is given nor the BLQ sample carries its
#' limit as `DV`.
#'
#' **Terminal slope.** With `lz_method = "best-fit"` the regression of
#' `log(concentration)` on time is fitted to the last 3, 4, 5, ... points after
#' Tmax (from `lz_min_points` to `lz_max_points`, 0 meaning no limit). Among the
#' fits with a negative slope the one with the highest adjusted R-squared is
#' chosen; fits within `lz_tolerance` of the best are considered equal and the
#' one with more points wins. A fit with an adjusted R-squared below
#' `lz_min_adj_r2` is rejected and no lambda z is reported. `"last-n"` uses
#' exactly the last `lz_points` points and `"time-range"` the points between
#' `lz_from` and `lz_to`; `"none"` skips the terminal phase. A terminal phase
#' shorter than `lz_min_span_ratio` half-lives is flagged (or refused with
#' `lz_require_span = TRUE`). Subject-specific choices go in `lz_overrides`, a
#' list of lists with elements among `id`, `occasion`, `from`, `to`, `points`,
#' `exclude` (times to leave out) and `skip` (no terminal phase).
#'
#' **Extrapolated area.** When the area beyond the last sample exceeds
#' `lz_max_extrap` percent of the area to infinity the profile is flagged
#' (`lz_extrap_action = "flag"`) or its infinity-based parameters are left
#' empty (`"omit"`).
#'
#' **Steady state (`steady_state`).** `"auto"` uses the `SS` flag and `II` of
#' the data, `"none"` treats every interval as a single dose, `"last"` treats
#' the last interval of each subject as steady state and `"all"` every
#' interval. The dosing interval `tau` comes from `II` (with the `SS` flag),
#' from the `tau` argument or from the gap to the next dose when the doses are
#' evenly spaced. A steady-state interval reports the measures over the
#' dosing interval (`AUCTAU`, `CAVG` and so on) and no infinity-based
#' parameters.
#'
#' **Partial areas.** `partial = list(c(0, 4), c(0, 12))` adds the parameters
#' `AUC_0_4` and `AUC_0_12`. A partial area beyond the last sample is
#' extrapolated with lambda z when `partial_extrapolate = TRUE`.
#'
#' **Units.** `dose_unit`, `time_unit` and `conc_unit` are labels used in
#' summary tables; nothing is converted. `dose_factor` multiplies the dose in
#' clearance and volume, for example 1000 to turn mg and ng/mL into L/h and L.
#' The unit labels are not changed by it, so set them to match.
#'
#' @param auc_method `"linear"`, `"linear-log"` or `"linear-up-log-down"`.
#' @param blq_rule `"zero"`, `"missing"`, `"lloq"`, `"half-lloq"`,
#'   `"tmax-zero-missing"` or `"first-zero-mid-half-last-missing"`.
#' @param lloq Limit of quantification, or `NULL` when none is known.
#' @param blq_flagged `TRUE` when a BLQ flag column is used, so that flagged
#'   samples are treated as BLQ ([hp_nca()] sets it when `blq_flag` is given).
#' @param intervals Dosing intervals to analyse: `"all"`, `"first"` or `"last"`.
#' @param steady_state `"auto"`, `"none"`, `"last"` or `"all"`.
#' @param tau Dosing interval to use at steady state when the data have none.
#' @param min_observations Fewest usable samples a profile must have.
#' @param max_time Last time after the dose that is analysed, or `NULL`.
#' @param duplicates What to do with several samples at one time: `"mean"`,
#'   `"first"`, `"last"`, `"max"` or `"error"`.
#' @param c0 Concentration at the dose time of an intravenous bolus:
#'   `"extrapolate"` (the observed value when there is a positive one at the
#'   dose time, otherwise back-extrapolated from the first two quantifiable
#'   samples), `"observed"` (the observed value only; the area then starts at
#'   the first sample when there is none) or `"first"` (the first quantifiable
#'   sample stands for it).
#' @param zero_at_dose Add a concentration of zero at the time of the first
#'   dose of an extravascular profile that has no sample there (otherwise the
#'   area starts at the first sample).
#' @param partial A list of `c(from, to)` pairs for partial areas, or `NULL`.
#' @param partial_extrapolate Extrapolate a partial area that extends beyond
#'   the last sample with lambda z.
#' @param lz_method `"best-fit"`, `"last-n"`, `"time-range"` or `"none"`.
#' @param lz_min_points,lz_max_points Smallest and largest number of points
#'   tried by `"best-fit"` (`lz_max_points = 0` means no limit).
#' @param lz_points Number of points used by `"last-n"`.
#' @param lz_from,lz_to Time range used by `"time-range"`.
#' @param lz_include_tmax Allow the point at Tmax in the terminal regression.
#' @param lz_min_adj_r2 Smallest adjusted R-squared accepted for a terminal
#'   regression.
#' @param lz_tolerance Fits whose adjusted R-squared is within this distance
#'   of the best count as equal, and the one with most points is taken.
#' @param lz_weighting `"uniform"`, `"1/y"` or `"1/y2"` weights of the
#'   regression.
#' @param lz_min_span_ratio Smallest span of the terminal phase, in
#'   half-lives, before a flag is raised.
#' @param lz_require_span `TRUE` refuses a terminal phase shorter than
#'   `lz_min_span_ratio`.
#' @param lz_max_extrap Percent of the area to infinity beyond the last sample
#'   that triggers the extrapolation flag.
#' @param lz_extrap_action `"flag"` or `"omit"`.
#' @param lz_overrides A list of lists with subject-specific terminal-slope
#'   choices, see Details.
#' @param dose_unit,time_unit,conc_unit Unit labels for the output tables.
#' @param dose_factor Multiplier of the dose in clearance and volume.
#' @return A list of class `hp_nca_options`, to be passed as `options` to
#'   [hp_nca()].
#' @seealso [hp_nca()], [hp_nca_parameters()]
#' @family non-compartmental analysis
#' @examples
#' # the defaults
#' o <- hp_nca_options()
#' o$auc_method
#'
#' # an analysis plan: linear trapezoids, half the limit for BLQ, a partial area
#' o <- hp_nca_options(auc_method = "linear", blq_rule = "half-lloq", lloq = 0.5,
#'                     partial = list(c(0, 12)), lz_min_points = 4)
#' o$lz$min_points
#'
#' # fix the terminal phase of subject 7 to the times 12 to 48 hours
#' o <- hp_nca_options(lz_overrides = list(list(id = "7", from = 12, to = 48)))
#'
#' # units shown in the tables, and mg / (ng/mL) -> L conversion of CL and V
#' o <- hp_nca_options(dose_unit = "mg", conc_unit = "ng/mL", dose_factor = 1000)
#' hp_nca_unit("CL_F_OBS", o)
#' @export
hp_nca_options <- function(auc_method = "linear-up-log-down", blq_rule = "tmax-zero-missing",
                           lloq = NULL, blq_flagged = FALSE, intervals = "all",
                           steady_state = "auto", tau = NULL, min_observations = 3,
                           max_time = NULL, duplicates = "mean", c0 = "extrapolate",
                           zero_at_dose = TRUE, partial = NULL, partial_extrapolate = TRUE,
                           lz_method = "best-fit", lz_min_points = 3, lz_max_points = 0,
                           lz_points = 3, lz_from = NULL, lz_to = NULL, lz_include_tmax = FALSE,
                           lz_min_adj_r2 = 0.7, lz_tolerance = 1e-4, lz_weighting = "uniform",
                           lz_min_span_ratio = 2, lz_require_span = FALSE, lz_max_extrap = 20,
                           lz_extrap_action = "flag", lz_overrides = NULL,
                           dose_unit = "mg", time_unit = "h", conc_unit = "ng/mL",
                           dose_factor = 1) {
  o <- list(
    auc_method = hp_choice(auc_method, list(
      "linear" = c("linear", "lin", "linear-trapezoidal"),
      "linear-log" = c("linear-log", "lin-log", "linear-trapezoidal-linear-log", "log"),
      "linear-up-log-down" = c("linear-up-log-down", "lin-up-log-down", "linear-up/log-down",
                               "linuplogdown", "linear-up-logarithmic-down")), "AUC method"),
    blq_rule = hp_choice(blq_rule, list(
      "zero" = c("zero", "0", "set-to-zero"),
      "missing" = c("missing", "omit", "drop", "exclude", "na"),
      "lloq" = c("lloq", "loq", "set-to-lloq"),
      "half-lloq" = c("half-lloq", "lloq/2", "half-loq", "lloq-2"),
      "tmax-zero-missing" = c("tmax-zero-missing", "winnonlin", "phoenix", "before-tmax-zero"),
      "first-zero-mid-half-last-missing" = c("first-zero-mid-half-last-missing", "monolix", "mid-half")),
      "BLQ rule"),
    lloq = lloq, blq_flagged = isTRUE(blq_flagged),
    intervals = hp_choice(intervals, list(all = c("all", "every"), first = c("first", "single"),
                                          last = "last"), "intervals"),
    steady_state = hp_choice(steady_state, list(auto = c("auto", "flag"),
                                                none = c("none", "no", "off", "never", "false"),
                                                last = "last", all = c("all", "yes", "true")),
                             "steady state"),
    tau = tau, min_observations = min_observations, max_time = max_time,
    duplicates = hp_choice(duplicates, list(mean = c("mean", "average"), first = "first",
                                            last = "last", max = c("max", "maximum"),
                                            error = c("error", "stop", "fail")), "duplicates"),
    c0 = hp_choice(c0, list(extrapolate = c("extrapolate", "back-extrapolate", "log-linear", "auto"),
                            observed = c("observed", "data", "as-is"),
                            first = c("first", "first-observation")), "C0"),
    zero_at_dose = isTRUE(zero_at_dose), partial = partial,
    partial_extrapolate = isTRUE(partial_extrapolate),
    lz = list(
      method = hp_choice(lz_method, list("best-fit" = c("best-fit", "bestfit", "best", "auto"),
                                         "last-n" = c("last-n", "lastn", "last-points", "points"),
                                         "time-range" = c("time-range", "range", "times"),
                                         none = c("none", "off", "no")), "lambda-z method"),
      min_points = lz_min_points, max_points = lz_max_points, points = lz_points,
      from = lz_from, to = lz_to, include_tmax = isTRUE(lz_include_tmax),
      min_adj_r2 = lz_min_adj_r2, tolerance = lz_tolerance,
      weighting = hp_choice(lz_weighting, list(uniform = c("uniform", "none", "1", "equal"),
                                               "1/y" = c("1/y", "1/yhat", "inverse-y", "1/c"),
                                               "1/y2" = c("1/y2", "1/y^2", "1/yhat2", "inverse-y2",
                                                          "1/c2", "1/c^2")), "weighting"),
      min_span_ratio = lz_min_span_ratio, require_span = isTRUE(lz_require_span),
      max_extrap = lz_max_extrap,
      extrap_action = hp_choice(lz_extrap_action, list(flag = c("flag", "warn"),
                                                       omit = c("omit", "drop", "blank")),
                                "extrapolation action"),
      overrides = lz_overrides),
    units = list(dose = dose_unit, time = time_unit, conc = conc_unit, dose_factor = dose_factor))
  class(o) <- c("hp_nca_options", "list")
  o
}

hp_method_code <- function(m) c("linear" = 0L, "linear-log" = 1L, "linear-up-log-down" = 2L)[[m]]

# ---------------------------------------------------------------------------
# parameter registry
# ---------------------------------------------------------------------------

#' The list of NCA parameters, with labels and units
#'
#' Returns the registry of every parameter that [hp_nca()] can compute, in the
#' order of the output tables, with a plain-language description and the kind of
#' quantity (which decides its unit). Use it as the data dictionary of an
#' analysis report.
#'
#' @details
#' Parameter names follow the WinNonlin convention: `CMAX`, `TMAX`, `AUCLAST`,
#' `AUCINF_OBS`, `LAMZ`, `THALF`, `CL_OBS` (intravenous) or `CL_F_OBS`
#' (extravascular), `VZ_OBS`, `VSS_OBS` and so on, with `_PRED` for the version
#' that uses the concentration predicted by the terminal regression at the last
#' time, `_D` for dose-normalised values, and `_SS` for steady-state values.
#' Partial areas are named by their limits, for example `AUC_0_4` (a decimal
#' point is written `p` and a minus sign `m`), and are not in this registry;
#' their unit is that of an area.
#'
#' @return A data frame with the columns `NAME`, `KIND` and `LABEL`, one row
#'   per parameter.
#' @seealso [hp_nca_unit()], [hp_nca()]
#' @family non-compartmental analysis
#' @examples
#' p <- hp_nca_parameters()
#' nrow(p)
#' head(p, 8)
#'
#' # the steady-state parameters
#' p[grepl("_SS$|TAU", p$NAME), c("NAME", "LABEL")]
#' @export
hp_nca_parameters <- function() {
  if (!is.null(.hp_cache$params)) return(.hp_cache$params)
  r <- function(n, k, l) data.frame(NAME = n, KIND = k, LABEL = l, stringsAsFactors = FALSE)
  out <- do.call(rbind, list(
    r("DOSE", "dose", "dose amount"), r("N_SAMPLES", "count", "samples used"),
    r("INF_DURATION", "time", "infusion duration"), r("TAU", "time", "dosing interval"),
    r("TLAG", "time", "lag time before the first measurable concentration (extravascular)"),
    r("C0", "conc", "concentration at the dose time (IV bolus: observed or back-extrapolated)"),
    r("CMAX", "conc", "maximum observed concentration"), r("TMAX", "time", "time of the maximum"),
    r("CMAX_D", "conc_per_dose", "CMAX / dose"),
    r("CMIN", "conc", "minimum observed concentration in the dosing interval"),
    r("TMIN", "time", "time of the minimum"), r("CLAST", "conc", "last measurable concentration"),
    r("TLAST", "time", "time of the last measurable concentration"),
    r("CLAST_PRED", "conc", "concentration at TLAST predicted by the terminal regression"),
    r("LAMZ", "rate", "terminal rate constant lambda-z"), r("THALF", "time", "terminal half-life ln2/lambda-z"),
    r("LAMZ_N", "count", "points in the terminal regression"),
    r("LAMZ_START", "time", "first time of the terminal regression"),
    r("LAMZ_END", "time", "last time of the terminal regression"),
    r("LAMZ_SPAN", "ratio", "span of the terminal regression in half-lives"),
    r("R2", "ratio", "R-squared of the terminal regression"),
    r("R2ADJ", "ratio", "adjusted R-squared of the terminal regression"),
    r("AUCLAST", "auc", "area under the curve from the dose to TLAST"),
    r("AUCLAST_D", "auc_per_dose", "AUCLAST / dose"),
    r("AUCINF_OBS", "auc", "area to infinity, extrapolated from the observed CLAST"),
    r("AUCINF_PRED", "auc", "area to infinity, extrapolated from the predicted CLAST"),
    r("AUCINF_D", "auc_per_dose", "AUCINF_OBS / dose"),
    r("AUC_PCT_EXTRAP_OBS", "pct", "percent of AUCINF_OBS beyond TLAST"),
    r("AUC_PCT_EXTRAP_PRED", "pct", "percent of AUCINF_PRED beyond TLAST"),
    r("AUMCLAST", "aumc", "area under the first-moment curve to TLAST"),
    r("AUMCINF_OBS", "aumc", "first-moment area to infinity (observed CLAST)"),
    r("AUMCINF_PRED", "aumc", "first-moment area to infinity (predicted CLAST)"),
    r("MRTLAST", "time", "mean residence time to TLAST"),
    r("MRTINF_OBS", "time", "mean residence time to infinity (observed CLAST)"),
    r("MRTINF_PRED", "time", "mean residence time to infinity (predicted CLAST)"),
    r("CL_OBS", "clearance", "clearance, dose / AUCINF_OBS (IV)"),
    r("CL_PRED", "clearance", "clearance, dose / AUCINF_PRED (IV)"),
    r("VZ_OBS", "volume", "volume of distribution in the terminal phase, dose / (lambda-z * AUCINF_OBS) (IV)"),
    r("VZ_PRED", "volume", "volume of distribution in the terminal phase (predicted CLAST) (IV)"),
    r("VSS_OBS", "volume", "volume of distribution at steady state, CL * MRT (IV)"),
    r("VSS_PRED", "volume", "volume of distribution at steady state (predicted CLAST) (IV)"),
    r("CL_F_OBS", "clearance", "apparent clearance CL/F, dose / AUCINF_OBS (extravascular)"),
    r("CL_F_PRED", "clearance", "apparent clearance CL/F (predicted CLAST) (extravascular)"),
    r("VZ_F_OBS", "volume", "apparent volume Vz/F (extravascular)"),
    r("VZ_F_PRED", "volume", "apparent volume Vz/F (predicted CLAST) (extravascular)"),
    r("AUCTAU", "auc", "area over the dosing interval"),
    r("AUMCTAU", "aumc", "first-moment area over the dosing interval"),
    r("CTAU", "conc", "concentration at the end of the dosing interval"),
    r("CAVG", "conc", "average concentration AUCTAU / tau"),
    r("FLUCT_PCT", "pct", "fluctuation 100 * (CMAX - CMIN) / CAVG"),
    r("SWING", "ratio", "swing (CMAX - CMIN) / CMIN"), r("PTR", "ratio", "peak-to-trough ratio CMAX / CMIN"),
    r("CL_SS", "clearance", "clearance at steady state, dose / AUCTAU (IV)"),
    r("CL_F_SS", "clearance", "apparent clearance at steady state (extravascular)"),
    r("VZ_SS", "volume", "terminal volume at steady state, dose / (lambda-z * AUCTAU) (IV)"),
    r("VZ_F_SS", "volume", "apparent terminal volume at steady state (extravascular)"),
    r("MRT_SS", "time", "mean residence time at steady state"),
    r("VSS_SS", "volume", "volume of distribution at steady state, CL_SS * MRT_SS (IV)"),
    r("RACC_PRED", "ratio", "predicted accumulation 1 / (1 - exp(-lambda-z * tau))"),
    r("RACC_AUC", "ratio", "accumulation AUCTAU(steady state) / AUCTAU or AUCINF of the first dose"),
    r("RACC_CMAX", "ratio", "accumulation CMAX(steady state) / CMAX of the first dose")))
  .hp_cache$params <- out
  out
}

#' Unit of an NCA parameter
#'
#' Builds the unit text of a parameter from the unit labels in the options: a
#' time for `TMAX`, a concentration for `CMAX`, concentration times time for
#' areas, dose per (concentration times time) for clearance, and so on.
#'
#' @param name Parameter name, for example `"AUCLAST"` or a partial area such
#'   as `"AUC_0_4"`.
#' @param options Settings from [hp_nca_options()], which carry the unit
#'   labels.
#' @return A character string, empty for dimensionless parameters.
#' @seealso [hp_nca_parameters()]
#' @family non-compartmental analysis
#' @examples
#' o <- hp_nca_options(dose_unit = "mg", time_unit = "h", conc_unit = "ng/mL")
#' hp_nca_unit("CMAX", o)
#' hp_nca_unit("AUCINF_OBS", o)
#' hp_nca_unit("CL_F_OBS", o)
#' hp_nca_unit("LAMZ", o)
#' hp_nca_unit("AUC_0_4", o)
#' @export
hp_nca_unit <- function(name, options = hp_nca_options()) {
  reg <- hp_nca_parameters()
  k <- reg$KIND[match(name, reg$NAME)]
  if (is.na(k)) k <- if (startsWith(name, "AUC_")) "auc" else "ratio"
  u <- options$units; d <- u$dose; t <- u$time; c <- u$conc
  switch(k, dose = d, time = t, conc = c, auc = paste0(c, "*", t), aumc = paste0(c, "*", t, "^2"),
         rate = paste0("1/", t), pct = "%", ratio = "", count = "",
         clearance = paste0(d, "/(", c, "*", t, ")"), volume = paste0(d, "/", c),
         conc_per_dose = paste0(c, "/", d), auc_per_dose = paste0(c, "*", t, "/", d), "")
}

hp_partial_name <- function(a, b) {
  f <- function(x) gsub("-", "m", gsub(".", "p", as.character(x), fixed = TRUE), fixed = TRUE)
  paste0("AUC_", f(a), "_", f(b))
}

# ---------------------------------------------------------------------------
# the dataset: subjects, doses, observations
# ---------------------------------------------------------------------------

# number formatting in notes and labels (shortest form)
hp_short <- function(x) format(x, digits = 15, trim = TRUE, scientific = FALSE, drop0trailing = TRUE)

#' Read a concentration-time dataset into subjects, doses and observations
#'
#' The first step of an NCA. It turns the records of a NONMEM-style dataset into
#' one object per subject, holding the doses (with additional doses `ADDL` and
#' `II` expanded and steady-state flags kept) and the observations, with the
#' route of every dose decided and the BLQ status of every sample set.
#'
#' @details
#' **Records.** A dose is a record with `EVID` 1 or 4 and `AMT > 0` (when there
#' is no `EVID` column, `AMT > 0`). An observation is a record with `EVID` 0,
#' `MDV` 0 and a numeric `DV`. A subject is a run of rows with the same `ID`.
#'
#' **Additional doses.** A dose record with `ADDL` and `II` (with `II` above 0)
#' is expanded into the repeated doses. With the `SS` flag and `II` above 0 the
#' first of them starts at steady state.
#'
#' **Route.** `route = "auto"` decides per dose: a record with `RATE` above 0 (or
#' -1 or -2) is an infusion; a dose into the compartment that holds most of the
#' observations is an intravenous bolus; any other dose is extravascular. When no
#' dose has a rate and the pooled profiles rise after the dose, intravenous
#' doses are re-classified as extravascular. Because this is a guess, state the
#' route explicitly when you know it: `"oral"` (any extravascular route),
#' `"iv"` (bolus, or infusion when there is a rate), `"bolus"`, `"infusion"` or
#' `"mixed"` (a `CMT` column tells oral and intravenous doses apart).
#'
#' **BLQ.** A sample is BLQ when the `blq_flag` column is non-zero, when it is
#' below `lloq`, or when it is not positive.
#'
#' **Filters.** `filter` and `obs_filter` are R expressions (see
#' [hp_filter()]); `obs_filter` removes observations only, never doses.
#'
#' @param data A data frame, usually from [hp_read_data()].
#' @param route `"auto"`, `"oral"`, `"iv"`, `"bolus"`, `"infusion"` or
#'   `"mixed"`.
#' @param id,time,dv,amt,evid,mdv,cmt,rate,addl,ii,ss Names of the columns that
#'   hold the subject, time, concentration, dose amount, event type, missing-DV
#'   flag, compartment, infusion rate, additional doses, interdose interval and
#'   steady-state flag.
#' @param groups Names of columns that define groups, for example
#'   `"DOSEGRP"`. A group is labelled with the column and the value, such as
#'   `DOSEGRP=100`.
#' @param lloq Limit of quantification; observations below it are BLQ.
#' @param blq_flag Name of a column that flags BLQ observations (non-zero).
#' @param dv_log `TRUE` when `DV` is the natural logarithm of the
#'   concentration.
#' @param filter R expression of the rows to keep, applied to every row.
#' @param obs_filter R expression of the observations to keep.
#' @return A list with `subjects` (one list per subject with its `id`,
#'   `group`, `doses` (a data frame with `time`, `amt`, `rate`, `kind`,
#'   `ss_ii` and `group`) and `obs` (a data frame with `time`, `conc`, `blq` and
#'   `row`)), `route` (the route setting used) and `notes` (messages about
#'   filters).
#' @seealso [hp_nca()], which calls it, [hp_nca_profiles()]
#' @family non-compartmental analysis
#' @examples
#' d <- hp_read_data(hp_example_file())
#' x <- hp_nca_read(d, route = "oral", lloq = 0.5)
#' length(x$subjects)
#' x$subjects[[1]]$doses
#' head(x$subjects[[1]]$obs)
#'
#' # a group column is read at the dose record
#' x <- hp_nca_read(d, route = "oral", groups = "DOSEGRP")
#' unique(vapply(x$subjects, function(s) s$group, ""))
#' @export
hp_nca_read <- function(data, route = "auto", id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                        evid = "EVID", mdv = "MDV", cmt = "CMT", rate = "RATE", addl = "ADDL",
                        ii = "II", ss = "SS", groups = NULL, lloq = NULL, blq_flag = NULL,
                        dv_log = FALSE, filter = NULL, obs_filter = NULL) {
  names(data) <- toupper(names(data))
  up <- toupper
  data$.ROW <- seq_len(nrow(data)) + 1L
  notes <- character()
  need <- function(n, what) if (!up(n) %in% names(data)) hp_stop("the ", what, " column `", n, "` is not in the dataset")
  need(id, "id"); need(time, "time"); need(dv, "dv"); need(amt, "amt")
  if (!is.null(filter)) {
    n0 <- nrow(data); data <- hp_filter(data, filter)
    notes <- c(notes, sprintf("the filter `%s` left out %d of %d rows", filter, n0 - nrow(data), n0))
  }
  g <- function(n, default = NA_real_) if (up(n) %in% names(data)) hp_num(data[[up(n)]]) else rep(default, nrow(data))
  tm <- g(time); a <- g(amt); a[!is.finite(a)] <- 0; y <- g(dv)
  ev <- g(evid); ev <- ifelse(is.finite(ev), as.integer(ev), ifelse(a > 0, 1L, 0L))
  md <- g(mdv); md <- ifelse(is.finite(md), as.integer(md), ifelse(ev != 0 | !is.finite(y), 1L, 0L))
  if (any(!is.finite(tm))) hp_stop("TIME is missing or not a number")
  rt <- g(rate); rt[!is.finite(rt)] <- 0
  ad <- g(addl); ad[!is.finite(ad)] <- 0
  iv_ <- g(ii); iv_[!is.finite(iv_)] <- 0
  sv <- g(ss); sv[!is.finite(sv)] <- 0
  cm <- g(cmt)
  ids <- as.character(data[[up(id)]])
  is_dose <- (ev == 1L | ev == 4L) & a > 0
  is_obs <- ev == 0L & md == 0L & is.finite(y)
  if (!is.null(obs_filter)) {
    keep <- eval(parse(text = obs_filter)[[1L]], envir = hp_eval_env(data)); keep[is.na(keep)] <- FALSE
    dropped <- is_obs & !keep
    is_obs <- is_obs & keep
    if (any(dropped)) notes <- c(notes, sprintf("the observation filter left out %d observation(s)", sum(dropped)))
  }
  obs_cmt <- if (any(is_obs & is.finite(cm))) {
    tb <- table(cm[ev == 0L & md == 0L & is.finite(cm)]); as.numeric(names(tb)[which.max(tb)])
  } else NA_real_
  route0 <- hp_choice(route, list(auto = "auto", oral = c("oral", "po", "extravascular", "ev", "sc", "im"),
                                  iv = c("iv", "intravenous"), bolus = c("bolus", "iv-bolus"),
                                  infusion = c("infusion", "iv-infusion", "inf"),
                                  mixed = c("mixed", "both", "oral+iv", "iv+oral")), "route")
  if (route0 == "mixed" && (is.na(obs_cmt))) hp_stop("route mixed needs a CMT column to tell oral from intravenous doses")
  has_rate <- rt > 0 | rt == -1 | rt == -2
  by_cmt <- ifelse(!is.na(obs_cmt) & is.finite(cm) & cm == obs_cmt, "bolus", "oral")
  kind <- ifelse(has_rate, "infusion", by_cmt)
  if (route0 == "bolus") kind[] <- "bolus"
  if (route0 == "infusion") kind[] <- "infusion"
  if (route0 == "oral") kind <- ifelse(has_rate, "infusion", "oral")
  if (route0 == "iv") kind <- ifelse(has_rate, "infusion", "bolus")
  gcols <- up(groups %||% character())
  for (gc in gcols) need(gc, "group")
  glabel <- function(rows) {
    if (!length(gcols)) return(rep("all", length(rows)))
    do.call(paste, c(lapply(gcols, function(gc) {
      v <- data[[gc]][rows]; vn <- hp_num(v)
      paste0(gc, "=", ifelse(is.finite(vn), vapply(vn, hp_short, ""), trimws(as.character(v))))
    }), sep = ", "))
  }
  flag <- if (!is.null(blq_flag)) { need(blq_flag, "blq"); hp_num(data[[up(blq_flag)]]) } else rep(NA_real_, nrow(data))
  starts <- which(c(TRUE, ids[-1L] != ids[-length(ids)]))
  ends <- c(starts[-1L] - 1L, length(ids))
  subjects <- vector("list", length(starts))
  for (s in seq_along(starts)) {
    rows <- starts[s]:ends[s]
    d_rows <- rows[is_dose[rows]]
    dl <- list()
    for (r in d_rows) {
      nrep <- if (iv_[r] > 0) max(ad[r], 0) else 0
      for (j in 0:nrep) dl[[length(dl) + 1L]] <- data.frame(
        time = tm[r] + j * iv_[r], amt = a[r], rate = if (rt[r] > 0) rt[r] else 0, kind = kind[r],
        ss_ii = if (j == 0 && sv[r] > 0 && iv_[r] > 0) iv_[r] else 0, group = glabel(r),
        stringsAsFactors = FALSE)
    }
    doses <- if (length(dl)) do.call(rbind, dl) else data.frame(time = numeric(), amt = numeric(), rate = numeric(),
                                                                  kind = character(), ss_ii = numeric(), group = character())
    if (nrow(doses)) doses <- doses[order(doses$time), , drop = FALSE]
    o_rows <- rows[is_obs[rows]]
    conc <- if (dv_log) exp(y[o_rows]) else y[o_rows]
    flagged <- is.finite(flag[o_rows]) & flag[o_rows] != 0
    blq <- flagged | (!is.null(lloq) & conc < (lloq %||% -Inf)) | !(conc > 0)
    obs <- data.frame(time = tm[o_rows], conc = conc, blq = blq, row = data$.ROW[o_rows], stringsAsFactors = FALSE)
    if (nrow(obs)) obs <- obs[order(obs$time), , drop = FALSE]
    subjects[[s]] <- list(id = ids[rows[1L]], group = glabel(rows[1L]), doses = doses, obs = obs)
  }
  if (!length(subjects)) hp_stop("the dataset has no subjects")
  if (!sum(vapply(subjects, function(s) nrow(s$obs), 1L))) hp_stop("no observations are left after the filters")
  if (!sum(vapply(subjects, function(s) nrow(s$doses), 1L))) hp_stop("the dataset has no dose records (AMT > 0)")
  # route "auto": the shape of the pooled profiles decides between oral and intravenous
  if (route0 == "auto" && !any(has_rate & is_dose)) {
    rises <- vapply(subjects, function(s) {
      if (nrow(s$doses) == 0L || nrow(s$obs) < 3L) return(NA)
      d <- s$doses$time[1L]; nx <- if (nrow(s$doses) > 1L) s$doses$time[2L] else Inf
      o <- s$obs[!s$obs$blq & s$obs$time > d + 1e-9 & s$obs$time < nx - 1e-9, ]
      if (nrow(o) < 3L || !(max(o$conc) > 0)) return(NA)
      im <- which.max(o$conc)
      im > 1L && o$conc[1L] / max(o$conc) < 0.85
    }, NA)
    if (any(!is.na(rises)) && mean(rises, na.rm = TRUE) >= 0.5) {
      for (s in seq_along(subjects)) subjects[[s]]$doses$kind[subjects[[s]]$doses$kind == "bolus"] <- "oral"
    }
  }
  list(subjects = subjects, route = route0, notes = notes)
}

# ---------------------------------------------------------------------------
# BLQ rules and duplicate samples
# ---------------------------------------------------------------------------

hp_status_label <- c(observed = "observed", blq_zero = "blq set to zero", blq_lloq = "blq set to lloq",
                     blq_half = "blq set to lloq/2", blq_missing = "blq left out",
                     negative_missing = "negative left out", boundary = "next-dose sample (end of interval)",
                     inserted_zero = "added zero at dose", inserted_c0 = "added C0",
                     predose_replaced = "zero at dose replaced by C0")
hp_status_used <- function(s) !(s %in% c("blq_missing", "negative_missing", "predose_replaced"))
hp_status_sample <- function(s) s %in% c("observed", "blq_zero", "blq_lloq", "blq_half")

# pts: data frame t, c, censored, row (sorted by time)
hp_nca_combine <- function(pts, how, id) {
  if (!nrow(pts)) return(list(pts = pts, n = 0L))
  grp <- cumsum(c(TRUE, abs(diff(pts$t)) > 1e-9 * (1 + abs(pts$t[-1L]))))
  if (!anyDuplicated(grp)) return(list(pts = pts, n = 0L))
  out <- list(); n <- 0L
  for (gi in unique(grp)) {
    g <- pts[grp == gi, , drop = FALSE]
    if (nrow(g) == 1L) { out[[length(out) + 1L]] <- g; next }
    n <- n + nrow(g) - 1L
    if (how == "error") hp_stop("subject ", id, ": ", nrow(g), " samples at time ", hp_short(g$t[1L]),
                                "; give the samples distinct times, or set duplicates to mean, first, last or max")
    live <- g[!g$censored, , drop = FALSE]
    if (!nrow(live)) { out[[length(out) + 1L]] <- g[1L, , drop = FALSE]; next }
    v <- live$c
    cc <- switch(how, mean = mean(v), first = v[1L], last = v[length(v)], max = max(v), mean(v))
    out[[length(out) + 1L]] <- data.frame(t = live$t[1L], c = cc, censored = FALSE, row = live$row[1L])
  }
  list(pts = do.call(rbind, out), n = n)
}

hp_nca_rule <- function(pts, rule, lloq) {
  n <- nrow(pts)
  if (!n) return(data.frame(t = numeric(), c = numeric(), raw = numeric(), censored = logical(),
                            status = character(), row = integer(), stringsAsFactors = FALSE))
  q <- pts[!pts$censored & pts$c > 0 & is.finite(pts$c), , drop = FALSE]
  tmax <- NA_real_
  if (nrow(q)) { best <- -Inf; for (i in seq_len(nrow(q))) if (q$c[i] > best) { best <- q$c[i]; tmax <- q$t[i] } }
  first_q <- if (nrow(q)) q$t[1L] else Inf
  last_q <- if (nrow(q)) q$t[nrow(q)] else -Inf
  limit <- function(i) {
    l <- lloq %||% pts$c[i]
    if (is.finite(l) && l > 0) l else hp_stop("the BLQ rule `", rule, "` needs the limit of quantification, but the sample at line ",
                                             pts$row[i], " (time ", hp_short(pts$t[i]), ") has none; give lloq, or code BLQ samples with DV = LLOQ")
  }
  cc <- pts$c; st <- rep("observed", n)
  for (i in seq_len(n)) {
    if (!pts$censored[i]) { if (pts$c[i] < 0) { cc[i] <- NA_real_; st[i] <- "negative_missing" }; next }
    set <- function(v, s) { cc[i] <<- v; st[i] <<- s }
    switch(rule,
      zero = set(0, "blq_zero"), missing = set(NA_real_, "blq_missing"),
      lloq = set(limit(i), "blq_lloq"), "half-lloq" = set(limit(i) / 2, "blq_half"),
      "tmax-zero-missing" = if (!is.finite(tmax) || pts$t[i] < tmax) set(0, "blq_zero") else set(NA_real_, "blq_missing"),
      "first-zero-mid-half-last-missing" = if (pts$t[i] < first_q) set(0, "blq_zero")
        else if (pts$t[i] < last_q) set(limit(i) / 2, "blq_half") else set(NA_real_, "blq_missing"))
  }
  data.frame(t = pts$t, c = cc, raw = pts$c, censored = pts$censored, status = st, row = pts$row,
             stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# profiles: one per dosing interval of each subject
# ---------------------------------------------------------------------------

#' Cut subjects into analysis profiles
#'
#' Turns the subjects of [hp_nca_read()] into profiles: one per dosing interval
#' of each subject, with times counted from the dose, BLQ and duplicate rules
#' applied, and the dosing interval and steady state decided. This is the step
#' where the choices of [hp_nca_options()] about BLQ data, repeated samples and
#' steady state take effect.
#'
#' @details
#' Each dose opens an interval that lasts until the next dose (or the last
#' sample). At steady state the interval is cut at `tau`, and a sample taken at the
#' end of the interval (the time of the next dose) is kept as the
#' end-of-interval value. Samples before the first dose
#' are not analysed. A profile with fewer than `min_observations` usable
#' samples, or with no quantifiable concentration, is kept with a reason in its
#' `skip` element and appears in the results as skipped.
#'
#' @param x The result of [hp_nca_read()].
#' @param options Settings from [hp_nca_options()].
#' @return A list of profiles. Each is a list with `id`, `group`, `occasion`
#'   (the number of the dose), `dose_time`, `dose`, `rate`, `route`
#'   (`"extravascular"`, `"iv bolus"` or `"iv infusion"`), `tau`,
#'   `steady_state`, `ss_without_tau`, `has_next_dose`, `points` (a data frame
#'   with `t`, `c`, `raw`, `censored`, `status` and `row`), `boundary` and
#'   `skip`. The attribute `"notes"` holds messages about what was done with the
#'   data.
#' @seealso [hp_nca_analyse()], [hp_nca()]
#' @family non-compartmental analysis
#' @examples
#' d <- hp_read_data(hp_example_file())
#' x <- hp_nca_read(d, route = "oral", lloq = 0.5)
#' p <- hp_nca_profiles(x, hp_nca_options(lloq = 0.5))
#' length(p)
#' p[[1]]$points
#'
#' # the rule applied to BLQ samples is recorded for every point
#' table(unlist(lapply(p, function(z) z$points$status)))
#' @export
hp_nca_profiles <- function(x, options = hp_nca_options()) {
  s <- options
  out <- list(); notes <- character(); n_dup <- 0L; n_pre <- 0L; n_neg <- 0L; unresolved <- 0L
  tol <- function(a, b) abs(a - b) <= 1e-9 * (1 + abs(a))
  for (sub in x$subjects) {
    dd <- sub$doses
    if (!nrow(dd)) next
    # merge records at the same time
    mg <- list()
    for (i in seq_len(nrow(dd))) {
      d <- dd[i, ]
      if (length(mg) && tol(mg[[length(mg)]]$time, d$time)) {
        m <- mg[[length(mg)]]; m$amt <- m$amt + d$amt; m$rate <- m$rate + d$rate
        if (d$ss_ii > 0 && m$ss_ii == 0) m$ss_ii <- d$ss_ii
        mg[[length(mg)]] <- m
      } else mg[[length(mg) + 1L]] <- as.list(d)
    }
    m <- length(mg)
    times <- vapply(mg, function(z) z$time, 1)
    n_pre <- n_pre + sum(sub$obs$time < times[1L] - 1e-9)
    picked <- switch(s$intervals, all = seq_len(m), first = 1L, last = m)
    regular <- m >= 3L && { g0 <- times[2L] - times[1L]; all(abs(diff(times) - g0) <= 0.01 * abs(g0) + 1e-9) }
    declared <- s$steady_state %in% c("last", "all")
    for (k in picked) {
      d <- mg[[k]]; start <- d$time
      nxt <- if (k < m) times[k + 1L] - start else NULL
      tau_flag <- if (d$ss_ii > 0) d$ss_ii else NULL
      gap <- if (regular || declared) nxt else NULL
      tau <- tau_flag %||% s$tau %||% gap
      ss_no_tau <- FALSE
      ss <- switch(s$steady_state, none = FALSE, auto = !is.null(tau_flag), last = k == m, all = TRUE)
      if (ss && is.null(tau)) {
        if (k > 1L) tau <- start - times[k - 1L] else { ss <- FALSE; ss_no_tau <- TRUE; unresolved <- unresolved + 1L }
      }
      end <- Inf
      if (!is.null(nxt)) end <- min(end, nxt)
      if (ss && !is.null(tau)) end <- min(end, tau + 1e-9)
      if (!is.null(s$max_time)) end <- min(end, s$max_time + 1e-9)
      cens <- function(b) b && (!is.null(s$lloq) || s$blq_flagged)
      boundary <- NULL
      o_t <- sub$obs$time; o_c <- sub$obs$conc; o_b <- sub$obs$blq; o_r <- sub$obs$row
      rel_all <- o_t - start
      take <- rep(FALSE, length(o_t))
      for (i in seq_along(o_t)) {
        rel <- rel_all[i]
        if (rel < -1e-9) next
        if (!is.null(nxt)) {
          if (abs(rel - nxt) <= 1e-9 * (1 + abs(nxt))) {
            boundary <- data.frame(t = nxt, c = o_c[i], censored = cens(o_b[i]), row = o_r[i]); next
          }
          if (rel > nxt) next
        }
        if (rel > end) next
        take[i] <- TRUE
      }
      cens_all <- if (!is.null(s$lloq) || s$blq_flagged) o_b else rep(FALSE, length(o_b))
      raw <- data.frame(t = pmax(rel_all[take], 0), c = o_c[take], censored = cens_all[take], row = o_r[take])
      cmb <- hp_nca_combine(raw, s$duplicates, sub$id)
      n_dup <- n_dup + cmb$n
      pts <- hp_nca_rule(cmb$pts, s$blq_rule, s$lloq)
      n_neg <- n_neg + sum(pts$status == "negative_missing")
      if (!is.null(boundary)) {
        b <- hp_nca_rule(boundary, s$blq_rule, s$lloq)
        boundary <- if (hp_status_used(b$status[1L])) { b$status[1L] <- "boundary"; b } else NULL
      }
      route <- switch(d$kind, oral = "extravascular", bolus = "iv bolus",
                      infusion = if (d$rate > 0) "iv infusion" else "iv bolus")
      prof <- list(id = sub$id, group = d$group,
                   occasion = k, dose_time = start, dose = d$amt, rate = if (d$kind == "infusion") d$rate else 0,
                   route = route, tau = tau, steady_state = ss, ss_without_tau = ss_no_tau,
                   has_next_dose = !is.null(nxt), points = pts, boundary = boundary, skip = NULL)
      usable <- sum(hp_status_used(pts$status)); quant <- sum(hp_status_used(pts$status) & pts$c > 0, na.rm = TRUE)
      if (usable < max(s$min_observations, 1)) {
        prof$skip <- sprintf("%d usable sample(s) (at least %d are needed)", usable, s$min_observations)
      } else if (quant == 0) prof$skip <- "no quantifiable concentration"
      out[[length(out) + 1L]] <- prof
    }
  }
  if (n_pre > 0) notes <- c(notes, sprintf("%d sample(s) taken before a subject's first dose are not analysed", n_pre))
  if (n_dup > 0) notes <- c(notes, sprintf("%d sample(s) at a time that already had a sample were combined (%s)", n_dup, s$duplicates))
  if (n_neg > 0) notes <- c(notes, sprintf("%d negative concentration(s) were left out", n_neg))
  if (unresolved > 0) notes <- c(notes, sprintf("%d steady-state interval(s) have no dosing interval (give II with the SS flag, or tau) and are analysed as single doses", unresolved))
  attr(out, "notes") <- notes
  out
}

# ---------------------------------------------------------------------------
# the terminal rate constant
# ---------------------------------------------------------------------------

hp_nca_lz <- function(obs, tmax, lz, ov = NULL) {
  none <- function(how, flags = "NO_LZ", cand = NULL) list(fit = NULL, candidates = cand, how = how, flags = flags, times = numeric())
  skip <- isTRUE(ov$skip)
  if (lz$method == "none" || skip)
    return(none(if (skip) "no terminal phase (override)" else "not estimated (method none)"))
  pos <- obs[is.finite(obs$c) & obs$c > 0, , drop = FALSE]
  after <- if (lz$include_tmax) pos$t >= tmax - 1e-12 else pos$t > tmax + 1e-12
  in_range <- function(df, from, to) (is.null(from) | df$t >= (from %||% -Inf) - 1e-9) & (is.null(to) | df$t <= (to %||% Inf) + 1e-9)
  excluded <- function(df) if (length(ov$exclude)) vapply(df$t, function(t) any(abs(ov$exclude - t) < 1e-9), NA) else rep(FALSE, nrow(df))
  w_code <- c(uniform = 0L, "1/y" = 1L, "1/y2" = 2L)[[lz$weighting]]
  manual_range <- !is.null(ov$from) || !is.null(ov$to)
  rng_text <- function(f, t) paste(if (is.null(f)) "start" else hp_short(f), "to", if (is.null(t)) "end" else hp_short(t))
  if (manual_range) {
    pool <- pos[in_range(pos, ov$from, ov$to) & !excluded(pos), , drop = FALSE]
    sizes <- nrow(pool); how <- paste0("set for this subject: ", rng_text(ov$from, ov$to)); manual <- TRUE
  } else if (!is.null(ov$points)) {
    pool <- pos[!excluded(pos), , drop = FALSE]; sizes <- ov$points
    how <- sprintf("last %d points (set for this subject)", ov$points); manual <- TRUE
  } else {
    pool <- pos[after & !excluded(pos), , drop = FALSE]
    if (lz$method == "last-n") { sizes <- lz$points; how <- sprintf("last %d points", lz$points); manual <- TRUE }
    else if (lz$method == "time-range") {
      pool <- pool[in_range(pool, lz$from, lz$to), , drop = FALSE]; sizes <- nrow(pool)
      how <- paste("time range", rng_text(lz$from, lz$to)); manual <- TRUE
    } else {
      maxn <- if (lz$max_points == 0) nrow(pool) else min(lz$max_points, nrow(pool))
      sizes <- if (maxn >= max(lz$min_points, 2)) seq.int(max(lz$min_points, 2), maxn) else integer()
      how <- "best fit"; manual <- FALSE
    }
  }
  sizes <- as.integer(sizes)
  if (manual && length(sizes) == 1L) sizes <- min(sizes, nrow(pool))
  cand <- hp_lz_windows(pool$t, pool$c, sizes, w_code)
  if (!nrow(cand)) {
    r <- none(if (all(sizes < 2L) || !length(sizes)) "too few points after the maximum" else "the regression failed")
    if (manual) r$how <- paste0(how, ": ", r$how)
    return(r)
  }
  cd <- as.data.frame(cand)
  cd$half_life <- ifelse(cd$lambda > 0, log(2) / cd$lambda, NA_real_)
  cd$span_ratio <- ifelse(is.finite(cd$half_life), (cd$end - cd$start) / cd$half_life, NA_real_)
  cd$selected <- FALSE
  chosen <- if (manual) {
    which(cd$lambda > 0)[1L]
  } else {
    el <- which(cd$lambda > 0 & is.finite(cd$adj_r2))
    if (length(el)) { best <- max(cd$adj_r2[el]); ok <- el[cd$adj_r2[el] >= best - lz$tolerance]; ok[which.max(cd$n[ok])] } else NA_integer_
  }
  if (is.na(chosen)) return(none(if (manual) "the chosen points do not decline (the slope is not negative)"
                                  else "no regression with a negative slope and an adjusted R\u00b2", cand = cd))
  f <- cd[chosen, ]
  flags <- character()
  if (!manual && f$adj_r2 < lz$min_adj_r2)
    return(none(sprintf("best adjusted R\u00b2 %.4f is below the minimum %s", f$adj_r2, hp_short(lz$min_adj_r2)),
                flags = c("NO_LZ", "LZ_R2ADJ_LOW"), cand = cd))
  if (manual && is.finite(f$adj_r2) && f$adj_r2 < lz$min_adj_r2) flags <- c(flags, "LZ_R2ADJ_LOW")
  if (is.finite(f$span_ratio) && f$span_ratio < lz$min_span_ratio) {
    if (lz$require_span)
      return(none(sprintf("the terminal phase spans %.2f half-lives, fewer than %s", f$span_ratio, hp_short(lz$min_span_ratio)),
                  flags = c("NO_LZ", "LZ_SPAN_SHORT"), cand = cd))
    flags <- c(flags, "LZ_SPAN_SHORT")
  }
  if (f$n < 3) flags <- c(flags, "LZ_FEW_POINTS")
  cd$selected[chosen] <- TRUE
  win <- pool[pool$t >= f$start - 1e-12 & pool$t <= f$end + 1e-12, , drop = FALSE]
  list(fit = as.list(cd[chosen, ]), candidates = cd, how = how, flags = flags, times = win$t)
}

# ---------------------------------------------------------------------------
# the parameters of one profile
# ---------------------------------------------------------------------------

hp_tail <- function(t0, c0, t1, lambda) {
  c1 <- c0 * exp(-lambda * (t1 - t0))
  c(c0 * (-expm1(-lambda * (t1 - t0))) / lambda, hp_seg_aumc(t0, c0, t1, c1, TRUE))
}

#' Compute the NCA parameters of one profile
#'
#' Analyses a single profile from [hp_nca_profiles()]: peak and trough
#' measures, the terminal slope, areas and moments, mean residence time,
#' clearance and volume, and the steady-state measures when the profile is at
#' steady state. Use it to look closely at one subject; [hp_nca()] applies it to
#' all.
#'
#' @details
#' Areas to the last measurable concentration use the trapezoid method of the
#' options. The area to infinity adds the extrapolated tail, computed twice:
#' from the observed last concentration (`AUCINF_OBS`) and from the
#' concentration predicted by the terminal regression (`AUCINF_PRED`). Clearance
#' and volumes use the dose times `dose_factor`. At steady state the
#' parameters over the dosing interval are computed instead of the
#' infinity-based ones. Flags that explain the quality of the result are
#' collected in the `flags` element:
#'
#' - `NO_LZ`: no terminal slope could be estimated.
#' - `LZ_R2ADJ_LOW`, `LZ_SPAN_SHORT`, `LZ_FEW_POINTS`: the adjusted R-squared of
#'   the terminal regression is below the minimum, its span is shorter than the
#'   minimum, or it is based on fewer than three points.
#' - `AUC_EXTRAP_GT_<n>`: more than n percent of the area to infinity is
#'   extrapolated.
#' - `INF_OMITTED`: infinity-based parameters left out for that reason.
#' - `INF_AFTER_EARLIER_DOSES`: an area to infinity after earlier doses, which
#'   includes the tail of those doses.
#' - `C0_BACKEXTRAP`, `C0_FIRST_SAMPLE`: how C0 of an intravenous bolus was
#'   obtained.
#' - `AUC_FROM_FIRST_SAMPLE`: the area starts at the first sample because
#'   there is none at the dose time.
#' - `SS_NO_TAU`, `TAU_NO_EXTRAP`: steady-state interval without a dosing
#'   interval, or one that cannot be completed.
#'
#' @param p A profile from [hp_nca_profiles()].
#' @param options Settings from [hp_nca_options()].
#' @return A list with `id`, `group`, `occasion`, `route`, `dose`, `tau`,
#'   `steady_state`, `params` (named list of the parameter values), `lz` (the
#'   terminal-slope search: `fit`, the chosen regression; `candidates`, every
#'   regression tried; `how`, how it was chosen; `flags`; `times`, the times
#'   used), `flags`, `points` (every point, including the zero or C0 added at
#'   the dose) and `skip` (the reason when the profile was not analysed).
#' @seealso [hp_nca()], [hp_nca_parameters()]
#' @family non-compartmental analysis
#' @examples
#' d <- hp_read_data(hp_example_file())
#' o <- hp_nca_options(lloq = 0.5)
#' p <- hp_nca_profiles(hp_nca_read(d, route = "oral", lloq = 0.5), o)
#' res <- hp_nca_analyse(p[[1]], o)
#' unlist(res$params[c("CMAX", "TMAX", "AUCLAST", "LAMZ", "THALF", "AUCINF_OBS")])
#' res$flags
#' res$lz$how
#' @export
hp_nca_analyse <- function(p, options = hp_nca_options()) {
  s <- options
  res <- list(id = p$id, group = p$group, occasion = p$occasion, route = p$route, dose = p$dose,
              tau = p$tau, steady_state = p$steady_state, params = list(), lz = NULL, flags = character(),
              points = p$points, skip = p$skip)
  if (!is.null(p$skip)) return(res)
  dose <- p$dose; dose_eff <- dose * s$units$dose_factor
  iv <- p$route != "extravascular"; ss <- p$steady_state; tau <- p$tau
  t_inf <- if (p$route == "iv infusion" && p$rate > 0) dose / p$rate else 0
  flags <- character()
  if (isTRUE(p$ss_without_tau)) flags <- c(flags, "SS_NO_TAU")
  pt <- p$points
  smp <- pt[hp_status_used(pt$status) & hp_status_sample(pt$status), , drop = FALSE]
  cmax <- -Inf; tmax <- 0
  for (i in seq_len(nrow(smp))) if (smp$c[i] > cmax) { cmax <- smp$c[i]; tmax <- smp$t[i] }
  tlast <- clast <- NA_real_
  for (i in seq_len(nrow(smp))) if (smp$c[i] > 0) { tlast <- smp$t[i]; clast <- smp$c[i] }
  par <- list(); set <- function(n, v) par[[n]] <<- v
  set("DOSE", dose); set("N_SAMPLES", nrow(smp))
  if (t_inf > 0) set("INF_DURATION", t_inf)
  if (!is.null(tau)) set("TAU", tau)
  cp <- data.frame(t = smp$t, c = smp$c)
  added <- NULL
  first_t <- if (nrow(cp)) cp$t[1L] else 0
  at_zero <- abs(first_t) <= 1e-9
  c0_value <- NA_real_; curve_from_first <- FALSE
  tau_sample <- if (!is.null(tau)) {
    if (!is.null(p$boundary) && abs(p$boundary$t - tau) <= 1e-9 * (1 + abs(tau))) p$boundary$c
    else { i <- which(abs(cp$t - tau) <= 1e-9 * (1 + abs(tau))); if (length(i)) cp$c[i[1L]] else NULL }
  } else NULL
  positive <- cp[cp$t > 1e-9 & cp$c > 0, , drop = FALSE]
  push_added <- function(c, st) added <<- rbind(added, data.frame(t = 0, c = c, raw = NA_real_, censored = FALSE, status = st, row = 0L))
  if (p$route == "iv bolus") {
    observed_c0 <- if (at_zero) cp$c[1L] else NULL
    usable_obs <- !is.null(observed_c0) && observed_c0 > 0
    if (s$c0 == "observed" || (usable_obs && s$c0 != "first")) {
      if (usable_obs) c0_value <- observed_c0 else curve_from_first <- !at_zero
    } else {
      est <- NULL
      if (s$c0 == "first") est <- if (nrow(positive)) positive$c[1L] else NULL
      else if (nrow(positive) >= 2L && positive$c[2L] < positive$c[1L]) {
        lam <- log(positive$c[1L] / positive$c[2L]) / (positive$t[2L] - positive$t[1L])
        est <- positive$c[1L] * exp(lam * positive$t[1L])
      } else { flags <- c(flags, "C0_FIRST_SAMPLE"); est <- if (nrow(positive)) positive$c[1L] else NULL }
      if (!is.null(est) && is.finite(est) && est > 0) {
        if (at_zero) {
          old <- cp[1L, ]; cp <- cp[-1L, , drop = FALSE]
          hit <- res$points$t == 0 & hp_status_sample(res$points$status) & res$points$c == old$c
          hit[is.na(hit)] <- FALSE
          if (any(hit)) res$points$status[which(hit)[1L]] <- "predose_replaced"
        }
        cp <- rbind(data.frame(t = 0, c = est), cp)
        push_added(est, "inserted_c0"); c0_value <- est
        if (s$c0 == "extrapolate" && !"C0_FIRST_SAMPLE" %in% flags) flags <- c(flags, "C0_BACKEXTRAP")
      } else curve_from_first <- !at_zero
    }
  } else if (ss && !at_zero) {
    if (!is.null(tau_sample) && is.finite(tau_sample)) {
      cp <- rbind(data.frame(t = 0, c = tau_sample), cp); push_added(tau_sample, "inserted_c0"); c0_value <- tau_sample
    } else curve_from_first <- TRUE
  } else if (!ss && !at_zero) {
    if (s$zero_at_dose && p$occasion == 1L) { cp <- rbind(data.frame(t = 0, c = 0), cp); push_added(0, "inserted_zero") }
    else curve_from_first <- TRUE
  } else if (at_zero) c0_value <- cp$c[1L]
  if (curve_from_first) flags <- c(flags, "AUC_FROM_FIRST_SAMPLE")
  if (p$route == "iv bolus" && is.finite(c0_value)) set("C0", c0_value)
  if (!is.null(added)) res$points <- rbind(res$points, added)
  set("CMAX", cmax); set("TMAX", tmax); set("CMAX_D", cmax / dose)
  if (!is.null(tau)) {
    cmin <- Inf; tmin <- 0
    for (i in seq_len(nrow(smp))) if (smp$c[i] < cmin) { cmin <- smp$c[i]; tmin <- smp$t[i] }
    if (!is.null(p$boundary) && p$boundary$c < cmin) { cmin <- p$boundary$c; tmin <- p$boundary$t }
    set("CMIN", cmin); set("TMIN", tmin)
  }
  set("CLAST", clast); set("TLAST", tlast)
  if (p$route == "extravascular" && !ss) {
    i <- which(smp$c > 0)[1L]
    tlag <- if (is.na(i)) NA_real_ else if (i == 1L) {
      zero_first <- nrow(cp) > 0 && cp$c[1L] == 0 && cp$t[1L] == 0
      if (abs(smp$t[1L]) <= 1e-9 || zero_first) 0 else NA_real_
    } else smp$t[i - 1L]
    set("TLAG", tlag)
  }
  curve_tmax <- if (is.finite(c0_value) && c0_value > cmax && p$route == "iv bolus") 0 else tmax
  t_start <- if (nrow(cp)) cp$t[1L] else 0
  keep <- c(TRUE, diff(cp$t) > 1e-12); cp <- cp[keep, , drop = FALSE]
  mcode <- hp_method_code(s$auc_method)
  auc_last <- hp_curve_area(cp$t, cp$c, mcode, curve_tmax, t_start, tlast, FALSE)
  aumc_last <- hp_curve_area(cp$t, cp$c, mcode, curve_tmax, t_start, tlast, TRUE)
  set("AUCLAST", auc_last); set("AUCLAST_D", auc_last / dose); set("AUMCLAST", aumc_last)
  mrt_corr <- if (p$route == "iv infusion") t_inf / 2 else 0
  if (auc_last > 0) set("MRTLAST", aumc_last / auc_last - mrt_corr)
  obs_pos <- smp[smp$c > 0, c("t", "c"), drop = FALSE]
  ov <- NULL
  for (o in s$lz$overrides) if (identical(as.character(o$id), as.character(p$id)) && (is.null(o$occasion) || o$occasion == p$occasion)) { ov <- o; break }
  lz <- hp_nca_lz(obs_pos, tmax, s$lz, ov)
  flags <- c(flags, lz$flags)
  lam <- if (!is.null(lz$fit)) lz$fit$lambda else NULL
  clast_pred <- NA_real_
  if (!is.null(lz$fit)) {
    f <- lz$fit; clast_pred <- exp(f$intercept - f$lambda * tlast)
    set("LAMZ", f$lambda); set("THALF", f$half_life); set("LAMZ_N", f$n); set("LAMZ_START", f$start)
    set("LAMZ_END", f$end); set("LAMZ_SPAN", f$span_ratio); set("R2", f$r2); set("R2ADJ", f$adj_r2)
    set("CLAST_PRED", clast_pred)
  }
  if (!is.null(lam) && !ss) {
    l <- lam
    inf_obs <- auc_last + clast / l; inf_pred <- auc_last + clast_pred / l
    pct_obs <- 100 * (inf_obs - auc_last) / inf_obs; pct_pred <- 100 * (inf_pred - auc_last) / inf_pred
    set("AUC_PCT_EXTRAP_OBS", pct_obs); set("AUC_PCT_EXTRAP_PRED", pct_pred)
    omit <- pct_obs > s$lz$max_extrap && s$lz$extrap_action == "omit"
    if (pct_obs > s$lz$max_extrap) flags <- c(flags, paste0("AUC_EXTRAP_GT_", hp_short(s$lz$max_extrap)))
    if (p$occasion > 1L) flags <- c(flags, "INF_AFTER_EARLIER_DOSES")
    if (!omit) {
      m_obs <- aumc_last + tlast * clast / l + clast / (l * l)
      m_pred <- aumc_last + tlast * clast_pred / l + clast_pred / (l * l)
      set("AUCINF_OBS", inf_obs); set("AUCINF_PRED", inf_pred); set("AUCINF_D", inf_obs / dose)
      set("AUMCINF_OBS", m_obs); set("AUMCINF_PRED", m_pred)
      mrt_obs <- m_obs / inf_obs - mrt_corr; mrt_pred <- m_pred / inf_pred - mrt_corr
      set("MRTINF_OBS", mrt_obs); set("MRTINF_PRED", mrt_pred)
      cl_obs <- dose_eff / inf_obs; cl_pred <- dose_eff / inf_pred
      vz_obs <- dose_eff / (l * inf_obs); vz_pred <- dose_eff / (l * inf_pred)
      if (iv) {
        set("CL_OBS", cl_obs); set("CL_PRED", cl_pred); set("VZ_OBS", vz_obs); set("VZ_PRED", vz_pred)
        set("VSS_OBS", cl_obs * mrt_obs); set("VSS_PRED", cl_pred * mrt_pred)
      } else {
        set("CL_F_OBS", cl_obs); set("CL_F_PRED", cl_pred); set("VZ_F_OBS", vz_obs); set("VZ_F_PRED", vz_pred)
      }
    } else flags <- c(flags, "INF_OMITTED")
  }
  if (!is.null(tau)) {
    tp <- cp
    if (!is.null(p$boundary) && abs(p$boundary$t - tau) <= 1e-9 * (1 + abs(tau)))
      tp <- rbind(tp, data.frame(t = tau, c = p$boundary$c))
    tp <- tp[c(TRUE, diff(tp$t) > 1e-9), , drop = FALSE]
    last <- if (nrow(tp)) tp[nrow(tp), ] else data.frame(t = 0, c = 0)
    ar <- function(a, b, mo) hp_curve_area(tp$t, tp$c, mcode, curve_tmax, a, b, mo)
    auc_tau <- aumc_tau <- c_tau <- NA_real_
    if (last$t >= tau - 1e-9) {
      auc_tau <- ar(t_start, tau, FALSE); aumc_tau <- ar(t_start, tau, TRUE)
      c_tau <- hp_curve_conc_at(tp$t, tp$c, mcode, curve_tmax, tau)
    } else if (last$c <= 0) {
      auc_tau <- ar(t_start, last$t, FALSE); aumc_tau <- ar(t_start, last$t, TRUE); c_tau <- 0
    } else if (!is.null(lam)) {
      tl <- hp_tail(last$t, last$c, tau, lam)
      auc_tau <- ar(t_start, last$t, FALSE) + tl[1L]; aumc_tau <- ar(t_start, last$t, TRUE) + tl[2L]
      c_tau <- last$c * exp(-lam * (tau - last$t))
    } else flags <- c(flags, "TAU_NO_EXTRAP")
    if (is.finite(auc_tau)) {
      set("AUCTAU", auc_tau); set("AUMCTAU", aumc_tau); set("CTAU", c_tau)
      cavg <- auc_tau / tau; set("CAVG", cavg)
      cmin <- par[["CMIN"]]
      if (!is.null(cmin) && is.finite(cmin) && cavg > 0) set("FLUCT_PCT", 100 * (cmax - cmin) / cavg)
      if (!is.null(cmin) && is.finite(cmin) && cmin > 0) { set("SWING", (cmax - cmin) / cmin); set("PTR", cmax / cmin) }
      if (ss) {
        cl <- dose_eff / auc_tau
        set(if (iv) "CL_SS" else "CL_F_SS", cl)
        if (!is.null(lam)) {
          set(if (iv) "VZ_SS" else "VZ_F_SS", dose_eff / (lam * auc_tau))
          mrt <- (aumc_tau + tau * c_tau / lam) / auc_tau - mrt_corr
          set("MRT_SS", mrt)
          if (iv) set("VSS_SS", cl * mrt)
          set("RACC_PRED", 1 / (1 - exp(-lam * tau)))
        }
      } else if (!is.null(lam)) set("RACC_PRED", 1 / (1 - exp(-lam * tau)))
    }
  }
  for (ab in s$partial) {
    a <- ab[1L]; b <- ab[2L]; v <- NA_real_
    if (a >= t_start - 1e-9) {
      last_t <- cp$t[nrow(cp)]; last_c <- cp$c[nrow(cp)]
      fa <- function(x, y) hp_curve_area(cp$t, cp$c, mcode, curve_tmax, x, y, FALSE)
      if (b <= last_t + 1e-9) v <- fa(a, min(b, last_t))
      else if (last_c <= 0) v <- fa(min(a, last_t), last_t)
      else if (s$partial_extrapolate && !is.null(lam)) {
        base <- if (a < last_t) fa(a, last_t) else 0
        from <- max(a, last_t); c_from <- last_c * exp(-lam * (from - last_t))
        v <- base + hp_tail(from, c_from, b, lam)[1L]
      }
    }
    set(hp_partial_name(a, b), v)
  }
  res$params <- par; res$lz <- lz; res$flags <- flags
  res
}

# ---------------------------------------------------------------------------
# summaries
# ---------------------------------------------------------------------------

hp_nca_accumulation <- function(results) {
  first <- Filter(function(r) is.null(r$skip) && !r$steady_state && r$occasion == 1L, results)
  for (i in seq_along(results)) {
    r <- results[[i]]
    if (!is.null(r$skip) || !r$steady_state) next
    f <- Filter(function(z) identical(z$id, r$id) && identical(z$route, r$route), first)
    if (!length(f)) next
    f <- f[[1L]]
    a1 <- f$params[["AUCTAU"]] %||% f$params[["AUCINF_OBS"]] %||% NA_real_
    c1 <- f$params[["CMAX"]] %||% NA_real_
    a <- r$params[["AUCTAU"]] %||% NA_real_; cm <- r$params[["CMAX"]] %||% NA_real_
    if (is.finite(a) && is.finite(a1) && a1 > 0) results[[i]]$params[["RACC_AUC"]] <- a / a1
    if (is.finite(cm) && is.finite(c1) && c1 > 0) results[[i]]$params[["RACC_CMAX"]] <- cm / c1
  }
  results
}

#' Table of the parameters of all profiles
#'
#' Collects the parameters of the analysed profiles into one data frame with one
#' row per profile and one column per parameter, ready to be saved, listed or
#' summarised.
#'
#' @param results A list of results from [hp_nca_analyse()], for example the
#'   `results` element returned by [hp_nca()].
#' @return A data frame with the columns `ID`, `GROUP`, `OCCASION`, `ROUTE`,
#'   `STEADY_STATE` (0 or 1), then every parameter in the order of
#'   [hp_nca_parameters()] (followed by partial areas; parameters that do not
#'   apply are `NA`), then `FLAGS` (the flags of the profile, separated by
#'   spaces) and `SKIPPED` (the reason a profile was not analysed, otherwise
#'   empty).
#' @seealso [hp_nca()], [hp_nca_summary()], [hp_save_table()]
#' @family non-compartmental analysis
#' @examples
#' r <- hp_nca(hp_read_data(hp_example_file()), route = "oral", lloq = 0.5)
#' tab <- hp_nca_table(r$results)
#' identical(tab, r$parameters)
#' tab[1:5, c("ID", "GROUP", "CMAX", "TMAX", "AUCLAST", "THALF", "FLAGS")]
#' @export
hp_nca_table <- function(results) {
  reg <- hp_nca_parameters()$NAME
  extra <- unique(unlist(lapply(results, function(r) names(r$params))))
  cols <- c(reg, setdiff(extra, reg))
  n <- length(results)
  m <- matrix(NA_real_, n, length(cols), dimnames = list(NULL, cols))
  for (i in seq_len(n)) {
    r <- results[[i]]
    for (nm in names(r$params)) { v <- r$params[[nm]]; if (!is.null(v)) m[i, nm] <- as.numeric(v) }
    if ("DOSE" %in% cols) m[i, "DOSE"] <- r$dose
  }
  out <- data.frame(ID = vapply(results, function(r) as.character(r$id), ""),
                    GROUP = vapply(results, function(r) as.character(r$group), ""),
                    OCCASION = vapply(results, function(r) r$occasion, 1L),
                    ROUTE = vapply(results, function(r) r$route, ""),
                    STEADY_STATE = vapply(results, function(r) as.integer(r$steady_state), 1L),
                    stringsAsFactors = FALSE)
  for (nm in cols) out[[nm]] <- m[, nm]
  out$FLAGS <- vapply(results, function(r) paste(r$flags, collapse = " "), "")
  out$SKIPPED <- vapply(results, function(r) r$skip %||% "", "")
  out
}

#' Every terminal regression that was tried
#'
#' Documents how the terminal rate constant was chosen. For every profile it
#' lists each candidate regression (number of points, time window, lambda z,
#' half-life, R-squared, adjusted R-squared, span in half-lives) and marks the
#' chosen one. This is the table to attach to a report when a reviewer asks why
#' a particular half-life was reported.
#'
#' @param results A list from [hp_nca_analyse()], for example the `results`
#'   element returned by [hp_nca()].
#' @return A data frame with the columns `ID`, `GROUP`, `OCCASION`, `ROUTE`,
#'   `CHOSEN` (1 for the chosen regression), `POINTS`, `START`, `END`,
#'   `LAMBDA_Z`, `HALF_LIFE`, `R2`, `R2ADJ`, `SPAN_RATIO`, `INTERCEPT` and `HOW`
#'   (how the choice was made). A profile with no regression has one row with
#'   empty values and the reason in `HOW`.
#' @seealso [hp_nca()], [hp_nca_options()]
#' @family non-compartmental analysis
#' @examples
#' r <- hp_nca(hp_read_data(hp_example_file()), route = "oral", lloq = 0.5)
#' lz <- hp_nca_lambda_z_table(r$results)
#' lz[lz$ID == lz$ID[1], c("ID", "CHOSEN", "POINTS", "START", "END", "R2ADJ", "HALF_LIFE")]
#' @export
hp_nca_lambda_z_table <- function(results) {
  rows <- list()
  for (r in results) {
    if (!is.null(r$skip)) next
    base <- function(chosen) data.frame(ID = r$id, GROUP = r$group, OCCASION = r$occasion, ROUTE = r$route,
                                        CHOSEN = chosen, stringsAsFactors = FALSE)
    cd <- r$lz$candidates
    if (is.null(cd) || !nrow(cd)) {
      rows[[length(rows) + 1L]] <- cbind(base(0L), POINTS = NA, START = NA, END = NA, LAMBDA_Z = NA, HALF_LIFE = NA,
                                         R2 = NA, R2ADJ = NA, SPAN_RATIO = NA, INTERCEPT = NA, HOW = r$lz$how)
      next
    }
    rows[[length(rows) + 1L]] <- data.frame(
      ID = r$id, GROUP = r$group, OCCASION = r$occasion, ROUTE = r$route, CHOSEN = as.integer(cd$selected),
      POINTS = cd$n, START = cd$start, END = cd$end, LAMBDA_Z = cd$lambda, HALF_LIFE = cd$half_life,
      R2 = cd$r2, R2ADJ = cd$adj_r2, SPAN_RATIO = cd$span_ratio, INTERCEPT = cd$intercept,
      HOW = ifelse(cd$selected, r$lz$how, ""), stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

#' Every point of every profile and what was done with it
#'
#' Lists each sample of each profile with the concentration as measured, the
#' concentration used in the analysis after the BLQ rule, the status explaining
#' the difference and whether the sample entered the terminal regression. It
#' also lists the zero or C0 added at the dose. It is the audit trail of the
#' analysis.
#'
#' @param results A list from [hp_nca_analyse()], for example the `results`
#'   element returned by [hp_nca()].
#' @return A data frame with the columns `ID`, `GROUP`, `OCCASION`, `ROUTE`,
#'   `DOSE`, `TIME`, `CONC` (as measured), `CONC_USED`, `STATUS` (`observed`,
#'   `blq set to zero`, `blq set to lloq`, `blq set to lloq/2`,
#'   `blq left out`, `negative left out`, `next-dose sample (end of interval)`,
#'   `added zero at dose`, `added C0` or `zero at dose replaced by C0`),
#'   `IN_LAMBDA_Z` (1
#'   when used in the terminal regression) and `LINE` (the row of the
#'   dataset, counting the header as row 1).
#' @seealso [hp_nca()], [hp_nca_plots()]
#' @family non-compartmental analysis
#' @examples
#' r <- hp_nca(hp_read_data(hp_example_file()), route = "oral", lloq = 0.5)
#' pt <- hp_nca_profile_table(r$results)
#' head(pt[pt$ID == pt$ID[1], ], 8)
#' table(pt$STATUS)
#' @export
hp_nca_profile_table <- function(results) {
  rows <- lapply(results, function(r) {
    pt <- r$points
    if (is.null(pt) || !nrow(pt)) return(NULL)
    lzt <- r$lz$times
    in_lz <- vapply(seq_len(nrow(pt)), function(i) as.integer(length(lzt) > 0 && any(abs(lzt - pt$t[i]) < 1e-9) && hp_status_sample(pt$status[i])), 1L)
    data.frame(ID = r$id, GROUP = r$group, OCCASION = r$occasion, ROUTE = r$route, DOSE = r$dose, TIME = pt$t,
               CONC = pt$raw, CONC_USED = pt$c, STATUS = unname(hp_status_label[pt$status]), IN_LAMBDA_Z = in_lz,
               LINE = ifelse(pt$row > 0, pt$row, NA), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' Descriptive statistics of NCA parameters by stratum
#'
#' Summarises the parameters of the profiles, by group, route, dose, occasion
#' and steady state, with the statistics usual in a pharmacokinetic report:
#' number, mean, SD, CV, geometric mean, geometric CV, median, minimum and
#' maximum.
#'
#' @details
#' Only the stratifying variables in `by` that take more than one value are
#' used, so a single-group single-dose study gets one stratum named `all`;
#' otherwise a stratum is named after its levels, for example
#' `DOSEGRP=100, dose 100`. Skipped profiles are left out. Geometric statistics
#' use the positive values only. Bookkeeping parameters (dose, number of
#' samples, terminal-regression window) are not summarised unless you name them
#' in `parameters`.
#'
#' @param table The parameter table from [hp_nca_table()].
#' @param by Any of `"group"`, `"route"`, `"dose"`, `"steady_state"` and
#'   `"occasion"`.
#' @param statistics Statistics to report, among `n`, `mean`, `sd`, `cv`,
#'   `sem`, `geomean`, `geocv`, `median`, `min`, `q1`, `q3`, `max` and
#'   `harmonic` (harmonic mean, used for half-lives). The default is `n`,
#'   `mean`, `sd`, `cv`, `geomean`, `geocv`, `median`, `min` and `max`.
#' @param parameters Parameters to summarise; `NULL` summarises all that have a
#'   value.
#' @param options Settings from [hp_nca_options()], for the units.
#' @return A data frame with the columns `STRATUM`, `PARAMETER`, `UNIT`, then
#'   one column per requested statistic (`N`, `MEAN`, `SD`, `CV_PCT`, `SEM`,
#'   `GEOMEAN`, `GEOCV_PCT`, `MEDIAN`, `MIN`, `Q1`, `Q3`, `MAX`,
#'   `HARMONIC_MEAN`). An empty data frame when nothing can be summarised.
#' @seealso [hp_nca()], [hp_nca_table()], [hp_describe()]
#' @family non-compartmental analysis
#' @examples
#' r <- hp_nca(hp_read_data(hp_example_file()), route = "oral", lloq = 0.5)
#'
#' # the summary that hp_nca() already made
#' head(r$summary[, c("STRATUM", "PARAMETER", "UNIT", "N", "GEOMEAN", "GEOCV_PCT")])
#'
#' # a different one: only Cmax and AUC, by dose, with the median and range
#' s <- hp_nca_summary(r$parameters, by = "dose",
#'                     statistics = c("n", "median", "min", "max"),
#'                     parameters = c("CMAX", "AUCLAST"), options = r$options)
#' s
#' @export
hp_nca_summary <- function(table, by = c("group", "route", "dose", "steady_state", "occasion"),
                           statistics = c("n", "mean", "sd", "cv", "geomean", "geocv", "median", "min", "max"),
                           parameters = NULL, options = hp_nca_options()) {
  t <- table[table$SKIPPED == "", , drop = FALSE]
  if (!nrow(t)) return(data.frame())
  level <- function(v) switch(v, group = t$GROUP, route = t$ROUTE, dose = paste("dose", vapply(t$DOSE, hp_short, "")),
                              occasion = paste("occasion", t$OCCASION),
                              steady_state = ifelse(t$STEADY_STATE == 1L, "steady state", "single dose"))
  vars <- Filter(function(v) length(unique(level(v))) > 1L, by)
  key <- lapply(vars, function(v) switch(v, dose = t$DOSE, occasion = t$OCCASION, steady_state = t$STEADY_STATE, level(v)))
  label <- if (!length(vars)) rep("all", nrow(t)) else do.call(paste, c(lapply(vars, level), sep = ", "))
  ord <- if (length(key)) do.call(order, key) else seq_len(nrow(t))
  strata <- unique(label[ord])
  pnames <- setdiff(names(t), c("ID", "GROUP", "OCCASION", "ROUTE", "STEADY_STATE", "FLAGS", "SKIPPED"))
  book <- c("DOSE", "N_SAMPLES", "LAMZ_N", "LAMZ_START", "LAMZ_END")
  hdr <- c(n = "N", mean = "MEAN", sd = "SD", cv = "CV_PCT", sem = "SEM", geomean = "GEOMEAN", geocv = "GEOCV_PCT",
           median = "MEDIAN", min = "MIN", q1 = "Q1", q3 = "Q3", max = "MAX", harmonic = "HARMONIC_MEAN")
  rows <- list()
  for (st in strata) for (n in pnames) {
    if (!is.null(parameters)) { if (!n %in% parameters) next } else if (n %in% book) next
    v <- t[[n]][label == st]
    if (!any(is.finite(v))) next
    d <- hp_describe_nca(v)
    rows[[length(rows) + 1L]] <- data.frame(STRATUM = st, PARAMETER = n, UNIT = hp_nca_unit(n, options),
                                            as.list(stats::setNames(d[statistics], hdr[statistics])),
                                            check.names = FALSE, stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

#' Non-compartmental analysis of a concentration-time dataset
#'
#' The main function of the NCA part of the package. It reads the dataset, cuts
#' every subject into dosing intervals, analyses every interval and summarises
#' the parameters, all with settings you control through [hp_nca_options()].
#' Oral, intravenous bolus and infusion dosing are supported, as are single doses,
#' repeated doses and steady state.
#'
#' @details
#' The analysis follows WinNonlin practice. For each profile it reports:
#'
#' - peak and trough: `CMAX`, `TMAX`, `CLAST`, `TLAST`, `TLAG` (extravascular
#'   single doses), `C0` (intravenous bolus) and, for dosing intervals,
#'   `CMIN` and `TMIN`;
#' - terminal phase: `LAMZ`, `THALF`, the number of points, the window and the
#'   R-squared and adjusted R-squared of the regression, and `CLAST_PRED`;
#' - areas: `AUCLAST`, `AUCINF_OBS`, `AUCINF_PRED`, the percent extrapolated
#'   (`AUC_PCT_EXTRAP_OBS`, `AUC_PCT_EXTRAP_PRED`), `AUMCLAST`, `AUMCINF_OBS`,
#'   `AUMCINF_PRED`, `MRTLAST`, `MRTINF_OBS`, `MRTINF_PRED` and partial areas
#'   such as `AUC_0_12`;
#' - clearance and volume: `CL_OBS`, `VZ_OBS`, `VSS_OBS` (and `_PRED`) for
#'   intravenous doses, `CL_F_OBS` and `VZ_F_OBS` (and `_PRED`) for
#'   extravascular doses;
#' - over the dosing interval: `AUCTAU`, `AUMCTAU`, `CTAU`, `CAVG`,
#'   `FLUCT_PCT`, `SWING`, `PTR`, and at steady state `CL_SS` (`CL_F_SS`),
#'   `VZ_SS` (`VZ_F_SS`), `MRT_SS`, `VSS_SS`, and the accumulation ratios
#'   `RACC_PRED`, `RACC_AUC` and `RACC_CMAX`;
#' - dose-normalised values `CMAX_D`, `AUCLAST_D` and `AUCINF_D`.
#'
#' Run [hp_nca_parameters()] for the full dictionary. Quality flags in the
#' `FLAGS` column mark profiles whose terminal phase is poor or whose extrapolated
#' area is large; see [hp_nca_analyse()].
#'
#' The areas, moments and terminal regressions are computed in C++, so a
#' study of thousands of profiles runs in seconds.
#'
#' State the `route` explicitly whenever you know it; `"auto"` is a
#' convenience that guesses from the shape of the data (see [hp_nca_read()]).
#'
#' @param data A data frame, usually from [hp_read_data()], with at least `ID`,
#'   `TIME`, `DV` and `AMT` columns.
#' @param options Settings from [hp_nca_options()]. The default is the
#'   standard settings with the `lloq` argument of this function.
#' @param route,groups,lloq,blq_flag,dv_log,filter,obs_filter,id,time,dv,amt,evid,mdv,cmt,rate,addl,ii,ss
#'   As in [hp_nca_read()]. Giving `blq_flag` also makes the flagged samples
#'   BLQ in the options.
#' @param summary_by,summary_statistics,summary_parameters As `by`,
#'   `statistics` and `parameters` of [hp_nca_summary()].
#' @return A list of class `hp_nca` with the elements
#'   - `parameters`: one row per profile, see [hp_nca_table()];
#'   - `summary`: descriptive statistics by stratum, see [hp_nca_summary()];
#'   - `lambda_z`: every terminal regression tried, see
#'     [hp_nca_lambda_z_table()];
#'   - `profiles`: every point and its status, see [hp_nca_profile_table()];
#'   - `results`: the list of per-profile results from [hp_nca_analyse()];
#'   - `notes`: messages about the data;
#'   - `options` and `route`: the settings used.
#'   Printing the object shows the number of profiles analysed and the notes.
#' @seealso [hp_nca_options()], [hp_nca_plots()], [hp_nca_parameters()],
#'   [hp_save_table()]
#' @family non-compartmental analysis
#' @examples
#' d <- hp_read_data(hp_example_file())
#'
#' # a standard analysis of an oral study; concentrations below 0.5 are BLQ
#' r <- hp_nca(d, route = "oral", lloq = 0.5)
#' r
#' head(r$parameters[, c("ID", "CMAX", "TMAX", "AUCLAST", "AUCINF_OBS", "THALF")])
#'
#' # by dose group, with the linear trapezoid rule and a partial area
#' opt <- hp_nca_options(auc_method = "linear", lloq = 0.5, partial = list(c(0, 12)))
#' r2 <- hp_nca(d, options = opt, route = "oral", groups = "DOSEGRP")
#' r2$summary[r2$summary$PARAMETER %in% c("CMAX", "AUC_0_12"),
#'            c("STRATUM", "PARAMETER", "N", "GEOMEAN", "GEOCV_PCT")]
#'
#' # profiles that were not analysed, and why
#' r$parameters$SKIPPED[r$parameters$SKIPPED != ""]
#' @export
hp_nca <- function(data, options = hp_nca_options(lloq = lloq), route = "auto", groups = NULL, lloq = NULL,
                   blq_flag = NULL, dv_log = FALSE, filter = NULL, obs_filter = NULL,
                   id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID", mdv = "MDV",
                   cmt = "CMT", rate = "RATE", addl = "ADDL", ii = "II", ss = "SS",
                   summary_by = c("group", "route", "dose", "steady_state", "occasion"),
                   summary_statistics = c("n", "mean", "sd", "cv", "geomean", "geocv", "median", "min", "max"),
                   summary_parameters = NULL) {
  if (is.null(options$lloq) && !is.null(lloq)) options$lloq <- lloq
  if (!is.null(blq_flag)) options$blq_flagged <- TRUE
  x <- hp_nca_read(data, route, id, time, dv, amt, evid, mdv, cmt, rate, addl, ii, ss, groups,
                   options$lloq, blq_flag, dv_log, filter, obs_filter)
  profiles <- hp_nca_profiles(x, options)
  if (!length(profiles)) hp_stop("the dataset has no dosing interval to analyse")
  results <- lapply(profiles, hp_nca_analyse, options = options)
  results <- hp_nca_accumulation(results)
  ok <- sum(vapply(results, function(r) is.null(r$skip), NA))
  if (ok == 0L) hp_stop("none of the ", length(results), " profile(s) can be analysed: ",
                        results[[1L]]$skip %||% "")
  tab <- hp_nca_table(results)
  out <- list(parameters = tab,
              summary = hp_nca_summary(tab, summary_by, summary_statistics, summary_parameters, options),
              lambda_z = hp_nca_lambda_z_table(results), profiles = hp_nca_profile_table(results),
              results = results, notes = c(x$notes, attr(profiles, "notes")), options = options, route = x$route)
  class(out) <- c("hp_nca", "list")
  out
}

#' @export
print.hp_nca <- function(x, ...) {
  ok <- sum(x$parameters$SKIPPED == "")
  cat("hyperutils non-compartmental analysis: ", ok, " profile(s) analysed of ", nrow(x$parameters), "\n", sep = "")
  for (n in x$notes) cat("  ", n, "\n", sep = "")
  invisible(x)
}
