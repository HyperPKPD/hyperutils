# Exploratory data analysis: the data frames behind the plots and tables.
#
# hp_eda_prepare() builds a record frame (every column of the dataset plus
# derived columns) and a subject frame (one row per individual) exactly as
# described in the documentation of the function.

hp_num_text <- function(v) {
  ifelse(!is.finite(v), "",
         ifelse(v == round(v) & abs(v) < 1e15, formatC(v, format = "f", digits = 0),
                sub("\\.?0+$", "", formatC(v, format = "f", digits = 6))))
}

# the usual names of the optional columns
hp_first_name <- function(cands, nms) { f <- cands[cands %in% nms]; if (length(f)) f[1L] else NULL }

#' Prepare a dataset for exploratory analysis
#'
#' The first step of every exploratory analysis. It takes a NONMEM-style
#' dataset (one row per dose or observation) and builds the two data frames
#' that the exploratory plots and tables work on: a record frame with derived
#' time and dose columns, and a subject frame with one row per individual that
#' holds baseline covariates, counts and the observed exposure. Doing this once
#' makes all later plots and tables fast and consistent.
#'
#' @details
#' **Record frame (`$rec`).** Every column of the dataset (names are put in
#' upper case) plus
#'
#' - `TAD`: time after the latest dose (the dataset's own `TAD` column is kept
#'   as is and the computed one is `TAD_CALC`);
#' - `TAFD`: time after the first dose;
#' - `DOSE`: the latest dose amount, and `DVNORM`, the dose-normalised `DV`;
#' - `BLQ`: 1 for an observation below the limit of quantification or flagged;
#' - `EVENT`: `observation`, `observation_mdv`, `dose` or `other`;
#' - `PKOBS`: 1 for an observation of the PK `DVID`.
#'
#' **Subject frame (`$subj`).** The first-record value of every column, the
#' baseline value of every covariate, `NOBS`, `NDOSE`, `NBLQ`, `DOSE1`,
#' `TOTDOSE` and the observed exposure `CMAX`, `TMAX`, `CMIN`, `CLAST`,
#' `TLAST`, `AUCLAST` (linear-up/log-down from the dose, a quick descriptive
#' value and not a replacement for [hp_nca()]) and the dose-normalised `CMAXD`
#' and `AUCLASTD`.
#'
#' **Covariate detection.** When `continuous` and `categorical` are not given
#' they are found from the data: a numeric column with up to six distinct
#' values, or a text column, is categorical and other numeric columns are
#' continuous. Covariates that change within an individual are listed as time
#' varying.
#'
#' **DVID.** Datasets with several kinds of observations (for example
#' concentrations and a response) identify them with a column such as `DVID`.
#' Plots of concentrations then use the `pk_dvid` observations only.
#'
#' @param data A data frame, usually from [hp_read_data()].
#' @param id,time,dv,amt,evid,mdv Names of the columns that hold the subject,
#'   time, observed value, dose amount, event type and missing-value flag.
#' @param dvid Column that tells the observation types apart. Found from
#'   `DVID`, `YTYPE`, `TYPE` and `DVTYPE` when `NULL`.
#' @param nominal_time Column of the nominal times (found from `NTIM`,
#'   `NTIME` and similar names when `NULL`).
#' @param tad Column with the time after dose, when the dataset has one.
#' @param addl,ii Columns of additional doses and the dosing interval.
#' @param lloq Limit of quantification; observations below it are BLQ.
#' @param blq_flag Column flagging BLQ observations (found from `BLQ`, `BQL`,
#'   `BLOQ` or `CENS` when `NULL` and no `lloq` is given).
#' @param dv_log `TRUE` when `DV` is the natural logarithm of the
#'   concentration; the limit of quantification is then compared with
#'   `exp(DV)`.
#' @param pk_dvid,pd_dvid `DVID` values of the concentrations and of the
#'   response (by default the smallest value and the next one).
#' @param color,strata Columns that colour and facet the standard plots.
#' @param continuous,categorical,time_varying Covariate names, see Details.
#' @param exposure_first Compute the observed exposure over the first dosing
#'   interval (otherwise over the whole profile).
#' @param filtered_out Number of records already removed by earlier filtering,
#'   reported in the data summary.
#' @param filter R expression of the rows to keep, evaluated on every row (see
#'   [hp_filter()]).
#' @param obs_filter R expression of the observation rows to keep; dose
#'   records always stay.
#' @return A list of class `hp_eda_data` with the elements `rec` (record
#'   frame), `subj` (subject frame), `info` (counts, covariate lists, DVID
#'   values and settings) and `names` (the resolved column names).
#' @seealso [hp_eda_plot()], [hp_eda_table()], [hp_eda_subjects()],
#'   [hp_eda_data_summary()]
#' @family exploratory data analysis
#' @examples
#' d <- hp_read_data(hp_example_file())
#' e <- hp_eda_prepare(d, dvid = "DVID", lloq = 0.5, color = "DOSEGRP")
#' class(e)
#' head(e$subj[, c("ID", "DOSE1", "CMAX", "AUCLAST", "WT")])
#' table(e$rec$EVENT)
#'
#' # restrict the analysis to part of the data
#' e2 <- hp_eda_prepare(d, dvid = "DVID", filter = "DOSEGRP == 100")
#' e2$info$individuals
#'
#' # the covariates found from the data
#' e$info$continuous
#' e$info$categorical
#' @export
hp_eda_prepare <- function(data, id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
                           mdv = "MDV", dvid = NULL, nominal_time = NULL, tad = NULL, addl = "ADDL",
                           ii = "II", lloq = NULL, blq_flag = NULL, dv_log = FALSE, pk_dvid = NULL,
                           pd_dvid = NULL, color = NULL, strata = NULL, continuous = NULL,
                           categorical = NULL, time_varying = NULL, exposure_first = TRUE,
                           filtered_out = 0, filter = NULL, obs_filter = NULL) {
  names(data) <- toupper(names(data))
  nm <- names(data)
  if (!is.null(filter) || !is.null(obs_filter)) {
    keep <- rep(TRUE, nrow(data))
    ev_ <- toupper(evid); am_ <- toupper(amt)
    zero <- function(col) { v <- data[[col]]; is.na(v) | as.character(v) %in% c("", ".") | (!is.na(suppressWarnings(as.numeric(as.character(v)))) & suppressWarnings(as.numeric(as.character(v))) == 0) }
    is_obs <- if (ev_ %in% nm) zero(ev_) else if (am_ %in% nm) zero(am_) else rep(TRUE, nrow(data))
    n_f <- n_o <- 0L
    if (!is.null(filter)) {
      k <- eval(parse(text = filter)[[1L]], envir = hp_eval_env(data)); k[is.na(k)] <- FALSE
      n_f <- sum(!k); keep <- keep & k
    }
    if (!is.null(obs_filter)) {
      k <- eval(parse(text = obs_filter)[[1L]], envir = hp_eval_env(data)); k[is.na(k)] <- FALSE
      drop <- keep & is_obs & !k; n_o <- sum(drop); keep <- keep & !drop
    }
    data <- data[keep, , drop = FALSE]; rownames(data) <- NULL
    filtered_out <- filtered_out + n_f + n_o
  }
  id <- toupper(id); time <- toupper(time); dv <- toupper(dv); amt <- toupper(amt)
  for (r in c(id, time, dv)) if (!r %in% nm) hp_stop("the ", r, " column is not in the dataset (the columns are: ", paste(nm, collapse = ", "), ")")
  dvid <- if (!is.null(dvid)) toupper(dvid) else hp_first_name(c("DVID", "YTYPE", "TYPE", "DVTYPE"), nm)
  if (is.null(blq_flag) && is.null(lloq)) blq_flag <- hp_first_name(c("BLQ", "BQL", "BLOQ", "CENS"), nm)
  nominal_time <- if (!is.null(nominal_time)) toupper(nominal_time) else hp_first_name(c("NTIM", "NTIME", "NOMTIME", "NOMT", "TNOM", "NT"), nm)
  n <- nrow(data)
  tm <- hp_num(data[[time]]); y <- hp_num(data[[dv]])
  a <- if (amt %in% nm) hp_num(data[[amt]]) else rep(NA_real_, n)
  ev <- if (toupper(evid) %in% nm) hp_num(data[[toupper(evid)]]) else ifelse(is.finite(a) & a > 0, 1, 0)
  ev[!is.finite(ev)] <- 0
  md <- if (toupper(mdv) %in% nm) hp_num(data[[toupper(mdv)]]) else ifelse(ev != 0, 1, 0)
  md[!is.finite(md)] <- 0
  kind <- ifelse(ev == 0 & md == 0, "obs", ifelse(ev == 0, "obsmdv", ifelse(ev %in% c(1, 4), "dose", "other")))
  idtxt <- as.character(data[[id]])
  idtxt[is.na(idtxt)] <- "(blank)"
  blk <- cumsum(c(TRUE, idtxt[-1L] != idtxt[-n]))
  nb <- max(blk)
  starts <- which(c(TRUE, diff(blk) != 0)); ends <- c(starts[-1L] - 1L, n)
  ad <- if (toupper(addl) %in% nm) hp_num(data[[toupper(addl)]]) else rep(NA_real_, n)
  ivl <- if (toupper(ii) %in% nm) hp_num(data[[toupper(ii)]]) else rep(NA_real_, n)
  # dose events of each block (the repeats of ADDL included), sorted by time then row
  doses <- vector("list", nb)
  for (b in seq_len(nb)) {
    rows <- starts[b]:ends[b]
    dr <- rows[kind[rows] == "dose" & is.finite(tm[rows])]
    tt <- numeric(); rr <- integer()
    for (r in dr) {
      tt <- c(tt, tm[r]); rr <- c(rr, r)
      if (is.finite(ad[r]) && ad[r] >= 1 && is.finite(ivl[r]) && ivl[r] > 0) {
        k <- seq_len(min(as.integer(ad[r]), 10000L)); tt <- c(tt, tm[r] + k * ivl[r]); rr <- c(rr, rep(r, length(k)))
      }
    }
    o <- order(tt, rr); doses[[b]] <- list(time = tt[o], row = rr[o])
  }
  dvid_txt <- if (!is.null(dvid)) as.character(data[[dvid]]) else NULL
  if (!is.null(dvid)) { v <- hp_num(data[[dvid]]); dvid_txt <- ifelse(is.finite(v), hp_num_text(v), dvid_txt); dvid_txt[is.na(dvid_txt)] <- "" }
  anyobs <- kind %in% c("obs", "obsmdv")
  dvids <- character()
  if (!is.null(dvid)) {
    dvids <- unique(dvid_txt[anyobs])
    nn <- suppressWarnings(as.numeric(dvids))
    dvids <- if (all(!is.na(nn))) dvids[order(nn)] else sort(dvids)
  }
  pk_dvid <- if (!is.null(pk_dvid)) as.character(pk_dvid) else if (length(dvids)) dvids[1L] else NULL
  pd_dvid <- if (!is.null(pd_dvid)) as.character(pd_dvid) else { o <- setdiff(dvids, pk_dvid); if (length(o)) o[1L] else NULL }
  is_pk <- if (!is.null(dvid) && !is.null(pk_dvid)) dvid_txt == pk_dvid else rep(TRUE, n)
  tad_v <- rep(NA_real_, n); tafd <- rep(NA_real_, n); dose <- rep(NA_real_, n); blq <- rep(NA_real_, n)
  event <- rep("", n); pkobs <- rep(0, n)
  lloq_v <- rep(lloq %||% NA_real_, n)
  t0 <- rep(NA_real_, nb)
  for (b in seq_len(nb)) {
    d <- doses[[b]]; t0[b] <- if (length(d$time)) d$time[1L] else NA_real_
    for (r in starts[b]:ends[b]) {
      t <- tm[r]
      tafd[r] <- if (is.finite(t0[b])) t - t0[b] else t
      if (length(d$time)) { k <- which(d$time <= t + 1e-9); if (length(k)) dose[r] <- a[d$row[max(k)]] }
      event[r] <- switch(kind[r], obs = "observation", obsmdv = "observation_mdv", dose = "dose", "other")
      if (anyobs[r]) {
        if (length(d$time)) {
          k <- sum(d$time < t | (d$time == t & d$row < r))
          if (k > 0L) tad_v[r] <- t - d$time[k]
        }
        flag <- if (!is.null(blq_flag) && toupper(blq_flag) %in% nm) { v <- hp_num(data[[toupper(blq_flag)]][r]); is.finite(v) && v != 0 } else FALSE
        below <- FALSE
        if (is.finite(lloq_v[r]) && is.finite(y[r])) { conc <- if (dv_log) exp(y[r]) else y[r]; below <- conc < lloq_v[r] * (1 - 1e-6) }
        blq[r] <- as.numeric(below || flag)
        if (kind[r] == "obs" && is_pk[r]) pkobs[r] <- 1
      } else if (kind[r] == "dose") tad_v[r] <- 0
    }
  }
  dvnorm <- ifelse(is.finite(dose) & dose > 0, y / dose, NA_real_)
  rec <- data
  if (!is.null(tad) && toupper(tad) %in% nm) rec$TAD_CALC <- tad_v else if ("TAD" %in% nm) rec$TAD_CALC <- tad_v else rec$TAD <- tad_v
  rec$TAFD <- tafd
  if (!"DOSE" %in% nm) rec$DOSE <- dose
  rec$DVNORM <- dvnorm
  if (!"BLQ" %in% nm) rec$BLQ <- blq
  rec$EVENT <- event
  rec$PKOBS <- pkobs
  blq_use <- blq
  # covariates
  std <- c(id, time, dv, amt, "EVID", "MDV", "CMT", "RATE", "ADDL", "II", "SS", "DVID", "YTYPE", "TYPE", "DVTYPE",
           "TAD", "TAD_CALC", "TAFD", "DOSE", "DVNORM", "BLQ", "EVENT", "PKOBS", "LLOQ", "NTIM", "NTIME", "BQL", "BLOQ",
           "CENS", "OCC", "FLAG", "PCMT", "MDV", "ADDL", "II", "SS", toupper(c(dvid, nominal_time, tad, blq_flag)),
           "TALD", "TSLD", "TSFD", "NTAD", "NOMTIME", "NOMT", "TNOM", "NT", "DAY", "VISIT", "LNDV", "LDV", "LOGDV",
           "LNCONC", "LOQ", "YTYPE", "TYPE", "DVTYPE", "DOSEGRP", "GROUP", "GRP", "ARM", "TRT", "TRTA", "COHORT",
           "STUDY", "PERIOD", "PART")
  cand <- setdiff(names(data), std)
  base_txt <- function(col) vapply(seq_len(nb), function(b) {
    v <- data[[col]][starts[b]:ends[b]]
    vn <- hp_num(v)
    f <- which(is.finite(vn) & vn != -99)
    if (length(f)) return(hp_num_text(vn[f[1L]]))
    tx <- as.character(v); tx <- tx[!is.na(tx) & nzchar(tx) & is.na(suppressWarnings(as.numeric(tx)))]
    if (length(tx)) tx[1L] else ""
  }, "")
  if (is.null(continuous) && is.null(categorical)) {
    continuous <- character(); categorical <- character(); time_varying <- time_varying %||% character()
    for (cv in cand) {
      v <- data[[cv]]
      if (all(is.na(v))) next
      varies <- any(tapply(as.character(v), blk, function(z) length(unique(z)) > 1L))
      bt <- base_txt(cv); bn <- suppressWarnings(as.numeric(bt))
      isnum <- is.numeric(v) && all(bt == "" | !is.na(bn))
      dist <- sort(unique(bn[is.finite(bn)]))
      cat_ <- !isnum || length(dist) == 2L || (all(dist == round(dist)) && length(dist) <= 5L)
      if (cat_) categorical <- c(categorical, cv) else continuous <- c(continuous, cv)
      if (varies && !cv %in% time_varying) time_varying <- c(time_varying, cv)
    }
  }
  continuous <- toupper(continuous %||% character()); categorical <- toupper(categorical %||% character())
  time_varying <- toupper(time_varying %||% character())
  # the subject frame: the first record of every individual
  subj <- rec[starts, , drop = FALSE]
  for (cv in unique(c(continuous, categorical, time_varying))) {
    if (!cv %in% names(data)) next
    bt <- base_txt(cv)
    subj[[cv]] <- if (all(bt == "" | !is.na(suppressWarnings(as.numeric(bt))))) suppressWarnings(as.numeric(bt)) else bt
  }
  nobs <- ndose <- nblq <- dose1 <- totdose <- numeric(nb)
  ex <- matrix(NA_real_, nb, 7, dimnames = list(NULL, c("cmax", "tmax", "cmin", "clast", "tlast", "auclast", "nobs")))
  for (b in seq_len(nb)) {
    d <- doses[[b]]; rows <- starts[b]:ends[b]
    dose1[b] <- if (length(d$time)) a[d$row[1L]] else NA_real_
    totdose[b] <- sum(a[d$row][is.finite(a[d$row])]); ndose[b] <- length(d$time)
    t_end <- if (exposure_first) (if (length(d$time) >= 2L) d$time[2L] else Inf) else Inf
    nobs[b] <- sum(kind[rows] == "obs"); nblq[b] <- sum(blq[rows] == 1, na.rm = TRUE)
    use <- rows[pkobs[rows] == 1 & !(blq[rows] %in% 1) & tm[rows] < t_end - 1e-9]
    tt <- if (is.finite(t0[b])) tm[use] - t0[b] else tm[use]
    ok <- tt >= -1e-9
    ex[b, ] <- if (any(ok)) hp_exposure(tt[ok], y[use][ok]) else hp_exposure(numeric(), numeric())
  }
  subj$NOBS <- nobs; subj$NDOSE <- ndose; subj$NBLQ <- nblq; subj$DOSE1 <- dose1; subj$TOTDOSE <- totdose
  subj$CMAX <- ex[, "cmax"]; subj$TMAX <- ex[, "tmax"]; subj$CMIN <- ex[, "cmin"]; subj$CLAST <- ex[, "clast"]
  subj$TLAST <- ex[, "tlast"]; subj$AUCLAST <- ex[, "auclast"]
  nrm <- function(v) ifelse(is.finite(dose1) & dose1 > 0, v / dose1, NA_real_)
  subj$CMAXD <- nrm(ex[, "cmax"]); subj$AUCLASTD <- nrm(ex[, "auclast"])
  rownames(subj) <- NULL
  lev <- sort(unique(dose1[is.finite(dose1)]))
  varying <- names(rec)[vapply(names(rec), function(cn) {
    v <- as.character(rec[[cn]]); any(tapply(v, blk, function(z) any(z != z[1L], na.rm = TRUE)))
  }, NA)]
  info <- list(records = n, individuals = nb, observations = sum(kind == "obs"), pk_observations = sum(pkobs == 1),
               doses = sum(kind == "dose"), blq = sum(blq == 1, na.rm = TRUE), dvid_values = dvids, dose_levels = lev,
               continuous = continuous, categorical = categorical, time_varying = time_varying,
               color = color, strata = strata, pk_dvid = pk_dvid, pd_dvid = pd_dvid, filtered_out = filtered_out)
  out <- list(rec = rec, subj = subj, info = info, varying = varying,
              names = list(id = id, time = time, dv = dv, amt = amt, dvid = dvid, nominal = nominal_time))
  class(out) <- c("hp_eda_data", "list")
  out
}

#' @export
print.hp_eda_data <- function(x, ...) {
  i <- x$info
  cat(sprintf("hyperutils exploratory data: %d individuals, %d records (%d observations, %d PK, %d doses, %d BLQ)\n",
              i$individuals, i$records, i$observations, i$pk_observations, i$doses, i$blq))
  if (length(i$continuous)) cat("  continuous covariates:", paste(i$continuous, collapse = ", "), "\n")
  if (length(i$categorical)) cat("  categorical covariates:", paste(i$categorical, collapse = ", "), "\n")
  invisible(x)
}

#' One row per individual
#'
#' Extracts the subject frame of a prepared dataset without its internal helper
#' columns (`EVENT`, `PKOBS`, `DVNORM`, `TAD_CALC`, `TAD`, `TAFD` and `BLQ`).
#' It is the table to start from when you want to model the exposure or the
#' covariates of the individuals, or to list the study population.
#'
#' @param d An `hp_eda_data` object from [hp_eda_prepare()].
#' @return A data frame with one row per individual: baseline covariates,
#'   counts (`NOBS`, `NDOSE`, `NBLQ`), doses and the observed exposure
#'   (`CMAX`, `TMAX`, `AUCLAST` and others).
#' @seealso [hp_eda_prepare()]
#' @family exploratory data analysis
#' @examples
#' e <- hp_eda_prepare(hp_read_data(hp_example_file()), dvid = "DVID")
#' s <- hp_eda_subjects(e)
#' dim(s)
#' # exposure by dose group
#' aggregate(CMAX ~ DOSEGRP, data = s, FUN = mean)
#' @export
hp_eda_subjects <- function(d) {
  skip <- c("EVENT", "PKOBS", "DVNORM", "TAD_CALC", "TAD", "TAFD", "BLQ")
  d$subj[, setdiff(names(d$subj), skip), drop = FALSE]
}
