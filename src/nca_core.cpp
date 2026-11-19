// hyperutils: C++ core of the non-compartmental analysis and of the exploratory
// exposure measures. The numerical rules are those of common NCA software:
// trapezoids and log-trapezoids written so that near-equal concentrations
// lose no accuracy, and the (weighted) straight-line regression of the
// terminal phase.

#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

// ln(c1 / c2) of two positive concentrations through log1p of the relative
// difference
static inline double log_ratio(double c1, double c2) {
  double d = (c1 - c2) / c2;
  if (std::isfinite(d)) return std::log1p(d);
  return std::log(c1) - std::log(c2);
}

// (e^l - 1 - l) / l^2, smooth at l = 0 (where it is 1/2)
static inline double exp_second_order(double l) {
  if (std::fabs(l) < 0.1) {
    double term = 0.5, sum = 0.5;
    for (int k = 1; k < 14; ++k) {
      term *= l / (k + 2.0);
      sum += term;
    }
    return sum;
  }
  return (std::expm1(l) - l) / (l * l);
}

static inline double seg_auc(double t1, double c1, double t2, double c2, bool lg) {
  double dt = t2 - t1;
  if (lg && c1 > 0.0 && c2 > 0.0) {
    double l = log_ratio(c1, c2);
    if (l == 0.0) return c2 * dt;
    return (c1 - c2) * dt / l;
  }
  return (c1 + c2) / 2.0 * dt;
}

static inline double seg_aumc(double t1, double c1, double t2, double c2, bool lg) {
  double dt = t2 - t1;
  if (lg && c1 > 0.0 && c2 > 0.0) {
    double l = log_ratio(c1, c2);
    return t1 * seg_auc(t1, c1, t2, c2, true) + dt * dt * c2 * exp_second_order(l);
  }
  return (t1 * c1 + t2 * c2) / 2.0 * dt;
}

// method: 0 linear, 1 linear-log, 2 linear-up/log-down
static inline bool log_segment(int method, double t1, double c1, double c2, double tmax) {
  if (!(c1 > 0.0 && c2 > 0.0) || c1 == c2) return false;
  switch (method) {
    case 0: return false;
    case 1: return t1 >= tmax - 1e-12;
    default: return c2 < c1;
  }
}

//' Area of one trapezoid under a concentration-time curve
//'
//' Computes the area under the curve between two consecutive sampling points,
//' either with the linear trapezoidal rule or with the logarithmic trapezoidal
//' rule. This is the building block of every area calculation in the package
//' and is implemented in C++ for speed.
//'
//' @details
//' The linear rule gives `(c1 + c2) / 2 * (t2 - t1)`. The logarithmic rule is
//' the exact integral of the exponential curve that passes through the two
//' points, `(c1 - c2) * (t2 - t1) / log(c1 / c2)`. It is used only when both
//' concentrations are positive; otherwise (a zero or a negative value) the
//' function falls back to the linear rule, as pharmacokinetic software does.
//' When the two concentrations are equal the logarithmic rule reduces to the
//' rectangle `c * (t2 - t1)`.
//'
//' The logarithm of the ratio is evaluated through `log1p` of the relative
//' difference, so nearly equal concentrations lose no precision.
//'
//' @param t1,c1 Time and concentration of the first point.
//' @param t2,c2 Time and concentration of the second point (`t2 > t1`).
//' @param log_rule `TRUE` for the logarithmic trapezoid, `FALSE` for the linear
//'   trapezoid.
//' @return A single number, the area of the segment, in concentration times
//'   time units.
//' @seealso [hp_seg_aumc()] for the first-moment area of a segment,
//'   [hp_curve_area()] for a whole curve.
//' @family area functions
//' @examples
//' # linear trapezoid: (2 + 4) / 2 * 1
//' hp_seg_auc(0, 2, 1, 4, log_rule = FALSE)
//'
//' # a monoexponential decline C = 10 * exp(-0.5 * t) between t = 0 and t = 2:
//' # the log trapezoid is exact, the linear trapezoid overestimates
//' c2 <- 10 * exp(-0.5 * 2)
//' exact <- (10 - c2) / 0.5
//' c(exact = exact,
//'   log = hp_seg_auc(0, 10, 2, c2, TRUE),
//'   linear = hp_seg_auc(0, 10, 2, c2, FALSE))
//' @export
// [[Rcpp::export]]
double hp_seg_auc(double t1, double c1, double t2, double c2, bool log_rule) {
  return seg_auc(t1, c1, t2, c2, log_rule);
}

//' Area of one trapezoid under the first-moment curve (time times concentration)
//'
//' The first-moment curve is `t * C(t)`. Its area is needed for the mean
//' residence time (AUMC divided by AUC). This function gives the area of one
//' segment between two sampling points with either the linear or the
//' logarithmic rule.
//'
//' @details
//' For the linear rule the result is `(t1 * c1 + t2 * c2) / 2 * (t2 - t1)`. For
//' the logarithmic rule, with both concentrations positive, the result is the
//' exact integral of `t * C(t)` over an exponential segment. A zero or negative
//' concentration falls back to the linear rule.
//'
//' @inheritParams hp_seg_auc
//' @return A single number, the first-moment area of the segment, in
//'   concentration times time squared.
//' @seealso [hp_seg_auc()], [hp_curve_area()] (with `moment = TRUE`).
//' @family area functions
//' @examples
//' hp_seg_aumc(0, 2, 1, 4, log_rule = FALSE)
//'
//' # area of t * C(t) over an exponential decline, exact versus linear
//' c2 <- 10 * exp(-0.5 * 2)
//' c(log = hp_seg_aumc(0, 10, 2, c2, TRUE),
//'   linear = hp_seg_aumc(0, 10, 2, c2, FALSE))
//' @export
// [[Rcpp::export]]
double hp_seg_aumc(double t1, double c1, double t2, double c2, bool log_rule) {
  return seg_aumc(t1, c1, t2, c2, log_rule);
}

//' Area under a whole concentration-time curve between two times
//'
//' Integrates a profile, point by point, from time `a` to time `b` with the
//' linear, linear-log or linear-up/log-down trapezoidal method. The limits do
//' not have to be sampling times: the first and the last segment are cut at `a`
//' and `b` using the interpolation of the chosen method, which is what a
//' partial area such as AUC from 0 to 4 hours needs.
//'
//' @details
//' The three methods are the ones of WinNonlin and similar software:
//'
//' - `method = 0`, **linear**: linear trapezoids everywhere.
//' - `method = 1`, **linear-log**: linear trapezoids up to the time of the
//'   maximum `tmax` and logarithmic trapezoids after it.
//' - `method = 2`, **linear-up/log-down**: linear trapezoids while the
//'   concentration rises or stays equal and logarithmic trapezoids while it
//'   falls (both concentrations positive).
//'
//' Segments with a zero or negative concentration are always linear. The area
//' is limited to the span of the points: nothing is added before the first or
//' after the last time. Extrapolation beyond the last point is the job of the
//' terminal slope, see [hp_nca()].
//'
//' The function is written in C++ and is the routine behind every area of
//' [hp_nca()].
//'
//' @param t Numeric vector of times, strictly increasing.
//' @param c Numeric vector of concentrations, same length as `t`.
//' @param method Integer code of the rule: 0 linear, 1 linear-log,
//'   2 linear-up/log-down.
//' @param tmax Time of the maximum concentration. Only used by `method = 1`,
//'   where the logarithmic rule starts at this time.
//' @param a,b Limits of the area (`a < b`).
//' @param moment `TRUE` returns the area under the first-moment curve
//'   `t * C(t)` instead of the area under `C(t)`.
//' @return A single number: the area, or `0` when the limits do not overlap the
//'   curve.
//' @seealso [hp_seg_auc()], [hp_curve_conc_at()], [hp_nca_options()] for the
//'   methods by name.
//' @family area functions
//' @examples
//' t <- c(0, 1, 2, 4, 8, 12)
//' c <- c(0, 4, 5, 3, 1, 0.4)
//'
//' # AUC from 0 to the last time, three methods
//' sapply(0:2, function(m) hp_curve_area(t, c, m, tmax = 2, a = 0, b = 12))
//'
//' # a partial area, from 0 to 6 hours, cutting the segment between 4 and 8
//' hp_curve_area(t, c, method = 2, tmax = 2, a = 0, b = 6)
//'
//' # the area under the first-moment curve and the mean residence time to 12 h
//' aumc <- hp_curve_area(t, c, 2, 2, 0, 12, moment = TRUE)
//' auc  <- hp_curve_area(t, c, 2, 2, 0, 12)
//' aumc / auc
//' @export
// [[Rcpp::export]]
double hp_curve_area(NumericVector t, NumericVector c, int method, double tmax,
                     double a, double b, bool moment = false) {
  int n = t.size();
  double total = 0.0;
  for (int i = 0; i + 1 < n; ++i) {
    double t1 = t[i], c1 = c[i], t2 = t[i + 1], c2 = c[i + 1];
    double lo = std::max(a, t1), hi = std::min(b, t2);
    if (!(hi > lo)) continue;
    bool lg = log_segment(method, t1, c1, c2, tmax);
    auto at = [&](double x) {
      double f = (x - t1) / (t2 - t1);
      return lg ? c1 * std::pow(c2 / c1, f) : c1 + (c2 - c1) * f;
    };
    double clo = at(lo), chi = at(hi);
    total += moment ? seg_aumc(lo, clo, hi, chi, lg) : seg_auc(lo, clo, hi, chi, lg);
  }
  return total;
}

//' Concentration at a given time, interpolated with the trapezoid rule
//'
//' Returns the concentration of a profile at time `x` using the same
//' interpolation as the area rules of [hp_curve_area()]: a straight line
//' between two points for a linear segment and an exponential curve for a
//' logarithmic segment. It is used, for example, to read the concentration at
//' the end of a partial-area interval.
//'
//' @details
//' Times outside the sampled range (beyond a small tolerance) give `NA`. A time
//' that falls on a sampling time returns that observation.
//'
//' @inheritParams hp_curve_area
//' @param x The time at which the concentration is wanted.
//' @return A single number, the interpolated concentration, or `NA` outside
//'   the sampled range.
//' @seealso [hp_curve_area()]
//' @family area functions
//' @examples
//' t <- c(0, 1, 2, 4, 8)
//' c <- c(0, 4, 5, 3, 1)
//' hp_curve_conc_at(t, c, method = 0, tmax = 2, x = 3)   # straight line: 4
//' hp_curve_conc_at(t, c, method = 2, tmax = 2, x = 3)   # exponential: 3.87
//' hp_curve_conc_at(t, c, method = 2, tmax = 2, x = 20)  # outside: NA
//' @export
// [[Rcpp::export]]
double hp_curve_conc_at(NumericVector t, NumericVector c, int method, double tmax, double x) {
  int n = t.size();
  if (n == 0 || x < t[0] - 1e-12 || x > t[n - 1] + 1e-12) return NA_REAL;
  for (int i = 0; i + 1 < n; ++i) {
    double t1 = t[i], c1 = c[i], t2 = t[i + 1], c2 = c[i + 1];
    if (x >= t1 - 1e-12 && x <= t2 + 1e-12) {
      double f = t2 > t1 ? (x - t1) / (t2 - t1) : 0.0;
      return log_segment(method, t1, c1, c2, tmax) ? c1 * std::pow(c2 / c1, f)
                                                   : c1 + (c2 - c1) * f;
    }
  }
  return c[n - 1];
}

// weighted least-squares line; returns intercept, slope, r2, adj_r2, se_slope, df
static NumericVector fit_line(const std::vector<double>& x, const std::vector<double>& y,
                              const std::vector<double>* w) {
  NumericVector bad(6, NA_REAL);
  size_t n = x.size();
  if (n < 2 || y.size() != n) return bad;
  auto wt = [&](size_t i) { return w ? (*w)[i] : 1.0; };
  double sw = 0.0, xm = 0.0, ym = 0.0;
  for (size_t i = 0; i < n; ++i) { sw += wt(i); xm += wt(i) * x[i]; ym += wt(i) * y[i]; }
  if (!(sw > 0.0)) return bad;
  xm /= sw; ym /= sw;
  double sxx = 0.0, sxy = 0.0, syy = 0.0;
  for (size_t i = 0; i < n; ++i) {
    sxx += wt(i) * (x[i] - xm) * (x[i] - xm);
    sxy += wt(i) * (x[i] - xm) * (y[i] - ym);
    syy += wt(i) * (y[i] - ym) * (y[i] - ym);
  }
  if (!(sxx > 0.0)) return bad;
  double slope = sxy / sxx, intercept = ym - slope * xm, sse = 0.0;
  for (size_t i = 0; i < n; ++i) {
    double r = y[i] - intercept - slope * x[i];
    sse += wt(i) * r * r;
  }
  double r2 = syy > 0.0 ? 1.0 - sse / syy : 1.0;
  double adj = n >= 3 ? 1.0 - (1.0 - r2) * (n - 1.0) / (n - 2.0) : NA_REAL;
  double se = n >= 3 ? std::sqrt(sse / (n - 2.0) / sxx) : NA_REAL;
  NumericVector out(6);
  out[0] = intercept; out[1] = slope; out[2] = r2; out[3] = adj; out[4] = se;
  out[5] = (double)(n > 2 ? n - 2 : 0);
  return out;
}

//' Weighted straight-line regression
//'
//' Fits `y = intercept + slope * x` by (weighted) least squares and returns the
//' statistics that the terminal phase of a pharmacokinetic profile needs. The
//' function is a fast C++ alternative to `lm()` for the many small regressions
//' done when the terminal slope is searched.
//'
//' @details
//' Weights are optional; with `w = NULL` all points count equally. The
//' adjusted R-squared and the standard error of the slope need at least three
//' points and are `NA` with two. A line cannot be fitted with fewer than two
//' points, with no variation in `x`, or when the weights sum to zero; in these
//' cases every element of the result is `NA`.
//'
//' @param x,y Numeric vectors of the same length.
//' @param w Optional numeric vector of non-negative weights.
//' @return A numeric vector of length six: intercept, slope, R-squared,
//'   adjusted R-squared, standard error of the slope and the residual degrees
//'   of freedom.
//' @seealso [hp_lz_windows()] for the regression over the last n points of a
//'   profile.
//' @family terminal phase
//' @examples
//' x <- c(4, 8, 12, 24)
//' y <- log(c(3, 1.8, 1.1, 0.25))
//' fit <- hp_fit_line(x, y)
//' names(fit) <- c("intercept", "slope", "r2", "adj_r2", "se_slope", "df")
//' fit
//'
//' # agrees with lm()
//' coef(lm(y ~ x))
//' @export
// [[Rcpp::export]]
NumericVector hp_fit_line(NumericVector x, NumericVector y, Nullable<NumericVector> w = R_NilValue) {
  std::vector<double> xv(x.begin(), x.end()), yv(y.begin(), y.end());
  if (w.isNotNull()) {
    NumericVector ww(w);
    std::vector<double> wv(ww.begin(), ww.end());
    return fit_line(xv, yv, &wv);
  }
  return fit_line(xv, yv, nullptr);
}

//' Terminal regressions over the last n points of a profile
//'
//' For every window size in `sizes`, regresses the natural logarithm of the
//' concentration on time over the last `n` points of the profile and returns the
//' slope (as the elimination rate constant), the intercept and the fit
//' statistics. The window with the best adjusted R-squared is the usual choice
//' for the terminal slope, lambda z, and this function computes all candidates
//' in one fast C++ call.
//'
//' @details
//' `t` and `c` are the points that can take part in the terminal phase, in time
//' order, with strictly positive concentrations. The regression is of
//' `log(c)` on `t`, so `lambda = -slope`. A window that cannot be fitted is left
//' out of the result.
//'
//' @param t,c The pool: times and concentrations of the points allowed in the
//'   terminal phase, in time order. Concentrations must be positive.
//' @param sizes Integer vector of the numbers of points to try, for example
//'   `3:6`.
//' @param weighting 0 for uniform weights, 1 for weights `1/y`, 2 for weights
//'   `1/y^2`, where `y` is the concentration.
//' @return A numeric matrix with one row per window that could be fitted and
//'   the columns: number of points, first time, last time, lambda, intercept,
//'   R-squared and adjusted R-squared.
//' @seealso [hp_fit_line()], [hp_nca_options()] for the rules that choose among
//'   the windows.
//' @family terminal phase
//' @examples
//' t <- c(2, 4, 6, 8, 12, 24)
//' c <- c(5, 3, 2, 1.4, 0.7, 0.1)
//' w <- hp_lz_windows(t, c, sizes = 3:6)
//' colnames(w) <- c("n", "start", "end", "lambda", "intercept", "r2", "adj_r2")
//' round(w, 4)
//'
//' # the window with the highest adjusted R-squared
//' w[which.max(w[, "adj_r2"]), ]
//' @export
// [[Rcpp::export]]
NumericMatrix hp_lz_windows(NumericVector t, NumericVector c, IntegerVector sizes, int weighting = 0) {
  int m = t.size();
  std::vector<std::vector<double>> rows;
  for (int k = 0; k < sizes.size(); ++k) {
    int n = sizes[k];
    if (n < 2 || n > m) continue;
    std::vector<double> x, y, w;
    for (int i = m - n; i < m; ++i) {
      x.push_back(t[i]);
      y.push_back(std::log(c[i]));
      if (weighting == 1) w.push_back(1.0 / c[i]);
      else if (weighting == 2) w.push_back(1.0 / (c[i] * c[i]));
    }
    NumericVector f = fit_line(x, y, weighting == 0 ? nullptr : &w);
    if (NumericVector::is_na(f[1])) continue;
    rows.push_back({(double)n, x.front(), x.back(), -f[1], f[0], f[2], f[3]});
  }
  NumericMatrix out(rows.size(), 7);
  for (size_t i = 0; i < rows.size(); ++i)
    for (int j = 0; j < 7; ++j) out(i, j) = rows[i][j];
  colnames(out) = CharacterVector::create("n", "start", "end", "lambda", "intercept", "r2", "adj_r2");
  return out;
}

//' Observed exposure measures of one profile
//'
//' Computes the quick, model-independent exposure measures of a single
//' concentration-time profile: the maximum, its time, the minimum, the last
//' measurable concentration with its time, and the area under the curve to the
//' last measurable concentration. It is the fast C++ routine behind the
//' subject-level exposure columns of the exploratory analysis.
//'
//' @details
//' The area uses linear-up/log-down trapezoids starting at the dose. When the
//' profile has no sample at time zero a concentration of zero at time zero is
//' assumed, as for extravascular dosing. Missing and non-finite values are
//' ignored. For a complete non-compartmental analysis (terminal slope, area to
//' infinity, clearance and so on) use [hp_nca()].
//'
//' @param t,c Numeric vectors: times after the dose and concentrations.
//' @return A named numeric vector with `cmax`, `tmax`, `cmin`, `clast`,
//'   `tlast`, `auclast` and `nobs` (the number of observations used).
//' @seealso [hp_nca()], [hp_eda_prepare()] (which adds these measures to the
//'   subject frame).
//' @family area functions
//' @examples
//' t <- c(0.5, 1, 2, 4, 8, 12, 24)
//' c <- c(1.2, 2.8, 3.5, 2.6, 1.3, 0.7, 0.1)
//' hp_exposure(t, c)
//' @export
// [[Rcpp::export]]
NumericVector hp_exposure(NumericVector t, NumericVector c) {
  std::vector<std::pair<double, double>> p;
  for (int i = 0; i < t.size(); ++i)
    if (std::isfinite(t[i]) && std::isfinite(c[i])) p.push_back({t[i], c[i]});
  std::stable_sort(p.begin(), p.end(),
                   [](const std::pair<double, double>& a, const std::pair<double, double>& b) { return a.first < b.first; });
  NumericVector out = NumericVector::create(NA_REAL, NA_REAL, NA_REAL, NA_REAL, NA_REAL, NA_REAL, 0.0);
  out.names() = CharacterVector::create("cmax", "tmax", "cmin", "clast", "tlast", "auclast", "nobs");
  if (p.empty()) return out;
  double cmax = -INFINITY, tmax = NA_REAL, cmin = INFINITY;
  for (auto& q : p) {
    if (q.second > cmax) { cmax = q.second; tmax = q.first; }
    cmin = std::min(cmin, q.second);
  }
  double tlast = NA_REAL, clast = NA_REAL;
  for (int i = (int)p.size() - 1; i >= 0; --i)
    if (p[i].second > 0.0) { tlast = p[i].first; clast = p[i].second; break; }
  std::vector<std::pair<double, double>> q = p;
  if (q[0].first > 0.0) q.insert(q.begin(), {0.0, 0.0});
  double auc = 0.0;
  for (size_t i = 0; i + 1 < q.size(); ++i) {
    if (std::isnan(tlast) || q[i + 1].first > tlast + 1e-12) break;
    auc += seg_auc(q[i].first, q[i].second, q[i + 1].first, q[i + 1].second, q[i + 1].second < q[i].second);
  }
  out[0] = cmax; out[1] = tmax; out[2] = cmin; out[3] = clast; out[4] = tlast;
  out[5] = std::isnan(tlast) ? NA_REAL : auc; out[6] = (double)p.size();
  return out;
}
