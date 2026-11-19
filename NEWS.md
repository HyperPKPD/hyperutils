# hyperutils 0.1.0

First release.

## Data preparation

* `hp_read_data()` reads comma, tab, semicolon and white-space delimited files
  and NONMEM table files, with `.` as a missing value and upper-cased names.
* `hp_mutate()`, `hp_filter()`, `hp_case_when()` and `hp_custom_prep()` run data
  steps written as R expressions in a fixed, documented order.
* `hp_describe()` and `hp_centre_spread()` give the statistics of a
  clinical pharmacology summary.

## Exploratory data analysis

* `hp_eda_prepare()` builds a record frame and a subject frame with the time
  after dose, dose-normalised concentrations, BLQ flags and the observed
  exposure of each individual.
* `hp_eda_spec()` and `hp_eda_plot()` draw scatter, line, summary, box,
  histogram, bar, correlation, pairs and PK-PD figures, split over pages when
  large.
* `hp_eda_tspec()` and `hp_eda_table()` build summary, count and listing
  tables; `hp_eda_data_summary()`, `hp_eda_missing()`, `hp_eda_correlation()`
  and `hp_eda_subjects()` give the standard data-quality tables.

## Non-compartmental analysis

* `hp_nca()` analyses every dosing interval of every subject: oral,
  intravenous bolus and infusion dosing, single and repeated doses, steady
  state.
* Three trapezoidal methods, six BLQ rules, best-fit, last-n and time-range
  selection of the terminal slope with per-subject overrides, partial areas,
  duplicate-sample handling and quality flags.
* `hp_nca_summary()`, `hp_nca_lambda_z_table()` and `hp_nca_profile_table()`
  document the statistics, the terminal regressions and the handling of every
  sample.
* `hp_nca_plots()` draws individual and mean profiles, box plots and
  individual fits.
* The step functions `hp_nca_read()`, `hp_nca_profiles()` and
  `hp_nca_analyse()` are exported for work on single profiles.

## Custom outputs

* `hp_custom_plot()` builds 'ggplot2' figures from a list of settings with
  twenty layers, facets, scales, reference lines, themes and palettes.
* `hp_custom_table()` builds listing, summary, demographics, frequency, pivot
  and correlation tables.

## Output

* `hp_save_plot()`, `hp_save_plots()` and `hp_save_table()` write PNG, SVG and
  CSV files only to a folder given by the user.
* `hp_theme()` is the plot theme of the package.

## Performance

* The areas, moments, concentration interpolation and terminal-phase
  regressions are written in C++ with 'Rcpp' and exported as
  `hp_seg_auc()`, `hp_seg_aumc()`, `hp_curve_area()`, `hp_curve_conc_at()`,
  `hp_fit_line()`, `hp_lz_windows()` and `hp_exposure()`.

## Documentation

* Every function has a help page with details, use cases and examples.
* Six vignettes: getting started, data preparation, exploratory analysis,
  non-compartmental analysis step by step, custom plots and tables, and a
  reproducible reporting workflow.
