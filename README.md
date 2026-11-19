# hyperutils

**Fast exploratory analysis, non-compartmental analysis and custom outputs for
clinical pharmacology and pharmacometrics.**

`hyperutils` is a standalone R package for the everyday data work of
pharmacokinetic (PK) and pharmacodynamic (PD) studies. It reads NONMEM-style
datasets, explores them with ready-made figures and tables, runs
non-compartmental analysis (NCA) of every dosing interval of every subject,
and builds custom figures and tables from plain lists of settings. The
numerical core is written in C++ with 'Rcpp'.

* **Exploratory data analysis (EDA)** in a few lines: individual and typical
  concentration-time profiles, dose-normalised plots, covariate and exposure
  plots, distributions, correlations, data summaries and missing-value tables.
* **Non-compartmental analysis** with the methods of standard NCA software:
  three trapezoidal rules, six rules for data below the limit of
  quantification, best-fit selection of the terminal slope, partial areas,
  oral, intravenous bolus and infusion dosing, repeated doses and steady state.
  Every choice is documented in tables you can attach to a report.
* **Custom plots and tables** for reports: a plot or a table is a plain list
  of settings, so it can be stored, reviewed and reused. Plots are ordinary
  'ggplot2' objects.
* **Reproducible by design**: nothing is written to disk unless you give a
  folder, every setting is an argument, and every function has a help page with
  examples.

## Installation

From CRAN:

```r
install.packages("hyperutils")
```

The development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("hyperpkpd/hyperutils")
```

### Compilation

The package contains C++ code (written with 'Rcpp'), so it must be compiled
when it is installed from source. Binary packages from CRAN need nothing
extra. To install from source you need a C++ compiler:

| System | What to install |
|--------|-----------------|
| Windows | [Rtools](https://cran.r-project.org/bin/windows/Rtools/) matching your R version |
| macOS | the Xcode command line tools: `xcode-select --install` |
| Linux | a compiler toolchain, for example `build-essential` on Debian and Ubuntu |

The package depends on 'ggplot2' and 'Rcpp', which are installed
automatically.

## Quick start

```r
library(hyperutils)

d <- hp_read_data(hp_example_file())        # a NONMEM-style dataset

# 1. Non-compartmental analysis
res <- hp_nca(d, route = "oral", lloq = 0.5, groups = "DOSEGRP")
res$parameters[1:5, c("ID", "CMAX", "TMAX", "AUCLAST", "AUCINF_OBS", "THALF")]
res$summary                                  # statistics by dose group

# 2. Exploratory figures and tables
e <- hp_eda_prepare(d, dvid = "DVID", lloq = 0.5, color = "DOSEGRP")
spec <- hp_eda_spec(kind = "line", x = "TAFD", y = "DV", color = "DOSEGRP",
                    pk_only = TRUE, log_y = TRUE)
hp_eda_plot(e, spec)[[1]]
hp_eda_table(e, hp_eda_tspec(kind = "summary", scope = "subjects",
                             vars = c("WT", "AGE"), by = "DOSEGRP"))

# 3. A custom figure and a custom table
subj <- hp_custom_prep(d, distinct = "ID")
hp_custom_plot(subj, list(geoms = c("boxplot", "jitter"), x = "factor(DOSEGRP)",
                          y = "WT", fill = "factor(DOSEGRP)", legend = "none"))
hp_custom_table(subj, list(kind = "demographics", vars = c("WT", "AGE", "SEX"),
                           by = "DOSEGRP", total = TRUE))

# 4. Save (every writing function takes the folder)
hp_save_table(res$parameters, "nca_parameters.csv", tempdir())
```

## What it does

### Reading and preparing data

| Function | Purpose |
|----------|---------|
| `hp_read_data()` | read comma, tab, semicolon or space delimited files and NONMEM table files |
| `hp_mutate()`, `hp_filter()`, `hp_case_when()` | data steps written as R expressions |
| `hp_custom_prep()` | mutate, filter, distinct, sort, head and select in one call |
| `hp_describe()` | n, mean, SD, CV, geometric mean, percentiles and confidence interval |

### Exploratory data analysis

`hp_eda_prepare()` derives the time after dose, dose-normalised
concentrations, BLQ flags and the observed exposure of each subject, and
returns a record frame and a subject frame. Figures are described by
`hp_eda_spec()` and drawn by `hp_eda_plot()`: scatter, line, summary, box,
histogram, bar, correlation, pairs and PK-PD plots. Tables are described by
`hp_eda_tspec()` and built by `hp_eda_table()`: summary statistics, counts and
listings. `hp_eda_data_summary()`, `hp_eda_missing()` and `hp_eda_correlation()`
give the standard data-quality tables.

### Non-compartmental analysis

`hp_nca()` returns, for every dosing interval:

* peak and trough: `CMAX`, `TMAX`, `CMIN`, `TLAG`, `C0`, `CLAST`, `TLAST`
* terminal phase: `LAMZ`, `THALF`, number of points, R-squared
* areas and moments: `AUCLAST`, `AUCINF_OBS`, `AUCINF_PRED`, `AUMC*`, `MRT*`,
  percent extrapolated, partial areas
* clearance and volume: `CL`, `VZ`, `VSS` (intravenous) and `CL_F`, `VZ_F`
  (extravascular)
* dosing-interval measures: `AUCTAU`, `CAVG`, `CTAU`, fluctuation, swing,
  peak-to-trough ratio, steady-state clearance, accumulation ratios
* dose-normalised values and quality flags

| Setting | Choices |
|---------|---------|
| Area method | `linear`, `linear-log`, `linear-up-log-down` |
| BLQ rule | `zero`, `missing`, `lloq`, `half-lloq`, `tmax-zero-missing`, `first-zero-mid-half-last-missing` |
| Terminal slope | `best-fit` (adjusted R-squared), `last-n`, `time-range`, `none`, per-subject overrides |
| Route | `oral`, `iv`, `bolus`, `infusion`, `mixed`, `auto` |
| Steady state | from the data (`SS`, `II`), `last`, `all`, `none` |
| Duplicates | `mean`, `first`, `last`, `max`, `error` |

The result also contains the summary statistics by stratum
(`hp_nca_summary()`), the table of every terminal regression tried
(`hp_nca_lambda_z_table()`), and the table of every sample and what was done
with it (`hp_nca_profile_table()`). `hp_nca_plots()` draws individual and mean
profiles, box plots of parameters and pages of individual fits.

### Custom plots and tables

`hp_custom_plot()` draws points, lines, smoothers, box and violin plots,
bars, histograms, densities, error bars, heat maps and more, with facets,
log axes, reference lines, themes and palettes. `hp_custom_table()` builds
listings, summaries, demographics tables, frequency tables, pivot tables and
correlation tables.

### Output

`hp_save_plot()`, `hp_save_plots()` and `hp_save_table()` write PNG, SVG and
CSV files to the folder you name. `hp_theme()` is the plot theme of the
package, available for your own plots.

### C++ kernels

The building blocks are exported for those who want to write their own
analyses: `hp_seg_auc()`, `hp_seg_aumc()`, `hp_curve_area()`,
`hp_curve_conc_at()`, `hp_fit_line()`, `hp_lz_windows()` and `hp_exposure()`.

## Learn more

The vignettes walk through each part with examples:

| Vignette | Content |
|----------|---------|
| `vignette("hyperutils")` | getting started |
| `vignette("data-preparation")` | reading files, deriving columns, filtering |
| `vignette("exploratory-analysis")` | EDA step by step |
| `vignette("nca-step-by-step")` | NCA step by step: BLQ, terminal slope, IV, steady state |
| `vignette("custom-outputs")` | custom plots and tables for reports |
| `vignette("reporting-workflow")` | a reproducible script, verification and scaling |

Every function has a help page: `?hp_nca`, `?hp_eda_plot`, `?hp_custom_table`.
`help(package = "hyperutils")` lists them all.

## Performance

The areas, moments and terminal-phase regressions are computed in C++, and
the surrounding bookkeeping is vectorised. An NCA of about a thousand profiles
takes a few seconds on an ordinary laptop, and the time grows in proportion
to the number of profiles. Exploratory preparation of a dataset of that size
takes a fraction of a second.

## Compatibility with nlmx

`hyperutils` is highly compatible with 'nlmx', a high-throughput,
high-efficiency pharmacometrics modelling software. The file layout, the
column conventions and the derivations agree, and the reproducible R scripts
that 'nlmx' can write for its exploratory, non-compartmental and custom outputs
use this package. 'nlmx' is not required: everything in `hyperutils` works on
its own.

## Citation

To cite the package, use `citation("hyperutils")`, or:

> Obianom O (2026). *hyperutils: Fast Exploratory Analysis, Non-Compartmental
> Analysis and Custom Outputs for Clinical Pharmacology*. R package version
> 0.1.0. <https://github.com/hyperpkpd/hyperutils>

## Getting help and contributing

Report problems and suggest features at
<https://github.com/hyperpkpd/hyperutils/issues>. A minimal reproducible
example (a few rows of data and the call) makes a report much easier to act
on.

## Author and licence

Written and maintained by Obinna Obianom,
<https://github.com/hyperpkpd>. Released under the MIT licence; see the
`LICENSE` file.
