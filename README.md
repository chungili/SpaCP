# SpaCP
Conformal prediction intervals for spatial frailty models with survival outcomes

## Installation

You can install the development version of SpaCP from GitHub using the `devtools` package:

```r
# Install devtools if not already installed
install.packages("devtools")

# Install SpaCP from GitHub
devtools::install_github("chungili/SpaCP")
```

## Dependencies

SpaCP requires the following packages to be installed:

```r
install.packages(c(
  "Rcpp", "RcppEigen", "Matrix", "survival",
  "nleqslv", "doParallel", "foreach",
  "kernlab", "doSNOW", "snow"
))
```

> **Note for Windows users:** [Rtools](https://cran.r-project.org/bin/windows/Rtools/) must be installed to compile the C++ code.

## Usage

### Conformal Prediction Intervals

Use `Predict.Sp()` to compute conformal prediction intervals for spatial survival outcomes:

```r
library(SpaCP)
# Load built-in datasets
data(train.dt)
data(test.dt)

# Compute 95% conformal prediction intervals
result <- Predict.Sp(
  train_dt = train.dt,
  test_dt  = test.dt,
  B        = 10,
  alpha    = 0.05
)

# View results
head(result)
```