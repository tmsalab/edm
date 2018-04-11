
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ecdm

[![Project Status: Active - The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Travis-CI Build
Status](https://travis-ci.org/tmsalab/ecdm.svg?branch=master)](https://travis-ci.org/tmsalab/ecdm)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/ecdm)](https://cran.r-project.org/package=ecdm)
[![CRAN RStudio mirror
downloads](http://cranlogs.r-pkg.org/badges/ecdm)](http://www.r-pkg.org/pkg/ecdm)
[![Coverage
Status](https://img.shields.io/codecov/c/github/tmsalab/ecdm/master.svg)](https://codecov.io/github/tmsalab/ecdm?branch=master)

The goal of ecdm is to provide a modeling framework for exploratory
cognitive diagnostic models and classical cognitive diagnostic models.

## Installation

The `ecdm` package is currently only available via GitHub. To install
`ecdm`, your computer will need to have a compiler. The following guides
are avaliable:

  - [Windows:
    Rtools](http://thecoatlessprofessor.com/programming/installing-rtools-for-compiled-code-via-rcpp/)
  - [macOS:
    Rtools](http://thecoatlessprofessor.com/programming/r-compiler-tools-for-rcpp-on-macos/)

From there, please use `devtools` to retrieve the latest development
version.

``` r
# install.packages("devtools")
devtools::install_github("tmsalab/ecdm")
```

## Usage

Load the `ecdm` package into *R*:

``` r
library(ecdm)
```

Exploratory CDM models can be estimated with:

``` r
edina_model = edina(<data>, <k>)
errum_model = errum(<data>, <k>)
```

Classical CDMs can be estimated using:

``` r
dina_model = dina(<data>, <q>)
rrum_model = rrum(<data>, <q>)
```

These classical CDMs are implemented in separate packages: `dina` and
`rrum`.

## Details

The `ecdm` package is designed to act more as a “virtual” package. The
main functionalities of `ecdm` are split across multiple packages. The
rationale for this is many areas of psychometrics have overlap in terms
of computational code used. By dividing the underlying source of the
`ecdm` package, we are enabling fellow psychometricians to be able to
incorporate established routines into their own code. In addition, we
are lowering the amount of redundancies, or copy and pasted code, within
the CDM framework we are building.

Specifically, the `ecdm` package imports:

  - `dina`: Estimating the Deterministic Input, Noisy “And” Gate (DINA)
    cognitive diagnostic model parameters using a Gibbs sampler.
  - `rrum`: Estimating the reduced Reparametrized Unified Model (rRUM)
    with a Gibbs sampler.
  - `shinyecdm`: User Interface for Modeling with Exploratory Models
  - `simcdm`: Simulate responses underneath a DINA or rRUM model.
  - `rgen`: Simulate Multivariate Probability Distributions

# License

GPL (\>= 2)
