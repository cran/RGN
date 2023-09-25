
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RGN

<!-- badges: start -->
<!-- badges: end -->

This repository contains the R package for the Robust Gauss-Newton (RGN)
algorithm, which is designed for solving optimization problems with a
sum of least squares objective function.

This R implementation is developed by David McInerney and Michael
Leonard, and is based on original RGN Fortran code developed by Youwei
Qin, Dmitri Kavetski and George Kuczera
(<https://github.com/eachonly/Robust-Gauss-Newton-Algorithm/>).

When using RGN please cite the following articles:

Qin Y, Kavetski D, Kuczera G (2018) A robust Gauss-Newton algorithm for
the optimization of hydrological models: From standard Gauss-Newton to
robust Gauss-Newton. Water Resources Research, 54.
<https://doi.org/10.1029/2017WR022488>

Qin Y, Kavetski D, Kuczera G (2018) A robust Gauss-Newton algorithm for
the optimization of hydrological models: Benchmarking against
industry-standard algorithms. Water Resources Research, 54.
<https://doi.org/10.1029/2017WR022489>

## Installation

You can install the development version of RGN from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ClimateAnalytics/RGN")
```

## Example

The first example is optimisation of a 2D Rosenbrock function:

``` r
library(RGN)
# Example 1: Rosenbrock
simFunc_rosenbrock=function(x) c(1.0-x[1],10.0*(x[2]-x[1]**2))
rgnOut = rgn(simFunc=simFunc_rosenbrock,
             par=c(-1.0,  0.0), lower=c(-1.5, -1.0), upper=c( 1.5,  3.0),
             simTarget=c(0,0))
rgnOut$par #optimal parameters
#> [1] 1 1
rgnOut$value #optimal objective function value
#> [1] 0
```

The second example is calibration of the 5 parameter hydrological model
HYMOD:

``` r
library(RGN)
# Example 2: Hymod
data("BassRiver") # load Bass River hydrological data
rgnOut = rgn(simFunc=simFunc_hymod,
             par=c(400.,0.5,0.1,0.2,0.1),
             lower=c(1.,0.1,0.05,0.000001,0.000001),
             upper=c(1000.,2.,0.95,0.99999,0.99999),
             simTarget=BassRiverData$Runoff.mm.day[365:length(BassRiverData$Date)],
             stateVal=c(100.0,30.0,27.0,25.0,30.0,0.0,0.0,0.0), # initial states for hymod
             nWarmUp=365,                                       # warmup period
             rain=BassRiverData$Rain.mm,                        # precip input
             pet=BassRiverData$ET.mm)                           # PET input
rgnOut$par #optimal parameters
#> [1] 146.7563960   0.3635988   0.1895957   0.9999900   0.7430698
rgnOut$value #optimal objective function value
#> [1] 6840.165
```
