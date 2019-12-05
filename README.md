Bayesian linear regression
================
Maria Izabel Cavassim Alves
2019-12-05

-   [Bayesian Linear Mixed Models](#bayesian-linear-mixed-models)
-   [How to install](#how-to-install)
-   [Functions and models associated with BayZ](#functions-and-models-associated-with-bayz)
-   [Examples](#examples)

Bayesian Linear Mixed Models
============================

[![Build Status](https://travis-ci.org/MarniTausen/BayzR.svg?branch=master)](https://travis-ci.org/MarniTausen/BayzR)[![Coverage Status coveralls](https://coveralls.io/repos/github/MarniTausen/BayzR/badge.svg?branch=master&kill_cache=1)](https://coveralls.io/github/MarniTausen/BayzR?branch=master)

Here write down the applicability of the software. For more information access the website: <http://www.bayz.biz/>.

How to install
==============

Bayz can be used in two ways: on a unix/linux/windows prompt, or from R (available for Mac and Windows). Use from R is limited to some standard models and is now mostly suited for using smaller data sets. See the Bayz in R section. Use from a prompt is more versatile but requires writing a bayz script, examples are in the Bayz2.5 Examples section.

To install this current development version, you first need to install devtools.

``` r
install.packages("devtools")
```

To install the R package on your system you can then use:

``` r
devtools::install_github("MarniTausen/BayzR")
```

Functions and models associated with BayZ
=========================================

Examples
========

``` r
library(BayzR)

x <- rnorm(500, 10, 5)

y <- 5 * x + 4 + rnorm(500, 0, 4)

example_data = data.frame(y = y, x = x)

fit <- bayz(y ~ ranf(x), data=example_data, chain=c(10000, 100, 10))
```

    ## cycle var.resid var.ranf.x mean.y 
    ## 2000 1.22833 805.1 50.6616 
    ## 4000 0.940035 772.305 48.9582 
    ## 6000 1.20571 783.969 48.2655 
    ## 8000 0.902765 996.047 47.3541 
    ## 10000 0.95025 753.195 47.1943

``` r
summary(fit)
```

    ## Summary statistics of Bayz object
    ## 
    ##    model formula: y ~ ranf(x) 
    ## 
    ## Random variable estimates:
    ##              postMean      postSD
    ## var.resid    1.010605   0.0934522
    ## var.ranf.x 763.445986 167.6535818
    ## 
    ## Variance explained (Heritability):
    ##   Variable Heritability
    ## 1    resid  0.001321991
    ## 2        x  0.998678009
