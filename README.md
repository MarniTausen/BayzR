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
example_data = data.frame(y = 1:20, x = 20:1)

fit <- bayz(y ~ ranf(x), data=example_data)
```

    ## Warning in bayz(y ~ ranf(x), data = example_data): running the default chain of 1100 cycles, this may be too short for many analyses

    ## 220 1.00931 71.6876 0.84113 
    ## 440 0.444949 1.29795 2.35149 
    ## 660 2.62573 24.901 -0.0419425 
    ## 880 1.03553 61.0391 -3.84536 
    ## 1100 0.648931 164.79 0.39281

``` r
summary(fit)
```

    ## Summary statistics of Bayz object
    ## 
    ##    model formula: y ~ ranf(x) 
    ## 
    ## Random variable estimates:
    ##             postMean     postSD
    ## var.resid   1.277064  0.5715017
    ## var.ranf.x 71.093073 80.7595373
    ## 
    ## Variance explained (Heritability):
    ##   Variable Heritability
    ## 1    resid   0.01764628
    ## 2        x   0.98235372
