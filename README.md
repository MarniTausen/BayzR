Bayesian Linear Mixed Models
============================

[![Build Status](https://travis-ci.org/MarniTausen/BayzR.svg?branch=master)](https://travis-ci.org/MarniTausen/BayzR)[![Coverage Status coveralls](https://coveralls.io/repos/github/MarniTausen/BayzR/badge.svg?branch=master)](https://coveralls.io/github/MarniTausen/BayzR?branch=master)

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

    ## cycle var.resid var.ranf.x mean.y 
    ## 220 1.02433 108.936 13.7026 
    ## 440 1.48455 872.084 13.8791 
    ## 660 0.609921 681.277 10.4569 
    ## 880 0.935697 4077.21 10.6522 
    ## 1100 0.753366 2481.22 12.9702

``` r
summary(fit)
```

    ## Summary statistics of Bayz object
    ## 
    ##    model formula: y ~ ranf(x) 
    ## 
    ## Random variable estimates:
    ##               postMean       postSD
    ## var.resid     1.226162    0.8235952
    ## var.ranf.x 1328.083949 1258.4308083
    ## 
    ## Variance explained (Heritability):
    ##   Variable Heritability
    ## 1    resid 0.0009224052
    ## 2        x 0.9990775948
