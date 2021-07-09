
Bayesian Linear Mixed Models
============================

[![Build Status](https://travis-ci.org/MarniTausen/BayzR.svg?branch=master)](https://travis-ci.org/MarniTausen/BayzR)[![Coverage Status coveralls](https://coveralls.io/repos/github/MarniTausen/BayzR/badge.svg?branch=master)](https://coveralls.io/github/MarniTausen/BayzR?branch=master)

-   [How to install](#how-to-install)
-   [Functions and models associated with Bayz](#functions-and-models-associated-with-bayz)
-   [Examples](#examples)

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

-   fx (Fixed effect as factor)
-   rn (Random effect factor)
-   rg (Fixed continuous variables)
-   rr (random regression on covariate table/matrix)

Examples
========

``` r
library(BayzR)

x <- rnorm(500, 10, 5)

y <- 5 * x + 4 + rnorm(500, 0, 4)

example_data = data.frame(y = y, x = x)

fit <- bayz(y ~ ranf(x), data=example_data, chain=c(10000, 100, 10))
```

    ## 2000 411.74 211.907 55.5413 
    ## 4000 526.394 82.3809 55.3901 
    ## 6000 341.156 275.283 52.6522 
    ## 8000 226.035 347.902 55.747 
    ## 10000 518.322 115.43 54.8652

``` r
summary(fit)
```

    ## Summary statistics of Bayz object
    ## 
    ##    model formula: y ~ ranf(x) 
    ## 
    ## Random variable estimates:
    ##            postMean   postSD
    ## var.resid  420.8849 176.6630
    ## var.ranf.x 190.5818 174.6422
    ## 
    ## Variance explained (Heritability):
    ##   Variable Heritability
    ## 1    resid    0.6883202
    ## 2        x    0.3116798
