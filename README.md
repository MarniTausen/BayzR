Bayesian Linear Mixed Models
============================

[![Build Status](https://travis-ci.org/MarniTausen/BayzR.svg?branch=master)](https://travis-ci.org/MarniTausen/BayzR)[![Coverage Status](https://coveralls.io/repos/github/MarniTausen/BayzR/badge.svg?branch=master)](https://coveralls.io/github/MarniTausen/BayzR?branch=master)

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

    ## 220 1.13587 146.52 11.7294 
    ## 440 0.905233 160.56 8.15912 
    ## 660 2.2259 214.22 5.71739 
    ## 880 2.06531 225.052 3.43328 
    ## 1100 0.582729 180.863 5.30218

``` r
summary(fit)
```

    ## Summary statistics of Bayz object
    ## 
    ##     model formula: y ~ ranf(x) 
    ## 
    ## Random variable estimates:
    ##              postMean      postSD
    ## var.resid    1.166754   0.5427654
    ## var.ranf.x 182.671712 114.4047595
    ## 
    ## Broad sense heritability:
    ## resid:        0.6346629 %
    ## x:        99.36534 %
