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

    ## 220 0.995131 396.876 9.12154 
    ## 440 1.35793 641.307 -6.45424 
    ## 660 0.795633 915.299 3.60962 
    ## 880 0.698552 682.898 7.90645 
    ## 1100 2.25099 255.516 7.85459

``` r
summary(fit)
```

    ## Summary statistics of Bayz object
    ## 
    ##     model formula: y ~ ranf(x) 
    ## 
    ## Random variable estimates:
    ##              postMean      postSD
    ## var.resid    1.269095   0.6153912
    ## var.ranf.x 853.870200 623.5912124
    ## 
    ## Broad sense heritability:
    ## resid:        0.148408 %
    ## x:        99.85159 %
