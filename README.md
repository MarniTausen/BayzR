# Bayesian mixed models, shrinkage, sparse and interaction kernel regression

[![Build
Status](https://travis-ci.org/MarniTausen/BayzR.svg?branch=master)](https://travis-ci.org/MarniTausen/BayzR)[![Coverage
Status
coveralls](https://coveralls.io/repos/github/MarniTausen/BayzR/badge.svg?branch=master)](https://coveralls.io/github/MarniTausen/BayzR?branch=master)

-   [How to install](#how-to-install)
-   [Short manual for using the bayz() main function](#Short-Manual)
-   [Example](#examples)
-   [Summarizing and using output](#README5.html)

The Rbayz packages has a single main function ‘bayz’ that fits
mixed-linear, Bayesian shrinkage, (sparse) kernel-regression, (kernel)
interaction and multi-trait models with complex covariance structures
using an extended R-formula syntax. The R help ?bayz contains basic help
information; full details are available on \<ljanss.github.io/Rbayz\>.

# How to install

Most users would like to install our binary packages, that we have
available for Windows 10 and MacOS in the bayz.biz website
<http://www.bayz.biz/>.

If you need or want to use the source code for installation, we have two
repositories on github: in ljanss/BayzR, and in MarniTausen/BayzR.
Installation from the ljanss/BayzR is recommended as first option, this
repository contains regularly ‘frozen’ version that are tested on
different data sets and compile-tested on Windows and Mac-OS. The
version in MarniTausen/BayzR is our development version and may only be
of interest to those that want to co-develop and submit changes to our
source code. You may run into bugs or compilation errors when using the
MarniTausen/BayzR version because we may just be in the middle of
developing / testing something.

To install from the github repositories, the following is needed:

1). download and install Rtools. This is not an R-package but a separate
part of the R system and can be found on r-project.org.

2). install the devtools package.

``` r
install.packages("devtools")
```

3). Install dependencies lme4, coda and Rcpp.

``` r
install.packages("lme4")
install.packages("coda")
install.packages("Rcpp")
```

4). Download and install/compile the BayzR package using:

``` r
library(devtools)
devtools::install_github("ljanss/BayzR")
```

You may be asked to install several other packages. If you are sure you
want the development-version, replace in step 3 to install from
MarniTausen/BayzR. If the compilation succeeds, you are ready to use the
main function bayz() after loading the package with library(BayzR).

# Short Manual

The bayz function fits various mixed-linear and Bayesian shrinkage
models with complex covariance structures using an extended R-formula
syntax.

The model formula in bayz has the basic syntax of an R formula but with
all explanatory (right-hand-side) terms wrapped by a function to specify
how to fit the explanatory variables in the model. This may look like
Yield ~ fx(Year) + rn(Variety) to fit Yield with Year as a fixed factor
and Variety as a random factor. The equivalent lme4 model would be Yield
~ factor(Year) + (1\|Variety). The list of model functions currently
available is: fx() : fixed factors (with interactions) rn() : random
factors (with interactions) rg() : fixed regressions (with interactions
or nested in a factor) rr() : random regressions The model functions
allow to specify interactions of variables, hierarchies, use of matrices
as input data, complex covariance structures, and prior distributions
can be changed to modify standard mixed model in Bayesian shrinkage
models.

Interactions between fators are specified using the colon, for instance
by writing for a fixed Year-Location interaction fx(Year:Location). Bayz
does NOT support automatic expansion with main effects by using the
‘star’ (Year\*Location) or ‘forward slash’ (Year/Location) syntaxes,
hence bayz requires to manually add the desired main effects in the
model (but note that bayz uses the ‘forward slash’ to specify
hierarchical models). Interactions between factors can be specified to
any degree. The ‘star’ syntax can be used to indicate interaction
between covariates, like rg(TempSum\*Precip), where it is simply
interpreted as multiplication. Interaction between a covariate and a
factor is specified using the ‘pipe’ character as in lme4 models. It can
be used to specify fixed nested regressions as
rg(TempSum\|Year:Location) to specify regressions on TempSum within each
Year-Location, or to speficy random slope models in rr(), for instance
rr(TempSum\|Variety).

The bayz call has a data argument to specify an input data-frame (the
“main data”) to contain model variables, but bayz will also search the R
environment if variables are not found in the main data. Typically,
input in the form of matrices such as large sets of covariates,
proportional variance-covariance / correlation / kernel / similarity
matrices used in variance models, as well as data for hierarchical
models, are not in the main data. Matrices / kernels should be prepared
with row-names to match a variable in the data. For matrices / kernels
used in variance models, the link is straightforward, for instance in
rn(Variety:Location, V=KG\*KE), KG and KE can be a genetic and
environmental kinship / kernel matrix, respectively, and KG should have
row-names matching Variety levels, while KE should have row-names
matching Location levels. To fit a large set of covariates, the match is
specified using a hierachical specification rr(Variety/Metabolites),
where Metabolites is then a matrix of covariates which must have
row-names matching Variety levels. In both cases, such kernel or
covariate matrices must have unique levels, but the main data may have
repeated levels and in different order. If one has repeated metabolite
data on each Variety, for instance at multiple time-points, then
consider that Variety is not the appropriate link to the data, but
Variety:Time is.

For random effects a variance-covariance structure can be specified
using a V= option within the model-term function, for example
rn(Variety, V=KG). When the fit is for an interaction of factors, the
variance specification should expand to include one term for each
variable, separated by stars (which should be read as Kronecker
products). The variance structure is then built up from a combination of
given (proportional) variance-covariance / correlation / kernel matrices
and predefined acronyms IDEN, DIAG and VCOV that indicate parameterized
matrices (with parameters to be estimated from the data). Alternatively,
the variance structure can be specified as a linear model using V=~,
which is interpreted as a use of a log-linear model for the variances.

# Examples

``` r
library(BayzR)
# A dummy test for bayz
y = rnorm(1000)
A = as.factor(rep(1:20,50))
B = as.factor(rep(1:10,each=100))
dat1 = data.frame(y,A,B)
library(BayzR)
fit1 = bayz(y~fx(A)+fx(B)+rn(A:B),data=dat1,chain=c(2000,100,20))
summary(fit1)
```
