
Bayesian Linear Mixed Models
============================

[![Build Status](https://travis-ci.org/MarniTausen/BayzR.svg?branch=master)](https://travis-ci.org/MarniTausen/BayzR)[![Coverage Status coveralls](https://coveralls.io/repos/github/MarniTausen/BayzR/badge.svg?branch=master)](https://coveralls.io/github/MarniTausen/BayzR?branch=master)

-   [How to install](#how-to-install)
-   [Short manual for using the bayz() main function](#Short-Manual)
-   [Examples](#examples)


How to install
==============

To install this current development version, you first need to install devtools.

``` r
install.packages("devtools")
```

To install the R package on your system you can then use:

``` r
devtools::install_github("MarniTausen/BayzR")
```

Short Manual
============

The bayz function fits various mixed-linear and Bayesian shrinkage models with complex
covariance structures using an extended R-formula syntax.

The model formula in bayz has the basic syntax of an R formula but with all
explanatory (right-hand-side) terms wrapped by a function to specify how to fit
the explanatory variables in the model. This may look like Yield ~ fx(Year) + rn(Variety)
to fit Yield with Year as a fixed factor and Variety as a random factor. 
The equivalent lme4 model would be Yield ~ factor(Year) + (1|Variety). 
The list of model functions currently available is:
   fx() : fixed factors (with interactions)
   rn() : random factors (with interactions)
   rg() : fixed regressions (with interactions or nested in a factor)
   rr() : random regressions 
The model functions allow to specify interactions of variables, hierarchies, use of
matrices as input data, complex covariance structures, and prior distributions can 
be changed to modify standard mixed model in Bayesian shrinkage models.

Interactions between fators are specified using the colon, for instance by writing for a 
fixed Year-Location interaction fx(Year:Location). Bayz does
NOT support automatic expansion with main effects by using the 'star' (Year*Location) or
'forward slash' (Year/Location) syntaxes, hence bayz requires to manually add the desired
main effects in the model (but 
note that bayz uses the 'forward slash' to specify hierarchical models).
Interactions between factors can be specified to any degree.
The 'star' syntax can be used to indicate interaction
between covariates, like rg(TempSum*Precip), where it is simply interpreted as
multiplication. 
Interaction between a covariate and a factor is specified using the 'pipe'
character as in lme4 models. It can be used to specify fixed nested regressions as 
rg(TempSum|Year:Location) to specify regressions on TempSum within each Year-Location, or
to speficy random slope models in rr(), for instance rr(TempSum|Variety). 

The bayz call has a data argument to specify an input data-frame (the "main data")
to contain model variables,
but bayz will also search the R environment if variables are not found in the main data.
Typically, input in the form of matrices such as large sets of covariates, proportional
variance-covariance / correlation / kernel / similarity matrices used in variance models,
as well as data for hierarchical models, are not in the main data.
Matrices / kernels should be prepared with row-names to match a variable
in the data. For matrices / kernels used in variance models, the link is straightforward, 
for instance in rn(Variety:Location, V=KG*KE), KG and KE can be a genetic and environmental
kinship / kernel matrix, respectively, and KG should have row-names matching Variety levels,
while KE should have row-names matching Location levels. 
To fit a large set of covariates, the match is specified using a hierachical specification
rr(Variety/Metabolites), where Metabolites is then a matrix of covariates which must have
row-names matching Variety levels. In both cases, such kernel or covariate matrices
must have unique levels, but the main data may have repeated levels and in different order.
If one has repeated metabolite data on each Variety, for instance at multiple time-points,
then consider that Variety is not the appropriate link to the data, but Variety:Time is.

For random effects a variance-covariance structure can be specified using a
V= option within the model-term function, for example rn(Variety, V=KG).
When the fit is for an interaction of factors, the variance specification should
expand to include one term for each variable, separated by stars (which should
be read as Kronecker products). The variance structure
is then built up from a combination of given (proportional) variance-covariance / correlation /
kernel matrices and predefined acronyms IDEN, DIAG and VCOV that indicate parameterized
matrices (with parameters to be estimated from the data). 
Alternatively, the variance structure can be specified as a linear model using V=~, 
which is interpreted as a use of a log-linear model for the variances. 

Examples
========

``` r
library(BayzR)
# A test data included in the package with Year, Location and Year-Location interaction effects
testdat1 = read.table("./tests/testdat1.txt", header=TRUE)
fit1 <- bayz(y~fx(YR)+fx(LC)+rn(YR:LC),data=testdat1,chain=c(5000,100,10))
summary(fit1)
# > produces output
#  Summary of bayz model fit
#  
#  model formula: y ~ fx(YR) + fx(LC) + rn(YR:LC) 
#  
#  Estimates for model coefficients(*):
#            postMean    postSD
#  mean     2.9314963 0.5929937
#  YR%yr1   0.0000000 0.0000000
#  YR%yr2  -1.9824071 0.6473416
#  YR%yr3  -0.4097616 0.5978790
#  YR%yr4  -1.1269487 0.5850317
#  YR%yr5  -2.2186498 0.5925470
#  YR%yr6  -2.7432130 0.5683635
#  LC%loc1  0.0000000 0.0000000
#  LC%loc2 -2.8570509 0.4706350
#  LC%loc3 -2.3987149 0.5011672
#  LC%loc4 -0.4855762 0.5704899
#  * Some estimates are not shown because they have more than maxLevel levels:
#    YR:LC
#  
#  Estimates and HPD intervals for hyper-parameters and other 'logged' parameters:
#             postMean     postSD   HPDleft HPDright
#  var.y     0.9390555 0.05569042 0.8333453 1.055571
#  mean      2.9314963 0.59359969 1.8454898 4.003390
#  var.YR:LC 0.6584807 0.31621030 0.1525434 1.317420
#  
#  Convergence diagnostics on 'logged' parameters:
#              effSize   GewekeZ        MCSE     MCCV%
#  var.y     448.30317 0.2232456 0.002630235 0.2800937
#  mean       14.95606 0.4236064 0.153491780 5.2359534
#  var.YR:LC 117.84447 2.2537794 0.029128721 4.4236258
#  
#  Estimates for explained variances (variances as proportions of total):
#             postMean    postSD   HPDleft  HPDright
#  var.y     0.6065921 0.1016172 0.3992172 0.8023483
#  var.YR:LC 0.3934079 0.1016172 0.1976517 0.6007828
```
