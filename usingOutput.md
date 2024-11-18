# R/bayz methods for summarizing, extracting and using output

-   [summary](#Summary%20of%20output)
-   [coef, fixef and ranef](#Coefficients%20estimates)
-   [vcomp](#Variance%20components%20and%20PVEs)
-   [Computing contrasts](#Computing%20constrasts)

# Summary of output

Use of summary() on the bayz output object produces a summary of
parameter estimates from the fitted model including convergence
diagnostics and Highest Posterior Density (HPD) regions.

The summary() method only lists a limited number of the so-called
“traced” parameters - these are model-parameters for which all MCMC
samples are saved in the output, allowing to compute convergence, HPD
regions, and to plot traces and densities (using plot()). The traced
parameters by default include: all scalar variance parameters, estimated
variance-covariances up to dimension 4x4, the model mean, scalar
regression coefficients, coefficient estimates from fixed effects with
up to 4 levels, and nested regressions with up to 4 levels.

# Coefficients estimates

Apart from using these functions, it is also quite straightforward to
extract estimates directly from the bayz output object. All estimates
are stored in the output in a list called
<output>

$Estimates, which has named elements according to parameter names.

# Variance components and PVEs

The vcomp() method extracts variance estimates from the output.
