#' Bayesian mixed linear and shrinkage models
#'
#' The bayz function fits various mixed-linear and Bayesian shrinkage models with complex
#' covariance structures using an extended R-formula syntax.
#'
#' The model formula in bayz has the basic syntax of an R formula but with all
#' explanatory (right-hand-side) terms wrapped by a function to specify how to fit
#' the explanatory variables in the model. This may look like Yield ~ fx(Year) + rn(Variety)
#' to fit Yield with Year as a fixed factor and Variety as a random factor. 
#' The list of model functions currently available is:
#' fx() (fixed factors), rn() (random factors), rg()
#' (fixed regressions), and rr() (random regressions). 
#' The model functions allow to specify interactions of variables, hierarchies, use of
#' matrices as input data, complex covariance structures, and prior distributions can 
#' be changed to modify standard mixed model in Bayesian shrinkage models.
#'
#' Interactions between fators are specified using the colon, for instance by writing for a 
#' fixed Year-Location interaction fx(Year:Location). Bayz does
#' NOT support automatic expansion with main effects by using the 'star' (Year*Location) or
#' 'forward slash' (Year/Location) syntaxes, hence bayz requires to manually add the desired
#' main effects in the model (but 
#' note that bayz uses the 'forward slash' to specify hierarchical models).
#' Interactions between factors can be specified to any degree.
#' The 'star' syntax can be used to indicate interaction
#' between covariates, like rg(TempSum*Precip), where it is simply interpreted as
#' multiplication. 
#' Interaction between a covariate and a factor is specified using the 'pipe'
#' character as in lme4 models. It can be used to specify fixed nested regressions as 
#' rg(TempSum|Year:Location) to specify regressions on TempSum within each Year-Location, or
#' to speficy random slope models in rr(), for instance rr(TempSum|Variety). 
#'
#' The bayz call has a data argument to specify an input data-frame (the "main data")
#' to contain model variables,
#' but bayz will also search the R environment if variables are not found in the main data.
#' Typically, input in the form of matrices such as large sets of covariates, proportional
#' variance-covariance / correlation / kernel / similarity matrices used in variance models,
#' as well as data for hierarchical models, are not in the main data.
#' Matrices / kernels should be prepared with row-names to match a variable
#' in the data. For matrices / kernels used in variance models, the link is straightforward, 
#' for instance in rn(Variety:Location, V=KG*KE), KG and KE can be a genetic and environmental
#' kinship / kernel matrix, respectively, and KG should have row-names matching Variety levels,
#' while KE should have row-names matching Location levels. 
#' To fit a large set of covariates, the match is specified using a hierachical specification
#' rr(Variety/Metabolites), where Metabolites is then a matrix of covariates which must have
#' row-names matching Variety levels. In both cases, such kernel or covariate matrices
#' must have unique levels, but the main data may have repeated levels and in different order.
#' If one has repeated metabolite data on each Variety, for instance at multiple time-points,
#' then Variety is not the appropriate link to the data, but Variety:Time is.
#'
#' For random effects a variance-covariance structure can be specified using a
#' V= option within the model-term function, for example rn(Variety, V=myGmat).
#' When the fit is for an interaction of factors, the variance specification should
#' expand to include one term for each variable, separated by stars (which should
#' be read as Kronecker products). The variance structure
#' is then built up from a combination of given (proportional) variance-covariance / correlation /
#' kernel matrices and predefined acronyms IDEN, DIAG and VCOV that indicate parameterized
#' matrices (with parameters to be estimated from the data). 
#' Alternatively, the variance structure can be specified as a linear model using V=~, 
#' which is interpreted as use of a log-linear model for the variances. 
#'
#'
#' @param model   A formula describing the model to be fitted
#' @param data    Data frame to collect data from
#' @param VE      Model for the residual variance
#' @param chain   Vector with length, burn-in and skip for the chain to run.
#' @param method  String to indicate analysis method: "Bayes" (full Bayesian, default), "BLUPMC"
#'                (BLUE/BLUP solutions with Monte Carlo to get SD/SE), "BLUP" (BLUE/BLUP solutions)
#' @param init    Initialisation/starting values: supply output from a previous bayz run and bayz
#'                will pick up estimates from a previous run to start a new chain.
#' @param verbose Integer to regulate printing to R console: 0(quiet), 1(some), >2(more)
#' @param ...     Additional parameters passed onto the Model function.
#'
#' @return A fitted bayz model
#' @import stats
#' @export
#'
#' @useDynLib BayzR, .registration = TRUE
#' @importFrom Rcpp sourceCpp
bayz <- function(model, VE="", data=NULL, chain=c(0,0,0), method="", verbose=1, init=NULL, ...){
	if(!inherits(model,"formula")) {
		stop("The first argument is not a valid formula")
	}
    if( !(class(VE)=="character" | class(VE)=="formula") ) {
        stop("VE must be given as a string or formula")
    }
    if( class(method)!="character" ) {
        stop("method must be given as a string")
    }
    if (is.null(data)){
        stop("The data= argument is missing")
    }
    if(class(VE)=="formula") {
        VE=deparse(VE)
    }
    if(method=="") {
        method="Bayes"
    }
    chain <- as.integer(chain)
    result <- rbayz_cpp(model, VE, data, chain, method, verbose, init)
    class(result) <- "bayz"
    #result[['modelname']] <- fct()
    #result[['modelfunction']] <- deparse(substitute(fct))
    #result[['terms']] <- terms(model)
    return(result)
}
