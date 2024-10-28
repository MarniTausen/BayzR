context("bayz.models")

# Because of the monte carlo sampling it is difficult to test if estimates are correct.
# Most tests are limited to checking if model runs or not and if it has right kind/sizes of output.

test_that("One factor as factor input", {
    my_data <- data.frame(x=as.factor(rep(1:2,50)), y=rnorm(100))
    expect_no_error(bayz(y~fx(x),data=my_data,chain=c(50,5,1), verbose=0))
}
)

test_that("One factor as numeric input", {
    my_data <- data.frame(x=rep(1:2,50), y=rnorm(100))
    expect_no_error(bayz(y~fx(x),data=my_data,chain=c(50,5,1), verbose=0))
}
)

test_that("One factor as character input", {
    my_data <- data.frame(x=rep(c("A","B"),50), y=rnorm(100))
    expect_no_error(bayz(y~fx(x),data=my_data,chain=c(50,5,1), verbose=0))
}
)

test_that("One factor 4 levels", {
    my_data <- data.frame(x=rep(c("A","B","C","D"),25), y=rnorm(100))
    expect_no_error(bayz(y~fx(x),data=my_data,chain=c(50,5,1), verbose=0))
}
)

test_that("Two factors 2 and 5 levels", {
    my_data <- data.frame(x1=rep(c("A","B"),50), x1=rep(c("A","B","C","D","E"),20), y=rnorm(100))
    expect_no_error(bayz(y~fx(x1)+fx(x2),data=my_data,chain=c(50,5,1), verbose=0))
}
)

test_that("Two factors and two covariates", {
    my_data <- data.frame(x1=rep(c("A","B"),50), x1=rep(c("A","B","C","D","E"),20),
                    x3=rnorm(100), x4=rnorm(100), y=rnorm(100))
    expect_no_error(bayz(y~fx(x1)+fx(x2)+rg(x3)+rg(x4),data=my_data,chain=c(50,5,1), verbose=0))
}
)

test_that("Factors and covariates in different orders", {
    my_data <- data.frame(x1=rep(c("A","B"),50), x1=rep(c("A","B","C","D","E"),20),
                    x3=rnorm(100), x4=rnorm(100), y=rnorm(100))
    expect_no_error(bayz(y~+rg(x3)+fx(x2)+fx(x1)+rg(x4),data=my_data,chain=c(50,5,1), verbose=0))
}
)


#test_that("Two fixed effects", {
#    testdat1 = read.table("../testdat1.txt",header=TRUE)
#    testdat1$YR = as.factor(testdat1$YR)
#    testdat1$LC = as.factor(testdat1$LC)
#    fit = bayz(y~fx(YR)+fx(LC),data=testdat1,chain=c(5000,100,10), silent=TRUE)
#    YRestimates = fit$Estimates[fit$Parameters["YR","EstStart"]:fit$Parameters["YR","EstEnd"],]
#    LCestimates = fit$Estimates[fit$Parameters["LC","EstStart"]:fit$Parameters["LC","EstEnd"],]
#    lmresults = summary(lm(y~YR+LC, data=testdat1))$coefficients
#    differences = c( (lmresults[2:6,"Estimate"] - YRestimates[2:6,"postMean"]), (lmresults[7:9,"Estimate"] - LCestimates[2:4,"postMean"]) )
#    expect_lte(sum(abs(differences)), 0.1)
#}
#)

#test_that("Fixed two-way interaction (no main effects)", {
#    testdat1 = read.table("../testdat1.txt",header=TRUE)
#    load("../testdat1.RData")
#    fit = bayz(y~fx(YR:LC),data=testdat1,chain=c(5000,100,10), silent=TRUE)
#    estimates = fit$Estimates[fit$Parameters["YR:LC","EstStart"]:fit$Parameters["YR:LC","EstEnd"],]
#    bayzmean = fit$Estimates[fit$Parameters["mean","EstStart"],"postMean"]
#    lmresults = summary(lm(y~YR:LC))$coefficients
#    differences = (estimates[c(1,17,10,3,19,12),"postMean"] + bayzmean) - (lmresults[c(2,6,10,14,18,22),"Estimate"] + lmresults[1,"Estimate"])
#    expect_lte(sum(abs(differences)), 0.1)
#}
#)

