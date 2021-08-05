context("bayz.models")

test_that("Two fixed effects", {
    testdat1 = read.table("../testdat1.txt",header=TRUE)
    testdat1$YR = as.factor(testdat1$YR)
    testdat1$LC = as.factor(testdat1$LC)
    fit = bayz(y~fx(YR)+fx(LC),data=testdat1,chain=c(5000,100,10), silent=TRUE)
    YRestimates = fit$Estimates[fit$Parameters["YR","EstStart"]:fit$Parameters["YR","EstEnd"],]
    LCestimates = fit$Estimates[fit$Parameters["LC","EstStart"]:fit$Parameters["LC","EstEnd"],]
    lmresults = summary(lm(y~YR+LC, data=testdat1))$coefficients
    differences = c( (lmresults[2:6,"Estimate"] - YRestimates[2:6,"postMean"]), (lmresults[7:9,"Estimate"] - LCestimates[2:4,"postMean"]) )
    expect_lte(sum(abs(differences)), 0.1)
}
)

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

