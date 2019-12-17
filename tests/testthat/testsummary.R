context("summary.bayz")

test_that("Validatation", {
    example_data = data.frame(y = 1:20, x = 20:1)
    fit <- bayz(y ~ ranf(x), data=example_data, chain=c(10000, 100, 10))

    summary_object <- summary(fit)

    expect_gte(summary_object$Heritability$Heritability[2], 0.90)
    expect_lte(summary_object$Heritability$Heritability[1], 0.1)

    capture.output(print(summary_object))

})
