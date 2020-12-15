context("bayz")

test_that("Invalid formula", {

    my_data <- data.frame(x=1:20, y=20:1)

    expect_error(bayz("noformula", data=my_data, silent=TRUE))

})

#test_that("Chain missing warning", {
#
#    my_data <- data.frame(x=1:20, y=20:1)
#
#    expect_warning(bayz(y ~ fx(x), data=my_data, silent=TRUE))
#    #expect_error(bayz(x ~ y, chain=c(10000, 100, 10)))
#
#})

#test_that("Chain missing settings", {
#
#    my_data <- data.frame(x=as.factor(1:20), y=20:1)
#
#    expect_error(bayz(y ~ fx(x), data=my_data, chain=c(10,2), silent=TRUE))
#
#})

#test_that("Chain wrong settings", {
#
#    my_data <- data.frame(x=as.factor(1:20), y=20:1)
#
#    expect_error(bayz(y ~ fx(x), data=my_data, chain=c(10,-1,-1), silent=TRUE))
#
#})


#test_that("Missing data input", {
#
#    my_data <- data.frame(x=as.factor(1:20), y=20:1)
#
#    expect_error(bayz(y ~ fx(x), chain=c(10000, 100, 10), silent=TRUE))
#
#})

#test_that("Missing model-function", {
#
#    my_data <- data.frame(x=1:20, y=20:1)
#
#    expect_error(bayz(y ~ x, data=my_data, chain=c(10000, 100, 10), silent=TRUE))
#
#})

#test_that("Wrong input", {
#
#    my_data <- data.frame(x=1:20, y=20:1)
#    expect_error(bayz(y ~ x, data=my_data, chain=c(10000,100,10)))
#
#    fit <- bayz(y ~ fx(x), data=my_data, chain=c(10000,100,10), silent=TRUE)
#    ## Run print
#    capture.output(print(fit))
#})

# test below gives error: variable is not a factor, when removing as.factor() for x

test_that("Working run with fx()", {

    my_data <- data.frame(x=as.factor(1:20), y=20:1)

    expect_equal(bayz(y ~ fx(x), data=my_data, chain=c(100,10,5), silent=TRUE)$nError, 0)

    capture.output(print(bayz(y ~ fx(x), data=my_data, chain=c(100,10,5), silent=TRUE)))

})

#test_that("Plotting", {
#
#     my_data <- data.frame(x=1:20, y=20:1)
#     fit <- bayz(y ~ fx(x), data=my_data, chain=c(10000,100,10), silent=TRUE)
#
#     plot(fit)
#
#     expect_equal(TRUE, TRUE)
#
#})
