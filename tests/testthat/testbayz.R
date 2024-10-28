context("bayz.interface")

test_that("Missing or invalid formula", {

    my_data <- data.frame(x=rep(1:2,10), y=20:1)

    expect_error(bayz("noformula", data=my_data, verbose=0))

})

test_that("Missing data argument", {

    my_data <- data.frame(x=rep(1:2,10), y=20:1)

    expect_error(bayz(y ~ fx(x), chain=c(10, 1, 1), verbose=0))

})


test_that("Chain missing warning", {

    my_data <- data.frame(x=rep(1:2,10), y=20:1)

    expect_warning(bayz(y ~ fx(x), data=my_data, verbose=0))

})

test_that("Chain setting too few elements", {

    my_data <- data.frame(x=rep(1:2,10), y=20:1)

    expect_error(bayz(y ~ fx(x), data=my_data, chain=c(10,2), verbose=0))

})

test_that("Chain wrong settings", {

    my_data <- data.frame(x=rep(1:2,10), y=20:1)

    expect_error(bayz(y ~ fx(x), data=my_data, chain=c(10,-1,-1), verbose=0))

})

test_that("Chain no output", {

    my_data <- data.frame(x=rep(1:2,10), y=20:1)

    expect_error(bayz(y ~ fx(x), data=my_data, chain=c(10,20,1), verbose=0))

})
test_that("Unknown function on response", {

    my_data <- data.frame(x=rep(1:2,10), y=20:1)

    expect_error(bayz(log(y) ~ fx(x), data=my_data, chain=c(10, 1, 1), verbose=0))

})

test_that("Missing model-function", {

    my_data <- data.frame(x=rep(1:2,10), y=20:1)

    expect_error(bayz(y ~ x, data=my_data, chain=c(10, 1, 1), verbose=0))

})

test_that("Unknown model-function", {

    my_data <- data.frame(x=rep(1:2,10), y=20:1)

    expect_error(bayz(y ~ bla(x), data=my_data, chain=c(10, 1, 1), verbose=0))

})

# There can be more tests on different data inputs (Integer, Numeric, matrix) for different model functions

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
