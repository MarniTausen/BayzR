context("bayz")


test_that("Chain warning", {

    my_data <- data.frame(x=1:20, y=20:1)

    expect_warning(bayz(y ~ fixf(x), data=my_data, silent=TRUE))
    #expect_error(bayz(x ~ y, chain=c(10000, 100, 10)))

})

test_that("Missing inputs", {

    my_data <- data.frame(x=1:20, y=20:1)

    expect_error(bayz())
    #expect_error(bayz(x ~ y, chain=c(10000, 100, 10)))

})

test_that("Error message", {

    my_data <- data.frame(x=1:20, y=20:1)

    expect_error(bayz(y ~ x, data=my_data))
    #expect_silent(bayz(y ~ x, data=my_data, chain=c(10000, 100, 10)))

})

test_that("Wrong input", {

    my_data <- data.frame(x=1:20, y=20:1)
    expect_error(bayz(y ~ x, data=my_data, chain=c(10000,100,10)))

    fit <- bayz(y ~ fixf(x), data=my_data, chain=c(10000,100,10), silent=TRUE)
    ## Run print
    capture.output(print(fit))
})

test_that("Working run", {

    my_data <- data.frame(x=1:20, y=20:1)

    expect_equal(bayz(y ~ fixf(x), data=my_data, chain=c(10000,100,10), silent=TRUE)$nError, 0)
    expect_equal(bayz(y ~ ranf(x), data=my_data, chain=c(10000,100,10), silent=TRUE)$nError, 0)

    capture.output(print(bayz(y ~ fixf(x), data=my_data, chain=c(10000,100,10), silent=TRUE)))

})

test_that("Plotting", {

     my_data <- data.frame(x=1:20, y=20:1)
     fit <- bayz(y ~ fixf(x), data=my_data, chain=c(10000,100,10), silent=TRUE)

     plot(fit)

     expect_equal(TRUE, TRUE)

})
