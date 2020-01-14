context("bayz")

test_that("Missing inputs", {

    my_data <- data.frame(x=1:20, y=20:1)

    expect_error(bayz())
    #expect_error(bayz(x ~ y, chain=c(10000, 100, 10)))

})

test_that("Warning messages", {

    my_data <- data.frame(x=1:20, y=20:1)

    expect_warning(bayz(y ~ x, data=my_data))
    #expect_silent(bayz(y ~ x, data=my_data, chain=c(10000, 100, 10)))

})

test_that("Wrong input", {

    my_data <- data.frame(x=1:20, y=20:1)
    fit <- bayz(y ~ x, data=my_data, chain=c(10000,100,10))

    expect_match(fit$Errors[1],
                 "Unknown wrapper function on data column")
    expect_match(fit$Errors[2],
                 "Bayz terminates after model building")
    expect_equal(fit$nError, 2)

    capture.output(print(fit))
})


test_that("Working run", {

    my_data <- data.frame(x=1:20, y=20:1)

    expect_equal(bayz(y ~ fixf(x), data=my_data, chain=c(10000,100,10))$nError, 0)
    expect_equal(bayz(y ~ ranf(x), data=my_data, chain=c(10000,100,10))$nError, 0)

    capture.output(print(bayz(y ~ fixf(x), data=my_data, chain=c(10000,100,10))))

})

test_that("Plotting", {

     my_data <- data.frame(x=1:20, y=20:1)
     fit <- bayz(y ~ fixf(x), data=my_data, chain=c(10000,100,10))

     plot(fit)

     expect_equal(TRUE, TRUE)

})
