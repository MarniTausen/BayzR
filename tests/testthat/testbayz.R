context("bayz")

test_that("Missing input", {

    my_data <- data.frame(x=1:20, y=20:1)

    expect_error(bayz())
    expect_error(bayz(x ~ y))

})
