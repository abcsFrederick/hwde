context("test-hwexact")

library(magrittr)

test_that("Basic hwexact functionality", {
    expect_equal(
    {
        hwexact(68, 28, 4) %>%
            round(3)
    }, 0.515)
})