context("test-test_diff")

test_that("parse_contrast works", {

  expect_equal(parse_contrast(A, levels = LETTERS[1:5]),
               c(A=1, B=0, C=0, D=0, E=0))

  expect_equal(parse_contrast(A - B, levels = LETTERS[1:5]),
               c(A=1, B=-1, C=0, D=0, E=0))

  expect_equal(parse_contrast("A - B", levels = LETTERS[1:5]),
               c(A=1, B=-1, C=0, D=0, E=0))

  expect_equal(parse_contrast(A * 3, levels = LETTERS[1:5]),
               c(A=3, B=0, C=0, D=0, E=0))

  expect_equal(parse_contrast(B, levels = c("Intercept", LETTERS[2:5]), reference_level = "A"),
               c(Intercept = 0, B=1, C=0, D=0, E=0))

  expect_error(parse_contrast(B - A, levels = c("Intercept", LETTERS[2:5]), reference_level = "A"))

  expect_equal(parse_contrast("Intercept", levels = c("Intercept")),
               c(Intercept = 1))

})
