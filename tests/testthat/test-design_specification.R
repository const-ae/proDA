context("test-design_specification")

test_that("Character to model matrix", {
  design_vec <- rep(c("A", "B"), each = 5)
  data <- matrix(0, nrow=100, ncol=10)

  mm <- convert_chr_vec_to_model_matrix(design_vec, reference_class = NULL)
  expect_equal(c(mm), rep(c(1,0,0,1), each=5))
  expect_equal(colnames(mm), c("A", "B"))

  design_vec <- rep(c("A", "B", "C"), each = 3)
  data <- matrix(0, nrow=100, ncol=9)

  mm <- convert_chr_vec_to_model_matrix(design_vec, reference_class = "C")
  expect_equal(c(mm), c(rep(1,times=9), rep(c(1,0,0), each=3), rep(c(0,1,0), each=3)))
  expect_equal(colnames(mm), c("Intercept","A_vs_C", "B_vs_C"))

})


test_that("Formula to model_matrix", {
  col_data <- data.frame(f1 = factor(rep(LETTERS[1:5],each=2), levels = LETTERS[1:7]),
                         f2 = factor(c("Good", "Neutral", "Neutral", "Bad", "Bad",
                                       "Bad", "Good", "Bad", "Neutral", "Bad"),
                                     levels = c("Bad", "Neutral", "Good"), ordered=TRUE),
                         f3 = factor(rep(c("hello", "world"),times=5)),
                         f4 = factor(rep(c("hello", "world", "foo", "bar", "foobar"),times=2)),
                         c1 = sample(rep(c("ABC", "xyz"), each=5)),
                         num = rnorm(10),
                         num2 = rnorm(10),
                         stringsAsFactors = FALSE)
  data <- matrix(0, nrow=100, ncol=10)

  mm <- convert_formula_to_model_matrix(~ f1 , col_data)
  expect_equal(colnames(mm), c("Intercept", paste0("f1", LETTERS[2:7])))
  expect_equal(unname(colSums(mm)), c(10, 2,2,2,2,0,0))


  mm <- convert_formula_to_model_matrix(~ f3 + f4, col_data, reference_class = "world")
  expect_equal(colnames(mm), c("Intercept", "f3hello", "f4bar", "f4foo", "f4foobar", "f4hello"))

  f3_mod <- relevel(col_data$f3, ref = "world")
  f4_mod <- relevel(col_data$f4, ref = "world")
  expect_equal(c(mm), c(model.matrix(~ f3_mod + f4_mod)))

})



