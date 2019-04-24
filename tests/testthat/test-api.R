context("test-api")

test_that("API of main functions feels natural", {
  dataset <- data("LFQ_Dataset")
  parameters <- fit_parameters(dataset, design = ~ Timepoint + MSRun + BioReplicate)
  result <- test_diff(parameters, contrast = Timepoint_0h - Timepoint_5h)

  result[order(result$pval), ]
})
