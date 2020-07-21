


#' Identify differentially abundant proteins
#'
#' The `test_diff()` function is used to test coefficients of a 'proDAFit'
#' object. It provides a Wald test to test individual
#' coefficients and a likelihood ratio F-test to compare the
#' original model with a reduced model. The \code{result_names}
#' method provides a quick overview which coefficients are
#' available for testing.
#'
#' To test if coefficient is different from zero with a Wald
#' test use the \code{contrast} function argument. To test if two
#' models differ with an F-test use the \code{reduced_model}
#' argument. Depending on the test that is conducted, the functions
#' returns slightly different data.frames.
#'
#' The function is designed to follow the principles of the
#' base R test functions (ie. \code{\link[stats]{t.test}} and
#' \code{\link[stats]{wilcox.test}}) and the functions designed
#' for collecting the results of high-throughput  testing
#' (ie. \code{limma::topTable} and \code{DESeq2::results}).
#'
#' @param fit an object of class 'proDAFit'. Usually, this is
#'   produced by calling \code{proDA()}
#' @param contrast an expression or a string specifying which
#'   contrast is tested. It can be a single coefficient (to see
#'   the available options use \code{result_names(fit)}) or any
#'   linear combination of them. The contrast is always compared
#'   against zero. Thus, to find out if two coefficients differ
#'   use \code{coef1 - coef2}. Remember if the coefficient is not
#'   a valid identifier in R, to escape it using back ticks. For 
#'   example if you test the interaction of A and B use 
#'   \code{`A:B`}.
#' @param reduced_model If you don't want to test an individual
#'   coefficient, you can can specify a reduced model and compare
#'   it with the original model using an F-test. This is useful
#'   to find out how a set of parameters affect the goodness of
#'   the fit. If neither a \code{contrast}, nor
#'   a \code{reduced_model} is specified, by default a comparison
#'   with an intercept model (ie. just the average across conditions)
#'   is done. Default: \code{~ 1}.
#' @param alternative a string that decides how the
#'   hypothesis test is done. This parameter is only relevant for
#'   the Wald-test specified using the `contrast` argument.
#'   Default: \code{"two.sided"}
#' @param pval_adjust_method a string the indicates the method
#'   that is used to adjust the p-value for the multiple testing.
#'   It must match the options in \code{\link[stats]{p.adjust}}.
#'   Default: \code{"BH"}
#' @param sort_by a string that specifies the column that is used
#'   to sort the resulting data.frame. Default: \code{NULL} which
#'   means the result is sorted by the order of the input matrix.
#' @param decreasing a boolean to indicate if the order is reversed.
#'   Default: \code{FALSE}
#' @param n_max the maximum number of rows returned by the method.
#'   Default: \code{Inf}
#' @param verbose boolean that signals if the method prints informative
#'   messages. Default: \code{FALSE}.
#'
#' @return
#'   The `result_names()` function returns a character vector.
#'
#'   The `test_diff()` function returns a \code{data.frame} with one row per protein
#'   with the key parameters of the statistical test. Depending what kind of test
#'   (Wald or F test) the content of the `data.frame` differs.
#'
#'   The Wald test, which can considered equivalent to a t-test, returns
#'   a `data.frame` with the following columns:
#'   \describe{
#'     \item{name}{the name of the protein, extracted from the rowname of
#'       the input matrix}
#'     \item{pval}{the p-value of the statistical test}
#'     \item{adj_pval}{the multiple testing adjusted p-value}
#'     \item{diff}{the difference that particular coefficient makes. In
#'       differential expression analysis this value is also called
#'       log fold change, which is equivalent to the difference on the
#'       log scale.}
#'     \item{t_statistic}{the \code{diff} divided by the standard
#'       error \code{se}}
#'     \item{se}{the standard error associated with the \code{diff}}
#'     \item{df}{the degrees of freedom, which describe the amount
#'       of available information for estimating the \code{se}. They
#'       are the sum of the number of samples the protein was observed
#'       in, the amount of information contained in the missing values,
#'       and the degrees of freedom of the variance prior.}
#'     \item{avg_abundance}{the estimate of the average abundance of
#'       the protein across all samples.}
#'     \item{n_approx}{the approximated information available for estimating
#'       the protein features, expressed as multiple of the information
#'       contained in one observed value.}
#'     \item{n_obs}{the number of samples a protein was observed in}
#'   }
#'
#'
#'   The F-test returns a `data.frame` with the following columns
#'   \describe{
#'     \item{name}{the name of the protein, extracted from the rowname of
#'       the input matrix}
#'     \item{pval}{the p-value of the statistical test}
#'     \item{adj_pval}{the multiple testing adjusted p-value}
#'     \item{f_statistic}{the ratio of difference of normalized deviances
#'      from original model and the reduced model, divided by the
#'      standard deviation.}
#'     \item{df1}{the difference of the number of coefficients in the
#'       original model and the number of coefficients in the reduced
#'       model}
#'     \item{df2}{the degrees of freedom, which describe the amount
#'       of available information for estimating the \code{se}. They
#'       are the sum of the number of samples the protein was observed
#'       in, the amount of information contained in the missing values,
#'       and the degrees of freedom of the variance prior.}
#'     \item{avg_abundance}{the estimate of the average abundance of
#'       the protein across all samples.}
#'     \item{n_approx}{the information available for estimating
#'       the protein features, expressed as multiple of the information
#'       contained in one observed value.}
#'     \item{n_obs}{the number of samples a protein was observed in}
#'   }
#'
#'
#' @seealso The contrast argument is inspired by
#'   \code{limma::makeContrasts}.
#'
#' @examples
#'   # "t-test"
#'   syn_data <- generate_synthetic_data(n_proteins = 10)
#'   fit <- proDA(syn_data$Y, design = syn_data$groups)
#'   result_names(fit)
#'   test_diff(fit, Condition_1 - Condition_2)
#'
#'
#'   suppressPackageStartupMessages(library(SummarizedExperiment))
#'   se <- generate_synthetic_data(n_proteins = 10,
#'                                 n_conditions = 3,
#'                                 return_summarized_experiment = TRUE)
#'   colData(se)$age <- rnorm(9, mean=45, sd=5)
#'   colData(se)
#'   fit <- proDA(se, design = ~ group + age)
#'   result_names(fit)
#'   test_diff(fit, "groupCondition_2",
#'             n_max = 3, sort_by = "pval")
#'
#'   # F-test
#'   test_diff(fit, reduced_model = ~ group)
#'
#'
#' @export
test_diff <- function(fit, contrast,
                      reduced_model = ~ 1,
                      alternative = c("two.sided", "greater", "less"),
                      pval_adjust_method = "BH",
                      sort_by = NULL,
                      decreasing = FALSE,
                      n_max = Inf,
                      verbose = FALSE){

  alternative <- match.arg(alternative, c("two.sided", "greater", "less"))

  if(missing(contrast)){
    # Do F test against reduced model
    # Handle the design parameter
    if(is.matrix(reduced_model)){
      red_model_matrix <- reduced_model
    }else if((is.vector(reduced_model) || is.factor(reduced_model)) && length(reduced_model) == ncol(fit)){
      red_model_matrix <- convert_chr_vec_to_model_matrix(reduced_model, fit$reference_level)
    }else if(inherits(reduced_model,"formula")){
      red_model_matrix <- convert_formula_to_model_matrix(reduced_model, colData(fit), fit$reference_level)
    }else{
      stop(paste0("red_model_matrix argment of class ", class(reduced_model), " is not supported. Please ",
                  "specify a `model_matrix`, a `character vector`, or a `formula`."))
    }
    if(verbose){
      message("F test: comparing the reduced model with the original model.")
    }
    rownames(red_model_matrix) <- colnames(fit)
    check_valid_model_matrix(red_model_matrix, fit)
    if(ncol(red_model_matrix) >= ncol(design(fit))){
      stop("Reduced model is larger or has the same size as the reference model.")
    }

    f_res <- run_nested_model_comparison(fit, red_model_matrix, verbose)



    res <- data.frame(name = rownames(fit),
                      pval = f_res[, "pval"],
                      adj_pval = p.adjust(f_res[, "pval"], method = pval_adjust_method),
                      # deviance_full = f_res[, "deviance_full"],
                      # deviance_reduced = f_res[, "deviance_reduced"],
                      f_statistic = f_res[, "f_statistic"],
                      df1 = f_res[, "df1"],
                      df2 = f_res[, "df2"],
                      avg_abundance = rowMeans(predict(fit, type="response")),
                      n_approx = feature_parameters(fit)$n_approx,
                      n_obs = feature_parameters(fit)$n_obs,
                      stringsAsFactors = FALSE)
  }else{
    if(any(reduced_model != formula(~ 1))){
      stop("You can only specify the contrast argument or the reduced_model argument.")
    }
    cntrst <- parse_contrast(contrast, levels =  colnames(design(fit)),
                   reference_level = reference_level(fit),
                   direct_call = FALSE)
    if(verbose){
      message("Wald test: comparing ", deparse(substitute(contrast)), " with zero.")
    }

    wald_res <- run_wald_parameter_test(fit, cntrst, alternative)

    res <- data.frame(name = rownames(fit),
               pval = wald_res[, "pval"],
               adj_pval = p.adjust(wald_res[, "pval"], method = pval_adjust_method),
               diff = wald_res[, "diff"],
               t_statistic = wald_res[, "t_statistic"],
               se = wald_res[, "se"],
               df = wald_res[, "df"],
               avg_abundance = rowMeans(predict(fit, type="response")),
               n_approx = feature_parameters(fit)$n_approx,
               n_obs = feature_parameters(fit)$n_obs,
               stringsAsFactors = FALSE)

  }

  rownames(res) <- NULL
  res <- if(is.null(sort_by)){
    res
  }else{
    res[order(res[[sort_by]], decreasing = decreasing), ]
  }
  res[seq_len(min(nrow(res), n_max)), ]

}


run_wald_parameter_test <- function(fit, contrast, alternative){
  dm <- design(fit)
  diff <- coefficients(fit) %*% contrast
  diff_var <- vapply(coefficient_variance_matrices(fit), function(mat) c(t(contrast) %*% mat %*% contrast),
                     FUN.VALUE = 0.0)
  diff_df <- nrow(dm) - ncol(dm)

  t_stat <- diff / sqrt(diff_var)
  pval <- pt(t_stat, df=diff_df)
  if(alternative == "two.sided"){
    pval <- 2 * pmin(pval, 1-pval)
  }else if(alternative == "greater"){
    pval <- 1 - pval
  }
  res <- cbind(diff,
        sqrt(diff_var),
        diff_df,
        t_stat,
        pval)
  colnames(res) <- c("diff", "se", "df",
                     "t_statistic", "pval")
  res

}



run_nested_model_comparison <- function(fit, red_model, verbose=FALSE){

  location_prior_mean <- fit$hyper_parameters$location_prior_mean
  location_prior_scale <- fit$hyper_parameters$location_prior_scale
  variance_prior_df <- fit$hyper_parameters$variance_prior_df
  variance_prior_scale <- fit$hyper_parameters$variance_prior_scale
  location_prior_df <- fit$hyper_parameters$location_prior_df

  moderate_location <- ! is.na(location_prior_mean)
  moderate_variance <- ! is.na(variance_prior_scale)

  # Fit the loglikelihood with the fixed \phi
  # using the reduced model
  res_mat <- mply_dbl(seq_len(nrow(fit)), function(idx){

    full_model <- design(fit)
    y <- fit$abundances[idx, ]
    yo <- y[! is.na(y)]

    rho <- fit$hyper_parameters$dropout_curve_position[is.na(y)]
    zeta <- fit$hyper_parameters$dropout_curve_scale[is.na(y)]

    sigma2 <- fit$feature_parameters$s2[idx]
    zetastar <- zeta * sqrt(1 + sigma2/zeta^2)
    df_full <- fit$feature_parameters$df[idx]
    if(is.na(sigma2)){
      return(c(f_statistic = NA, pval = NA,
               df1 = ncol(full_model) - ncol(red_model), df2 = df_full,
               deviance_full = NA, deviance_reduced = NA))
    }

    nl_res_full <- nlminb(start = fit$coefficients[idx, ], objective = function(beta){
      - objective_fnc(y, yo,
                      X = full_model,
                      Xm = full_model[is.na(y), , drop=FALSE],
                      Xo = full_model[!is.na(y), , drop=FALSE],
                      beta, sigma2, rho, zetastar,
                      location_prior_mean, location_prior_scale,
                      variance_prior_df, variance_prior_scale,
                      location_prior_df, moderate_location, moderate_variance)
    })
    if(nl_res_full$convergence >= 2){
      if(verbose){
        warning("Model didn't not properly converge\n")
        warning(nl_res_full$message, "\n")
        warning(y,"\n")
      }
      lik_full <- NA
    }else{
      lik_full <- nl_res_full$objective
    }


    beta_init <- vapply(colnames(red_model), function(n){
      if(is.na(nl_res_full$par[n])){
        0
      }else{
        nl_res_full$par[n]
      }
    }, FUN.VALUE = 0.0)
    red_res <- nlminb(start = beta_init, objective = function(beta){
       - objective_fnc(y, yo,
                       X = red_model,
                       Xm = red_model[is.na(y), , drop=FALSE],
                       Xo = red_model[!is.na(y), , drop=FALSE],
                       beta, sigma2, rho, zetastar,
                       location_prior_mean, location_prior_scale,
                       variance_prior_df, variance_prior_scale,
                       location_prior_df, moderate_location, moderate_variance)
    })
    if(red_res$convergence >= 2){
      if(verbose){
        warning("Model didn't not properly converge\n")
        warning(red_res$message, "\n")
        warning(y,"\n")
      }
      lik_red <- NA
    }else{
      lik_red <- red_res$objective
    }

    deviance_red <- lik_red * 2 * sigma2
    deviance_full <- lik_full * 2 * sigma2

    LR_stat <- (deviance_red - deviance_full) / (ncol(full_model) - ncol(red_model)) / sigma2
    pval <- pf(LR_stat, df1 = ncol(full_model) - ncol(red_model), df2 = df_full,
               lower.tail = FALSE)
    c(f_statistic = LR_stat, pval = pval,
      df1 = ncol(full_model) - ncol(red_model), df2 = df_full,
      deviance_full = deviance_full, deviance_reduced = deviance_red)
  }, ncol = 6)

}



parse_contrast <- function(contrast, levels, reference_level = NULL, direct_call=TRUE){

  env <- if(direct_call){
    environment()
  }else{
    parent.frame()
  }

  if(missing(contrast)){
    stop("No contrast argument was provided! The option is any linear combination of:\n",
         paste0(levels, collapse = ", "))
  }
  cnt_capture <- substitute(contrast, env = env)


  stopifnot(! is.null(levels))
  if(is.factor(levels)){
    levels <- levels(levels)
  }else if(! is.character(levels)){
    stop("levels must be either a character vector or a factor")
  }

  indicators <- diag(nrow=length(levels))
  rownames(indicators) <- levels
  colnames(indicators) <- levels

  level_environment <- new.env(parent = parent.frame(n = 1 + (!direct_call)))

  for(lvl in levels){
    ind <- indicators[, lvl]
    names(ind) <- levels
    assign(lvl, ind, level_environment)
  }
  if(! is.null(reference_level)){
    assign(reference_level, NA_real_, level_environment)
  }
  res <- eval(cnt_capture, envir= level_environment)
  if(! is.numeric(res)){
    if(is.character(res)){
      # If contrast was a string, eval will just spit it out the same way
      res <- eval(parse(text = res), envir= level_environment)
    }
  }
  if(any(is.na(res))){
    stop("Error the reference_level ", reference_level, " was used in the contrast argument. \n",
         "But it is already absorbed into the intercept, so it is not necessary to subtract it here.")
  }
  res
}


#' @rdname test_diff
#' @export
setMethod("result_names", signature = "proDAFit", function(fit){
  names <- colnames(coefficients(fit))
  # Check if names contains any illegal characters
  # this regex is approx. correct, but false positives are not so problematic
  valid <- grepl("^[_.]?[[:alpha:]]+[[:alnum:]_.]*$", names)
  ifelse(valid, names, paste0("`", names, "`"))
})



