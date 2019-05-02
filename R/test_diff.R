


#' Main proDA function to identify significant coefficients
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
  diff_var <- feature_parameters(fit)$s2 * c(t(contrast) %*% solve(t(dm) %*% dm) %*% contrast)
  diff_df <- feature_parameters(fit)$df

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

    # if(idx == 1 || idx == 4 || idx == 15 || idx == 35) browser()
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

    beta_init <- c(mean(design(fit) %*% fit$coefficients[idx, ]), rep(0, times=ncol(red_model) - 1))
    red_res <- nlminb(start = c(beta_init), objective = function(beta){
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

    full_model <- design(fit)
    lik_full <- - objective_fnc(y, yo,
                                X = full_model,
                                Xm = full_model[is.na(y), , drop=FALSE],
                                Xo = full_model[!is.na(y), , drop=FALSE],
                                beta = fit$coefficients[idx, ],
                                sigma2, rho, zetastar,
                                location_prior_mean, location_prior_scale,
                                variance_prior_df, variance_prior_scale,
                                location_prior_df, moderate_location, moderate_variance)

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

  level_environment <- new.env()

  for(lvl in levels){
    ind <- indicators[, lvl]
    names(ind) <- levels
    assign(lvl, ind, level_environment)
  }
  if(! is.null(reference_level)){
    assign(reference_level, NA_real_, level_environment)
  }

  res <- eval(cnt_capture, envir= level_environment, parent.frame(n = 1 + direct_call))
  if(! is.numeric(res)){
    if(is.character(res)){
      # If contrast was a string, eval will just spit it out the same way
      res <- eval(parse(text = cnt_capture), envir= level_environment, parent.frame(n = 1 + direct_call))
    }
  }
  if(any(is.na(res))){
    stop("Error the reference_level ", reference_level, " was used in the contrast argument. \n",
         "But it is already absorbed into the intercept, so it is not necessary to subtract it here.")
  }
  res
}


