


#' proDA: Identify differentially abundant proteins in label-free mass spectrometry
#'
#'
#'
#'
#' @docType  package
#' @name proDA_package
NULL

#' @import methods BiocGenerics SummarizedExperiment
NULL













#' Main proDA function to determine the hyper and protein parameters
#'
#'
#' @export
proDA <- function(data, design=~ 1,
                  col_data = NULL,
                  reference_level = NULL,
                  data_is_log_transformed = TRUE,
                  location_prior_df = 3,
                  moderate_location = TRUE,
                  moderate_variance = TRUE,
                  max_iter = 20,
                  epsilon = 1e-3,
                  verbose=FALSE, ...){

  # Validate Data
  stopifnot(is.matrix(data) || is(data, "SummarizedExperiment"))
  n_samples <- ncol(data)
  n_rows <- nrow(data)


  # Handle the design parameter
  if(is.matrix(design)){
    model_matrix <- design
    design_formula <- NULL
  }else if((is.vector(design) || is.factor(design)) && length(design) == n_samples){
    model_matrix <- convert_chr_vec_to_model_matrix(design, reference_level)
    design_formula <- NULL
  }else if(inherits(design,"formula")){
    if(design == formula(~ 1) && is.null(col_data)){
      col_data <- as.data.frame(matrix(numeric(0), nrow=n_samples))
    }
    compl_col_data <- if(is(data, "SummarizedExperiment")){
      if(is.null(col_data)) colData(data)
      else cbind(col_data, colData(data))
    }else{
      col_data
    }
    model_matrix <- convert_formula_to_model_matrix(design, compl_col_data, reference_level)
    design_formula <- design
  }else{
    stop(paste0("design argment of class ", class(design), " is not supported. Please ",
                "specify a `model_matrix`, a `character vector`, or a `formula`."))
  }
  rownames(model_matrix) <- colnames(data)
  check_valid_model_matrix(model_matrix, data)


  # Extract the raw data matrix
  if(is.matrix(data)){
    if(! data_is_log_transformed){
      data <- log2(data)
    }
    data_mat <- data
  }else if(inherits(data, "SummarizedExperiment")){
    if(! data_is_log_transformed){
      assay(data) <- log2(assay(data))
    }
    data_mat <- assay(data)

  }else{
    stop("data of tye ", class(data), " is not supported.")
  }

  fit_result <- fit_parameters_loop(data_mat, model_matrix,
                                    location_prior_df = location_prior_df,
                                    moderate_location = moderate_location,
                                    moderate_variance = moderate_variance,
                                    max_iter = max_iter,
                                    epsilon = epsilon,
                                    verbose = verbose)

  feat_df <- as.data.frame(mply_dbl(fit_result$feature_parameters, function(f){
    unlist(f[-1])
  }, ncol = 5))
  feat_df$rss <- fit_result$rss
  coef_mat <- mply_dbl(fit_result$feature_parameters, function(f){
    f$coefficients
  }, ncol=ncol(model_matrix))
  colnames(coef_mat) <- names(fit_result$feature_parameters[[1]]$coefficients)

  proDAFit(data, col_data,
           dropout_curve_position = fit_result$hyper_parameters$dropout_curve_position,
           dropout_curve_scale = fit_result$hyper_parameters$dropout_curve_scale,
           feature_parameters = feat_df,
           coefficients = coef_mat,
           design_matrix = model_matrix,
           design_formula = design_formula,
           reference_level = reference_level,
           location_prior_mean = fit_result$hyper_parameters$location_prior_mean,
           location_prior_scale = fit_result$hyper_parameters$location_prior_scale,
           location_prior_df = location_prior_df,
           variance_prior_scale = fit_result$hyper_parameters$variance_prior_scale,
           variance_prior_df = fit_result$hyper_parameters$variance_prior_df,
           convergence = fit_result$convergence, ...)
}






fit_parameters_loop <- function(Y, model_matrix, location_prior_df,
                                moderate_location, moderate_variance,
                                max_iter, epsilon, verbose=FALSE){
  if(verbose){
    message("Fitting the hyper-parameters for the probabilistic dropout model.")
  }
  # Initialization
  n_samples <- ncol(Y)


  Y_compl <- Y
  Y_compl[is.na(Y)] <- rnorm(sum(is.na(Y)), mean=quantile(Y, probs=0.1, na.rm=TRUE), sd=sd(Y,na.rm=TRUE)/5)
  res_init <- lapply(seq_len(nrow(Y)), function(i){
    pd_lm.fit(Y_compl[i, ], model_matrix,
              dropout_curve_position = rep(NA, n_samples),
              dropout_curve_scale =rep(NA, n_samples))
  })
  Pred_init <- msply_dbl(res_init, function(x) x$coefficients) %*% t(model_matrix)
  s2_init <-  vapply(res_init, function(x) x[["s2"]], 0.0)
  df_init <- vapply(res_init, function(x) x[["df"]], 0.0)
  if(moderate_location){
    lp <- location_prior(model_matrix, Pred = Pred_init, s2 = s2_init)
    mu0 <- lp$mu0
    sigma20 <- lp$sigma20
  }else{
    mu0 <- NA_real_
    sigma20 <- NA_real_
  }
  dc <- dropout_curves(Y, model_matrix, Pred_init, s2_init)
  rho <- dc$rho
  zetainv <- dc$zetainv
  if(moderate_variance){
    vp <- variance_prior(s2_init, df_init)
    tau20 <- vp$tau20
    df0_inv <- vp$df0_inv
  }else{
    tau20 <- NA_real_
    df0_inv <- NA_real_
  }

  last_round_params <- list(mu0, sigma20, rho, zetainv, tau20, df0_inv)
  converged <- FALSE
  iter <- 1
  while(! converged && iter <= max_iter){
    if(verbose){
      message(paste0("Starting iter: ", iter))
    }

    res_unreg <- lapply(seq_len(nrow(Y)), function(i){
      pd_lm.fit(Y[i, ], model_matrix,
                dropout_curve_position = rho, dropout_curve_scale = 1/zetainv)
    })
    if(moderate_location || moderate_variance){
      res_reg <- lapply(seq_len(nrow(Y)), function(i){
        pd_lm.fit(Y[i, ], model_matrix,
                  dropout_curve_position = rho, dropout_curve_scale = 1/zetainv,
                  location_prior_mean = mu0, location_prior_scale = sigma20,
                  variance_prior_scale = tau20, variance_prior_df = 1/df0_inv,
                  location_prior_df = location_prior_df)
      })
    }else{
      res_reg <- res_unreg
    }

    Pred_unreg <- msply_dbl(res_unreg, function(x) x$coefficients) %*% t(model_matrix)
    Pred_reg <- msply_dbl(res_reg, function(x) x$coefficients) %*% t(model_matrix)
    s2_unreg <-  vapply(res_unreg, function(x) x[["s2"]], 0.0)
    s2_reg <-  vapply(res_reg, function(x) x[["s2"]], 0.0)
    df_unreg <-vapply(res_unreg, function(x) x[["df"]], 0.0)

    if(moderate_location){
      lp <- location_prior(model_matrix, Pred = Pred_unreg,
                           # mu0 = mean( Pred_reg, na.rm=TRUE),
                           s2 = s2_unreg)
      mu0 <- lp$mu0
      sigma20 <- lp$sigma20
    }
    dc <- dropout_curves(Y, model_matrix, Pred_reg, s2_reg)
    rho <- dc$rho
    zetainv <- dc$zetainv
    if(moderate_variance){
      vp <- variance_prior(s2_unreg, df_unreg)
      tau20 <- vp$tau20
      df0_inv <- vp$df0_inv
    }

    error <- sum(mapply(function(new, old) {
      sum(new - old, na.rm=TRUE)/length(new)
    }, list(mu0, sigma20, rho, zetainv, tau20, df0_inv), last_round_params)^2)
    if (error < epsilon) {
      if(verbose){
        message("Converged!")
      }
      converged <- TRUE
    }

    last_round_params <- list(mu0, sigma20, rho, zetainv, tau20, df0_inv)
    if(verbose){
      print(last_round_params)
      print(paste0("Error: ", sprintf("%.2g", error)))
    }
    iter <- iter + 1
  }


  convergence <- list(successful = converged, iterations = iter-1, error = error)
  names(last_round_params) <- c("location_prior_mean", "location_prior_scale",
                                "dropout_curve_position", "dropout_curve_scale",
                                "variance_prior_scale", "variance_prior_df")
  last_round_params[["dropout_curve_scale"]] <- 1/last_round_params[["dropout_curve_scale"]]
  last_round_params[["variance_prior_df"]] <- 1/last_round_params[["variance_prior_df"]]

  rss_unreg <- vapply(res_unreg, function(x) x[["rss"]], 0.0)
  list(hyper_parameters = last_round_params,
       convergence = convergence,
       feature_parameters = res_reg,
       rss = rss_unreg)

}




variance_prior <- function(s2, df){
  stopifnot(length(s2) == length(df) || length(df) == 1)
  if(any(df <= 0, na.rm=TRUE)){
    stop(paste0("All df must be positive. ", paste0(which(df < 0), collapse=", "), " are not."))
  }
  if(any(s2 <= 0, na.rm=TRUE)){
    stop(paste0("All s2 must be positive. ", paste0(which(s2 < 0), collapse=", "), " are not."))
  }
  opt_res <- optim(par=c(tau=1, dfinv=1), function(par){
    if(par[1] < 0 || par[2] < 0 ) return(Inf)
    -sum(df(s2/par[1], df1=df, df2=1/par[2], log=TRUE) - log(par[1]), na.rm=TRUE)
  })
  if(opt_res$convergence != 0){
    warning("Model didn't not properly converge\n")
    print(opt_res)
  }
  list(tau20 = unname(opt_res$par[1]), df0_inv=unname(opt_res$par[2]))
}



location_prior <- function(X, Pred, s2,
                           mu0 = mean(Pred, na.rm=TRUE),
                           min_var = 0, max_var = 1e3){
  if(any(s2 <= 0, na.rm=TRUE)){
    stop(paste0("All s2 must be positive. ", paste0(which(s2 < 0), collapse=", "), " are not."))
  }

  pred <- c(Pred)
  larger_than_mu0 <- which(pred > mu0)
  pred <- (pred - mu0)[larger_than_mu0]
  pred_var <- c(s2 %*% t(sapply(seq_len(nrow(X)), function(i) t(X[i,]) %*% solve(t(X) %*% X) %*% X[i,])))[larger_than_mu0]

  objective_fun <- function(A){
    sum((pred^2 - pred_var) / (2 * (A + pred_var)^2), na.rm=TRUE) / sum(1/(2 * (A + pred_var)^2), na.rm=TRUE) - A
  }

  if(sign(objective_fun(min_var)) == sign(objective_fun(max_var))){
    root <- list(root = sum(pred^2) / length(pred))
  }else{
    root <- uniroot(objective_fun, lower=min_var, upper=max_var)
  }

  list(mu0 = mu0, sigma20 = root$root)
}





dropout_curves <- function(Y, X, Pred, s2){
  if(any(s2 <= 0, na.rm=TRUE)){
    stop(paste0("All s2 must be positive. ", paste0(which(s2 < 0), collapse=", "), " are not."))
  }
  n_samples <- nrow(X)
  Pred_var <- s2 %*% t(sapply(seq_len(nrow(X)), function(i) t(X[i,]) %*% solve(t(X) %*% X) %*% X[i,]))

  mu0 <- median(Pred, na.rm=TRUE)
  sigma20 <- median((mu0 - Pred)^2, na.rm=TRUE)

  rho <- rep(NA, n_samples)
  zetainv <- rep(NA, n_samples)
  for(colidx in seq_len(n_samples)){
    y <- Y[, colidx]
    yo <- y[! is.na(y)]
    predm <- Pred[is.na(y), colidx]
    pred_var_m <- Pred_var[is.na(y), colidx]
    if(any(is.na(y))){
      opt_res <- optim(par=c(rho=0, zetainv=-1/5), function(par){
        if(par[2]  >= 0) return(Inf)
        val <- 0 +
          dnorm(par[1], mu0, sd=sqrt(sigma20), log=TRUE) +
          sum(invprobit(yo, par[1], 1/par[2], log=TRUE, oneminus = TRUE), na.rm=TRUE) +
          sum(invprobit(predm, par[1], sign(par[2]) * sqrt(1/par[2]^2 +  pred_var_m), log=TRUE), na.rm=TRUE)
        -val
      })
      if(opt_res$convergence != 0){
        warning("Dropout curve estimation did not properly converge")
      }
      rho[colidx] <- opt_res$par[1]
      zetainv[colidx] <- opt_res$par[2]
    }else{
      rho[colidx] <- NA_real_
      zetainv[colidx] <- NA_real_
    }
  }

  list(rho=rho, zetainv=zetainv)
}







check_valid_model_matrix <- function(matrix, data){
  stopifnot(is.matrix(matrix))
  stopifnot(nrow(matrix) == ncol(data))
}



convert_chr_vec_to_model_matrix <- function(design, reference_level){
  if(! is.factor(design)){
    design_fct <- as.factor(design)
  }else{
    design_fct <- design
  }
  if(is.null(reference_level)){
    helper_df <- data.frame(x_ = design_fct)
    mm <- model.matrix(~ x_ - 1, helper_df)
    colnames(mm) <- sub("^x_", "", colnames(mm))
  }else{
    design_fct <- relevel(design_fct, ref = reference_level)
    helper_df <- data.frame(x_ = design_fct)
    mm <- model.matrix(~ x_ + 1, helper_df)
    colnames(mm)[-1] <- paste0(sub("^x_", "", colnames(mm)[-1]), "_vs_", reference_level)
  }
  colnames(mm)[colnames(mm) == "(Intercept)"] <- "Intercept"
  mm
}


convert_formula_to_model_matrix <- function(formula, col_data, reference_level=NULL){
  if(! is.null(reference_level)){
    has_ref_level <- vapply(col_data, function(x){
      any(x == reference_level)
    }, FUN.VALUE = FALSE)
    if(all(has_ref_level == FALSE)){
      stop("None of the columns contains the specified reference_level.")
    }
    col_data[has_ref_level] <- lapply(col_data[has_ref_level], relevel, ref = reference_level)
  }
  mm <- model.matrix(formula, col_data)
  colnames(mm)[colnames(mm) == "(Intercept)"] <- "Intercept"
  mm
}










