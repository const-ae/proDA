


#' proDA: Identify differentially abundant proteins in label-free mass spectrometry
#'
#'
#'
#'
#' @docType  package
#' @name proDA_package
NULL

#' @import stats methods
#' @importFrom SummarizedExperiment SummarizedExperiment rowData rowData<-
#'   colData colData<- mcols mcols<- assay assay<- assayNames assayNames<-
#'   assays assays<- rbind cbind
#' @importFrom BiocGenerics design
#' @importFrom utils .DollarNames
NULL


## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("Condition1", "Condition2"))











#' Main function to determine the hyper and protein parameters
#'
#' The function fits a linear probabilistic dropout model and infers
#' the hyper-parameters for the location prior, the variance prior,
#' and the dropout curves. In addition it infers for each protein
#' the coefficients that best explain the observed data and the
#' associated uncertainty.
#'
#' By default, the method is moderating the locations and the variance
#' of each protein estimate. The variance moderation is fairly standard
#' in high-throughput experiments and can boost the power to detect
#' differentially abundant proteins. The location moderation is important
#' to handle extreme cases where in one conditio a protein is not observed
#' in any sample. In addition it can help to get more precise estimates
#' of the difference between conditions. Unlike 'DESeq2', which moderates
#' the coefficient estimates (ie. the "betas") to be centered around zero,
#' 'proDA' penalizes predicted intensities that strain far from the other
#' observed intensities.
#'
#' @param data a matrix like object (\code{matrix()} or
#'   \code{SummarizedExperiment()}) with the one column per
#'   sample and one row per protein. Missing values should be
#'   coded \code{NA}.
#' @param design a specification of the experimental design that
#'   is used to fit the linear model. It can be a \code{model.matrix()}
#'   with one row for each sample and one column for each
#'   coefficient. It can also be a formula with the entries refering
#'   to global objects, columns in the \code{col_data} argument or
#'   columns in the \code{colData(data)} if data is a
#'   \code{SummarizedExperiment}. Thirdly, it can be a vector that
#'   for each sample specifies the condition of that sample.
#'   Default: \code{~ 1}, which means that all samples are treated
#'   as if they are in the same condition.
#' @param col_data a data.frame with one row for each sample in
#'   \code{data}. Default: \code{NULL}
#' @param reference_level a string that specifies which level in a
#'   factor coefficient is used for the intercept.  Default:
#'   \code{NULL}
#' @param data_is_log_transformed the raw intensities from mass
#'   spectrometry experiments have a linear mean-variance relation.
#'   This is undesirable and can be removed by working on the log
#'   scale. The easiest way to find out if the data is already log-
#'   transformed is to see if the intensities are in the range of
#'   0 to 100 in which case they are transformed or if they rather
#'   are between 1e5 to 1e12, in which case they need to be
#'   transformed. Default: \code{TRUE}
#' @param moderate_location,moderate_variance boolean values
#'   to indicate if the location and the variances are
#'   moderated. Default: \code{TRUE}
#' @param location_prior_df the number of degrees of freedom used
#'   for the location prior. A large number (> 30) means that the
#'   prior is approximately Normal. Default: \code{3}
#' @param max_iter the maximum of iterations \code{proDA()} tries
#'   to converge to the hyper-parameter estimates. Default:
#'   \code{20}
#' @param epsilon if the remaining error is smaller than \code{epsilon}
#'   the model has converged. Default: \code{1e-3}
#' @param verbose boolean that signals if the method prints informative
#'   messages. Default: \code{FALSE}
#' @param ... additional parameters for the construction of the
#'   'proDAFit' object
#'
#' @return
#'   An object of class 'proDAFit'. The object contains information
#'   on the hyper-parameters and feature parameters, the convergence,
#'   the experimental design etc. Internally, it is a sub-class of
#'   \code{SummarizedExperiment} which means the object is subsettable.
#'
#'
#' @examples
#'
#' # Quick start
#' set.seed(1)
#' library(proDA)
#' syn_data <- generate_synthetic_data(n_proteins = 10)
#' fit <- proDA(syn_data$Y, design = syn_data$groups)
#' fit
#' result_names(fit)
#' test_diff(fit, Condition_1 - Condition_2)
#'
#' # SummarizedExperiment
#' se <- generate_synthetic_data(n_proteins = 10,
#'                      return_summarized_experiment = TRUE)
#' se
#' proDA(se, design = ~ group)
#'
#' # Design using model.matrix()
#' data_mat <- matrix(rnorm(5 * 10), nrow=10)
#' colnames(data_mat) <- paste0("sample", 1:5)
#' annotation_df <- data.frame(names = paste0("sample", 1:5),
#'                      condition = c("A", "A", "A", "B", "B"),
#'                      age = rnorm(5, mean=40, sd=10))
#'
#' design_mat <- model.matrix(~ condition + age,
#'                            data=annotation_df)
#' design_mat
#' proDA(data_mat, design_mat, col_data = annotation_df)
#'
#'
#' @export
proDA <- function(data, design=~ 1,
                  col_data = NULL,
                  reference_level = NULL,
                  data_is_log_transformed = TRUE,
                  moderate_location = TRUE,
                  moderate_variance = TRUE,
                  location_prior_df = 3,
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
    if(any(! is.na(data) & data == 0)){
      warning(paste0("The data contains", sum(! is.na(data) & data == 0) ," exact zeros. ",
                     "Replacing them with 'NA's."))
      data[! is.na(data) & data == 0] <- NA
    }
    if(! data_is_log_transformed){
      data <- log2(data)
    }
    data_mat <- data
  }else if(inherits(data, "SummarizedExperiment")){
    if(any(! is.na(assay(data)) & assay(data) == 0)){
      warning(paste0("The data contains", sum(! is.na(assay(data)) & assay(data) == 0) ," exact zeros. ",
                     "Replacing them with 'NA's."))
      assay(data)[! is.na(assay(data)) & assay(data) == 0] <- NA
    }
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
  }, ncol = 4))
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
  error <- NA
  res_reg <- res_init
  res_unreg <- res_init
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
                           mu0 = median( Pred_reg, na.rm=TRUE),
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
      log_parameters(last_round_params)
      message(paste0("Error: ", sprintf("%.2g", error)), "\n")
    }
    iter <- iter + 1
  }


  convergence <- list(successful = converged, iterations = iter-1, error = error)
  names(last_round_params) <- c("location_prior_mean", "location_prior_scale",
                                "dropout_curve_position", "dropout_curve_scale",
                                "variance_prior_scale", "variance_prior_df")
  last_round_params[["dropout_curve_scale"]] <- 1/last_round_params[["dropout_curve_scale"]]
  last_round_params[["variance_prior_df"]] <- 1/last_round_params[["variance_prior_df"]]

  list(hyper_parameters = last_round_params,
       convergence = convergence,
       feature_parameters = res_reg)

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
                           mu0 = median(Pred, na.rm=TRUE),
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
  if(sigma20 == 0){
    sigma20 <- 5 # Not ideal. But what else can I do...
  }


  rho <- rep(NA, n_samples)
  zetainv <- rep(NA, n_samples)
  for(colidx in seq_len(n_samples)){
    y <- Y[, colidx]
    yo <- y[! is.na(y)]
    predm <- Pred[is.na(y), colidx]
    pred_var_m <- Pred_var[is.na(y), colidx]
    if(any(is.na(y))){
      opt_res <- optim(par=c(rho=mu0, zetainv=-1/sqrt(sigma20)), function(par){
        if(par[2]  >= 0) return(Inf)
        val <- 0 +
          dnorm(par[1], mu0, sd=sqrt(sigma20), log=TRUE) +
          min(log(abs(par[2])), log(1e4)) +
          sum(invprobit(yo, par[1], 1/par[2], log=TRUE, oneminus = TRUE), na.rm=TRUE) +
          sum(invprobit(predm, par[1], sign(par[2]) * sqrt(1/par[2]^2 +  pred_var_m), log=TRUE), na.rm=TRUE)
        -val
      })
      if(opt_res$convergence != 0){
        warning("Dropout curve estimation did not properly converge")
      }
      rho[colidx] <- opt_res$par[1]
      zetainv[colidx] <- if(abs(opt_res$par[2]) > 1e4){
        -1e4
      }else{
        opt_res$par[2]
      }
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

  if(length(levels(design_fct)) == 1){
    # All entries are identical build an intercept only model
    mm <- matrix(1, nrow=length(design_fct), ncol=1)
    colnames(mm) <- levels(design_fct)
  }else if(is.null(reference_level)){
    helper_df <- data.frame(x_ = design_fct)
    mm <- stats::model.matrix.default(~ x_ - 1, helper_df)
    colnames(mm) <- sub("^x_", "", colnames(mm))
  }else{
    design_fct <- stats::relevel(design_fct, ref = reference_level)
    helper_df <- data.frame(x_ = design_fct)
    mm <- stats::model.matrix.default(~ x_ + 1, helper_df)
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
    col_data[has_ref_level] <- lapply(col_data[has_ref_level], stats::relevel, ref = reference_level)
  }
  mm <- stats::model.matrix.default(formula, col_data)
  colnames(mm)[colnames(mm) == "(Intercept)"] <- "Intercept"
  mm
}





log_parameters <- function(hp){
  names(hp) <- c("location_prior_mean", "location_prior_scale",
                             "dropout_curve_position", "dropout_curve_scale",
                             "variance_prior_scale", "variance_prior_df")
  hyper_para_txt <- paste0("The inferred parameters are:\n",
                           paste0(vapply(seq_along(hp), function(idx){
                             pretty_num <- if(names(hp)[idx] == "dropout_curve_scale"){
                               scales <- hp[[idx]]
                               ifelse(is.na(scales) | 1/scales > -100,
                                      formatC(1/scales, digits=3, width=1, format="g"),
                                      "< -100")
                             }else if(names(hp)[idx] == "variance_prior_df"){
                               if(is.na(hp[[idx]]) || 1/hp[[idx]] < 100){
                                 formatC(1/hp[[idx]], digits=3, width=1, format="g")
                               }else{
                                 "> 100"
                               }
                             }else{
                               formatC(hp[[idx]], digits=3, width=1, format="g")
                             }
                             paste0(names(hp)[idx], ":",
                                    paste0(rep(" ", times=24-nchar(names(hp)[idx])), collapse=""),
                                    paste0(pretty_num, collapse=", "))
                           }, FUN.VALUE = ""), collapse = "\n"))
  message(hyper_para_txt)
}




