

#' Fit a single linear probabilistic dropout model
#'
#' The function works similar to the classical \code{\link[stats]{lm}}
#' but with special handling of \code{NA}'s. Whereas \code{lm} usually
#' just ignores response value that are missing, \code{pd_lm} applies
#' a probabilistic dropout model, that assumes that missing values
#' occur because of the dropout curve. The dropout curve describes for
#' each position the chance that that a value is missed. A negative
#' \code{dropout_curve_scale} means that the lower the intensity was,
#' the more likely it is to miss the value.
#'
#' @param formula a formula that specifies a linear model
#' @param data an optional data.frame whose columns can be used to
#'   specify the \code{formula}
#' @param subset an optional selection vector for data to subset it
#' @param dropout_curve_position the value where the chance to
#'   observe a value is 50\%. Can either be a single value that is
#'   repeated for each row or a vector with one element for each
#'   row. Not optional.
#' @param dropout_curve_scale the width of the dropout curve. Smaller
#'   values mean that the sigmoidal curve is steeper.
#'   Can either be a single value that is
#'   repeated for each row or a vector with one element for each
#'   row. Not optional.
#' @param location_prior_mean,location_prior_scale the optional mean
#'   and variance of the prior around
#'   which the predictions are supposed to scatter. If no value is
#'   provided no location regularization is applied.
#' @param variance_prior_scale,variance_prior_df the optional scale and degrees of
#'   freedom of the variance prior.
#'   If no value is provided no variance regularization is applied.
#' @param location_prior_df The degrees of freedom for the t-distribution
#'   of the location prior. If it is large (> 30) the prior is approximately
#'   Normal. Default: 3
#' @param method one of 'analytic_hessian', 'analytic_gradient', or
#'   'numeric'. If 'analytic_hessian' the \code{\link[stats]{nlminb}}
#'   optimization routine is used, with the hand derived first and
#'   second derivative. Otherwise, \code{\link[stats]{optim}} either
#'   with or without the first derivative is used.
#' @param verbose boolean that signals if the method prints informative
#'   messages. Default: \code{FALSE}.
#'
#'
#'
#' @return a list with the following entries
#'   \describe{
#'     \item{coefficients}{a named vector with the fitted values}
#'     \item{coef_variance_matrix}{a \code{p*p} matrix with the variance associated
#'       with each coefficient estimate}
#'     \item{n_approx}{the estimated "size" of the data set (n_hat - variance_prior_df)}
#'     \item{df}{the estimated degrees of freedom (n_hat - p)}
#'     \item{s2}{the estimated unbiased variance}
#'     \item{n_obs}{the number of response values that were not `NA`}
#'   }
#'
#' @examples
#'   # Without missing values
#'   y <- rnorm(5, mean=20)
#'   lm(y ~ 1)
#'   pd_lm(y ~ 1,
#'         dropout_curve_position = NA,
#'         dropout_curve_scale = NA)
#'
#'   # With some missing values
#'   y <- c(23, 21.4, NA)
#'   lm(y ~ 1)
#'   pd_lm(y ~ 1,
#'         dropout_curve_position = 19,
#'         dropout_curve_scale = -1)
#'
#'
#'   # With only missing values
#'   y <- c(NA, NA, NA)
#'   # lm(y ~ 1)  # Fails
#'   pd_lm(y ~ 1,
#'         dropout_curve_position = 19,
#'         dropout_curve_scale = -1,
#'         location_prior_mean = 21,
#'         location_prior_scale = 3,
#'         variance_prior_scale = 0.1,
#'         variance_prior_df = 2)
#'
#'
#' @export
pd_lm <- function(formula, data = NULL, subset = NULL,
                  dropout_curve_position,
                  dropout_curve_scale,
                  location_prior_mean = NULL,
                  location_prior_scale = NULL,
                  variance_prior_scale = NULL,
                  variance_prior_df = NULL,
                  location_prior_df = 3,
                  method =  c("analytic_hessian", "analytic_grad", "numeric"),
                  verbose = FALSE){

  method <- match.arg(method, c("analytic_hessian", "analytic_grad", "numeric"))
  if(! is.null(subset) && ! is.null(data)){
    data <- data[subset, ,drop=FALSE]
  }else if(is.null(data)){
    data <- environment(formula)
  }

  current.na.action <- options('na.action')
  options(na.action = "na.pass")
  X <- stats::model.matrix.default(formula, data=data)
  options('na.action' = unname(unlist(current.na.action)))
  colnames(X)[colnames(X) == "(Intercept)"] <- "Intercept"

  y <- eval(formula[[2]], data, environment(formula))
  if(is.null(y)){
    stop("Formula ", formula, " is one sided.")
  }

  if(missing(dropout_curve_position) || missing(dropout_curve_scale)){
    stop("Neither dropout_curve_position nor dropout_curve_scale must be missing.")
  }
  if(length(dropout_curve_position) == 1){
    dropout_curve_position <- rep(dropout_curve_position, nrow(X))
  }
  if(length(dropout_curve_scale) == 1){
    dropout_curve_scale <- rep(dropout_curve_scale, nrow(X))
  }

  res <- pd_lm.fit(y, X,
            dropout_curve_position=dropout_curve_position,
            dropout_curve_scale=dropout_curve_scale,
            location_prior_mean = location_prior_mean,
            location_prior_scale = location_prior_scale,
            variance_prior_scale = variance_prior_scale,
            variance_prior_df = variance_prior_df,
            location_prior_df = location_prior_df,
            method = method,
            verbose = verbose)
  res[c("coefficients", "coef_variance_matrix", "n_approx", "df", "s2", "n_obs")]
}


#' The work horse for fitting the probabilistic dropout model
#'
#' If there is no location and variance moderation and no missing values,
#' the model is fitted with `lm`.
#'
#' @return a list with the following entries
#'   \describe{
#'     \item{coefficients}{a named vector with the fitted values}
#'     \item{n_approx}{the estimated "size" of the data set (n_hat - variance_prior_df)}
#'     \item{df}{the estimated degrees of freedom (n_hat - p)}
#'     \item{s2}{the estimated unbiased variance}
#'     \item{n_obs}{the number of response values that were not `NA`}
#'   }
#' @keywords internal
pd_lm.fit <- function(y, X,
  dropout_curve_position,
  dropout_curve_scale,
  location_prior_mean = NULL,
  location_prior_scale = NULL,
  variance_prior_scale = NULL,
  variance_prior_df = NULL,
  location_prior_df = 3,
  method = c("analytic_hessian", "analytic_grad", "numeric"),
  verbose = FALSE){
  
  method <- match.arg(method, c("analytic_hessian", "analytic_grad", "numeric"))
    
  moderate_location <- !missing(location_prior_mean) && 
    !is.null(location_prior_mean) && !is.na(location_prior_mean)
  moderate_variance <- !missing(variance_prior_scale) && 
    !is.null(variance_prior_scale)  && !is.na(variance_prior_scale)
    
  if (!moderate_location && !moderate_variance && nrow(X) < ncol(X) + 1) {
    stop("Underdetermined system. There are more parameters to estimate than available rows.")
  }
    
  y_na <- is.na(y)
  Xo <- X[!y_na, , drop = FALSE]
  Xm <- X[y_na, , drop = FALSE]
  yo <- y[!y_na]
  p <- ncol(X)
  n <- nrow(X)
    
  all_observed <- all(!y_na)
  all_missing <- !all_observed ##all(y_na)
    
    rho <- dropout_curve_position[y_na]
    zeta <- dropout_curve_scale[y_na]
    
    ## create initial value for beta
    if (moderate_location) {
        beta_init <- c(location_prior_mean, rep(0, times = p - 1))
    } else if (length(yo) == 0) {
        beta_init <- rep(0, times = p)
    } else {
        if (has_intercept(X)) {
            beta_init <- c(mean(yo), rep(0, times = p - 1))
        } else {
            beta_init <- rep(mean(yo), times = p)
        }
    }
    
    ## create initial value for sigma2
    if (moderate_variance) {
        sigma2_init <- variance_prior_df * variance_prior_scale / (variance_prior_df + 2)
    } else {
        sigma2_init <- 1
    }
    
    beta_sel <- seq_len(p)
    
    fit_beta <- rep(NA, p)
    names(fit_beta) <- colnames(X)
    
    failed_result <- list(coefficients = rep(NA, p),
        coef_variance_matrix = matrix(NA, nrow = p, ncol = p),
        correction_factor = matrix(NA, nrow = p, ncol = p),
        n_approx = NA, df = NA, s2 = NA, n_obs = length(yo))
    
    if (all_observed && !moderate_variance && !moderate_location) {
        ## Run lm if there are no missing values
        lm_res <- lm(yo ~ Xo - 1)
        fit_beta <- coefficients(lm_res)
        fit_sigma2 <- summary(lm_res)$sigma ^ 2 * (n - p) / n
        coef_hessian <-  1 / fit_sigma2 * (t(X) %*% X)
        fit_sigma2_var <- 2 * fit_sigma2 ^ 2 / n
      
    } else if (all_missing && !moderate_variance) {
        return(failed_result)
      
    } else if (method == "numeric") {
        opt_res <- stats::optim(par = c(beta_init, sigma2_init), function(par) {
            beta <- par[beta_sel]
            sigma2 <- par[p + 1]
            if (sigma2 <= 0) return (10000)
            zetastar <- zeta * sqrt(1 + sigma2 / zeta ^ 2)
            -objective_fnc(y, yo, X, Xm, Xo,
                beta, sigma2, rho, zetastar,
                location_prior_mean, location_prior_scale,
                variance_prior_df, variance_prior_scale,
                location_prior_df, moderate_location, moderate_variance)
        },
        method = "Nelder-Mead", hessian = TRUE)
        
        if (opt_res$convergence != 0) {
            return(failed_result)
        }
        
        fit_beta <- opt_res$par[beta_sel]
        coef_hessian <- opt_res$hessian[beta_sel, beta_sel, drop = FALSE]
        fit_sigma2 <- opt_res$par[p + 1]
        fit_sigma2_var <- 1 / opt_res$hessian[p + 1, p + 1]
      
    } else if (method == "analytic_grad") {
        # Run optim
        opt_res <- stats::optim(par = c(beta_init, sigma2_init), function(par) {
            beta <- par[beta_sel]
            sigma2 <- par[p + 1]
            if (sigma2 <= 0) return(10000)
            zetastar <- zeta * sqrt(1 + sigma2 / zeta ^ 2)
            -objective_fnc(y, yo, X, Xm, Xo,
                beta, sigma2, rho, zetastar,
                location_prior_mean, location_prior_scale,
                variance_prior_df, variance_prior_scale,
                location_prior_df, moderate_location, moderate_variance)
            },
            gr = function(par) {
                beta <- par[beta_sel]
                sigma2 <- par[p + 1]
                if (sigma2 <= 0) return(10000)
                zetastar <- zeta * sqrt(1 + sigma2 / zeta ^ 2)
                -grad_fnc(y, yo, X, Xm, Xo,
                    beta, sigma2, rho, zetastar,
                    location_prior_mean, location_prior_scale,
                    variance_prior_df, variance_prior_scale,
                    location_prior_df, moderate_location, moderate_variance)
            }, method = "BFGS", hessian=TRUE)
            
            if (opt_res$convergence != 0) {
                return(failed_result)
            }
          
            fit_beta <- opt_res$par[beta_sel]
            coef_hessian <- opt_res$hessian[beta_sel, beta_sel, drop = FALSE]
            fit_sigma2 <- opt_res$par[p + 1]
            fit_sigma2_var <- 1 / opt_res$hessian[p + 1, p + 1]
      
    } else if (method == "analytic_hessian") {
        ## Run nlminb
        nl_res <- nlminb(start = c(beta_init, sigma2_init),
            objective = function(par) {
                beta <- par[beta_sel]
                sigma2 <- par[p + 1]
                zetastar <- zeta * sqrt(1 + sigma2 / zeta ^ 2)
                -objective_fnc(y, yo, X, Xm, Xo,
                    beta, sigma2, rho, zetastar,
                    location_prior_mean, location_prior_scale,
                    variance_prior_df, variance_prior_scale,
                    location_prior_df, moderate_location, moderate_variance)
            },
            gradient = function(par) {
                beta <- par[beta_sel]
                sigma2 <- par[p + 1]
                zetastar <- zeta * sqrt(1 + sigma2 / zeta ^ 2)
                -grad_fnc(y, yo, X, Xm, Xo,
                    beta, sigma2, rho, zetastar,
                    location_prior_mean, location_prior_scale,
                    variance_prior_df, variance_prior_scale,
                    location_prior_df, moderate_location, moderate_variance)
            },
            hessian = function(par) {
                beta <- par[beta_sel]
                sigma2 <- par[p + 1]
                zetastar <- zeta * sqrt(1 + sigma2 / zeta ^ 2)
                -hess_fnc(y, yo, X, Xm, Xo,
                    beta, sigma2, rho, zetastar,
                    location_prior_mean, location_prior_scale,
                    variance_prior_df, variance_prior_scale,
                    location_prior_df, moderate_location, moderate_variance,
                    beta_sel, p)
        }, lower= c(rep(-Inf, length(beta_init)), 0))
        
        if (nl_res$convergence != 0) {
            return(failed_result)
        }
        
        fit_beta <- nl_res$par[beta_sel]
        fit_sigma2 <- nl_res$par[p + 1]
        zetastar <- zeta * sqrt(1 + fit_sigma2 / zeta ^ 2)
        hessian <- -hess_fnc(y, yo, X, Xm, Xo,
            fit_beta, fit_sigma2, rho, zetastar,
            location_prior_mean, location_prior_scale,
            variance_prior_df, variance_prior_scale,
            location_prior_df, moderate_location, moderate_variance,
            beta_sel, p)
        fit_sigma2_var <- 1 / hessian[p + 1, p + 1]
        coef_hessian <- hessian[beta_sel, beta_sel, drop = FALSE]
    }
    
    if (fit_sigma2_var < 0){
      return(failed_result)
    }
    
    Var_coef <- invert_hessian_matrix(coef_hessian, p, verbose)
    
    # Use fitted sigma2 and associated uncertainty to estimate df and unbiased sigma2
    sigma2_params <- calculate_sigma2_parameters(fit_sigma2, fit_sigma2_var,
        variance_prior_scale, variance_prior_df, moderate_variance, n, p)
    n_approx <- sigma2_params$n_approx
    df_approx <- sigma2_params$df_approx
    s2_approx <- sigma2_params$s2_approx
    
    # Calculate correction factor for skew
    # the skew means that the fit is bad on both sides. We only care about the
    # right side, so we will calculate a factor that improves that one
    zetastar <- zeta * sqrt(1 + fit_sigma2 / zeta ^ 2)
    Correction_Factor <- calculate_skew_correction_factors(y, yo, X, Xm, Xo,
        fit_beta, fit_sigma2, Var_coef, rho, zetastar, location_prior_mean, 
        location_prior_scale, variance_prior_df, variance_prior_scale,
        location_prior_df, moderate_location, moderate_variance, out_factor = 8)
    
    # Correct Var_coef to make it unbiased
    # Plugging unbiased s2_approx into Hessian calculation
    hessian <- - hess_fnc(y, yo, X, Xm, Xo, fit_beta, s2_approx, rho, zetastar,
        location_prior_mean, location_prior_scale, variance_prior_df, 
        variance_prior_scale, location_prior_df, moderate_location, 
        moderate_variance, beta_sel, p)
    coef_hessian <- hessian[beta_sel, beta_sel, drop = FALSE]
    
    ## Set hessian to zero, so that the variances become infinite when one
    ## diagonal element is negative (this is bad)
    if (any(diag(coef_hessian) < 0)) {
        coef_hessian <- matrix(0, nrow = p, ncol = p)
    }
    Var_coef_unbiased <- invert_hessian_matrix(coef_hessian, p, verbose)
    
    ## Apply correction factor to the unbiased Var_coef
    Var_coef_unbiased <- Correction_Factor %*% Var_coef_unbiased %*% Correction_Factor
    
    ## Set estimates to NA if no reasonable inference
    coef_should_be_inf <- diag(Var_coef) > 1e6
    for (idx in which(coef_should_be_inf)) {
        Var_coef_unbiased[idx, idx] <- Inf
    }
    
    fit_beta[coef_should_be_inf] <- NA
    if (all(coef_should_be_inf)) {
        n_approx <- df_approx <- s2_approx <- NA
    }
    
    # Make everything pretty and return results
    names(fit_beta) <- colnames(Var_coef_unbiased) <- 
        rownames(Var_coef_unbiased) <- colnames(X)
    
    list(coefficients = fit_beta,
        coef_variance_matrix = Var_coef_unbiased,
        correction_factor = Correction_Factor,
        n_approx = n_approx, df = df_approx,
        s2 = s2_approx,
        n_obs = length(yo))
}

objective_fnc <- function(y, yo, X, Xm, Xo, beta, sigma2, rho, zetastar, mu0, 
    sigma20, df0, tau20, location_prior_df, moderate_location, 
    moderate_variance) {
    val <- 0
    if (moderate_location) {
        val <- sum(dt.scaled(X %*% beta, df = location_prior_df, 
            mean = mu0, sd = sqrt(sigma20), log = TRUE))
    }
    if (moderate_variance) {
        val <- val + extraDistr::dinvchisq(sigma2, df0, tau20, log = TRUE) + 
            log(sigma2)
    }
    val + sum(dnorm(Xo %*% beta, yo, sd = sqrt(sigma2), log = TRUE)) + 
        sum(invprobit(Xm %*% beta, rho, zetastar, log = TRUE))
}

grad_fnc <- function (y, yo, X, Xm, Xo, beta, sigma2, rho, zetastar, mu0, 
    sigma20, df0, tau20, location_prior_df, moderate_location, 
    moderate_variance) {
  
    imr <- proDA:::inv_mills_ratio(Xm %*% beta, rho, zetastar)
  
    if (moderate_location) {
        Xbmu <- X %*% beta - mu0
        dbeta_p <- -(location_prior_df + 1) * t(X) %*% 
            (Xbmu / (location_prior_df * sigma20 + Xbmu^2))
    } else {
        dbeta_p <- 0
    }
  
    Xo_t <- t(Xo)
    dbeta_o <- -(Xo_t %*% Xo %*% beta - Xo_t %*% yo) / sigma2
    dbeta_m <- t(Xm) %*% imr
  
  
    if (moderate_variance) {
        dsig2_p <- -(1 + df0 / 2) / sigma2 + df0 * tau20 / (2 * sigma2 ^ 2) + 
            1 / sigma2
    } else {
        dsig2_p <- 0
    }
    
    dsig2_o <- sum(((Xo %*% beta - yo) ^ 2 - sigma2) / (2 * sigma2 ^ 2))
    dsig2_m <- -sum((Xm %*% beta - rho)/(2 * zetastar ^ 2) * imr)
    c(dbeta_p + dbeta_o + dbeta_m, dsig2_p + dsig2_o + dsig2_m)
}


hess_fnc <- function(y, yo, X, Xm, Xo, beta, sigma2, rho, zetastar, mu0,
    sigma20, df0, tau20, location_prior_df, moderate_location,
    moderate_variance, beta_sel, p){
    
    imr <- inv_mills_ratio(Xm %*% beta, rho, zetastar)
    
    ## precalculate values
    zetastar_2 <- zetastar ^ 2
    zetastar_4 <- zetastar_2 ^ 2
    Xmbr <- Xm %*% beta - rho
    X0by0 <- Xo %*% beta - yo
    q <- p + 1
    
    if (moderate_location) {
        Xbm <- (X %*% beta - mu0) ^ 2
        t_prior_fact <- (location_prior_df * sigma20 - Xbm) / (location_prior_df * sigma20 + Xbm) ^ 2
        dbb_p <- -(location_prior_df + 1) * t(X) %*% diag(t_prior_fact, nrow = nrow(X)) %*% X
    } else {
        dbb_p <- 0
    }
    dbb_o <- -2 * t(Xo) %*% Xo / (2 * sigma2)
    dbb_m <- - t(Xm) %*% diag(c((imr ^ 2 + Xmbr / zetastar_2 * imr)), nrow(Xm)) %*% Xm

    if (moderate_variance) {
        dss_p <- (1 + df0 / 2)/ (sigma2 ^ 2) - df0 * tau20 / (sigma2 ^ 3) - 1 / sigma2 ^ 2
    } else {
        dss_p <- 0
    }
    dss_o <- sum((sigma2 - 2 * Xoby0 ^ 2) / (2 * sigma2 ^ 3))
    dss_m <- sum(Xmbr / (4 * zetastar_4) * imr *
        (3 - Xmbr * imr - Xmbr ^ 2 / zetastar_2))

    dbs_o <- t(Xo) %*% X0by0 / sigma2 ^ 2
    dbs_m <- t(Xm) %*% (Xmbr / (2 * zetastar_2) * imr ^ 2 -
        (zetastar_2 - Xmbr^2) / (2 * zetastar_4) * imr)

    res <- matrix(NA, nrow = q, ncol = q)
    res[beta_sel, beta_sel] <- dbb_p + dbb_o + dbb_m
    res[q, q] <- dss_p + dss_o + dss_m
    res[q, beta_sel] <- res[beta_sel, q] <- dbs_o + dbs_m
    res
}


has_intercept <- function(X){

  any(vapply(seq_len(ncol(X)), function(idx){
    all(X[, idx] == 1)
  }, FUN.VALUE = FALSE))

}



calculate_skew_correction_factors <- function(y, yo, X, Xm, Xo, fit_beta, fit_sigma2, Var_coef, rho, zetastar,
                                             mu0, sigma20, df0, tau20, location_prior_df,
                                             moderate_location, moderate_variance, out_factor = 2){
  p <- length(fit_beta)
  res <- vapply(seq_len(p), function(idx){
    if(any(is.na(fit_beta))){
      1
    }else{
      offset <- objective_fnc(y = y,
                              yo = yo,
                              X = X,
                              Xm = Xm,
                              Xo = Xo,
                              beta = fit_beta,
                              sigma2 = fit_sigma2,
                              rho = rho,
                              zetastar = zetastar,
                              mu0 = mu0,
                              sigma20 = sigma20,
                              df0 = df0,
                              tau20 = tau20,
                              location_prior_df = location_prior_df,
                              moderate_location = moderate_location,
                              moderate_variance = moderate_variance)

      beta_shift <- rep(0, length(fit_beta))
      # Need to adapt the variance for the conditioning
      Mat_remain_inv <- tryCatch(solve(Var_coef[-idx, -idx,drop=FALSE]),
                                 error = function(err){
                                   matrix(NA, ncol=p-1, nrow=p-1)
                                 })
      var_coef_idx_corrected <- Var_coef[idx, idx ,drop=FALSE] -
        Var_coef[idx, -idx,drop=FALSE] %*% Mat_remain_inv %*% Var_coef[-idx, idx,drop=FALSE]
      if(is.na(var_coef_idx_corrected) || var_coef_idx_corrected < 0){
        beta_shift[idx] <- NA
      }else{
        beta_shift[idx] <- sqrt(out_factor * var_coef_idx_corrected)
      }

      diff <- objective_fnc(y = y,
                            yo = yo,
                            X = X,
                            Xm = Xm,
                            Xo = Xo,
                            beta = fit_beta + beta_shift,
                            sigma2 = fit_sigma2,
                            rho = rho,
                            zetastar = zetastar,
                            mu0 = mu0,
                            sigma20 = sigma20,
                            df0 = df0,
                            tau20 = tau20,
                            location_prior_df = location_prior_df,
                            moderate_location = moderate_location,
                            moderate_variance = moderate_variance)

      (abs(diff - offset) / (out_factor / 2))^(-1)
    }
  }, FUN.VALUE = 0.0)

  diag(sqrt(res), nrow=length(fit_beta))
}



invert_hessian_matrix <- function(coef_hessian, p, verbose = FALSE){
  Var_coef <- diag(Inf, nrow=p)
  # Make hessian robust for inversion!
  very_small_entry <- which(diag(coef_hessian)  < 1e-10)
  if(length(very_small_entry) == p){
    # All entries of matrix are practically zero
    Var_coef <- diag(Inf, nrow=p)
  }else{
    for(idx in very_small_entry){
      coef_hessian[idx, idx] <- 1e-10
    }
    # Reduce all very large numbers in the table
    coef_hessian[coef_hessian > 1e10] <- 1e10
    tryCatch({
      Var_coef <- solve(coef_hessian)
    }, error = function(err){
      if(verbose){
        warning(err)
      }
    })
  }
  Var_coef
}

calculate_sigma2_parameters <- function(fit_sigma2, fit_sigma2_var,
                                        variance_prior_scale, variance_prior_df,
                                        moderate_variance, n, p){

  n_approx <- 2 * fit_sigma2^2 / fit_sigma2_var
  rss_approx <- 2 * fit_sigma2^3 / fit_sigma2_var
  s2_approx <- rss_approx / (n_approx - p)

  if(s2_approx < 0 || n_approx <= p){
    df_approx <- 1e-3
    n_approx <- df_approx + p
    s2_approx <-  sqrt(fit_sigma2_var * (n_approx)^3/ (2 * df_approx^2))
  }else{
    if(moderate_variance){
      n_approx <- n_approx - variance_prior_df
      df_approx <- n_approx - p + variance_prior_df
      if(df_approx > 30 * n && df_approx > 100){
        df_approx <- Inf
        s2_approx <- variance_prior_scale   # Check if this is correct
      }
    }else{
      df_approx <- n_approx - p
    }
  }
  list(n_approx = n_approx, df_approx = df_approx, s2_approx = s2_approx)
}

