

#' Fit a single linear probabilistic dropout model
#'
#'
#'
#' @return a list with the following entries
#'   \describe{
#'     \item{coefficients}{a named vector with the fitted values}
#'     \item{n_approx}{the estimated "size" of the data set (n_hat - variance_prior_df)}
#'     \item{df}{the estimated degrees of freedom (n_hat - p)}
#'     \item{s2}{the estimated unbiased variance}
#'     \item{rss}{the estimated sum of the squared residuals. `NA` if the variance is moderated}
#'     \item{n_obs}{the number of response values that were not `NA`}
#'   }
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
  X <- model.matrix(formula, data=data)
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

  pd_lm.fit(y, X,
            dropout_curve_position=dropout_curve_position,
            dropout_curve_scale=dropout_curve_scale,
            location_prior_mean = location_prior_mean,
            location_prior_scale = location_prior_scale,
            variance_prior_scale = variance_prior_scale,
            variance_prior_df = variance_prior_df,
            location_prior_df = location_prior_df,
            method = method,
            verbose = verbose)
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
#'     \item{rss}{the estimated sum of the squared residuals. `NA` if the variance is moderated}
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

  moderate_location <- !missing(location_prior_mean) && ! is.null(location_prior_mean) && ! is.na(location_prior_mean)
  moderate_variance <- !missing(variance_prior_scale) && ! is.null(variance_prior_scale)  && ! is.na(variance_prior_scale)

  if(! moderate_location && ! moderate_variance && nrow(X) < ncol(X) + 1){
    stop("Underdetermined system. There are more parameters to estimate than available rows.")
  }

  Xo <- X[!is.na(y), , drop=FALSE]
  Xm <- X[is.na(y), , drop=FALSE]
  yo <- y[!is.na(y)]
  p <- ncol(X)
  n <- nrow(X)

  all_observed <- all(! is.na(y))
  all_missing <- all(is.na(y))

  rho <- dropout_curve_position[is.na(y)]
  zeta <- dropout_curve_scale[is.na(y)]

  if(moderate_location){
    beta_init <- c(location_prior_mean, rep(0, times=p-1))
  }else if(length(yo) == 0){
    beta_init <- rep(0, times=p)
  }else{
    if(has_intercept(X)){
      beta_init <- c(mean(yo), rep(0, times=p-1))
    }else{
      beta_init <- rep(mean(yo), times=p)
    }

  }
  if(moderate_variance){
    sigma2_init <- variance_prior_df * variance_prior_scale / (variance_prior_df + 2)
  }else{
    sigma2_init <- 1
  }
  beta_sel <- seq_len(p)

  fit_beta <- rep(NA, p)
  names(fit_beta) <- colnames(X)

  if(all_observed && ! moderate_variance && ! moderate_location){
    # Run lm
    lm_res <- lm(yo ~ Xo - 1)
    fit_beta <- coefficients(lm_res)
    fit_sigma2 <- summary(lm_res)$sigma^2 * (n-p) / n
    fit_sigma2_var <- 2 * fit_sigma2^2 / n
  }else if(all_missing && ! moderate_variance){
    return(list(coefficients=fit_beta, n_approx=NA, df=NA, s2=NA, rss=NA, n_obs = length(yo)))
  }else if(method == "numeric"){
    opt_res <- stats::optim(par = c(beta_init, sigma2_init), function(par){
      beta <- par[beta_sel]
      sigma2 <- par[p+1]
      if(sigma2 <= 0) return(10000)
      zetastar <- zeta * sqrt(1 + sigma2/zeta^2)
      - objective_fnc(y, yo, X, Xm, Xo,
                    beta, sigma2, rho, zetastar,
                    location_prior_mean, location_prior_scale,
                    variance_prior_df, variance_prior_scale,
                    location_prior_df, moderate_location, moderate_variance)
    },
       method = "Nelder-Mead", hessian = TRUE)
    if(opt_res$convergence != 0){
      if(verbose){
        warning("Model didn't not properly converge\n")
        warning(opt_res$message, "\n")
        warning(y,"\n")
      }
      return(list(coefficients=fit_beta, n_approx=NA, df=NA, s2=NA, rss=NA, n_obs = length(yo)))
    }

    fit_beta <- opt_res$par[beta_sel]
    fit_sigma2 <- opt_res$par[p+1]
    fit_sigma2_var <- 1/opt_res$hessian[p+1, p+1]
  }else if(method == "analytic_grad"){
    # Run optim
    opt_res <- stats::optim(par = c(beta_init, sigma2_init), function(par){
      beta <- par[beta_sel]
      sigma2 <- par[p+1]
      if(sigma2 <= 0) return(10000)
      zetastar <- zeta * sqrt(1 + sigma2/zeta^2)
      - objective_fnc(y, yo, X, Xm, Xo,
                      beta, sigma2, rho, zetastar,
                      location_prior_mean, location_prior_scale,
                      variance_prior_df, variance_prior_scale,
                      location_prior_df, moderate_location, moderate_variance)
    },
    gr = function(par){
      beta <- par[beta_sel]
      sigma2 <- par[p+1]
      if(sigma2 <= 0) return(10000)
      zetastar <- zeta * sqrt(1 + sigma2/zeta^2)
      - grad_fnc(y, yo, X, Xm, Xo,
                 beta, sigma2, rho, zetastar,
                 location_prior_mean, location_prior_scale,
                 variance_prior_df, variance_prior_scale,
                 location_prior_df, moderate_location, moderate_variance)
    }, method = "BFGS", hessian=TRUE)
    if(opt_res$convergence != 0){
      if(verbose){
        warning("Model didn't not properly converge\n")
        warning(opt_res$message, "\n")
        warning(y,"\n")
      }
      return(list(coefficients=fit_beta, n_approx=NA, df=NA, s2=NA, rss=NA, n_obs = length(yo)))
    }

    fit_beta <- opt_res$par[beta_sel]
    fit_sigma2 <- opt_res$par[p+1]
    fit_sigma2_var <- 1/opt_res$hessian[p+1, p+1]
  }else if(method == "analytic_hessian"){
    # Run nlminb
    nl_res <- nlminb(start = c(beta_init, sigma2_init),
       objective = function(par){
         beta <- par[beta_sel]
         sigma2 <- par[p+1]
         zetastar <- zeta * sqrt(1 + sigma2/zeta^2)
         - objective_fnc(y, yo, X, Xm, Xo,
                         beta, sigma2, rho, zetastar,
                         location_prior_mean, location_prior_scale,
                         variance_prior_df, variance_prior_scale,
                         location_prior_df, moderate_location, moderate_variance)
       },
       gradient = function(par){
         beta <- par[beta_sel]
         sigma2 <- par[p+1]
         zetastar <- zeta * sqrt(1 + sigma2/zeta^2)
         - grad_fnc(y, yo, X, Xm, Xo,
                    beta, sigma2, rho, zetastar,
                    location_prior_mean, location_prior_scale,
                    variance_prior_df, variance_prior_scale,
                    location_prior_df, moderate_location, moderate_variance)
       },
       hessian = function(par){
         beta <- par[beta_sel]
         sigma2 <- par[p+1]
         zetastar <- zeta * sqrt(1 + sigma2/zeta^2)
         - hess_fnc(y, yo, X, Xm, Xo,
                    beta, sigma2, rho, zetastar,
                    location_prior_mean, location_prior_scale,
                    variance_prior_df, variance_prior_scale,
                    location_prior_df, moderate_location, moderate_variance,
                    beta_sel, p)
       }, lower= c(rep(-Inf, length(beta_init)), 0))
    if(nl_res$convergence != 0){
      if(verbose){
        warning("Model didn't not properly converge\n")
        warning(nl_res$message, "\n")
        warning(y,"\n")
      }
      return(list(coefficients=fit_beta, n_approx=NA, df=NA, s2=NA, rss=NA, n_obs = length(yo)))
    }

    fit_beta <- nl_res$par[beta_sel]
    fit_sigma2 <- nl_res$par[p+1]
    zetastar <- zeta * sqrt(1 + fit_sigma2/zeta^2)
    hessian <- - hess_fnc(y, yo, X, Xm, Xo,
                          fit_beta, fit_sigma2, rho, zetastar,
                          location_prior_mean, location_prior_scale,
                          variance_prior_df, variance_prior_scale,
                          location_prior_df, moderate_location, moderate_variance,
                          beta_sel, p)
    fit_sigma2_var <- 1/hessian[p+1, p+1]
  }

  # Set estimates to NA if no reasonable inference
  # zetastar <- zeta * sqrt(1 + fit_sigma2/zeta^2)
  # coef_should_be_inf <- vapply(seq_len(ncol(Xo)), function(colidx){
  #   sel <- Xm[, colidx] != 0
  #   any(invprobit(fit_beta[colidx], rho[sel], zetastar[sel], oneminus = TRUE) <  1e-8)
  # }, FUN.VALUE = FALSE)
  # fit_beta[coef_should_be_inf] <- NA
  # if(all(coef_should_be_inf)){
  #   n_approx <- NA
  #   df_approx <- NA
  #   rss_approx <- NA
  #   s2_approx <- NA
  # }else{
    n_approx <- 2 * fit_sigma2^2 / fit_sigma2_var
    rss_approx <- 2 * fit_sigma2^3 / fit_sigma2_var
    s2_approx <- rss_approx / (n_approx - p)

    if(s2_approx < 0 || n_approx <= p){
      df_approx <- 1e-3
      s2_approx <- sqrt(fit_sigma2_var / df_approx^2 * (df_approx + 2)^2 * (df_approx + 2) / 2)
      n_approx <- df_approx + p
      rss_approx <- s2_approx * df_approx
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
  # }

  if(moderate_variance){
    rss_approx <- NA
  }

  names(fit_beta) <- colnames(X)

  list(coefficients=fit_beta,
       n_approx=n_approx, df=df_approx,
       s2=s2_approx, rss = rss_approx,
       n_obs = length(yo))

}


objective_fnc <- function(y, yo, X, Xm, Xo, beta, sigma2, rho, zetastar, mu0, sigma20, df0, tau20, location_prior_df, moderate_location, moderate_variance){
  val <- 0
  if(moderate_location){
    val <- val + mean(dt.scaled(X %*% beta, df=location_prior_df, mean= mu0, sd=sqrt(sigma20), log=TRUE))
  }
  if(moderate_variance){
    val <- val +
      extraDistr::dinvchisq(sigma2, df0, sqrt(tau20), log=TRUE) +
      log(sigma2)   # important to remove the implicit 1/sigma2 prior in the Inv-Chisq distr
  }
  val <- val +
    sum(dnorm(Xo %*% beta, yo, sd=sqrt(sigma2), log=TRUE)) +
    sum(invprobit(Xm %*% beta, rho, zetastar, log=TRUE))
  val
}


grad_fnc <- function(y, yo, X, Xm, Xo, beta, sigma2, rho, zetastar, mu0, sigma20, df0, tau20, location_prior_df, moderate_location, moderate_variance){
  imr <- inv_mills_ratio(Xm %*% beta, rho, zetastar)

  if(moderate_location){
    dbeta_p <-  -(location_prior_df + 1) * t(X) %*% ((X %*% beta - mu0) / (location_prior_df * sigma20 + (X %*% beta - mu0)^2)) / nrow(X)
  }else{
    dbeta_p <-  0
  }
  dbeta_o <- - (t(Xo) %*% Xo %*% beta - t(Xo) %*% yo) / sigma2
  dbeta_m <- t(Xm) %*% imr

  if(moderate_variance){
    dsig2_p <- -(1 + df0/2) / sigma2 + df0 * sqrt(tau20) / (2 * sigma2^2) + 1/sigma2
  }else{
    dsig2_p <- 0
  }
  dsig2_o <- sum(((Xo %*% beta - yo)^2 - sigma2) / (2 * sigma2^2))
  dsig2_m <- -sum((Xm %*% beta - rho) / (2 * zetastar^2) * imr)

  c(dbeta_p + dbeta_o + dbeta_m, dsig2_p + dsig2_o + dsig2_m)
}


hess_fnc <- function(y, yo, X, Xm, Xo, beta, sigma2, rho, zetastar, mu0, sigma20, df0, tau20, location_prior_df,
                     moderate_location, moderate_variance, beta_sel, p){
    imr <- inv_mills_ratio(Xm %*% beta, rho, zetastar)

    if(moderate_location){
      t_prior_fact <- c((location_prior_df * sigma20 - (X %*% beta - mu0)^2) / (location_prior_df * sigma20 + (X %*% beta - mu0)^2)^2)
      dbb_p <- -(location_prior_df + 1) * t(X) %*% diag(t_prior_fact, nrow=nrow(X)) %*% X / nrow(X)
    }else{
      dbb_p <- 0
    }
    dbb_o <- -2 * t(Xo) %*% Xo / (2 * sigma2)
    dbb_m <- - t(Xm) %*% diag(c((imr^2 + (Xm %*% beta - rho) / zetastar^2 * imr)), nrow(Xm)) %*% Xm

    if(moderate_variance){
      dss_p <- (1 + df0/2)/ (sigma2^2) - df0 * sqrt(tau20) / (sigma2^3) - 1/sigma2^2
    }else{
      dss_p <- 0
    }
    dss_o <- sum((sigma2 - 2 * (Xo %*% beta - yo)^2) / (2 * sigma2^3))
    dss_m <- sum((Xm %*% beta - rho) / (4 * zetastar^4) * imr *
                   (3 - (Xm %*% beta - rho) * imr - (Xm %*% beta - rho)^2 / zetastar^2))

    dbs_o <- t(Xo) %*% (Xo %*% beta - yo) / sigma2^2
    dbs_m <- t(Xm) %*% ((Xm %*% beta - rho) / (2 * zetastar^2) * imr^2 -
                          (zetastar^2 - (Xm %*% beta - rho)^2) / (2 * zetastar^4) * imr)

    res <- matrix(NA, nrow=p+1, ncol=p + 1)
    res[beta_sel, beta_sel] <- dbb_p + dbb_o + dbb_m
    res[p+1, p+1] <- dss_p + dss_o + dss_m
    res[p+1, beta_sel] <- c(dbs_o + dbs_m)
    res[beta_sel, p+1] <- c(dbs_o + dbs_m)
    res
}


has_intercept <- function(X){

  any(vapply(seq_len(ncol(X)), function(idx){
    all(X[, idx] == 1)
  }, FUN.VALUE = FALSE))

}

