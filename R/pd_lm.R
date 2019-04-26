
#' Fit linear probabilistic dropout model
#'
pd_lm <- function(y, X, rho, zeta, mu0, sigma20, df0, tau20,
                  location_prior_df = 3,
                  method = c("analytic_grad", "numeric_grad")){
  stopifnot(is.numeric(y))
  stopifnot(is.matrix(X))
  stopifnot(length(y) == nrow(X))
  stopifnot(length(y) == length(rho))
  stopifnot(length(y) == length(zeta))
  method <- match.arg(method,  c("analytic_grad", "numeric_grad"))

  p <- ncol(X)

  Xo <- X[!is.na(y), , drop=FALSE]
  Xm <- X[is.na(y), , drop=FALSE]
  yo <- y[!is.na(y)]

  rho <- rho[is.na(y)]
  zeta <- zeta[is.na(y)]

  beta_init <- c(mu0, rep(0, times=p-1))
  sigma2_init <- df0 * tau20 / (df0 + 2)


  if(method == "numeric_grad"){
    opt_res <- stats::optim(par = c(beta_init, sigma2_init), function(par){
      beta <- par[seq_len(ncol(X))]
      sigma2 <- par[p+1]
      if(sigma2 <= 0) return(10000)


      zetastar <- zeta * sqrt(1 + sigma2/zeta^2)
      val <- 0 +
        mean(dt.scaled(X %*% beta, df=3, mean= mu0, sd=sqrt(sigma20), log=TRUE)) +
        # mean(dnorm(X %*% beta, mean= mu0, sd=sqrt(sigma20), log=TRUE)) +
        extraDistr::dinvchisq(sigma2, df0, sqrt(tau20), log=TRUE) +
        log(sigma2) +   # important to remove the implicit 1/sigma2 prior in the Inv-Chisq distr
        sum(dnorm(Xo %*% beta, yo, sd=sqrt(sigma2), log=TRUE)) +
        sum(invprobit(Xm %*% beta, rho, zetastar, log=TRUE))
      -val
    }, hessian=TRUE)
    if(opt_res$convergence != 0){
      warning("Model didn't not properly converge\n")
      warning(opt_res$message, "\n")
      warning(y,"\n")
      return(list(beta=rep(NA, p), n=NA, df=NA, s2=NA, n_obs = length(yo)))
    }

    fit_beta <- opt_res$par[seq_len(p)]
    fit_sigma2 <- opt_res$par[p+1]
    fit_sigma2_var <- 1/opt_res$hessian[p+1, p+1]   # Get the last element on the diagonal
  }else if(method == "analytic_grad"){
    lik_res <- maxLik::maxNR(function(par){
      beta <- par[seq_len(p)]
      sigma2 <- par[p+1]
      if(sigma2 <= 0) return(-Inf)
      pred <- c(X %*% beta)

      zetastar <- zeta * sqrt(1 + sigma2/zeta^2)
      mean(dt.scaled(X %*% beta, df=3, mean= mu0, sd=sqrt(sigma20), log=TRUE)) +
        # mean(dnorm(X %*% beta, mean= mu0, sd=sqrt(sigma20), log=TRUE)) +
        extraDistr::dinvchisq(sigma2, df0, sqrt(tau20), log=TRUE) +
        log(sigma2) +
        sum(dnorm(pred, y, sd=sqrt(sigma2), log=TRUE), na.rm=TRUE) +
        sum(invprobit(pred[is.na(y)], rho, zetastar, log=TRUE))
    }, grad = function(par){
      beta <- par[seq_len(p)]
      sigma2 <- par[p+1]
      if(sigma2 <= 0) return(NA)
      zetastar <- zeta * sqrt(1 + sigma2/zeta^2)
      imr <- inv_mills_ratio(Xm %*% beta, rho, zetastar)

      dbeta_p <-  -(location_prior_df + 1) * t(X) %*% ((X %*% beta - mu0) / (location_prior_df * sigma20 + (X %*% beta - mu0)^2)) / nrow(X)
      # dbeta_p <- -(beta - mu0) / (sigma20)
      dbeta_o <- - (t(Xo) %*% Xo %*% beta - t(Xo) %*% yo) / sigma2
      dbeta_m <- t(Xm) %*% imr

      dsig2_p <- -(1 + df0/2) / sigma2 + df0 * sqrt(tau20) / (2 * sigma2^2) + 1/sigma2
      dsig2_o <- sum(((Xo %*% beta - yo)^2 - sigma2) / (2 * sigma2^2))
      dsig2_m <- -sum((Xm %*% beta - rho) / (2 * zetastar^2) * imr)

      c(dbeta_p + dbeta_o + dbeta_m, dsig2_p + dsig2_o + dsig2_m)
    }, hess = function(par){
      beta <- par[seq_len(p)]
      sigma2 <- par[p+1]
      if(sigma2 <= 0) return(NA)
      zetastar <- zeta * sqrt(1 + sigma2/zeta^2)
      imr <- inv_mills_ratio(Xm %*% beta, rho, zetastar)

      t_prior_fact <- c((location_prior_df * sigma20 - (X %*% beta - mu0)^2) / (location_prior_df * sigma20 + (X %*% beta - mu0)^2)^2)
      dbb_p <- -(location_prior_df + 1) * t(X) %*% diag(t_prior_fact, nrow=nrow(X)) %*% X / nrow(X)
      # dbb_p <- -1/sigma20 * diag(nrow=p)
      dbb_o <- -2 * t(Xo) %*% Xo / (2 * sigma2)
      dbb_m <- - t(Xm) %*% diag(c((imr^2 + (Xm %*% beta - rho) / zetastar^2 * imr)), nrow(Xm)) %*% Xm

      dss_p <- (1 + df0/2)/ (sigma2^2) - df0 * sqrt(tau20) / (sigma2^3) - 1/sigma2^2
      dss_o <- sum((sigma2 - 2 * (Xo %*% beta - yo)^2) / (2 * sigma2^3))
      dss_m <- sum((Xm %*% beta - rho) / (4 * zetastar^4) * imr *
                     (3 - (Xm %*% beta - rho) * imr - (Xm %*% beta - rho)^2 / zetastar^2))

      dbs_o <- t(Xo) %*% (Xo %*% beta - yo) / sigma2^2
      dbs_m <- t(Xm) %*% ((Xm %*% beta - rho) / (2 * zetastar^2) * imr^2 -
                            (zetastar^2 - (Xm %*% beta - rho)^2) / (2 * zetastar^4) * imr)

      res <- matrix(NA, nrow=length(par), ncol=length(par))
      res[seq_along(beta), seq_along(beta)] <- dbb_p + dbb_o + dbb_m
      res[length(par), length(par)] <- dss_p + dss_o + dss_m
      res[length(par), seq_along(beta)] <- c(dbs_o + dbs_m)
      res[seq_along(beta), length(par)] <- c(dbs_o + dbs_m)
      res
    }, start = c(beta_init, sigma2_init),
    finalHessian = TRUE)

    if(lik_res$code > 2){
      warning("Model didn't not properly converge\n")
      warning(lik_res$message)
      warning(y,"\n")
      return(list(beta=rep(NA, p), n=NA, df=NA, s2=NA, n_obs = length(yo)))
    }
    fit_beta <- lik_res$estimate[seq_len(p)]
    fit_sigma2 <- lik_res$estimate[p+1]
    fit_sigma2_var <- -1/lik_res$hessian[p+1, p+1]   # Get the last element on the diagonal
  }

  n_approx <- 2 * fit_sigma2^2 / fit_sigma2_var
  rss_approx <- 2 * fit_sigma2^3 / fit_sigma2_var
  s2_approx <- rss_approx / (n_approx - p)

  if(s2_approx < 0 || n_approx <= p){
    df_mod <- 1e-3
    s2_approx <- sqrt(fit_sigma2_var / df_mod^2 * (df_mod + 2)^2 * (df_mod + 2) / 2)
    n_approx <- df_mod + p
    rss_approx <- s2_approx * df_mod
  }

  names(fit_beta) <- colnames(X)

  list(coefficients=fit_beta,
       n_approx=n_approx - df0, df=n_approx - p,
       s2=s2_approx, rss = rss_approx,
       n_obs = length(yo))
}



pd_lm_unreg <- function(y, X, rho, zeta,
                  method = c("analytic_grad", "numeric_grad")){
  stopifnot(is.numeric(y))
  stopifnot(is.matrix(X))
  stopifnot(length(y) == nrow(X))
  stopifnot(length(y) == length(rho))
  stopifnot(length(y) == length(zeta))
  method <- match.arg(method,  c("analytic_grad", "numeric_grad"))

  p <- ncol(X)

  Xo <- X[!is.na(y), , drop=FALSE]
  Xm <- X[is.na(y), , drop=FALSE]
  yo <- y[!is.na(y)]

  rho <- rho[is.na(y)]
  zeta <- zeta[is.na(y)]

  beta_init <-  c(mean(y,na.rm=TRUE), rep(0, times=p-1))
  sigma2_init <- 1

  if(length(yo) < p+1){
    # There are less observation than parameters
    return(list(beta=rep(NA, p), n=NA, df=NA, s2=NA, n_obs = length(yo)))
  }


  if(method == "numeric_grad"){
    opt_res <- optim(par = c(beta_init, sigma2_init), function(par){
      beta <- par[seq_len(ncol(X))]
      sigma2 <- par[p+1]
      if(sigma2 <= 0) return(10000)


      zetastar <- zeta * sqrt(1 + sigma2/zeta^2)
      val <- 0 +
        sum(dnorm(Xo %*% beta, yo, sd=sqrt(sigma2), log=TRUE)) +
        sum(invprobit(Xm %*% beta, rho, zetastar, log=TRUE))
      -val
    }, hessian=TRUE)
    if(opt_res$convergence != 0){
      warning("Model didn't not properly converge\n")
      print(y)
      print(opt_res)
      return(list(beta=rep(NA, p), n=NA, df=NA, s2=NA, n_obs = length(yo)))
    }

    fit_beta <- opt_res$par[seq_len(p)]
    fit_sigma2 <- opt_res$par[p+1]
    fit_sigma2_var <- 1/opt_res$hessian[p+1, p+1]   # Get the last element on the diagonal
  }else if(method == "analytic_grad"){
    lik_res <- maxLik::maxNR(function(par){
      beta <- par[seq_len(p)]
      sigma2 <- par[p+1]
      if(sigma2 <= 0) return(-Inf)
      pred <- c(X %*% beta)

      zetastar <- zeta * sqrt(1 + sigma2/zeta^2)
      sum(dnorm(pred, y, sd=sqrt(sigma2), log=TRUE), na.rm=TRUE) +
        sum(invprobit(pred[is.na(y)], rho, zetastar, log=TRUE))
    }, grad = function(par){
      beta <- par[seq_len(p)]
      sigma2 <- par[p+1]
      if(sigma2 <= 0) return(-Inf)
      zetastar <- zeta * sqrt(1 + sigma2/zeta^2)
      imr <- inv_mills_ratio(Xm %*% beta, rho, zetastar)

      dbeta_o <- - (t(Xo) %*% Xo %*% beta - t(Xo) %*% yo) / sigma2
      dbeta_m <- t(Xm) %*% imr

      dsig2_o <- sum(((Xo %*% beta - yo)^2 - sigma2) / (2 * sigma2^2))
      dsig2_m <- -sum((Xm %*% beta - rho) / (2 * zetastar^2) * imr)

      c(dbeta_o + dbeta_m, dsig2_o + dsig2_m)
    }, hess = function(par){
      beta <- par[seq_len(p)]
      sigma2 <- par[p+1]
      if(sigma2 <= 0) return(NA)
      zetastar <- zeta * sqrt(1 + sigma2/zeta^2)
      imr <- inv_mills_ratio(Xm %*% beta, rho, zetastar)

      dbb_o <- -2 * t(Xo) %*% Xo / (2 * sigma2)
      dbb_m <- - t(Xm) %*% diag(c((imr^2 + (Xm %*% beta - rho) / zetastar^2 * imr)), nrow(Xm)) %*% Xm

      dss_o <- sum((sigma2 - 2 * (Xo %*% beta - yo)^2) / (2 * sigma2^3))
      dss_m <- sum((Xm %*% beta - rho) / (4 * zetastar^4) * imr *
                     (3 - (Xm %*% beta - rho) * imr - (Xm %*% beta - rho)^2 / zetastar^2))

      dbs_o <- t(Xo) %*% (Xo %*% beta - yo) / sigma2^2
      dbs_m <- t(Xm) %*% ((Xm %*% beta - rho) / (2 * zetastar^2) * imr^2 -
                            (zetastar^2 - (Xm %*% beta - rho)^2) / (2 * zetastar^4) * imr)

      res <- matrix(NA, nrow=length(par), ncol=length(par))
      res[seq_along(beta), seq_along(beta)] <- dbb_o + dbb_m
      res[length(par), length(par)] <- dss_o + dss_m
      res[length(par), seq_along(beta)] <- c(dbs_o + dbs_m)
      res[seq_along(beta), length(par)] <- c(dbs_o + dbs_m)
      res
    }, start = c(beta_init, sigma2_init),
    finalHessian = TRUE)

    if(lik_res$code > 2){
      warning("Model didn't not properly converge\n")
      print(y)
      print(lik_res)
      return(list(beta=rep(NA, p), n=NA, df=NA, s2=NA, n_obs = length(yo)))
    }
    fit_beta <- lik_res$estimate[seq_len(p)]
    fit_sigma2 <- lik_res$estimate[p+1]
    fit_sigma2_var <- -1/lik_res$hessian[p+1, p+1]   # Get the last element on the diagonal
  }else if(method == "lm"){
    lin_res <- lm(yo ~ Xo - 1)

    fit_beta <- unname(coef(lin_res))
    fit_sigma2 <- summary(lin_res)$sigma^2

    n_approx <- sum(summary(lin_res)$df[1:2])   # See help for summary.lm why
    s2_approx <- summary(lin_res)$sigma^2

    return(list(beta=fit_beta, n=n_approx, df=n_approx - p, s2=s2_approx, n_obs = length(yo)))
  }


  n_approx <- 2 * fit_sigma2^2 / fit_sigma2_var
  rss_approx <- 2 * fit_sigma2^3 / fit_sigma2_var
  s2_approx <- rss_approx / (n_approx - p)

  if(s2_approx < 0 || n_approx <= p){
    df_mod <- 1e-3
    s2_approx <- sqrt(fit_sigma2_var / df_mod^2 * (df_mod + 2)^2 * (df_mod + 2) / 2)
    n_approx <- df_mod + p
    rss_approx <- s2_approx * df_mod
  }
  names(fit_beta) <- colnames(X)

  list(coefficients=fit_beta,
       n_approx=n_approx, df=n_approx - p,
       s2=s2_approx, rss = rss_approx, n_obs = length(yo))
}


