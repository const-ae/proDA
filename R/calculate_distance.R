

#' Distance method for 'proDAFit' object
#'
#' The method calculates either the euclidean distance between samples or proteins
#' taking into account the missing values and the associated uncertainty.
#' Because with missing value no single deterministic distance can be
#' calculated two objects are returned: the mean and the associated
#' standard deviation of the distance estimates.
#'
#' @param object the 'proDAFit' object for which we calculate the distance or a matrix like
#'   object for which 'proDAFit' is created internally
#' @param by_sample a boolean that indicates if the distances is calculated
#'   between the samples (`by_sample = TRUE`) or between the proteins
#'   (`by_sample = FALSE`). Default: `TRUE`
#' @param blind fit an intercept model for the missing values to make
#'   sure that the results are not biased for the expected result.
#'   Default: `TRUE`
#' @param ... additional argument to \code{proDA()} in case object is a
#'   \code{SummarizedExperiment} or a \code{matrix}
#'
#'
#' @return a list with two elements: `mean` and `sd` both are formally
#'   of class "dist"
#'
#'
#' @name dist_approx_impl
#' @aliases dist_approx,proDAFit-method
#' @export
setMethod("dist_approx", signature = "proDAFit", function(object,
                                                          by_sample=TRUE,
                                                          blind = TRUE){

  if(blind){
    alt_fit <- predict(object, newdesign = ~ 1, type="feature_parameters")
  }else{
    alt_fit <- object
  }

  dist_approx_mean_var(Y = abundances(alt_fit),
                       Pred = predict(alt_fit),
                       X  = design(alt_fit),
                       coef_var = coefficient_variance_matrices(alt_fit),
                       by_sample = by_sample)

})


#' @rdname dist_approx_impl
#' @export
setMethod("dist_approx", signature = "SummarizedExperiment", function(object,
                                                                      by_sample=TRUE,
                                                                      blind = TRUE, ...){
  fit <- proDA(object, ...)
  dist_approx(fit, by_sample = by_sample, blind = blind)
})


#' @rdname dist_approx_impl
#' @export
setMethod("dist_approx", signature = "ANY", function(object,
                                                     by_sample=TRUE,
                                                     blind = TRUE, ...){
  fit <- proDA(object, ...)
  dist_approx(fit, by_sample = by_sample, blind = blind)
})




dist_approx_mean_var <- function(Y, Pred, X, coef_var, by_sample = FALSE){
  stopifnot(ncol(Y) == ncol(Pred), nrow(Y) == nrow(Pred), ncol(Y) == nrow(X),
            nrow(Y) == length(coef_var))
  Mu <- ifelse(! is.na(Y),
               Y,
               Pred)

  Pred_var <- mply_dbl(seq_len(nrow(Y)), function(i){
    sapply(seq_len(nrow(X)), function(j) t(X[j,]) %*% coef_var[[i]] %*% X[j,])
  }, ncol=ncol(Y))
  Mu_var <- ifelse(! is.na(Y),
                   0,
                   Pred_var)

  if(by_sample){
    Y <- t(Y)
    Mu <- t(Mu)
    Mu_var <- t(Mu_var)
  }

  mu_res <- matrix(0, nrow=nrow(Y), ncol=nrow(Y),
                   dimnames = list(rownames(Y), rownames(Y)))
  var_res <- matrix(0, nrow=nrow(Y), ncol=nrow(Y),
                    dimnames = list(rownames(Y), rownames(Y)))

  for(idx in seq_len(nrow(Y)-1)){
    for(idx2 in idx:nrow(Y)){
      ds <- distance_sq(Mu[idx, ], Mu_var[idx, ], Mu[idx2, ], Mu_var[idx2, ])
      mu_res[idx2, idx] <- sqrt(ds$mean)
      var_res[idx2, idx] <- ds$var  * 1/(4 * ds$mean)
    }
  }
  list(mean=as.dist(mu_res), sd=as.dist(sqrt(var_res)))
}



#' Square distance between two Gaussian distributions
#'
#' The function takes the mean and the diagonal of the covariance matrix
#' as vector and calculates the mean and variance of their distance distribution.
#' The formulas are based on [1] page 53.
#'
#' @return a list with elements `mean` and `var`
#'
#' 1. Mathai, A. & Provost, S. Quadratic Forms in Random Variables. (1992).
#' @keywords internal
distance_sq <- function(mu1, sigma1, mu2, sigma2){
  mu <- mu2 - mu1
  sigma <- sigma1 + sigma2

  analyt_mean <- sum(sigma) + sum(mu^2)
  analyt_var <- 2 *  sum(sigma^2) + 4 * sum(mu^2 * sigma)
  list(mean=analyt_mean, var=analyt_var)
}

