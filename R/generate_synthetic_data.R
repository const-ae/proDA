

generate_synthetic_data <- function(n_proteins, n_conditions = 2,
                                    n_replicates = 3,
                                    frac_changed = 0.1,
                                    dropout_curve_position = 18.5,
                                    dropout_curve_scale = -1.2,
                                    location_prior_mean = 20,
                                    location_prior_scale = 4,
                                    variance_prior_scale  = 0.05,
                                    variance_prior_df = 2) {

  if(length(n_replicates) == 1){
    n_replicates <- rep(n_replicates, n_conditions)
  }else if(length(n_replicates) != n_conditions){
    stop("The number of elements in n_replicates must match n_conditions")
  }

  groups <- unlist(lapply(seq_len(length(n_replicates)), function(idx){
    rep(paste0("Condition_", idx), times = n_replicates[idx])
  }))
  groups <- factor(groups, levels = unique(groups))

  if(length(dropout_curve_position) == 1){
    dropout_curve_position <- rep(dropout_curve_position, times = sum(n_replicates))
  }else if(length(dropout_curve_position) != sum(n_replicates)){
    stop("length(dropout_curve_position)=", length(dropout_curve_position)," must either be one or match the",
         "  number of samples (sum(n_replicates)=", sum(n_replicates), ")")
  }
  if(length(dropout_curve_scale) == 1){
    dropout_curve_scale <- rep(dropout_curve_scale, times = sum(n_replicates))
  }else if(length(dropout_curve_scale) != sum(n_replicates)){
    stop("length(dropout_curve_scale)=", length(dropout_curve_scale)," must either be one or match the",
         " number of samples (sum(n_replicates)=", sum(n_replicates), ")")
  }


  t_sigma2 <- rchisq(n_proteins, df=variance_prior_df) *
    variance_prior_df * variance_prior_scale


  n_changed <- round(n_proteins * frac_changed)
  t_mu <- matrix(rep(rnorm(n_proteins-n_changed, location_prior_mean,
                          sd=sqrt(location_prior_scale)), times=n_conditions),
                ncol=n_conditions)
  t_mu <- rbind(t_mu, t(mply_dbl(seq_len(n_conditions), ncol=n_changed, function(cond){
    rnorm(n_changed, location_prior_mean, sd=sqrt(location_prior_scale))
  })))


  Z <- do.call(cbind, lapply(as.numeric(groups), function(cond){
    rnorm(n_proteins, t_mu[, cond], sd=sqrt(t_sigma2))
  }))

  Y <- t(vapply(seq_len(n_proteins), function(idx){
    y <- Z[idx, ]
    dropout_prob <-  invprobit(y, dropout_curve_position[idx], dropout_curve_scale[idx])
    y[runif(length(y)) < dropout_prob] <- NA
    y
  }, FUN.VALUE = rep(0, length(groups))))

  changed <- c(rep(FALSE, n_proteins-n_changed), rep(TRUE, n_changed))


  colnames(Y) <- paste0(groups, "-", as_replicate(groups))
  colnames(Z) <- paste0(groups, "-", as_replicate(groups))
  colnames(t_mu) <- levels(groups)


  prot_names <- if(nrow(Y) == 0) character(0) else paste0("protein_", seq_len(nrow(Y)))
  rownames(Y) <- prot_names
  rownames(Z) <- prot_names
  rownames(t_mu) <- prot_names
  names(t_sigma2) <- prot_names
  names(changed) <- prot_names

  return(list(Y=Y, Z=Z, t_mu=t_mu, t_sigma2=t_sigma2,
              changed=changed,
              groups = groups))

}






#' Get numeric vector with the count of the replicate for each element
#'
#' For a vector with repeated values return a vector where each element
#' is the count how often the element was observed previously
#'
#' @param x a vector with repeated elements
#'
#' @return numeric vector
#'
#' @keywords internal
#'
#' @seealso \code{\link[base]{order}}, \code{\link[base]{rank}}
#' @examples
#'   x <- c("a", "b", "a", "b", "b", "d")
#'   all(proDD:::as_replicate(x) == c(1,1,2,2,3,1))
#'
as_replicate <- function(x){
  counter <- as.list(rep(1, times=length(unique(x))))
  names(counter) <- unique(x)

  replicate <- rep(NA, times=length(x))
  for(idx in seq_along(x)){
    replicate[idx] <- counter[[x[idx]]]
    counter[[x[idx]]] <- counter[[x[idx]]] + 1
  }
  replicate
}




