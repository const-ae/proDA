
#' Generate a dataset according to the probabilistic dropout model
#'
#'
#' @param n_proteins the number of rows in the dataset
#' @param n_conditions the number of conditions. Default: 2
#' @param n_replicates the number of replicates per condition.
#'   Can either be a single number or a vector with
#'   \code{length(n_replicates) == n_conditions}. Default: 3
#' @param frac_changed the fraction of proteins that actually
#'   differ between the conditions. Default: 0.1
#' @param dropout_curve_position the point where the chance
#'   to observe a value is 50%. Can be a single number or
#'   a vector of \code{length(dropout_curve_position) == n_conditions * n_replicates}.
#'   Default: 18.5
#' @param dropout_curve_scale The width of the dropout curve.
#'   Negative numbers mean that lower intensities are more likely
#'   to be missing.
#'   Can be a single number or a vector of
#'   \code{length(dropout_curve_position) == n_conditions * n_replicates}.
#'   Default: -1.2
#' @param location_prior_mean,location_prior_scale
#'   the position and the variance around which the individual
#'   condition means (\code{t_mu}) scatter. Default: 20 and 4
#' @param variance_prior_scale,variance_prior_df
#'    the scale and the degrees of freedom of the inverse
#'    Chi-squared distribution used as a prior for the
#'    variances. Default: 0.05 and 2
#' @param effect_size the standard deviation that used to draw
#'   different values for the \code{frac_changed} part of the
#'   proteins. Default: 2
#' @param return_summarized_experiment a boolean indicator if
#'   the method should return a \code{SummarizedExperiment}
#'   object instead of a list. Default: \code{FALSE}
#'
#' @return a list with the following elements
#'   \describe{
#'     \item{Y}{the intensity matrix including the missing values}
#'     \item{Z}{the intensity matrix before dropping out values}
#'     \item{t_mu}{a matrix with \code{n_proteins} rows and
#'        \code{n_conditions} columns that contains the underlying
#'        means for each protein.}
#'     \item{t_sigma2}{a vector with the true variance for each
#'        protein.}
#'     \item{changed}{a vector with boolean values if the
#'        protein is actually changed.}
#'     \item{group}{the group structure mapping samples to conditions}
#'   }
#'   if \code{return_summarized_experiment} is \code{FALSE}. Otherwise
#'   returns a \code{SummarizedExperiment} with the same information.
#'
#' @examples
#'   syn_data <- generate_synthetic_data(n_proteins = 10)
#'   names(syn_data)
#'   se <- generate_synthetic_data(n_proteins = 10, return_summarized_experiment = TRUE)
#'   se
#'
#' @export
generate_synthetic_data <- function(n_proteins, n_conditions = 2,
                                    n_replicates = 3,
                                    frac_changed = 0.1,
                                    dropout_curve_position = 18.5,
                                    dropout_curve_scale = -1.2,
                                    location_prior_mean = 20,
                                    location_prior_scale = 4,
                                    variance_prior_scale  = 0.05,
                                    variance_prior_df = 2,
                                    effect_size = 2,
                                    return_summarized_experiment = FALSE) {

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


  t_sigma2 <- extraDistr::rinvchisq(n_proteins, nu = variance_prior_df, tau = variance_prior_scale)


  n_changed <- round(n_proteins * frac_changed)
  t_mu <- matrix(rep(rnorm(n_proteins-n_changed, location_prior_mean,
                          sd=sqrt(location_prior_scale)), times=n_conditions),
                ncol=n_conditions)

  t_mu <- rbind(t_mu, mply_dbl(seq_len(n_changed), ncol = n_conditions, function(idx){
    avg <- rnorm(1, mean=location_prior_mean, sd=sqrt(location_prior_scale))
    eff <- rnorm(n_conditions, mean = 0, sd = effect_size)
    avg + eff
  }))


  Z <- do.call(cbind, lapply(as.numeric(groups), function(cond){
    rnorm(n_proteins, t_mu[, cond], sd=sqrt(t_sigma2))
  }))

  Y <- t(vapply(seq_len(n_proteins), function(idx){
    y <- Z[idx, ]
    dropout_prob <-  invprobit(y, dropout_curve_position, dropout_curve_scale)
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

  if(return_summarized_experiment){
    cdf <- data.frame(group = groups,
                      true_dropout_curve_position = dropout_curve_position,
                      true_dropout_curve_scale = dropout_curve_scale)
    rownames(cdf) <- colnames(Y)
    rdf <- data.frame(changed = changed,
                      true_s2 = t_sigma2)
    feat_df <- as.data.frame(t_mu)
    colnames(feat_df) <- paste0("true_", colnames(feat_df))
    rdf <- cbind(rdf, feat_df)
    rownames(rdf) <- rownames(Y)

    SummarizedExperiment(assays = list(abundances = Y,
                                       full_observations = Z),
                         colData = cdf,
                         rowData = rdf)
  }else{
    return(list(Y=Y, Z=Z, t_mu=t_mu, t_sigma2=t_sigma2,
                changed=changed,
                groups = groups))
  }
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
#'   all(proDA:::as_replicate(x) == c(1,1,2,2,3,1))
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




