





#' Column wise median normalization of the data matrix
#'
#' The method calculates for each sample the median change (i.e. the difference
#' between the observed value and the row average) and substracts it from
#' each row. Missing values are ignored in the procedure. The method is based
#' on the assumption that a majority of the rows did not change.
#'
#'
#' @param X a matrix or SummarizedExperiment of proteins and samples
#' @return the normalized matrix
#'
#' @examples
#'   syn_data <- generate_synthetic_data(n_proteins = 10)
#'   normalized_data <- median_normalization(syn_data$Y)
#'   normalized_data
#'
#' @export
median_normalization <- function(X){
  if(is.matrix(X)){
    stopifnot(length(dim(X)) == 2)
    Xnorm <- X
    for(idx in seq_len(ncol(X))){
      Xnorm[, idx] <- X[, idx, drop=FALSE] -
        median(X[, idx, drop=FALSE] - rowMeans(X, na.rm=TRUE), na.rm=TRUE )
    }
    Xnorm
  }else if(is(X, "SummarizedExperiment")){
    new_assay <- median_normalization(SummarizedExperiment::assay(X))
    SummarizedExperiment::assay(X) <- new_assay
    X
  }else if(canCoerce(X, "SummarizedExperiment")){
    se <- as(X, "SummarizedExperiment")
    se_norm <- median_normalization(se)
    as(se_norm, class(X)[1])
  }else{
    stop("Cannot handle argument X of type", class(X))
  }
}

