





#' Column wise median normalization of the data matrix
#'
#' The method calculates for each sample the median change (i.e. the difference
#' between the observed value and the row average) and subtracts it from
#' each row. Missing values are ignored in the procedure. The method is based
#' on the assumption that a majority of the rows did not change.
#'
#'
#' @param X a matrix or SummarizedExperiment of proteins and samples
#' @param spike_in_rows a numeric or boolean vector that is used to
#'   to normalize the intensities across samples. Default: \code{NULL}
#'   which means that all rows are used.
#' @return the normalized matrix
#'
#' @examples
#'   syn_data <- generate_synthetic_data(n_proteins = 10)
#'   normalized_data <- median_normalization(syn_data$Y)
#'   normalized_data
#'   
#'   # If we assume that the first 5 proteins are spike-ins
#'   normalized_data2 <- median_normalization(syn_data$Y, spike_in_rows = 1:5)
#'
#' @export
median_normalization <- function(X, spike_in_rows = NULL){
  if(is.matrix(X)){
    stopifnot(length(dim(X)) == 2)
    if(is.null(spike_in_rows)){
      # Use all rows for normalization
      spike_in_rows <- seq_len(nrow(X))
    }else if(! is.numeric(spike_in_rows) && ! is.logical(spike_in_rows)){
      stop("spike_in_rows must either be a numeric or a logical vector")
    }else if(is.logical(spike_in_rows)){
      if(length(spike_in_rows) != nrow(X)){
        stop("The spike_in_rows logical vector must have one entry for each row.")
      }else if(all(spike_in_rows == FALSE)){
        stop("Not all elements of the spike_in_rows vector can be FALSE.")
      }
    }


    Xnorm <- X
    for(idx in seq_len(ncol(X))){
      Xnorm[, idx] <- X[, idx, drop=FALSE] -
        median(X[spike_in_rows, idx, drop=FALSE] - rowMeans(X[spike_in_rows, , drop=FALSE], na.rm=TRUE), na.rm=TRUE )
    }
    Xnorm
  }else if(is(X, "SummarizedExperiment")){
    new_assay <- median_normalization(SummarizedExperiment::assay(X), spike_in_rows)
    SummarizedExperiment::assay(X) <- new_assay
    X
  }else if(canCoerce(X, "SummarizedExperiment")){
    se <- as(X, "SummarizedExperiment")
    se_norm <- median_normalization(se, spike_in_rows)
    as(se_norm, class(X)[1])
  }else{
    stop("Cannot handle argument X of type", class(X))
  }
}

