

#' apply function that always returns a numeric matrix
#'
#' The function is modelled after `vapply`, but always returns a matrix
#' with one row for each iteration. You need to provide the number
#' of elements each function call produces beforehand (i.e. the number of
#' resulting columns). For a more flexible version where you don't need to
#' provide the number of columns see \code{\link{msply_dbl}}
#' @describeIn mply_dbl apply function that always returns a numeric matrix
#'
#' @param x a vector that will be passed to `vapply` or a matrix that will be
#'   passed to apply with \code{MARGIN=1}.
#' @param FUN the function that returns a vector of length ncol
#' @param ncol the length of the vector returned by `FUN`.
#' @param ... additional arguments to FUN
#'
#' @return a matrix of size \code{length(x) x ncol}
mply_dbl <- function(x, FUN, ncol=1, ...){
  if(is.vector(x)){
    res <- vapply(x, FUN, FUN.VALUE=rep(0.0, times=ncol), ...)
  }else{
    res <- apply(x, 1, FUN, ...) * 1.0
    if(nrow(x) > 0 && length(res) == 0){
      # Empty result, make matrix
      res <- matrix(numeric(0),nrow=0, ncol=nrow(x))
    }else if(nrow(x) == 0){
      res <- matrix(numeric(0), nrow=ncol, ncol=0)
    }
    if((ncol == 1 && ! is.vector(res)) || (ncol > 1 && nrow(res) != ncol)){
      stop(paste0("values must be length ", ncol,
                  ", but result is length ", nrow(res)))
    }
  }

  if(ncol == 1){
    as.matrix(res, nrow=length(res), ncol=1)
  }else{
    t(res)
  }
}

#' @describeIn mply_dbl flexible version that automatically infers the number
#'   of columns
msply_dbl <- function(x, FUN, ...){
  if(is.vector(x)){
    res <- sapply(x, FUN, ...)
  }else{
    res <- apply(x, 1, FUN, ...)
  }

  if(is.list(res)){
    if(all(vapply(res, function(x) is.numeric(x) && length(x) == 0, FUN.VALUE = FALSE))){
      res <- matrix(numeric(0),nrow=0, ncol=length(res))
    }else{
      stop("Couldn't simplify result to a matrix")
    }
  }
  if(is.matrix(x) && length(res) == 0){
    # Empty result, make matrix
    res <- matrix(numeric(0),nrow=0, ncol=nrow(x))
  }

  if(is.numeric(res)){
    # Do nothing
  }else if(is.logical(res)){
    res <- res * 1.0
  }else{
    stop(paste0("Result is of type ", typeof(res), ". Cannot convert to numeric."))
  }

  if(is.matrix(res)){
    t(res)
  }else{
    as.matrix(res, nrow=length(res))
  }
}


#' Helper function that makes sure that NA * 0 = 0 in matrix multiply
#'
#' @keywords internal
`%zero_dom_mat_mult%` <- function(X, Y){
  X[is.infinite(X)] <- NA
  Y[is.infinite(Y)] <- NA
  X_cp <- X
  X_cp[is.na(X_cp)] <- 0
  Y_cp <- Y
  Y_cp[is.na(Y_cp)] <- 0

  res <- X_cp %*% Y_cp
  mask1 <- (is.na(X)) %*% (is.na(Y) | Y != 0)
  mask2 <- (is.na(X) | X != 0) %*% (is.na(Y))
  res[mask1 + mask2 != 0] <- NA
  res
}



