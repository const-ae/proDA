


#' Get the abundance matrix
#'
#' @param object the object to get from
#' @param ... additional arguments used by the concrete implementation
#'
#' @seealso \link{accessor_methods} for the implementation for a 'proDAFit' object
#' @export
setGeneric("abundances", function(object, ...) standardGeneric("abundances"))

#' Get the hyper parameters
#'
#' @param object the object to get from
#' @param ... additional arguments used by the concrete implementation
#'
#' @seealso \link{accessor_methods} for the implementation for a 'proDAFit' object
#' @export
setGeneric("hyper_parameters", function(object, ...) standardGeneric("hyper_parameters"))

#' Get the feature parameters
#'
#' @param object the object to get from
#' @param ... additional arguments used by the concrete implementation
#'
#' @seealso \link{accessor_methods} for the implementation for a 'proDAFit' object
#' @export
setGeneric("feature_parameters", function(object, ...) standardGeneric("feature_parameters"))

#' Get the coefficients
#'
#' @param object the object to get from
#' @param ... additional arguments used by the concrete implementation
#'
#' @seealso \link{accessor_methods} for the implementation for a 'proDAFit' object
#' @export
setGeneric("coefficients", function(object, ...) standardGeneric("coefficients"))

#' Get the reference level
#'
#' @param object the object to get from
#' @param ... additional arguments used by the concrete implementation
#' @seealso \link{accessor_methods} for the implementation for a 'proDAFit' object
#' @export
setGeneric("reference_level", function(object, ...) standardGeneric("reference_level"))


#' Get the convergence information
#'
#' @param object the object to get from
#' @param ... additional arguments used by the concrete implementation
#'
#' @seealso \link{accessor_methods} for the implementation for a 'proDAFit' object
#' @export
setGeneric("convergence", function(object, ...) standardGeneric("convergence"))


#' Calculate an approximate distance for 'object'
#'
#' @param object the object for which the distance is approximated
#' @param ... additional arguments used by the concrete implementation
#'
#' @seealso \code{\link[stats]{dist}} for the base R function and
#'   \code{\link[proDA:dist_approx,proDAFit-method]{dist_approx()}} concrete implementation
#'   for 'proDAFit' objects
#' @export
setGeneric("dist_approx", function(object, ...) standardGeneric("dist_approx"))


#' Get the result_names
#'
#' @param fit the fit to get the result_names from
#' @param ... additional arguments used by the concrete implementation
#'
#' @export
setGeneric("result_names", function(fit, ...) standardGeneric("result_names"))

