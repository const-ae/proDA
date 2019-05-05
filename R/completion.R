


.methods_to_suggest <- c("abundances", "hyper_parameters", "feature_parameters", "coefficients",
                         "convergence", "design", "reference_level", "colData", "rowData")

#' @rdname cash-proDAFit-method
#' @export
.DollarNames.proDAFit <- function(x, pattern = ""){
  grep(pattern, .methods_to_suggest, value = TRUE)
}

#' Fluent use of accessor methods
#'
#' The 'proDAFit' object overwrites the dollar function to make it easy
#' to call functions to access values inside the object. This has the
#' advantage that it is very easy to discover the relevant methods
#' but nonetheless have an isolated implementation. Unlike the
#' \code{`@`} operator which directly accesses the underlying implementation,
#' the \code{`$`} operator only exposes a limited set of functions
#' \itemize{
#'   \item abundances
#'   \item hyper_parameters
#'   \item feature_parameters
#'   \item coefficients
#'   \item convergence
#'   \item design
#'   \item reference_level
#'   \item colData
#'   \item rowData.
#' }
#'
#' @param x an object of class 'proDAFit' produced by \code{proDA()}
#' @param name one of the functions listed above
#' @param value \strong{Warning:} modifying the content of a 'proDAFit'
#'   object is not allowed
#' @param pattern the regex pattern that is provided by the IDE
#'
#' @return whatever the function called \code{name} returns.
#'
#' @seealso \link{accessor_methods} for more documentation on the
#'   accessor functions.
#' @aliases dollar_methods
setMethod("$", "proDAFit",
function(x, name){
  if(! name %in% .methods_to_suggest){
    stop("Illegal name after '$' sign: ", name)
  }
  getGeneric(name)(x)
})

#' @rdname cash-proDAFit-method
setReplaceMethod("$", "proDAFit",
function(x, name, value){
  # if(name %in% c("rowData", "feature_parameters")){
  #   getGeneric(paste0(name, "<-"))(x, value = value)
  # }else{
  #   stop("It is illegal to modify the content of ", name)
  # }
  stop("It is illegal to modify the content of proDAFit object")
})

