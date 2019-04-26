


.methods_to_suggest <- c("abundances", "hyper_parameters", "feature_parameters", "coefficients",
                         "convergence", "design", "reference_level", "colData", "rowData")


#' @export
.DollarNames.proDAFit <- function(x, pattern = ""){
  grep(pattern, .methods_to_suggest, value = TRUE)
}


setMethod("$", "proDAFit",
function(x, name){
  if(! name %in% .methods_to_suggest){
    stop("Illegal name after '$' sign: ", name)
  }
  getGeneric(name)(x)
})

setReplaceMethod("$", "proDAFit",
function(x, name, value){
  # if(name %in% c("rowData", "feature_parameters")){
  #   getGeneric(paste0(name, "<-"))(x, value = value)
  # }else{
  #   stop("It is illegal to modify the content of ", name)
  # }
  stop("It is illegal to modify the content of proDAFit object")
})

