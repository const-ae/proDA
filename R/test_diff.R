


#' Main proDA function to identify significant parameters
#'
#'
#' @export
test_diff <- function(parameters, contrast = `???`){

  cntrst <- parse_contrast(contrast, levels =  colnames(design(parameters)),
                 reference_level = reference_level(parameters),
                 direct_call = FALSE)


  # stop("Not yet implemented")

}



parse_contrast <- function(contrast, levels, reference_level = NULL, direct_call=TRUE){

  env <- if(direct_call){
    environment()
  }else{
    parent.frame()
  }

  cnt_capture <- substitute(contrast, env = env)


  stopifnot(! is.null(levels))
  if(is.factor(levels)){
    levels <- levels(levels)
  }else if(! is.character(levels)){
    stop("levels must be either a character vector or a factor")
  }

  indicators <- diag(nrow=length(levels))
  rownames(indicators) <- levels
  colnames(indicators) <- levels

  level_environment <- new.env()

  for(lvl in levels){
    assign(lvl, indicators[, lvl], level_environment)
  }
  if(! is.null(reference_level)){
    assign(reference_level, NA_real_, level_environment)
  }

  res <- eval(cnt_capture, envir= level_environment, parent.frame(n = 1 + direct_call))
  if(! is.numeric(res)){
    if(is.character(res)){
      # If contrast was a string, eval will just spit it out the same way
      res <- eval(parse(text = cnt_capture), envir= level_environment, parent.frame(n = 1 + direct_call))
    }
  }
  browser()
  if(any(is.na(res))){
    stop("Error the reference_level ", reference_level, " was used in the contrast argument. \n",
         "But it is already absorbed into the intercept, so it is not necessary to subtract it here.")
  }
  res[levels]
}


