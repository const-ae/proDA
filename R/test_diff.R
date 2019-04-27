


#' Main proDA function to identify significant parameters
#'
#'
#' @export
test_diff <- function(parameters, contrast,
                      alternative = c("two.sided", "greater", "less"),
                      pval_adjust_method = "BH"){

  alternative <- match.arg(alternative, c("two.sided", "greater", "less"))

  cntrst <- parse_contrast(contrast, levels =  colnames(design(parameters)),
                 reference_level = reference_level(parameters),
                 direct_call = FALSE)

  dm <- design(parameters)
  diff <- coefficients(parameters) %*% cntrst
  diff_var <- feature_parameters(parameters)$s2 * c(t(cntrst) %*% solve(t(dm) %*% dm) %*% cntrst)
  diff_df <- feature_parameters(parameters)$df

  t_stat <- diff / sqrt(diff_var)
  pval <- pt(t_stat, df=diff_df)
  if(alternative == "two.sided"){
    pval <- 2 * pmin(pval, 1-pval)
  }else if(alternative == "greater"){
    pval <- 1 - pval
  }

  avg_abundance <- rowMeans(coefficients(parameters) %*% t(design(parameters)))

  data.frame(name = rownames(parameters),
             pval = pval,
             adj_pval = p.adjust(pval, method = pval_adjust_method),
             diff = diff,
             avg_abundance = avg_abundance,
             t_statistic = t_stat,
             se = sqrt(diff_var),
             df = diff_df,
             n_approx = feature_parameters(parameters)$n_approx,
             n_obs = feature_parameters(parameters)$n_obs,
             stringsAsFactors = FALSE)

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
    ind <- indicators[, lvl]
    names(ind) <- levels
    assign(lvl, ind, level_environment)
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
  if(any(is.na(res))){
    stop("Error the reference_level ", reference_level, " was used in the contrast argument. \n",
         "But it is already absorbed into the intercept, so it is not necessary to subtract it here.")
  }
  res
}


