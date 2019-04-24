

#' Main proDA function to determine the hyper and protein parameters
#'
#'
#' @export
fit_parameters <- function(data, design=`???`,
                           col_data = NULL,
                           reference_class = NULL){

  # Validate Data
  stopifnot(is.matrix(data))
  n_samples <- ncol(data)
  n_rows <- nrow(data)


  model_matrix <- NULL
  # Handle the design parameter
  if(is.matrix(design)){
    model_matrix <- design
  }else if(is.vector(design) && length(design) == n_samples){
    model_matrix <- convert_chr_vec_to_model_matrix(design, reference_class)
  }else if(inherits(design,"formula")){


  }
  check_valid_model_matrix(design, data)

  stop("Not yet implemented")

}


check_valid_model_matrix <- function(matrix, data){
  stopifnot(is.matrix(matrix))
  stopifnot(nrow(matrix) == ncol(data))
}



convert_chr_vec_to_model_matrix <- function(design, reference_class){
  if(! is.factor(design)){
    design_fct <- as.factor(design)
  }else{
    design_fct <- design
  }
  if(is.null(reference_class)){
    helper_df <- data.frame(x_ = design_fct)
    mm <- model.matrix(~ x_ - 1, helper_df)
    colnames(mm) <- sub("^x_", "", colnames(mm))
  }else{
    design_fct <- relevel(design_fct, ref = reference_class)
    helper_df <- data.frame(x_ = design_fct)
    mm <- model.matrix(~ x_ + 1, helper_df)
    colnames(mm)[-1] <- paste0(sub("^x_", "", colnames(mm)[-1]), "_vs_", reference_class)
  }
  colnames(mm)[colnames(mm) == "(Intercept)"] <- "Intercept"
  mm
}


convert_formula_to_model_matrix <- function(formula, col_data){
  mm <- model.matrix(formula, col_data)
  colnames(mm)[colnames(mm) == "(Intercept)"] <- "Intercept"
  mm
}


