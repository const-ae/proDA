

#' Main proDA function to determine the hyper and protein parameters
#'
#'
#' @export
fit_parameters <- function(data, design=~ 1,
                           col_data = NULL,
                           reference_class = NULL){

  # Validate Data
  stopifnot(is.matrix(data))
  n_samples <- ncol(data)
  n_rows <- nrow(data)


  # Handle the design parameter
  if(is.matrix(design)){
    model_matrix <- design
  }else if(is.vector(design) && length(design) == n_samples){
    model_matrix <- convert_chr_vec_to_model_matrix(design, reference_class)
  }else if(inherits(design,"formula")){
    if(design == formula(~ 1) && is.null(col_data)){
      col_data <- as.data.frame(matrix(numeric(0), nrow=10))
    }
    model_matrix <- convert_formula_to_model_matrix(design, col_data, reference_class)
  }else{
    stop(paste0("design argment of class ", class(design), " is not supported. Please ",
                "specify a `model_matrix`, a `character vector`, or a `formula`."))
  }
  check_valid_model_matrix(model_matrix, data)

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


convert_formula_to_model_matrix <- function(formula, col_data, reference_class=NULL){
  if(! is.null(reference_class)){
    has_ref_level <- vapply(col_data, function(x){
      any(x == reference_class)
    }, FUN.VALUE = FALSE)
    if(all(has_ref_level == FALSE)){
      stop("None of the columns contains the specified reference_class.")
    }
    col_data[has_ref_level] <- lapply(col_data[has_ref_level], relevel, ref = reference_class)
  }
  mm <- model.matrix(formula, col_data)
  colnames(mm)[colnames(mm) == "(Intercept)"] <- "Intercept"
  mm
}


