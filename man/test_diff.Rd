% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test_diff.R
\name{test_diff}
\alias{test_diff}
\alias{result_names,proDAFit-method}
\title{Identify differentially abundant proteins}
\usage{
test_diff(
  fit,
  contrast,
  reduced_model = ~1,
  alternative = c("two.sided", "greater", "less"),
  pval_adjust_method = "BH",
  sort_by = NULL,
  decreasing = FALSE,
  n_max = Inf,
  verbose = FALSE
)

\S4method{result_names}{proDAFit}(fit)
}
\arguments{
\item{fit}{an object of class 'proDAFit'. Usually, this is
produced by calling \code{proDA()}}

\item{contrast}{an expression or a string specifying which
contrast is tested. It can be a single coefficient (to see
the available options use \code{result_names(fit)}) or any
linear combination of them. The contrast is always compared
against zero. Thus, to find out if two coefficients differ
use \code{coef1 - coef2}. Remember if the coefficient is not
a valid identifier in R, to escape it using back ticks. For 
example if you test the interaction of A and B use 
\code{`A:B`}.}

\item{reduced_model}{If you don't want to test an individual
coefficient, you can can specify a reduced model and compare
it with the original model using an F-test. This is useful
to find out how a set of parameters affect the goodness of
the fit. If neither a \code{contrast}, nor
a \code{reduced_model} is specified, by default a comparison
with an intercept model (ie. just the average across conditions)
is done. Default: \code{~ 1}.}

\item{alternative}{a string that decides how the
hypothesis test is done. This parameter is only relevant for
the Wald-test specified using the `contrast` argument.
Default: \code{"two.sided"}}

\item{pval_adjust_method}{a string the indicates the method
that is used to adjust the p-value for the multiple testing.
It must match the options in \code{\link[stats]{p.adjust}}.
Default: \code{"BH"}}

\item{sort_by}{a string that specifies the column that is used
to sort the resulting data.frame. Default: \code{NULL} which
means the result is sorted by the order of the input matrix.}

\item{decreasing}{a boolean to indicate if the order is reversed.
Default: \code{FALSE}}

\item{n_max}{the maximum number of rows returned by the method.
Default: \code{Inf}}

\item{verbose}{boolean that signals if the method prints informative
messages. Default: \code{FALSE}.}
}
\value{
The `result_names()` function returns a character vector.

  The `test_diff()` function returns a \code{data.frame} with one row per protein
  with the key parameters of the statistical test. Depending what kind of test
  (Wald or F test) the content of the `data.frame` differs.

  The Wald test, which can considered equivalent to a t-test, returns
  a `data.frame` with the following columns:
  \describe{
    \item{name}{the name of the protein, extracted from the rowname of
      the input matrix}
    \item{pval}{the p-value of the statistical test}
    \item{adj_pval}{the multiple testing adjusted p-value}
    \item{diff}{the difference that particular coefficient makes. In
      differential expression analysis this value is also called
      log fold change, which is equivalent to the difference on the
      log scale.}
    \item{t_statistic}{the \code{diff} divided by the standard
      error \code{se}}
    \item{se}{the standard error associated with the \code{diff}}
    \item{df}{the degrees of freedom, which describe the amount
      of available information for estimating the \code{se}. They
      are the sum of the number of samples the protein was observed
      in, the amount of information contained in the missing values,
      and the degrees of freedom of the variance prior.}
    \item{avg_abundance}{the estimate of the average abundance of
      the protein across all samples.}
    \item{n_approx}{the approximated information available for estimating
      the protein features, expressed as multiple of the information
      contained in one observed value.}
    \item{n_obs}{the number of samples a protein was observed in}
  }


  The F-test returns a `data.frame` with the following columns
  \describe{
    \item{name}{the name of the protein, extracted from the rowname of
      the input matrix}
    \item{pval}{the p-value of the statistical test}
    \item{adj_pval}{the multiple testing adjusted p-value}
    \item{f_statistic}{the ratio of difference of normalized deviances
     from original model and the reduced model, divided by the
     standard deviation.}
    \item{df1}{the difference of the number of coefficients in the
      original model and the number of coefficients in the reduced
      model}
    \item{df2}{the degrees of freedom, which describe the amount
      of available information for estimating the \code{se}. They
      are the sum of the number of samples the protein was observed
      in, the amount of information contained in the missing values,
      and the degrees of freedom of the variance prior.}
    \item{avg_abundance}{the estimate of the average abundance of
      the protein across all samples.}
    \item{n_approx}{the information available for estimating
      the protein features, expressed as multiple of the information
      contained in one observed value.}
    \item{n_obs}{the number of samples a protein was observed in}
  }
}
\description{
The `test_diff()` function is used to test coefficients of a 'proDAFit'
object. It provides a Wald test to test individual
coefficients and a likelihood ratio F-test to compare the
original model with a reduced model. The \code{result_names}
method provides a quick overview which coefficients are
available for testing.
}
\details{
To test if coefficient is different from zero with a Wald
test use the \code{contrast} function argument. To test if two
models differ with an F-test use the \code{reduced_model}
argument. Depending on the test that is conducted, the functions
returns slightly different data.frames.

The function is designed to follow the principles of the
base R test functions (ie. \code{\link[stats]{t.test}} and
\code{\link[stats]{wilcox.test}}) and the functions designed
for collecting the results of high-throughput  testing
(ie. \code{limma::topTable} and \code{DESeq2::results}).
}
\examples{
  # "t-test"
  syn_data <- generate_synthetic_data(n_proteins = 10)
  fit <- proDA(syn_data$Y, design = syn_data$groups)
  result_names(fit)
  test_diff(fit, Condition_1 - Condition_2)


  suppressPackageStartupMessages(library(SummarizedExperiment))
  se <- generate_synthetic_data(n_proteins = 10,
                                n_conditions = 3,
                                return_summarized_experiment = TRUE)
  colData(se)$age <- rnorm(9, mean=45, sd=5)
  colData(se)
  fit <- proDA(se, design = ~ group + age)
  result_names(fit)
  test_diff(fit, "groupCondition_2",
            n_max = 3, sort_by = "pval")

  # F-test
  test_diff(fit, reduced_model = ~ group)


}
\seealso{
The contrast argument is inspired by
  \code{limma::makeContrasts}.
}
