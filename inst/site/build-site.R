# This file modifies the internal pkgdown:::build_rmarkdown_format function
# to make sure the `df_print` is set to "paged"
# Copied from https://github.com/crew102/slowraker/blob/146f442085f652d824177e91fd1c38b29802b621/inst/site/build-site.R


swap_render_fun <- function() {

  # Alter pkgdown:::build_rmarkdown_format
  build_rmarkdown_format2 <- function(pkg = ".",
                                      name,
                                      depth = 1L,
                                      data = list(),
                                      toc = TRUE) {
    template <- pkgdown:::rmarkdown_template(pkg, name, depth = depth, data = data)

    out <- rmarkdown::html_document(
      toc = toc,
      toc_depth = rlang::`%||%`(pkg$meta$toc$depth, 2),
      self_contained = FALSE,
      theme = NULL,
      template = template$path,
      df_print = "tibble"
    )
    out$knitr$opts_chunk <- pkgdown:::fig_opts_chunk(pkg$figures, out$knitr$opts_chunk)

    attr(out, "__cleanup") <- template$cleanup

    out
  }
  assignInNamespace(
    "build_rmarkdown_format", build_rmarkdown_format2, ns = "pkgdown"
  )
}

build_site <- function(...) {

  # Change render function in pkgdown so it uses paged df printing
  swap_render_fun()

  # Build site
  pkgdown::build_site(..., new_process = FALSE)
}
