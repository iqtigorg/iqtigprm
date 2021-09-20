#' Apply precomputations to data
#'
#' Apply a specified set of precomputation rules to preprocessed patient survey data
#'
#' @param data A \code{data.frame} containing the preprocessed patient survey data.
#' @param precomputations_list A \code{list} of named precomputations. Each component is also a \code{list}
#' containing a computation rule \code{.$expr} and value labels \code{.$label}.
#'
#' @return This takes the patient survey data and performs the specified precomputations by adding new columns equal
#' to the number of precomputations to the data set. Each new column is named after one precomputation and filled with values
#' as obtained from the corresponding expression \code{.$expr} applied (rowwise) to the data set.
#' The extended survey data set is returned as a \code{data.frame}.
#'
#' @export
apply_precomputations <- function(data, precomputations_list) {
  Reduce(function(acc, el) {
    precomputation <- precomputations_list[[el]]
    expr <- precomputation$expr
    labels <- precomputation$label
    prototype_field_name <- as.character(precomputation$prototype)
    if (is.null(precomputation$prototype)) prototype_field_name <- NULL

    ifelse <- function(test, yes, no) {
      val <- dplyr::if_else(test, as.numeric(yes), as.numeric(no))
      if (is.null(prototype_field_name)) {
        val
      } else {
        to_haven(
          acc[[prototype_field_name]],
          val,
          labels = labels
        )
      }
    }
    acc[[el]] <- withr::with_package("iqtigfunctions", {
      rlang::eval_tidy(expr, data = acc)
    })

    # Questions for computed fields should not have the original question as label
    # Also we need to make sure item levels can be overwritten
    if (inherits(acc[[el]], "haven_labelled")) {
      if ("label" %in% names(attributes(acc[[el]]))) {
        attr(acc[[el]], "label") <- "Funktion"
      }
      if (!is.null(labels)) {
        attr(acc[[el]], "labels") <- labels
      } else {
        if (!is.null(labels)) {
          acc[[el]] <- labelled::labelled(acc[[el]], labels)
        }
      }
    }
    acc
  }, names(precomputations_list), data)
}

# bring in haven format
to_haven <- function(base_vec, val, labels = NULL) {
  stopifnot(is.numeric(val))
  attributes(val) <- attributes(base_vec)
  if (!is.null(labels)) {
    stopifnot(is.numeric(labels), !is.null(names(labels)))
    attr(val, "labels") <- labels
  }
  val
}
