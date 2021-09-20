utils::globalVariables(c("KN_ID", "ID", "meta_unit", "ID_LE"))

#' Apply category mappings
#'
#' @description Create a new data set from raw survey data by applying specified mappings of certain variables
#'
#' @param raw_data A \code{data.frame} containing the raw patient survey data.
#' @param mappings A \code{list} containing the mapping specifications, i.e. points value mappings and affected data fields.
#'
#' @return Takes the survey data \code{dat} and applies the point mappings as listed in \code{mappings$point_mappings} and
#' the replacement mappings in \code{mappings$ausweichkategorien}.
#' Additional \code{ID} and \code{meta_unit} columns are also added to the output.
#'
#' @export
preprocess_data <- function(raw_data, mappings = NULL) {
  # ensure correct format
  raw_data <- ensure_recent_haven_version(raw_data)

  # we add an ID column and a meta unit for now
  data <- dplyr::mutate(
    raw_data,
    ID = dplyr::row_number(),
    meta_unit = ID_LE
  )

  # Here we do the following things
  # 1) Map answers to values between 0 and 100
  # 2) Convert Ausweichkategorien (neutral responses) to NA, but keep track of them
  for (ausweichkategorie in mappings$ausweichkategorien) {
    data <- dplyr::mutate(data, dplyr::across(tidyselect::any_of(ausweichkategorie$felder), remap_values, ausweichkategorie$value, NA_real_))
  }

  for (mapping in mappings$point_mappings) {
    data <- dplyr::mutate(data, dplyr::across(tidyselect::any_of(mapping$felder), remap_category, as.numeric(names(mapping$mapping)), mapping$mapping))
  }
  data
}

## bring data in correct haven format
ensure_recent_haven_version <- function(data) {
  # the data was saved with an old haven format
  # so we need to convert it by saving it and reading it again
  tmp_file <- tempfile()
  on.exit(fs::file_delete(tmp_file))
  haven::write_sav(data, tmp_file)
  haven::read_sav(tmp_file)
}

## perform category mapping
remap_category <- function(values, from, to) {
  if (!inherits(values, "haven_labelled")) {
    return(values)
  }
  if (!(is.numeric(values) && length(unique(values)) < 7)) {
    return(values)
  }
  categories <- attr(values, "labels")
  categories <- categories[categories >= 0 & categories <= 100]
  if (!all(stats::na.omit(values) %in% categories)) {
    return(values)
  }

  categories <- stats::setNames(
    to,
    sapply(from, function(from_val) {
      names(categories[categories == from_val])
    })
  )
  stopifnot(length(categories) == length(to))
  vals <- remap_values(values,
    from = from,
    to = to
  )
  stopifnot(length(setdiff(unique(vals), categories)) > 0)
  attr(vals, "labels") <- categories[categories >= 0 & categories <= 100]
  vals
}

## perform value mapping
remap_values <- function(values, from, to) {
  stopifnot(length(from) == length(to))
  Reduce(function(acc, i) {
    acc[acc %in% from[[i]]] <- to[[i]]
    acc
  }, seq_along(from), init = values)
}
