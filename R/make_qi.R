#' Construct a QIDB data pair
#'
#' @description Generate a full set of computable QIs joined by corresponding data.
#'
#' @param qidb A \code{list} containing specification of multiple QIs. Each list component represents the raw information to one QI.
#' @param data A \code{data.frame} containing the pre-processed patient survey data the QI should be evaluated on.
#'
#' @return The patient survey data is joined to every qi from the provided \code{qidb}, resulting in a \code{list} of objects of class \code{qi}.
#' @export
connect_qidb2data <- function(qidb, data) {
  lapply(qidb, function(qi) {
    join_data2qi(qi, data)
  })
}

#' Make a QI data pair
#'
#' @description Generate a pair of QI specifications and corresponding data and compute some meta information.
#'
#' @param qi a list containing necessary QI specifications
#' @param data a data.frame containing the pre-processed patient survey data the QI should be evaluated on
#'
#' @return The function takes the QI specifications and data and generates a more extensive list object of class \code{qi} to be used for computation of the QI result.
#' For instance, the output contains some meta-information such as number of Merkmale or item categories. The output is ready to use for QI evaluation, i.e.
#' can serve as input for the function \code{evaluateQI}.
#'
#' @export
join_data2qi <- function(qi, data) {
  q <- list(Name = qi$Name, KN_ID = qi$KN_ID, GG = qi$GG, RefArt = qi$RefArt, RefVal = qi$RefVal, RefOp = qi$RefOp)

  ## KN_ID as character
  q$KN_ID <- if (is.numeric(q[["KN_ID"]])) sprintf("%0.2d", q[["KN_ID"]]) else q[["KN_ID"]]

  ## Merkmale

  q$Merkmale <- lapply(qi$Merkmale, function(m) {
    merkmal(m)
  })

  ## Number of Merkmale
  q$M <- length(q$Merkmale)
  ## Vector with the Merkmale
  q$L <- rlang::flatten_int(lapply(q$Merkmale, function(m) m[["L"]]))
  q$Lvec <- rep(1:q$M, times = q[["L"]])

  ## Names of the columns in the data.frame relevant for the QI
  q$colNames <- rlang::flatten_chr(lapply(q$Merkmale, function(m) m[["Items"]]))
  q$colNames <- stats::setNames(q[["colNames"]], q[["colNames"]])

  ## Extract number of categories for each item
  q$item_levels <- sapply(q[["colNames"]], function(colName) {
    vals <- stats::na.omit(unique(attr(data[, colName, drop = TRUE], "labels")))
    stopifnot(all(c(0, 100) %in% vals))
    length(vals)
  })

  q$observed_item_levels <- sapply(q[["colNames"]], function(colName) {
    unique(data[, colName, drop = TRUE])
  })

  ## compute cleaned data
  dat_clean <- data
  for (colName in q[["colNames"]]) {
    dat_clean <- cleanNA(dat_clean, colName = colName)
  }

  ## check if
  ## a) all items belonging to same Merkmal have equal number of categories
  ## b) points are rounded to zero digits
  ## c) points are equidistant and between 0 and 100
  check_points <- lapply(seq_len(q$M), function(i, data) {
    ## Extract response of relevant items from the data
    idx <- which(q$Lvec == i)
    nClass <- q$item_levels[idx]
    ## Check if all Items belonging to the same Merkmal have equal number of categories
    browser(expr = !all(nClass == nClass[1]))
    stopifnot(all(nClass == nClass[1]))
    nClass <- nClass[1]
    nItems <- length(idx)
    y <- unlist(dat_clean[, q$colNames[idx]])
    ## Base for points
    base <- seq(0, 100, length = nClass)
    ## The points themselves
    points <- round(base / 100, 2)
    # Extract non-missing data
    y_notNA <- y[!is.na(y)]
    ## Check if points are rounded to 0 digits, i.e., y is rounded to 2 digits
    stopifnot(isTRUE(all.equal(as.numeric(round(y_notNA, 2) - y_notNA), rep(0, length(y_notNA)))))
    ## Check if points are equidistant and between 0 and 100
    browser(expr = !all(y_notNA %in% points))
    stopifnot(all(y_notNA %in% points))
  }, data = dat_clean)

  ## Extract the question asked for each item
  q$item_fragen <- sapply(q[["colNames"]], function(colName) {
    attr(data[, colName, drop = TRUE], "label")
  })

  # Vector with the base weight of each item
  q$item_weight <- with(q, rep(1 / (L * M), times = L))

  # classificator for "statistisch auffaellig"
  q$Auff_Stat <- function(res) {
    with(
      q,
      if (RefOp == ">=") {
        as.numeric(res$interval["upper"]) <= RefVal
      } else {
        as.numeric(res$interval["lower"]) >= RefVal
      }
    )
  }

  ## Consistency checks
  stopifnot(length(q[["colNames"]]) == sum(q[["L"]]))

  # Set object class
  class(q) <- c(class(q), "qi")

  ## restrict data to relevant rows and columns
  q[["data"]] <- prepareData(q, dat_clean)

  return(q)
}

## set NA values to FALSE
na2FALSE <- function(x) {
  x[is.na(x)] <- FALSE
  x
}

## Generate equidistant point scale for item with K categories
cats <- function(K) 100 * seq(0, 1, length = K)

## Transform negative values to NAs and divide other values by 100
utils::globalVariables(':=')
cleanNA <- function(data, colName) {
  dplyr::mutate(
    data,
    !!colName := dplyr::if_else((!!rlang::sym(colName)) < 0, 
      NA_real_, as.numeric(!!rlang::sym(colName)) / 100)
  )
}

## Function to define a Merkmal
merkmal <- function(m) {
  m$L <- length(m[["Items"]])
  return(m)
}


#' Prepare data for the QI-computation.
#'
#' @description Data preparation includes restricting to the GG, defining the static
#' part of the weights, etc. We introduce this data preparation step
#' in order to be able to do subsequent bootstrap resampling faster.
#'
#' @param qi A qi (i.e. class \code{qi})
#' @param clean_data A \code{data.frame} containing the cleaned patient survey data
#' (all non-missing values are relevant for computation)
#'
#' @return Takes the cleaned data set, applies GG condition, restricts to relevant
#' columns and returns the resulting data set as a \code{data.frame}.
prepareData <- function(qi, clean_data) {

  # Extract from qi
  colNames <- qi[["colNames"]]
  item_weight <- qi[["item_weight"]]
  L <- qi[["L"]]
  M <- qi[["M"]]

  # check that all items have levels
  if (any(qi$item_levels == 0)) {
    stop("all items must have levels")
  }

  ## Restrict data to GG
  subset_GG <- function(data) {
    expr <- rlang::parse_expr(qi[["GG"]])
    withr::with_package("iqtigfunctions", {
      dplyr::filter(data, !!expr)
    })
  }

  relev_data <- subset_GG(clean_data)

  ## extract relevant columns from data - use Indikatorinformation for this
  relev_data <- dplyr::select(relev_data, ID, meta_unit, !!!colNames)

  ## Generate weights, i.e., assign each Item weight determined in qi()
  weightsItemList <- lapply(
    stats::setNames(colNames, paste0("w_", colNames)),
    function(colName) {
      idx <- grep(paste0("^", colName, "$"), colNames)
      item_weight[idx]
    }
  )

  ## Compute the weights
  relev_data <- dplyr::mutate(relev_data, !!!weightsItemList)

  return(relev_data)
}
