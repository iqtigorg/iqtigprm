#' Bayesian inference for theta for one Merkmal
#'
#' @description This is a helper function for the Bayesian QI computation function.
#' The function computes the parameters of the Beta posterior for the underlying
#' merkmals parameter theta as well as a credibility interval for theta.
#'
#' @param y Observed input of all Items of the Merkmal with points
#' @param nClass Number of Item categories (supposed do be the same for each Item)
#' @param a Prior shape parameter of Beta distribution
#' @param b Prior scale parameter of Beta distribution
#' @param conf_level Confidence level of the resulting uncertainty interval.
#'
#' @return A \code{list} containing the prior and posterior parameters of the Beta distribution,
#' the corresponding posterior mean value for theta and a credibility interval.
#'
#' @export
theta_bayes <- function(y, nClass, a, b, conf_level) {
  ## Translate conf_level to alpha
  alpha <- 1 - conf_level

  ## Base for points
  base <- seq(0, 1, length = nClass)
  ## The points themselves
  points <- round(100 * base) / 100

  ## Check that original responses are whole points, i.e., y is rounded to two digits
  y_notNA <- y[!is.na(y)]
  stopifnot(isTRUE(all.equal(
    as.numeric(round(y_notNA, 2) - y_notNA),
    rep(0, length(y_notNA))
  )))

  ## Check that points are equidistant by checking that no unvalid values occur
  stopifnot(all(y_notNA %in% points))

  ## Tabulate the responses (in points) into frequencies (note: NAs are discarded)
  ytab <- table(factor(y, levels = points))

  ## Update parameters to get the posterior
  astar <- a + sum(ytab * (0:(nClass - 1)))
  bstar <- b + sum(ytab * (nClass - 1):0)

  ## Point estimate (posterior mean)
  hat.post.mean <- astar / (astar + bstar)

  ## Uncertainty interval
  interval <- stats::qbeta(c(lower = alpha / 2, upper = 1 - alpha / 2), astar, bstar)
  ## Add zero and one to the interval if observation is at these boundaries.
  ## Inspired by binom::binom.bayes
  if (ytab[1] == sum(ytab)) interval[1] <- 0
  if (ytab[nClass] == sum(ytab)) interval[2] <- 1

  ## Done
  res <- list(a = a, b = b, astar = astar, bstar = bstar, hat.post.mean = hat.post.mean, interval = interval)
  return(res)
}


#' Computation of one QI with attached data and LE
#'
#' @description Apply Bayesian computation method to one QI with an attached data set for one LE.
#' Each individual merkmal of the QI is computed separately. Overall QI results are computed by MC-sampling
#' from the Merkmal-specific posteriors. The seed for MC-sampling is fixed such that repeated
#' evaluation of the QI yields the same results with respect to the credibility interval.
#'
#' @param qi The QI to compute (object of class \code{qi}).
#' @param meta_unit The meta_unit (i.e. LE) to compute results for. Default (\code{NULL}) means that all data are used for computation of overall results.
#' @param conf_level Confidence level of the resulting uncertainty interval
#' @param a first prior parameter of Beta distribution for each Merkmal
#' @param b second prior parameter of Beta distribution for each Merkmal
#' @param nMC Number of Monte Carlo Samples (needed for QIs with >1 Merkmale)
#'
#' @returns A list containing the point estimate and the uncertainty interval.
#' @field QI_hat Posterior mean estimate for QI value.
#' @field J Number of patients in the data set for the LE.
#' @field interval A vector of length two containing the uncertainty interval.
#'
#' @export
evaluateQI <- function(qi, meta_unit = NULL, conf_level = 0.95, a = 0.5, b = 0.5, nMC = 1e5) {
  alpha <- 1 - conf_level

  ## restrict data to provided meta_unit
  if (is.null(meta_unit)) {
    warning("No meta_unit provided. All data will be used for computation.")
    data <- qi$data
  } else {
    if (!meta_unit %in% qi$data$meta_unit) {
      warning("meta_unit not contained in data.")
    }
    data <- qi$data[qi$data$meta_unit == meta_unit, ]
  }

  J <- nrow(data)

  ## If there is nothing to compute, because data are empty, then stop here.
  if (J == 0) {
    res <- list(QI_hat = NA, J = J, Jbar = 0, interval = c(lower = NA, upper = NA))
    return(res)
  }

  # Merkmal results:
  ## Compute updated parameters of the beta-posterior for each theta (each Merkmal)
  merkmal_results <- lapply(seq_len(qi$M), function(i) {
    ## Extract response of relevant items from the data
    idx <- which(qi$Lvec == i)
    ## Prepare variables for use in the estimation function.
    nClass <- qi$item_levels[idx]
    stopifnot(all(nClass == nClass[1])) # All items of one Merkmal are supposed to have same number of categories
    nClass <- nClass[1]
    nItems <- length(idx)

    ## Extract data
    y <- unlist(data[, qi$colNames[idx]])
    ## Bayes inference
    theta_bayes(
      y = y, nClass = nClass, a = a, b = b,
      conf_level = conf_level
    )
  })

  merkmal_postparams <- dplyr::bind_rows(
    lapply(merkmal_results, function(m) {
      data.frame(astar = m$astar, bstar = m$bstar, postmean = m$hat.post.mean)
    }) 
  )

  # QI results:
  ## QI value
  Qi_hat <- mean(merkmal_postparams$postmean)

  ## QI interval
  ## if only one Merkmal
  if (qi[["M"]] == 1) {
    interval <- merkmal_results[[1]][["interval"]]
  } else {
    ## interval by mean of independent betas (for each Merkmal) by Monte Carlo simulation
    set.seed(123)
    theSum <- rep(0, nMC)
    for (i in seq_len(nrow(merkmal_postparams))) { # nrow(merkmal_postparams) = qi[["M"]]
      theSum <- theSum + stats::rbeta(nMC, merkmal_postparams[i, "astar"], merkmal_postparams[i, "bstar"])
    }
    theMean <- theSum / nrow(merkmal_postparams)
    ## Very simple quantile based (1-alpha)*100% tolerance interval
    interval <- stats::quantile(theMean, prob = c(alpha / 2, 1 - alpha / 2), type = 6)
    if (isTRUE(all.equal(merkmal_postparams[["astar"]], a))) {
      interval[1] <- 0
    }
    if (isTRUE(all.equal(merkmal_postparams[["bstar"]], b))) {
      interval[2] <- 1
    }
  }
  names(interval) <- c("lower", "upper")

  ## Consistency checks
  if (!is.na(Qi_hat)) {
    if ((Qi_hat < interval["lower"]) | (Qi_hat > interval["upper"])) {
      warning("Point estimate not in interval due to data priors.")
    }
  } else {
    interval[1] <- interval[2] <- NA
  }

  ## Make the result object
  res <- list(
    QI_hat = as.numeric(Qi_hat),
    J = J, interval = interval
  )

  res$Auff_stat <- qi$Auff_Stat(res)

  ## Done
  return(res)
}

#' Compute all QIs in a QIDB data pair.
#'
#' @description Applies the Bayesian computation method to all QIs in a QIDB-data-pair, separately to all LEs as
#' identified through the \code{meta_unit} column in the data.
#'
#' @param qidb_data_pair A list of QIs with linked data (of class \code{qi})
#' @param conf_level Confidence level of the resulting uncertainty interval
#' @param a first prior parameter of Beta distribution for each QI Merkmal
#' @param b second prior parameter of Beta distribution for each QI Merkmal
#' @param nMC Number of Monte Carlo Samples (needed for QIs with >1 Merkmale)
#'
#' @return A \code{data.frame} containing all computation results, one row per result
#' @field KN_ID the idenifier of the QI
#' @field meta_unit the identifier of the Leistungserbringer (LE)
#' @field QI_hat Posterior mean estimate for QI value
#' @field J overall number of patients in the data set.
#' @field lower lower border of the uncertainty interval.
#' @field upper upper border of the uncertainty interval.
#' @field Auff_stat indicates if result is "statistisch auffaellig"
#' @field RefVal reference value for computation of "statistisch auffaellig"
#'
#' @export
compute_results <- function(qidb_data_pair, conf_level = 0.95, a = 0.5, b = 0.5, nMC = 1e5) {
  ## run the computation for each Leistungserbringer (LE)
  res_all <- lapply(qidb_data_pair, function(qi) {
    ## Compute qi and uncertainty interval for each LE
    res <- lapply(
      unique(qi$data[["meta_unit"]]),
      function(theLE) {
        message(paste0("QI ", qi[["KN_ID"]], " for LE ", theLE))
        le_result <- tryCatch(evaluateQI(qi,
          meta_unit = theLE, conf_level = conf_level,
          nMC = nMC, a = a, b = b
        ))
        le_result$meta_unit <- theLE
        le_result$KN_ID <- qi[["KN_ID"]]
        le_result$RefOp <- qi[["RefOp"]]
        le_result$RefVal <- qi[["RefVal"]]
        purrr::flatten_df(le_result)
      }
    )

    ## Flatten the lists
    if (length(res) > 0) {
      res_flat <- dplyr::select(dplyr::bind_rows(res), KN_ID, meta_unit, dplyr::everything())
    } else {
      res_flat <- NULL
    }
    res_flat
  })

  res_all <- dplyr::bind_rows(res_all)

  ## Done
  return(res_all)
}
