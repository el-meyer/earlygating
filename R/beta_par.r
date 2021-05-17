#' Parameters of Beta distribution given Historical Data
#'
#' Function for calculating the parameters of the beta distribution used to
#' average the operating characteristics, given historical data.
#'
#' @param mu_cov A character vector containing the names of covariates in data
#'               that should be used to model the parameter \eqn{\mu} in the pdf.
#'
#' @param phi_cov A character vector containing the names of covariates in data that
#'                should be used to model the parameter \eqn{\phi} in the pdf.
#'                Default is NULL, so \eqn{\phi} will not be modelled with respect to the covariates.
#'
#' @param orr Character containing the name of the variable in data that represents the objective response rate.
#'
#' @param data Data frame containing all the covariates and the ORR.
#'
#' @param newdata Data frame containing a single value for each of the specified covariates that will
#'                be used to estimate the parameters of the Beta distribution.
#'
#' @param link Link function for \eqn{\mu}. Corresponds to \eqn{g}{g}. Default is NULL, which means the
#'             link function will be automatically chosen as the one yielding the highest log-likelihood
#'             for the given data and covariates.
#'
#' @param weights Weights that should be used for regression. Default is NULL, so no weights.
#'
#' @param plot Plots yes or no. Default is TRUE.
#'
#' @examples
#'
#' mu_cov <- c("date", "Phase")
#' orr <- "ORR"
#' newdata <- data.frame(
#'   "date" = 2017,
#'   "Phase" = factor(3)
#'   )
#' studs <- data.frame(
#'   "ORR"= c(0.693, 0.580, 0.693, 0.477, 0.609,
#'            0.727, 0.727, 0.591, 0.362, 0.593,
#'            0.792, 0.620, 0.550, 0.690, 0.776),
#'   "date" = c(2011, 2008.5, 2009, 1996, 2001,
#'              2003.5, 2002.5, 2008, 2000,
#'              2006, 2005, 2007.5, 2009.5,
#'              2010.5, 2010),
#'    "Phase" = factor(c(3, 2, 3, 3, 2, 2, 3, 3,
#'                       3, 3, 2, 3, 3, 3, 2)),
#'    "N" = c(293, 69, 336, 235, 92, 110, 131,
#'            208, 94, 123, 53, 182, 267, 239, 237)
#'  )
#'
#' beta_par(
#'   mu_cov = mu_cov,
#'   orr = orr,
#'   data = studs,
#'   newdata = newdata,
#'   weights = studs$N/mean(studs$N)
#' )
#'
#' @export
beta_par <- function(mu_cov, phi_cov=NULL, orr, data, newdata, link = NULL, weights = NULL, plot = T){ # names in cov and orr

#require(betareg)

  fmla <- paste(orr, " ~ ", paste(mu_cov, collapse= "+"))

  if(!is.null(phi_cov)){
    fmla <- paste(fmla, "|", paste(phi_cov, collapse= "+"))
  }
  fmla <- stats::as.formula(fmla)

  if(is.null(link)){
    if(is.null(weights)){
      modx <- betareg::betareg(fmla, data = data)
    }else{
      modx <- betareg::betareg(fmla, data = data, weights = weights)

    }
  a <- sapply(c("logit", "probit", "cloglog", "cauchit", "loglog"),
         function(x) stats::logLik(stats::update(modx, link = x)))
  link <- names(which.max(a))
  }

  if(is.null(weights)){
  mod1 <- betareg::betareg(fmla, data = data, link = link)
  }else{
    mod1 <- betareg::betareg(fmla, data = data, link = link, weights = weights)

  }

  alpha <- mod1$coefficients$precision *
                stats::predict(mod1, newdata = newdata)
  beta <- mod1$coefficients$precision - alpha

  if(plot){
    x <- NULL
     graphics::curve(stats::dbeta(x,alpha, beta), xlab = "Response Rate", ylab = "Density",
           main = paste("Probability density function of the Beta (", round(alpha,2),";", round(beta,2),") distribution", sep = ""))
    }

  ret <- c(alpha, beta)
  names(ret) <- c("Alpha", "Beta")
  ret

}

