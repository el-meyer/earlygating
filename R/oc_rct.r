#' RCT Operating Characteristics
#'
#' Function for calculating the operating characteristics of the RCT
#' Bayesian designs in setting 3 and 4 for early gating.
#'
#' @param N_c Sample Size in the control group. Can be either a single value or a vector,
#'            but needs to be the same length as N_e.
#'
#' @param N_e Sample Size in the experimental group. Can be either a single value or a vector,
#'            but needs to be the same length as N_c.
#'
#' @param delta Required superiority to make a "GO" decision. Corresponds to \eqn{\delta}.
#'
#' @param delta_power Superiority, at which decision power will be evaluated.
#'                    Corresponds to \eqn{\bar{\delta}}{\delta bar}.
#'
#' @param confidence Required confidence to make "GO" decision. Corresponds to \eqn{\gamma}{\gamma}.
#'
#' @param e_a Alpha parameter of Beta Prior Distribution for the experimental response rate.
#'            Corresponds to \eqn{\alpha_e}{\alpha e}. Default is \eqn{\frac{1}{2}}{1/2}.
#'
#' @param e_b Beta parameter of Beta Prior Distribution for the experimental response rate.
#'            Corresponds to \eqn{\beta_e}{\beta e}. Default is \eqn{\frac{1}{2}}{1/2}.
#'
#' @param c_a Alpha parameter of Beta Prior Distribution for the control response rate.
#'            Corresponds to \eqn{\alpha_c}{\alpha c}. Default is \eqn{\frac{1}{2}}{1/2}.
#'
#' @param c_b Beta parameter of Beta Prior Distribution for the control response rate.
#'            Corresponds to \eqn{\beta_c}{\beta c}. Default is \eqn{\frac{1}{2}}{1/2}.
#'
#' @param h_a Alpha parameter of Beta Prior Distribution for the historical control response rate.
#'            Corresponds to \eqn{\alpha_h}{\alpha h}. Only needs to be specified, if RR_h, N_h
#'            and w are also specified. Default is \eqn{\frac{1}{2}}{1/2}.
#'
#' @param h_b Beta parameter of Beta Prior Distribution for the historical control response rate.
#'            Corresponds to \eqn{\beta_h}{\beta h}. Only needs to be specified, if RR_h, N_h and w
#'            are also specified. Default is \eqn{\frac{1}{2}}{1/2}.
#'
#' @param RR_h Historical control response rate. Corresponds to \eqn{p_h}{p h}.
#'             If specified together with N_h and w, function will use setting 4 from pdf.
#'
#' @param N_h Historical control sample size. Corresponds to \eqn{n_h}{n h}.
#'            If specified together with RR_h and w, function will use setting 4 from pdf.
#'
#' @param w Level of dynmaic borrowing. Corresponds to \eqn{w}{w}.
#'
#' @param trues Sequence of true control response rates and experimental response rates, at which the
#'              Probability to Go will be computed. Default is seq(0,1,0.01) to ensure continuous
#'              plots and accurate results.
#'
#' @param plot Plots yes or no. Default is TRUE.
#'
#' @param legend Logical; whether or not to include legend in plot. Default is TRUE.
#'
#' @param legend.pos Position of legend. Default is "topleft".
#'
#' @return A matrix containing the decision power and decision alpha with respect to the true control response rate.
#'
#' @examples
#'
#' # Setting 3
#' oc_rct(
#'   N_c = 25, N_e = 25, delta = 0.08,
#'   delta_power = 0.13, confidence = 0.6
#' )
#'
#' # Setting 4
#' oc_rct(
#'   N_c = 25, N_e = 25, delta = 0.08,
#'   delta_power = 0.13, confidence = 0.6,
#'   RR_h = 0.5, N_h = 50, w = 0.3
#' )
#'
#' @export
oc_rct <- function(N_c, N_e, delta, delta_power, confidence, e_a=0.5, e_b=0.5, c_a=0.5,
                   c_b=0.5, h_a=0.5, h_b=0.5, RR_h=NULL, N_h=NULL, w=NULL,
                   trues=seq(0,1,0.01), plot=T, legend = T, legend.pos="topleft"){

  oc_rct_nodyn <- function(N_c, N_e, delta, delta_power, confidence, e_a, e_b, c_a, c_b,
                           trues = seq(0,1,0.01), plot = T, legend = T, legend.pos="topleft"){

    ocs <- matrix(nrow = length(trues), ncol = 2)
    resp_cl <- req_resp_rct(N_c, N_e, delta, confidence, e_a, e_b, c_a, c_b, plot=plot)[,2]
    resp <- resp_cl
    resp[is.na(resp)] <- N_e+1
    for(j in 1:length(trues)){
      alphas <- stats::pbinom(resp-1, N_e, trues[j], lower.tail = FALSE)
      weighted_alphas <- alphas * stats::dbinom(0:N_c, N_c, trues[j])

      powers <- suppressWarnings(stats::pbinom(resp-1, N_e, trues[j]+delta_power, lower.tail = FALSE))
      powers[is.nan(powers)] <- 0
      weighted_powers <- powers * stats::dbinom(0:N_c, N_c, trues[j])

      alpha <- sum(weighted_alphas)
      power <- sum(weighted_powers)

      ocs[j,1] <- alpha
      ocs[j,2] <- power
    }

    if(plot){

      plot(trues, ocs[,1], xlab = "True control response rate", ylab = "Probability", type = "l", ylim = c(0,max(ocs)), lty = 2,
           main = "Decision power and Decision type 1 error \n with respect to true control response rate")
      graphics::lines(trues, ocs[,2], type = "l")
      if(legend){
      legend(legend.pos, legend = c("Decision type 1 error", "Decision power"), bty = "n", lty = c(2,1))
      }
      }

    return(ocs)
  }


  oc_rct_dyn <- function(N_c, N_e, delta, delta_power, confidence, e_a, e_b, c_a, c_b, RR_h,
                         N_h, h_a, h_b, w, trues = seq(0,1,0.01), plot = T, legend = T, legend.pos="topleft"){

    ocs <- matrix(nrow = length(trues), ncol = 2)
    resp_cl <- req_resp_rct(N_c, N_e, delta, confidence, e_a, e_b, c_a, c_b, h_a, h_b, RR_h, N_h, w, plot)[,2]
    resp <- resp_cl
    resp[is.na(resp)] <- N_e+1
    for(j in 1:length(trues)){
      alphas <- stats::pbinom(resp-1, N_e, trues[j], lower.tail = FALSE)
      weighted_alphas <- alphas * stats::dbinom(0:N_c, N_c, trues[j])

      powers <- suppressWarnings(stats::pbinom(resp-1, N_e, trues[j]+delta_power, lower.tail = FALSE))
      powers[is.nan(powers)] <- 0
      weighted_powers <- powers * stats::dbinom(0:N_c, N_c, trues[j])

      alpha <- sum(weighted_alphas)
      power <- sum(weighted_powers)

      ocs[j,1] <- alpha
      ocs[j,2] <- power
    }

    if(plot){
      plot(trues, ocs[,1], xlab = "True control response rate", ylab = "Probability", type = "l", ylim = c(0,max(ocs)), lty = 2,
           main = "Decision power and Decision type 1 error \n with respect to true control response rate")
      graphics::lines(trues, ocs[,2], type = "l")
      if(legend){
      legend(legend.pos, legend = c("Decision type 1 error", "Decision power"), bty = "n", lty = c(2,1))
      }
    }

    return(ocs)
  }

  if(!is.null(RR_h)&!is.null(N_h)&!is.null(w)){
    ret <- oc_rct_dyn(N_c, N_e, delta, delta_power, confidence, e_a, e_b, c_a, c_b, RR_h,
                      N_h, h_a, h_b, w, trues, plot, legend, legend.pos)
  }else{
    ret <- oc_rct_nodyn(N_c, N_e, delta, delta_power, confidence, e_a, e_b, c_a, c_b, trues, plot, legend, legend.pos)
  }

  res <- cbind(ret, trues)
  colnames(res) <- c("Dec. Alpha", "Dec. Power", "True control RR")
  round(res,4)

}
