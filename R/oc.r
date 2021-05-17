#' Single Arm Operating Characteristics
#'
#' Function for calculating the operating characteristics of the single arm
#' Bayesian designs in setting 1 and 2 for early gating.
#'
#' @param N_e Sample Size in the experimental group. Can be either a single value or a vector.
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
#' @param h_a Alpha parameter of Beta Prior Distribution for the historical control response rate.
#'            Corresponds to \eqn{\alpha_h}{\alpha h}. Only needs to be specified, if RR_h and N_h
#'            are also specified. Default is \eqn{\frac{1}{2}}{1/2}.
#'
#' @param h_b Beta parameter of Beta Prior Distribution for the historical control response rate.
#'            Corresponds to \eqn{\beta_h}{\beta h}. Only needs to be specified, if RR_h and N_h
#'            are also specified. Default is \eqn{\frac{1}{2}}{1/2}.
#'
#' @param RR_h Historical control response rate. Corresponds to \eqn{p_h}{p h}.
#'             If specified together with N_h, function will use setting 2 from pdf.
#'
#' @param N_h Historical control sample size. Corresponds to \eqn{n_h}{n h}.
#'            If specified together with RR_h, function will use setting 2 from pdf.
#'
#' @param hist_RR_c Point estimate of historical control repsonse rate. Corresponds to \eqn{\hat{p_h}}{p h hat}.
#'                  If specified, while RR_h and N_h are not specified, function will use setting 1 from pdf.
#'
#' @param trues Sequence of true control response rates and experimental response rates, at which the
#'              Probability to Go will be computed. Default is seq(0,1,0.01) to ensure continuous
#'              plots and accurate results.
#'
#' @param adapt Level of adapting of experimental control rate to account for patient selection bias
#'              from phase II to phase III. Corresponds to \eqn{\xi}{\xi}. Default is 1, so no adapting.
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
#' # Setting 1
#' oc(
#'   N_e = 50, delta = 0.08, delta_power = 0.13,
#'   confidence = 0.6, hist_RR_c = 0.5
#' )
#'
#' # Setting 2
#' oc(
#'   N_e = 50, delta = 0.08, delta_power = 0.13,
#'   confidence = 0.6, RR_h = 0.5, N_h = 50
#' )
#'
#' @export
oc <- function(N_e, delta, delta_power, confidence, e_a=0.5, e_b=0.5, h_a=0.5, h_b=0.5,
               RR_h=NULL, N_h=NULL, hist_RR_c=NULL, trues = seq(0,1,0.01),
               adapt = 1, plot = T, legend = T, legend.pos="topleft"){

oc_dist <- function(N_e, true_RR_c, delta, delta_power, confidence, e_a=1, e_b=1, h_a=1, h_b=1,
                    RR_h, N_h, trues = seq(0,1,0.01), adapt = 1, plot = T, legend = T, legend.pos="topleft"){

  true_RR_e <- trues
  # firstly calculate how many "successes" are needed for "GO" Decision
  crit <- req_resp(N_e=N_e, delta=delta, confidence=confidence, e_a=e_a, e_b=e_b,
                   h_a=h_a, h_b=h_b, RR_h=RR_h, N_h=N_h, adapt=adapt)

  if(plot){
    P_go <- stats::pbinom(crit-1, N_e, true_RR_e, lower.tail = FALSE)
  }# depending on how many successes are needed for GO, what is probability
  #   with respect to different experimental response rates
  if(true_RR_c+delta_power<=1){# power at true_RR_c + delta power, unless true_RR_c + delta power > 1, in which case power = 1
    power_ach <- stats::pbinom(crit-1, N_e, true_RR_c+delta_power, lower.tail = FALSE)
  }else{power_ach <- 0}
  alpha_ach <- stats::pbinom(crit-1, N_e, true_RR_c, lower.tail = FALSE) # alpha at true_RR_c

  if(plot) {
    plot(true_RR_e, P_go, type = "l", xlab = "True experimental response rate",
         ylab = "Prob(GO)", main = "Prob(GO) w/r to True experimental response rate")
    graphics::abline(v = true_RR_e[which(round(true_RR_e,5) == round(true_RR_c+delta_power, 5))],
           h = P_go[which(round(true_RR_e,5) == round(true_RR_c+delta_power, 5))])
    graphics::abline(v = true_RR_e[which(round(true_RR_e,5) == round(true_RR_c, 5))],
           h = P_go[which(round(true_RR_e,5) == round(true_RR_c, 5))])
  }

 return(c(alpha_ach, power_ach, crit))
}


oc_true <- function(N_e, true_RR_c, hist_RR_c, delta, delta_power, confidence, e_a=1, e_b=1,
               trues = seq(0,1,0.01), adapt=1, plot = T, legend = T, legend.pos="topleft"){

  true_RR_e <- trues
  # firstly calculate how many "successes" are needed for "GO" Decision
  crit <- req_resp(N_e=N_e, hist_RR_c=hist_RR_c, delta=delta, confidence=confidence,
                   e_a = e_a, e_b=e_b, adapt=adapt)

  if(plot){
    P_go <- stats::pbinom(crit-1, N_e, true_RR_e, lower.tail = FALSE)
  }# depending on how many successes are needed for GO, what is probability
  #   with respect to different experimental response rates
  if(true_RR_c+delta_power<=1){ # power at true_RR_c + delta power, unless true_RR_c + delta power > 1, in which case power = 1
    power_ach <- stats::pbinom(crit-1, N_e, true_RR_c+delta_power, lower.tail = FALSE)
  }else{power_ach <- 0}
  alpha_ach <- stats::pbinom(crit-1, N_e, true_RR_c, lower.tail = FALSE) # alpha at true_RR_c
  if(plot) {
    plot(true_RR_e, P_go, type = "l", xlab = "True experimental response rate",
         ylab = "Prob(GO)", main = "Prob(GO) w/r to True experimental response rate")
    graphics::abline(v = true_RR_e[which(round(true_RR_e,5) == round(true_RR_c+delta_power, 5))],
           h = P_go[which(round(true_RR_e,5) == round(true_RR_c+delta_power, 5))])
    graphics::abline(v = true_RR_e[which(round(true_RR_e,5) == round(true_RR_c, 5))],
           h = P_go[which(round(true_RR_e,5) == round(true_RR_c, 5))])

  }
  return(c(alpha_ach, power_ach, crit))
  #return achieved power, alpha and critical value
}

if(!is.null(RR_h) & !is.null(N_h)){
  ret <- matrix(nrow = length(trues), ncol = 2)
  for(i in 1:length(trues)){
  ret[i,] <- oc_dist(N_e, trues[i], delta, delta_power, confidence, e_a, e_b, h_a, h_b,
                 RR_h, N_h, trues, adapt, plot=F)[1:2]
  }
}else{
  ret <- matrix(nrow = length(trues), ncol = 2)
  for(i in 1:length(trues)){
  ret[i,] <- oc_true(N_e, trues[i], hist_RR_c, delta, delta_power, confidence, e_a, e_b,
                 trues, adapt, plot=F)[1:2]
  }
}
res <- cbind(ret,trues)
colnames(res) <- c("Dec. Alpha", "Dec. Power", "True Control RR")
if(plot){
plot(trues,res[,1], type = "l", lty = 2, ylab = "Probability", xlab = "True control response rate", main = "Decision power and Decision type 1 error \n w/r to True Control response rate")
graphics::lines(trues,res[,2], type = "l", lty = 1)
 if(legend){
 legend(legend.pos, legend = c("Decision type 1 error", "Decision power"), bty = "n", lty = c(2,1))
 }
}
round(res,4)

}
