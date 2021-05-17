#' Required Responders for GO decision Single Arm
#'
#' Function for calculating the minimum required number of responders in the
#' experimental group to make a GO decision in Settings 1 and 2.
#'
#' @param N_e Sample Size in the experimental group.
#'
#' @param delta Required superiority to make a "GO" decision. Corresponds to \eqn{\delta}.
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
#' @param adapt Level of adapting of experimental control rate to account for patient selection bias
#'              from phase II to phase III. Corresponds to \eqn{\xi}{\xi}. Default is 1, so no adapting.
#'
#' @return Integer.
#'
#' @examples
#'
#' # Setting 1
#' req_resp(
#'   N_e = 50, delta = 0.08,
#'   confidence = 0.6, hist_RR_c = 0.5
#' )
#'
#' # Setting 2
#' req_resp(
#'   N_e = 50, delta = 0.08,
#'   confidence = 0.6, RR_h = 0.5, N_h = 50
#' )
#'
#' @export
req_resp <- function(N_e, delta, confidence, e_a=0.5, e_b=0.5, h_a=0.5, h_b=0.5,
                     RR_h=NULL, N_h=NULL, hist_RR_c=NULL, adapt = 1){

  req_resp_wr_true <- function(N_e, hist_RR_c, delta, confidence, e_a=1, e_b=1, adapt=1){
    crit <- 0
    e_a_temp <- e_a # starting values for Beta, i.e. 0 responders
    e_b_temp <- e_b + N_e
    a <- try(while(stats::qbeta(1-confidence, e_a_temp, e_b_temp) <= hist_RR_c+delta){ # iteratively solving the inequation
      crit <- crit + 1
      e_a_temp <- e_a + adapt*crit
      e_b_temp <- e_b - adapt*crit + N_e
      if(crit > N_e){ # if impossible to achieve number of responders with given sample size, then "NA"
        crit <- NA
        break
      }
    }, silent = T) # calculate critical value of "successes"
    if(class(a) == "try-error"){crit <- NA}
    return(crit)
  }

  req_resp_wr_hists <- function(N_e, delta, confidence, e_a=1, e_b=1, h_a=1,
                                h_b=1, RR_h, N_h, adapt=1){

    h_a_new <- h_a + RR_h*N_h
    h_b_new <- h_b - RR_h*N_h + N_h # it is assumed that "RR_h*N_h" responders were observed in historical control group
    crit <- 0
    e_a_temp <- e_a # starting values for Beta, i.e. 0 responders
    e_b_temp <- e_b + N_e
    a <- try(while(stats::integrate(function(y) { stats::dbeta(y, e_a_temp, e_b_temp)*
        sapply(y, function(y) {
          stats::pbeta(y-delta, h_a_new, h_b_new)
        })
    }, delta, 1)$value <= confidence){ # iteratively solving the inequation
      crit <- crit + 1
      e_a_temp <- e_a + adapt*crit
      e_b_temp <- e_b - adapt*crit + N_e
      if(crit>N_e){# if impossible to achieve number of responders with given sample size, then "NA"
        crit <- NA
        break
      }
    }, silent = T)

    if(class(a) == "try-error"){crit <- NA}
    return(crit)
  }

  if(!is.null(RR_h) & !is.null(N_h)){
    ret <- req_resp_wr_hists(N_e, delta, confidence, e_a, e_b, h_a,
                      h_b, RR_h, N_h, adapt)
  }else{
    ret <- req_resp_wr_true(N_e, hist_RR_c, delta, confidence, e_a, e_b, adapt)
  }
  ret
}
