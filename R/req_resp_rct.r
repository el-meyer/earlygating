#' Required Responders for GO decision RCT
#'
#' Function for calculating the minimum required number of responders in the
#' experimental group to make a GO decision in Settings 3 and 4.
#'
#' @param N_c Sample Size in the control group.
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
#' @param plot Plots yes or no. Default is TRUE.
#'
#' @return Matrix containing pairs of successes in control group and respective required successes in experimental group.
#'
#' @examples
#'
#' # Setting 3
#' req_resp_rct(
#'   N_c = 25, N_e = 25,
#'   delta = 0.08, confidence = 0.6
#' )
#'
#' # Setting 4
#' req_resp_rct(
#'   N_c = 25, N_e = 25,
#'   delta = 0.08, confidence = 0.6,
#'   RR_h = 0.5, N_h = 50, w = 0.3
#' )
#'
#' @export
req_resp_rct <- function(N_c, N_e, delta, confidence, e_a=0.5, e_b=0.5, c_a=0.5, c_b=0.5,
                         h_a=0.5, h_b=0.5, RR_h=NULL, N_h=NULL, w=NULL, plot=T){

  req_resp_rct_nodyn <- function(N_c, N_e, delta, confidence, e_a, e_b, c_a, c_b){
    resp <- rep(NA, N_c+1)
    for(i in 0:N_c){

      suc_c <- i
      c_a_new <- c_a + suc_c
      c_b_new <- c_b + N_c - suc_c

      crit <- 0
      e_a_temp <- e_a
      e_b_temp <- e_b + N_e

      a <- try(
        while(stats::integrate(function(y) {stats::dbeta(y, e_a_temp, e_b_temp)*
            sapply(y, function(y) {
              stats::pbeta(y-delta, c_a_new, c_b_new)
            })
        }, delta, 1)$value <= confidence){
          e_a_temp <- e_a + crit
          e_b_temp <- e_b - crit + N_e
          crit <- crit + 1
          if(crit>N_e){
            crit <- N_e + 1
            break
          }
        }, silent = T)


      if(class(a) == "try-error"){resp[i+1] <- NA}else{
        resp[i+1] <- crit
      }
    }
    resp_cl <- resp
    resp_cl[resp_cl == N_e+1] <- NA

    return(resp_cl)
  }

  req_resp_rct_dyn <- function(N_c, N_e, delta, confidence, e_a, e_b, c_a, c_b, RR_h, N_h, h_a, h_b, w){
    resp <- rep(NA, N_c+1)
    for(i in 0:N_c){

      suc_c <- i

      w1 <- w*beta(suc_c + N_h*RR_h + h_a, N_h + N_c - suc_c - N_h*RR_h + h_b)/beta(N_h*RR_h + h_a, N_h - N_h*RR_h + h_b)
      w2 <- (1-w)*beta(suc_c + c_a, N_c - suc_c + c_b)/beta(c_a, c_b)
      w1n <- w1/(w1+w2)
      w2n <- w2/(w1+w2)

      m_a <- w1n*(suc_c + N_h*RR_h + h_a) + w2n*(suc_c + c_a)
      m_b <- w1n*(N_h + N_c - suc_c - N_h*RR_h + h_b) + w2n*(N_c - suc_c + c_b)

      crit <- 0
      e_a_temp <- e_a
      e_b_temp <- e_b + N_e

      a <- try(
        while(stats::integrate(function(y) {stats::dbeta(y, e_a_temp, e_b_temp)*
            sapply(y, function(y) {
              stats::pbeta(y-delta, m_a, m_b)
            })
        }, delta, 1)$value <= confidence){
          e_a_temp <- e_a + crit
          e_b_temp <- e_b - crit + N_e
          crit <- crit + 1
          if(crit>N_e){
            crit <- N_e + 1
            break
          }
        }, silent = T)

      if(class(a) == "try-error"){resp[i+1] <- NA}else{
        resp[i+1] <- crit
      }
    }
    resp_cl <- resp
    resp_cl[resp_cl == N_e+1] <- NA

    return(resp_cl)
  }

  if(!is.null(RR_h)&!is.null(N_h)&!is.null(w)){
    ret <- req_resp_rct_dyn(N_c, N_e, delta, confidence, e_a, e_b, c_a, c_b, RR_h, N_h, h_a, h_b, w)
  }else{
    ret <- req_resp_rct_nodyn(N_c, N_e, delta, confidence, e_a, e_b, c_a, c_b)
  }
if(plot){
  plot(0:N_c,ret, xlab = "Successes in Control Group", ylab = "Successes in Experimental Group",
       main = "Minimum required number of successes in Exp Group to make \n  GO decision w/r to Successes in Con Group", type = "s")
  }
  res <- cbind(0:N_c,ret)
  colnames(res) <- c("Suc. in Con.", "Req. Suc. in Exp.")
  res
}
