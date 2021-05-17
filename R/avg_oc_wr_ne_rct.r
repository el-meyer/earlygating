#' RCT Average Operating Characteristics
#'
#' Function for calculating the average operating characteristics of
#' two RCT bayesian designs for early gating with respect to the sample size
#' in the experimental group, the sample size in the control group and possible historical data.
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
#' @param alpha_c Alpha parameter of Beta Distribution for the control response rate used to
#'                calculate average operating characteristics. Corresponds to \eqn{\alpha_c}{\alpha c}.
#'
#' @param beta_c Beta parameter of Beta Distribution for the control response rate used to calculate
#'               average operating characteristics. Corresponds to \eqn{\beta_c}{\beta c}.
#'
#' @param trues Sequence of true control response rates and experimental response rates, at which the
#'              Probability to Go will be computed. Default is seq(0,1,0.01) to ensure continuous
#'              plots and accurate results.
#'
#' @param plot Plots yes or no. Default is TRUE.
#'
#' @param coresnum Number of cores used for parallel computing, in case N_e is a vector.
#'                 Default is the number of total cores - 1.
#'
#' @param legend Logical; whether or not to include legend in plot. Default is TRUE.
#'
#' @param legend.pos Position of legend. Default is "topleft".
#'
#' @return Either a vector containing the average decision power and average alpha (if N_e has length 1) or
#'         a matrix containing the average decision power and average decision alpha (if N_e has length > 1),
#'         where every row corresponds to one value of N_e.
#'
#' @examples
#'
#' # Setting 3
#' avg_oc_wr_ne_rct(
#' N_c = 25, N_e = 25, delta = 0.08,
#' delta_power = 0.13, confidence = 0.6,
#' alpha_c = 15, beta_c = 13
#' )
#'
#'
#'
#' # Setting 4
#' avg_oc_wr_ne_rct(
#' N_c = 25, N_e = 25, delta = 0.08,
#' delta_power = 0.13, confidence = 0.6,
#' alpha_c = 15, beta_c = 13,
#' RR_h = 0.5, N_h = 100, w = 0.3
#' )
#'
#' @export
avg_oc_wr_ne_rct <- function(N_c, N_e, delta, delta_power, confidence, e_a = 0.5, e_b = 0.5,
                             c_a = 0.5, c_b = 0.5, h_a = 0.5, h_b = 0.5, N_h = NULL, RR_h = NULL,
                             w = NULL, alpha_c, beta_c, trues = seq(0,1,0.01), plot = T,
                             coresnum = NULL, legend = T, legend.pos = "topleft") {

  req_resp_rct <- function(N_c, N_e, delta, confidence, e_a, e_b, c_a, c_b){
    resp <- rep(NA, N_c+1)
    for(i in 0:N_c){

      suc_c <- i
      c_a_new <- c_a + suc_c
      c_b_new <- c_b + N_c - suc_c

      crit <- 0
      e_a_temp <- e_a
      e_b_temp <- e_b + N_e

      while(stats::integrate(function(y) { stats::dbeta(y, e_a_temp, e_b_temp)*
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
      }

      resp[i+1] <- crit
    }
    resp_cl <- resp
    resp_cl[resp_cl == N_e+1] <- NA

    return(resp_cl)
  }

  oc_rct <- function(N_c, N_e, delta, delta_power, confidence, e_a, e_b, c_a, c_b, trues = seq(0,1,0.01), plot = T, legend=T, legend.pos){

    ocs <- matrix(nrow = length(trues), ncol = 2)
    resp_cl <- req_resp_rct(N_c, N_e, delta, confidence, e_a, e_b, c_a, c_b)
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

      plot(0:N_c,resp_cl, xlab = "Successes in Control Group", ylab = "Successes in Experimental Group",
           main = "Minimum required number of successes in Exp Group to make \n  GO decision w/r to Successes in Con Group", type = "s")

      plot(trues, ocs[,1], xlab = "True control response rate", ylab = "Probability", type = "l", ylim = c(0,max(ocs)), lty = 2,
           main = "Decision power and Decision type 1 error \n with respect to true control response rate")
      graphics::lines(trues, ocs[,2], type = "l")
      if(legend){
        legend(legend.pos, legend = c("Decision type 1 error", "Decision power"), bty = "n", lty = c(2,1))
      }
    }

    return(ocs)
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

      while(stats::integrate(function(y) { stats::dbeta(y, e_a_temp, e_b_temp)*
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
      }

      resp[i+1] <- crit
    }
    resp_cl <- resp
    resp_cl[resp_cl == N_e+1] <- NA

    return(resp_cl)
  }


  oc_rct_dyn <- function(N_c, N_e, delta, delta_power, confidence, e_a, e_b, c_a, c_b, RR_h, N_h, h_a, h_b, w, trues = seq(0,1,0.01), plot = T, legend = T, legend.pos){

    ocs <- matrix(nrow = length(trues), ncol = 2)
    resp_cl <- req_resp_rct_dyn(N_c, N_e, delta, confidence, e_a, e_b, c_a, c_b, RR_h, N_h, h_a, h_b, w)
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
      plot(0:N_c,resp_cl, xlab = "Successes in Control Group", ylab = "Successes in Experimental Group",
           main = "Minimum required number of successes in Exp Group to make \n  GO decision w/r to Successes in Con Group", type = "s")

      plot(trues, ocs[,1], xlab = "True control response rate", ylab = "Probability", type = "l", ylim = c(0,max(ocs)), lty = 2,
           main = "Decision power and Decision type 1 error \n with respect to true control response rate")
      graphics::lines(trues, ocs[,2], type = "l")
      if(legend){
        legend(legend.pos, legend = c("Decision type 1 error", "Decision power"), bty = "n", lty = c(2,1))
      }
    }

    return(ocs)
  }


  if(!is.null(RR_h) & !is.null(N_h) & !is.null(h_a) & !is.null(h_b) & !is.null(w)){


    if(length(N_e)>1){

      av_power <- numeric(length(N_e))
      av_alpha <- numeric(length(N_e))

      if(is.null(coresnum)){coresnum <- parallel::detectCores()-1}
      cores <- coresnum
      cl <- parallel::makePSOCKcluster(cores)
      doParallel::registerDoParallel(cl)

      "%dopar%" <- foreach::"%dopar%"
      i <- NULL
      res <- foreach::foreach(i = 1:length(N_e), .combine = 'rbind') %dopar%{

        ocs <- oc_rct_dyn(N_c[i], N_e[i], delta, delta_power, confidence, e_a, e_b, c_a, c_b, RR_h, N_h, h_a, h_b, w, plot = F)

        av_power[i] <- mean(ocs[,2]*stats::dbeta(trues,alpha_c,beta_c))
        av_alpha[i] <- mean(ocs[,1]*stats::dbeta(trues,alpha_c,beta_c))

        return(c(av_alpha[i], av_power[i], N_c[i], N_e[i]))

      }

      # end parallel
      doParallel::stopImplicitCluster()
      # closeAllConnections()
      parallel::stopCluster(cl)

      if(plot){
        graphics::par(mar = c(5,6,7,5))
        plot(N_e, res[,1], type = "l", xlab = "Sample size in experimental group", ylab = "Probability", lty = 2, ylim = c(0,1))
        graphics::title("Average decision power and decision type 1 error with respect to \n sample size and uncertainty about the true control response rate", line = 4.7)
        graphics::lines(N_e, res[,2])
        graphics::abline(h = 0.10, lty = 3)
        graphics::abline(h = 0.60, lty = 3)
        if(legend){
          legend(legend.pos, legend = c("Avg. Dec. Alpha", "Avg. Dec. Power"), bty = "n", lty = c(2,1))
        }
        graphics::par(new = T)
        plot(N_c, rep(5, length(N_c)), xlab = "", ylab = "", xaxt = "n", yaxt = "n", pch = NA_integer_)
        graphics::axis(side=3)
        graphics::mtext("Sample Size in the control group", side = 3, line = 2.5)
      }

      colnames(res) <- c("Avg. Dec. Alpha", "Avg. Dec. Power", "Ne", "Nc")
      return(round(res,2))

    }else{

      res <- oc_rct_dyn(N_c, N_e, delta, delta_power, confidence, e_a, e_b, c_a, c_b, RR_h, N_h, h_a, h_b, w, plot = plot, legend.pos=legend.pos)
      av_power <- mean(res[,2]*stats::dbeta(trues,alpha_c,beta_c))
      av_alpha <- mean(res[,1]*stats::dbeta(trues,alpha_c,beta_c))

      if(plot){
        plot(trues, res[,1], xlab = "True control response rate", ylab = "Probability", type = "l", ylim = c(0,max(res)), lty = 2, main = paste("Average Decision alpha, Decision power = \n", round(av_alpha,2),",", round(av_power,2)))
        graphics::lines(trues, res[,2], type = "l")
        x <- NULL
        graphics::curve(1/(2*max(alpha_c, beta_c))*stats::dbeta(x,alpha_c,beta_c), from = 0, to = 1, n=101, add = T, lty = 3)
        if(legend){
          legend(legend.pos, legend = c("Decision type 1 error", "Decision power", "Beta density"), bty = "n", lty = c(2,1,3))
        }
      }

      ret <- round(c(av_alpha, av_power),2)
      names(ret) <- c("Avg. Dec. Alpha", "Avg. Dec. Power")
      return(ret)


    }


  }else{

    if(length(N_e)>1){

      av_power <- numeric(length(N_e))
      av_alpha <- numeric(length(N_e))

      if(is.null(coresnum)){coresnum <- parallel::detectCores()-1}
      cores <- coresnum
      cl <- parallel::makePSOCKcluster(cores)
      doParallel::registerDoParallel(cl)

      "%dopar%" <- foreach::"%dopar%"
      res <- foreach::foreach(i = 1:length(N_e), .combine = 'rbind') %dopar%{

        ocs <- oc_rct(N_c[i], N_e[i], delta, delta_power, confidence, e_a, e_b, c_a, c_b, plot = F)

        av_power[i] <- mean(ocs[,2]*stats::dbeta(trues,alpha_c,beta_c))
        av_alpha[i] <- mean(ocs[,1]*stats::dbeta(trues,alpha_c,beta_c))

        return(c(av_alpha[i], av_power[i], N_c[i], N_e[i]))

      }

      # end parallel
      doParallel::stopImplicitCluster()
      # closeAllConnections()
      parallel::stopCluster(cl)

      if(plot){
        graphics::par(mar = c(5,6,7,5))
        plot(N_e, res[,1], type = "l", xlab = "Sample size in experimental group", ylab = "Probability", lty = 2, ylim = c(0,1))
        graphics::title("Average decision power and decision type 1 error with respect to \n sample size and uncertainty about the true control response rate", line = 4.7)
        graphics::lines(N_e, res[,2])
        graphics::abline(h = 0.10, lty = 3)
        graphics::abline(h = 0.60, lty = 3)
        if(legend){
          legend(legend.pos, legend = c("Avg. Dec. Alpha", "Avg. Dec. Power"), bty = "n", lty = c(2,1))
        }
        graphics::par(new = T)
        plot(N_c, rep(5, length(N_c)), xlab = "", ylab = "", xaxt = "n", yaxt = "n", pch = NA_integer_)
        graphics::axis(side=3)
        graphics::mtext("Sample Size in the control group", side = 3, line = 2.5)
      }

      colnames(res) <- c("Avg. Dec. Alpha", "Avg. Dec. Power", "Ne", "Nc")
      return(round(res,2))

    }else{

      res <- oc_rct(N_c, N_e, delta, delta_power, confidence, e_a, e_b, c_a, c_b, plot = plot, legend.pos=legend.pos)
      av_power <- mean(res[,2]*stats::dbeta(trues,alpha_c,beta_c))
      av_alpha <- mean(res[,1]*stats::dbeta(trues,alpha_c,beta_c))

      if(plot){
        plot(trues, res[,1], xlab = "True control response rate", ylab = "Probability", type = "l", ylim = c(0,max(res)), lty = 2, main = paste("Average Decision alpha, Decision power = \n", round(av_alpha,2),",", round(av_power,2)))
        graphics::lines(trues, res[,2], type = "l")
        graphics::curve(1/(2*max(alpha_c, beta_c))*stats::dbeta(x,alpha_c,beta_c), from = 0, to = 1, n=101, add = T, lty = 3)
        if(legend){
          legend(legend.pos, legend = c("Decision type 1 error", "Decision power", "Beta density"), bty = "n", lty = c(2,1,3))
        }

      }

      ret <- round(c(av_alpha, av_power),2)
      names(ret) <- c("Avg. Dec. Alpha", "Avg. Dec. Power")
      return(ret)

    }
  }

}
