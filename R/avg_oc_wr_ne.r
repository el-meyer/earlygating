#' Single Arm Average Operating Characteristics
#'
#' Function for calculating the average operating characteristics of two
#' single arm bayesian designs for early gating with respect to the sample size
#' in the experimental group and possible historical data.
#'
#' @param N_e Sample Size in the experimental group. Can be either a single value or a vector.
#'
#' @param true_RR_c Default value is NULL. If specified, will be used in the generated plots,
#'                  indicating the true achieved decision power and decision type 1 error.
#'                  If not specified, will be set to either RR_h or hist_RR_c,
#'                  depending on which was specified by the user.
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
#' @param adapt Level of adapting of experimental control rate to account for patient selection bias
#'              from phase II to phase III. Corresponds to \eqn{\xi}{\xi}. Default is 1, so no adapting.
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
#' @return Either a vector containing the average decision power and average alpha (if N_e has length 1),
#'         or a matrix containing the average decision power and average decision alpha (if N_e has length > 1),
#'         where every row corresponds to one value of N_e.
#'
#' @examples
#'
#' # Setting 1
#' avg_oc_wr_ne(
#'   N_e = 50, delta = 0.08, delta_power = 0.13,
#'   confidence = 0.6, hist_RR_c = 0.5,
#'   alpha_c = 15, beta_c = 13
#' )
#'
#' # Setting 2
#' avg_oc_wr_ne(
#'   N_e = 50, delta = 0.08, delta_power = 0.13,
#'   confidence = 0.6, RR_h = 0.5, N_h = 50,
#'   alpha_c = 15, beta_c = 13
#' )
#'
#' @export
avg_oc_wr_ne <- function(N_e, true_RR_c = NULL, delta, delta_power,
                         confidence, e_a = 0.5, e_b = 0.5, h_a = 0.5,
                         h_b = 0.5, RR_h = NULL, N_h = NULL, hist_RR_c = NULL,
                         alpha_c, beta_c, trues = seq(0, 1, 0.001), adapt = 1,
                         plot = T, coresnum = NULL, legend = T, legend.pos = "topleft"){

  true_RR_e <- trues # Probability to GO will be evaluated along the same sequence as power and type 1 error

  if(is.null(true_RR_c)){ # if no true_RR_c was submitted, use either hist_RR_c or RR_h. Needed for evaluating power and type 1 error
    if(!is.null(hist_RR_c)){ # in the plot where \theta_c doesnt vary
      true_RR_c <- hist_RR_c
    }else{
      true_RR_c <- RR_h
    }
  }

  # Function to compute the minimum required number of responders w/r to N_e, hist_RR_c, delta, confidence, e_a, e_b and adapt
  # This is function for setting 1
  req_resp_wr_true <- function(N_e, hist_RR_c, delta, confidence, e_a, e_b, adapt){
    crit <- 0
    e_a_temp <- e_a # starting values for Beta, i.e. 0 responders
    e_b_temp <- e_b + N_e
    while(stats::qbeta(1-confidence, e_a_temp, e_b_temp) <= hist_RR_c+delta){ # iteratively solving the inequation
      crit <- crit + 1
      e_a_temp <- e_a + adapt*crit
      e_b_temp <- e_b - adapt*crit + N_e
      if(crit > N_e){ # if impossible to achieve number of responders with given sample size, then "NA"
        crit <- NA
        break
      }
    } # calculate critical value of "successes"
    return(crit)
  }

  # Function to compute the minimum required number of responders w/r to N_e, RR_h, N_h, delta, confidence, e_a, e_b and adapt
  # This is function for setting 2
  req_resp_wr_hists <- function(N_e, delta, confidence, e_a, e_b, h_a, h_b, RR_h, N_h, adapt){

    h_a_new <- h_a + RR_h*N_h
    h_b_new <- h_b - RR_h*N_h + N_h # it is assumed that "RR_h*N_h" responders were observed in historical control group
    crit <- 0
    e_a_temp <- e_a # starting values for Beta, i.e. 0 responders
    e_b_temp <- e_b + N_e
    while(stats::integrate(function(y) { stats::dbeta(y, e_a_temp, e_b_temp)*
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
    }
    return(crit)
  }

  # Function for Setting 1 in order to calculate power and type 1 error for a given setting and a fixed \theta_c
  # ret == should there be a return value
  # other parameters as above
  oc <- function(true_RR_e, N_e, true_RR_c, hist_RR_c, delta, confidence, e_a, e_b, delta_power,
                 adapt, plot = F, ret = T){

    # firstly calculate how many "successes" are needed for "GO" Decision
    crit <- req_resp_wr_true(N_e, hist_RR_c, delta, confidence, e_a, e_b, adapt)

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
    if(ret){return(c(alpha_ach, power_ach, crit))}
    #return achieved power, alpha and critical value
  }

  # Function for Setting 2 in order to calculate power and type 1 error for a given setting and a fixed \theta_c
  # ret == should there be a return value
  # other parameters as above
  oc_dist <- function(true_RR_e, N_e, true_RR_c, delta, confidence, e_a, e_b, h_a, h_b,
                      delta_power, RR_h, N_h, adapt, plot = F, ret = T){

    # firstly calculate how many "successes" are needed for "GO" Decision
    crit <- req_resp_wr_hists(N_e, delta, confidence, e_a, e_b, h_a, h_b, RR_h, N_h, adapt)

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

    if(ret){return(c(alpha_ach, power_ach, crit))}
  }

  # Function to compute average OCs in Setting 1. Additionally to function "oc" and "oc_dist", the user can now specify alpah_c and beta_c,
  # the parameters of the Beta distribution of the true control RR.
  # depending on whether (RR_h AND N_h) OR hist_RR_c were specified, the function uses either
  # the functions from Setting 1 or the functions from Setting 2
  average_oc <- function(true_RR_e, N_e, true_RR_c, delta, delta_power, confidence,
                         e_a, e_b, h_a, h_b, RR_h, N_h, hist_RR_c, adapt, alpha_c, beta_c, trues, plot = T, legend.pos){

    if(is.null(N_h) | is.null(RR_h)){ # find out whether Setting 1 or Setting 2

      n_resp_wr_trues <- numeric(length(trues))
      if(plot){
        for(k in 1:length(trues)){
          n_resp_wr_trues[k] <- try(req_resp_wr_true(N_e, trues[k], delta, confidence, e_a, e_b, adapt), silent = T)
          if(class(n_resp_wr_trues[k]) == "character"){n_resp_wr_trues[k] <- NA}
        } # firstly compute minimum required number of successes with respect to control RR
        if(plot){
          plot(trues, n_resp_wr_trues, ylab = "N Successes Req for GO", xlab = "Historical Control response rate", type = "s", main = "Required Number of successes with \n respect to the historical control response rate")
          graphics::abline(h = req_resp_wr_true(N_e, hist_RR_c, delta, confidence, e_a, e_b, adapt),
                 v = hist_RR_c)
        }
      }

      if(plot){
        # Now draw only the plots (== ret=F)
        oc(true_RR_e, N_e, true_RR_c, hist_RR_c, delta, confidence, e_a, e_b, delta_power, adapt, plot = plot, ret = F)
      }

      # Now compute OCs with respect to true_control_RR (i.e. \theta_c is now varying)
      oc_wr_trues <- matrix(nrow = length(trues), ncol = 2)
      crit <- req_resp_wr_true(N_e, hist_RR_c, delta, confidence, e_a, e_b, adapt)
      for(j in 1:length(trues)){
        if(trues[j]+delta_power<=1){# power at true_RR_c + delta power, unless true_RR_c + delta power > 1, in which case power = 1
          oc_wr_trues[j,2] <- stats::pbinom(crit-1, N_e, trues[j]+delta_power, lower.tail = FALSE)
        }else{oc_wr_trues[j,2] <- 0}
        oc_wr_trues[j,1] <- stats::pbinom(crit-1, N_e, trues[j], lower.tail = FALSE) # alpha at true_RR_c
      }

      if(plot){
        plot(trues,oc_wr_trues[,1], type = "l", lty = 2, ylab = "Probability", xlab = "True control response rate", main = "Decision power and Decision type 1 error \n w/r to True Control response rate")
        graphics::lines(trues,oc_wr_trues[,2], type = "l", lty = 1)
        if(legend){
        legend(legend.pos, legend = c("Decision type 1 error", "Decision power"), bty = "n", lty = c(2,1))
        }
      }

      # calculate average OCs as mean of OC*pdf of Beta for true control RR
      av_power <- mean(oc_wr_trues[,2]*stats::dbeta(trues,alpha_c,beta_c))
      av_alpha <- mean(oc_wr_trues[,1]*stats::dbeta(trues,alpha_c,beta_c))

      if(plot){
        plot(trues,oc_wr_trues[,1], type = "l", lty = 2, ylab = "Probability", xlab = "True control response rate",
             main = paste("Average Decision alpha, Decision power = \n ", round(av_alpha,2), ",", round(av_power,2)))
        graphics::lines(trues,oc_wr_trues[,2], type = "l", lty = 1)
        graphics::abline(v = hist_RR_c)
        x <- NULL
        graphics::curve(1/max(alpha_c, beta_c)*stats::dbeta(x,alpha_c,beta_c), from = 0, to = 1, n=101, add = T, lty = 3)
        if(legend){
        legend(legend.pos, legend = c("Decision type 1 error", "Decision power", "Beta density"), bty = "n", lty = c(2,1,3))
        }
      }
    }else{ # if in Setting 2, do everything analogously

      hists <- seq(0, 1, 0.01)
      if(plot){
        n_resp_wr_hists <- numeric(length(hists))
        for(k in 1:length(hists)){
          n_resp_wr_hists[k] <- try(req_resp_wr_hists(N_e, delta, confidence, e_a, e_b, h_a, h_b, hists[k], N_h, adapt), silent = T)
          if(class(n_resp_wr_hists[k]) == "character"){n_resp_wr_hists[k] <- NA}
        }
        plot(hists, n_resp_wr_hists, ylab = "N Successes Req for GO", xlab = "historical control response rate",type = "s", main = "Required Number of successes with \n respect to the historical control response rate")
        graphics::abline(h = req_resp_wr_hists(N_e, delta, confidence, e_a, e_b, h_a, h_b, RR_h, N_h, adapt), v = RR_h)
      }

      if(plot){
        oc_dist(true_RR_e, N_e, true_RR_c, delta, confidence, e_a, e_b, h_a, h_b, delta_power, RR_h, N_h, adapt, plot=plot, ret = F)
      }

      oc_wr_trues_dist <- matrix(nrow = length(trues), ncol = 2)
      crit <- req_resp_wr_hists(N_e, delta, confidence, e_a, e_b, h_a, h_b, RR_h, N_h, adapt)
      for(j in 1:length(trues)){
        if(trues[j]+delta_power<=1){# power at true_RR_c + delta power, unless true_RR_c + delta power > 1, in which case power = 1
          oc_wr_trues_dist[j,2] <- stats::pbinom(crit-1, N_e, trues[j]+delta_power, lower.tail = FALSE)
        }else{oc_wr_trues_dist[j,2] <- 0}
        oc_wr_trues_dist[j,1] <- stats::pbinom(crit-1, N_e, trues[j], lower.tail = FALSE) # alpha at true_RR_c
      }

      if(plot){
        plot(trues,oc_wr_trues_dist[,1], type = "l", lty = 2, ylab = "Probability",
             xlab = "True Control response rate",
             main = "Decision power and Decision type 1 error \n w/r to True Control response rate")
        graphics::lines(trues, oc_wr_trues_dist[,2], type = "l", lty = 1)
        if(legend){
        legend(legend.pos, legend = c("Decision type 1 error", "Decision power"), bty = "n", lty = c(2,1))
        }
      }

      av_power <- mean(oc_wr_trues_dist[,2]*stats::dbeta(trues,alpha_c,beta_c))
      av_alpha <- mean(oc_wr_trues_dist[,1]*stats::dbeta(trues,alpha_c,beta_c))

      if(plot){
        plot(trues,oc_wr_trues_dist[,1], type = "l", lty = 2, ylab = "Probability",
             xlab = "True control response rate",
             main = paste("Average Decision alpha, Decision power = \n ", round(av_alpha,2), ",", round(av_power,2)))
        graphics::lines(trues,oc_wr_trues_dist[,2], type = "l", lty = 1)
        graphics::abline(v = hist_RR_c)
        graphics::curve(1/max(alpha_c, beta_c)*stats::dbeta(x,alpha_c,beta_c), from = 0, to = 1, n=101, add = T, lty = 3)
        if(legend){
        legend(legend.pos, legend = c("Decision type 1 error", "Decision power", "Beta density"), bty = "n", lty = c(2,1,3))
        }
      }

    }
    ret <- c(av_alpha, av_power)
    names(ret) <- c("Avg. Dec. Alpha", "Avg. Dec. Power")
    return(ret)

  }


  if(length(N_e)==1){average_oc(true_RR_e=true_RR_e, N_e=N_e, true_RR_c=true_RR_c, delta=delta,
                                delta_power=delta_power, confidence=confidence, e_a=e_a, e_b=e_b, h_a=h_a,
                                h_b=h_b, RR_h=RR_h, N_h=N_h, hist_RR_c=hist_RR_c, adapt=adapt,
                                alpha_c=alpha_c, beta_c=beta_c, trues=trues, plot=plot, legend.pos=legend.pos)

  }else{

    avg_oc_wr_Ne <- matrix(nrow = length(N_e), ncol = 2)

    if(is.null(coresnum)){coresnum <- parallel::detectCores()-1}
    cores <- coresnum
    cl <- parallel::makePSOCKcluster(cores)
    doParallel::registerDoParallel(cl)
    "%dopar%" <- foreach::"%dopar%"
    z <- NULL
    avg_oc_wr_Ne <- foreach::foreach(z = 1:length(N_e), .combine = 'rbind') %dopar% {
      average_oc(true_RR_e=true_RR_e, N_e=N_e[z], true_RR_c=true_RR_c, delta=delta,
                 delta_power=delta_power, confidence=confidence, e_a=e_a, e_b=e_b,
                 h_a=h_a, h_b=h_b,
                 RR_h=RR_h, N_h=N_h, hist_RR_c=hist_RR_c, adapt=adapt,
                 alpha_c=alpha_c, beta_c=beta_c, trues=trues, plot = F)
    }

    if(plot){
      plot(N_e, avg_oc_wr_Ne[,1], type = "l", ylab = "Probability", xlab = "Sample Size Experimental",
           ylim = c(0,1), lty = 2, main = "Average OCs w/r to experimental group sample size")
      graphics::lines(N_e, avg_oc_wr_Ne[,2], type = "l", xlab = "Sample Size Experimental")
      if(legend){
        legend(legend.pos, legend = c("Avg. Dec. Alpha", "Avg. Dec. Power"), bty = "n", lty = c(2,1))

      }
    }

    # end parallel
    doParallel::stopImplicitCluster()
    # closeAllConnections()
    parallel::stopCluster(cl)

    colnames(avg_oc_wr_Ne) <- c("Avg. Dec. Alpha", "Avg. Dec. Power")
    rownames(avg_oc_wr_Ne) <- N_e
    return(avg_oc_wr_Ne)
  }
  }
