#' Internal function to compute the negative log likelihood
#'
#' @param parvect Vector with parameters
#' @param X.list The data
#' @param N The number of sates
#' @param p The number of covariates
#' @param sld The step-length distribution
#' @param tad The turn-angle distribution
#' @param stationary Is Markov Chain assumed to be stationary?
#' @name internal_functions

nlogLike <- function(parvect, X.list, N, p, sld,
                     tad, stationary){

  # get natural parameters
  lpn <- w2n(parvect, N, p, sld, tad, stationary)
  # initialise log-likelihood value
  lscale_all <- 0

  # loop across bursts
  X.list$allprobs <- purrr::map(X.list$data, ~ {
    xx <- cpr_Nt(as.matrix(.x), lpn$beta_clr)
    xx[is.na(xx) & !is.nan(xx)] <- 1 # Could go to the cpr_Nt function?
    xx
  })


  for(run_ID in unique(X.list$burst_)){
    # compute matrix including clr-probs for each state and each step, respectively
    # matrix dimension TxN (T=nSteps[run_ID])
    # missing values are replaced by 1
    tdat <- dplyr::filter(X.list, burst_ == run_ID) # Data of the current burtst
    if (nrow(tdat) < 2) {
      stop(paste0("Burst, ", run_ID, " has < 2 steps"))
    }

    # This is now implemented in cpp
    
    # foo = lpn$delta * tdat$allprobs[[1]] # t=1
    # sumfoo = sum(foo)
    # lscale = log(sumfoo) # current loglikelihood-value
    # foo = foo / sumfoo
    # for(ts in 2:nrow(tdat)) { # loop across all time-steps
    #   foo = foo %*% lpn$gamma * tdat$allprobs[[ts]]
    #   sumfoo = sum(foo)
    #   lscale = lscale + log(sumfoo) # updated loglikelihood-value at t=ts
    #   foo = foo / sumfoo
    # }
    # lscale<-nll_Rcpp(allprobs, lpn$gamma, lpn$delta, nSteps[run_ID])
    # lscale_all <- lscale_all + lscale # sum log-likelihood values of all bursts
    
    allprobs = do.call(rbind,tdat$allprobs) # allprobs muss eine Matrix sein (keine Liste wie bisher)! Dies könnte auch direct im Xlist- oder tdat-Objekt geändert werden  
    lscale <- ll_Rcpp(allprobs, lpn$gamma, lpn$delta, nrow(allprobs))
    lscale_all <- lscale_all + lscale # sum log-likelihood values of all bursts    

  }
  # return negative log-likelihood value (as nlm minimises functions)
  return(-lscale_all)
}

#' Transform parameters from the working scale to the natural scale
#' @rdname internal_functions

w2n <- function(parvect, N, p, sld, tad, stationary=FALSE) {
  # re-create beta-matrix - but now including movement kernel parameters (beta_clr)
  beta_clr <- matrix(parvect[1:(N*p)],p,N)

  # !note!: instead of shape and rate, shape-shape* and rate-rate* are used in clr (as in standard iSSA)
  # shape* and rate* depend on sampling procedure for sl of available steps
  # i) gamma-distr. sl: shape* and rate* are shape and rate values used for sampling
  # ii) uniform sampled sl: shape*=1 and rate*=0
  # !note!: shape and rate - related parameters are in last and second-last row
  beta_clr[p-1,] <- exp(beta_clr[p-1,]) - sld$params$shape # shape-shape*
  beta_clr[p,] <- exp(beta_clr[p,]) - 1/sld$params$scale # rate-rate*

  # adjust coefficient for kappa if ta_vonMises=TRUE
  # !note!: kappa-kappa* is used in clr instead of kappa (as in standard iSSA)
  # kappa* depend on sampling procedure for angle of available steps
  # i) von-Mises distr. ta: kappa* is kappa value used for sampling
  # ii) uniform sampled ta: kappa*=0
  # !note!: kappa - related parameters are in third-last row
  beta_clr[p - 2, ] <- exp(beta_clr[p - 2, ]) - tad$params$kappa # kappa-kappa*

  # re-create tpm gamma
  gamma<-diag(N)
  gamma[!gamma]<-exp(parvect[N*p +1:(N*(N-1))])
  gamma<-gamma/apply(gamma,1,sum)

  # compute initial distr. delta
  if(stationary){
    delta<-solve(t(diag(N)-gamma+1),rep (1,N))
  }else{
    foo<-c(1,exp(parvect[(N*p+N*(N-1))+1:(N-1)]))
    delta<-foo/sum(foo)
  }
  return(list(beta_clr=beta_clr,gamma=gamma,delta=delta))
}

cpr_Nt <- function(X,beta_clr){
  S<-exp(X%*%beta_clr) # exp of matrix product Mxp %*% pxN
  S<-S[1,]/apply(S,2,sum) # fraction with exp(X_used %>% beta_clr) in numerator, sum in denominator
  return(S)
}


n2w <- function(start.values, stationary){

  N <- nrow(start.values$gamma)

  # no transformation needed for beta
  tbeta <- start.values$beta

  # add log of kappa to tbeta matrix if kappa is not NULL
  # (log transformation assures kappa>0 in estimation process)
  tkappa <- if(!is.null(start.values$kappa)) log(start.values$kappa) else NULL

  # add log of shape and rate to tbeta matrix if shape is not NULL
  # (log transformation assures shape>0 and rate>0 in estimation process)
  tshape <- if(!is.null(start.values$shape)) log(start.values$shape) else NULL
  trate <- if (!is.null(start.values$rate)) log(start.values$rate) else NULL

  # mlogit transformation for tpm to assure row sum = 1 and tps in (0,1)
  foo <- with(start.values, log(gamma/diag(gamma)))
  tgamma <- as.vector(foo[!diag(N)])

  # no values for delta needed in case if stationarity
  if(stationary){
    tdelta <- NULL
  }else {
  # mlogit transformation for delta in case of non-stationarity
    tdelta <- with(start.values, log(delta[-1] / delta[1]))
  }

  c(as.vector(rbind( tbeta, tkappa, tshape, trate)), tgamma, tdelta)
}





