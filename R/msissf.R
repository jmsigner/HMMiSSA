
#' Fit a Markov-Switching integrated Step-Selection Function
#'
#' @param formula A list with three formulas.
#' @param N The number of states.
#' @param data An object from the `amt` package of class `random_steps`.
#' @param start.values Starting values, output from `generate_starting_values()`.
#' @param burst_ Column name of columns with the bursts.
#' @param stationary logical, TRUE if stationarity for underlying latent MC is assumed
#' @param iterlim maximum number of iterations for nlm, default is 1000
#' @param gradtol tolerance value for convergence in nlm, default is 1e-06 (default of nlm)
#' @param decode.states decode states
#' @param print.level print level for nlm, default is 2(all iterations printed)
#' @param stepmax maximum allowable scaled step length in nlm, default NULL corresponds to nlm default
#' @param conf.int return CI?
#' @param alpha alpha value for CI.
#'
#' @return
#' @export
#'



fit_msissf <- function(
  formula, data, N,
  start.values,
  burst_ = "burst_",
  stationary = FALSE,
  iterlim = 1000, gradtol = 1e-06,
  decode.states = FALSE,
  print.level = 2, stepmax=NULL, conf.int = TRUE, alpha = 0.05) {


  # 1) step IDs and covariates
  # place of stepID in formula
  stepID <- all.vars(formula$habitat.selection)[
    attr(stats::terms(formula$habitat.selection, "strata"), "specials")$strata]

  # vector of covariate names in formula
  covs <- c(formula$habitat.selection[[3]],
            formula$turn.angle[[2]],
            formula$step.length[[2]])

  covs <- covs |> lapply(\(x) x |> deparse() |> strsplit("\\+")) |> unlist()
  covs <- covs[!grepl("strata\\(.*\\)", covs)]
  covs <- trimws(covs)


  # number of covariates
  p <- length(covs)

  # name of used column
  used <- all.vars(formula$habitat.selection[[2]])

  # Make sure we have only complete observations
  data <- data[stats::complete.cases(data), ]

  # model matrix
  f <- stats::as.formula(paste("~ -1 +", paste0(covs, collapse = "+")))
  mm <- stats::model.matrix(f, data = data, na.action = stats::na.pass)

  # 2) bursts
  # vector including unique IDs
  vecIDs <- unique(data[[burst_]])
  used <- data[[used]]

  # number of IDs
  nID <- length(vecIDs)

  # vector for number of steps per ID
  nSteps <- numeric(nID)

  # X.list: double list including design matrix X for each step and each ID
  # first list index: ID
  # second list index: step
  X.list <- tidyr::nest(
    dplyr::mutate(
      as.data.frame(mm),
      burst_ = data[[burst_]], step_ = data[[stepID]]),
    data = -c(burst_, step_))

  # compute working parameters for initial values
  par.vect <- n2w(start.values, stationary)

  # set stepmax value for nlm
  if(is.null(stepmax)){
    stepmax <- max(1000 * sqrt(sum((par.vect)^2)), 1000)
  }

  sld <- amt::sl_distr(data)
  tad <- amt::ta_distr(data)

  # use nlm to find maximum likelihood estimates
  mod <- stats::nlm(nlogLike, par.vect, X.list, N, p, sld, tad,
             stationary,
             iterlim=iterlim,gradtol=gradtol,print.level=print.level,
             stepmax=stepmax, hessian=FALSE)

  # calculate natural parameters from maximum likelihood estimate vector
  pn <- w2n(mod$estimate, N, p, sld, tad, stationary)

  # beta without constraint movement parameters
  beta <- matrix(pn$beta_clr, p, N)
  beta.sel <- beta[1:(nrow(beta) - 3), ]


  # Compute the Hessian
  if (conf.int) {

    hes <- numDeriv::hessian(nlogLike, mod$estimate, X.list = X.list, N = N, p = p,
                             sld = sld, tad = tad,
                             stationary=stationary)

    # standard errors for beta_clr
    std.err <- try(sqrt(diag(MASS::ginv(hes))[c(1:(N*p))]))

    if(!methods::is(std.err, "try-error")) {
      # standard errors for beta
      se_beta <- matrix(std.err,ncol=N)[1:(p - 3),]
      # confidence limits
      CI_up <- beta.sel + stats::qnorm(1 - alpha/2) * se_beta
      CI_low <- beta.sel + stats::qnorm(alpha / 2) * se_beta
      p_beta <- abs(beta.sel / se_beta)
      p_beta <- (1 - stats::pnorm(p_beta))*2

      CI = list(CI_up = CI_up, CI_low = CI_low, p.value = p_beta)

    } else {
      se_beta<-CI_up<-CI_low<-p_beta<-matrix(NA,p_beta,N)
      CI = list()
    }
    mod$hessian <- hes
  } else {
    CI = list()
  }

  # add shape, rate and kappa parameters
  shape <- rate <- kappa <- NULL

  # compute shape parameters if sl_gamma=TRUE
  shape <- beta[p - 1, ] + sld$params$shape
  rate <- beta[p, ] + 1/sld$params$scale

  # compute kappa parameters if ta_vonMises=TRUE
  kappa <- beta[p - 2, ] + tad$params$kappa

  # Decode
  if (decode.states) {
    iv.list <- list()
    X.list$allprobs <- purrr::map(X.list$data, ~ {
      xx <- cpr_Nt(as.matrix(.x), pn$beta_clr)
      xx[is.na(xx) & !is.nan(xx)] <- 1 # Could go to the cpr_Nt function?
      xx
    })

    for(run_ID in unique(X.list$burst_)){
      tdat <- dplyr::filter(X.list, burst_ == run_ID)

      nSteps_run <- nrow(tdat)
      xi <- matrix(0,nSteps_run,N)

      foo <- pn$delta * tdat$allprobs[[1]]
      xi[1,] <- foo/sum(foo)

      for (i in 2:nSteps_run){
        foo <- apply(xi[i-1,] * pn$gamma, 2 , max) * tdat$allprobs[[i]]
        xi[i,] <- foo/sum(foo)
      }

      iv <- numeric(nSteps_run)
      iv[nSteps_run] <- which.max(xi[nSteps_run,])


      for (i in (nSteps_run - 1):1){
        iv[i] <- which.max(pn$gamma[, iv[i+1]] * xi[i, ])
      }
      iv.list <- c(iv.list, iv)
    }
    states <- do.call("c", iv.list)
  } else {
    states <- NULL
  }

  return(list(beta = beta.sel, shape = shape, rate = rate,
              kappa = kappa, gamma = pn$gamma,
              delta = pn$delta, beta_clr = pn$beta_clr,
              nllk = mod$minimum, estimate = mod$estimate, nit = mod$iterations,
              conv = mod$code, hessian = mod$hessian, CI = CI, decoded.states = states))
}
