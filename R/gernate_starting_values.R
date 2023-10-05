#' Generate Starting Values
#'
#' Helper function to automatically generate starting values given a formula and data.
#'
#' @inheritParams fit_msissf
#' @export
#'

generate_starting_values <- function(formula, N, data) {

  checkmate::assert_list(formula)
  checkmate::assert_number(N, lower = 2)
  checkmate::assert_class(data, "random_steps")

  stepID <- as.numeric(attr(stats::terms(formula$habitat.selection, "strata"), "specials")) - 1
  covs <- all.vars(formula$habitat.selection[[3]])[-stepID]

  beta.start <- matrix(sample(seq(-1, 1, len = N),N*length(covs), replace=TRUE), ncol = N, nrow = length(covs))

  # gamma distribution:
  # a) sample random mean values based on quantiles of observed step length
  quantiles_sl <- quantile(data[data$case_==TRUE,]$sl_,seq(0.1, 0.9, len = N+1)) |> as.numeric()
  mean_sl<-runif(N,quantiles_sl[1:N],quantiles_sl[2:(N+1)])

  # b) sample random standard deviations based on mean values
  sd_sl<-mean_sl*runif(N,min=0.25,max=2)

  # c) compute corresponding shape and rate parameter
  shape <- (mean_sl/sd_sl)^2
  rate<-mean_sl/sd_sl^2

  # von Mises: draw random concentration values:
  kappa <- amt::ta_distr(data)$params$kappa*runif(N,min=0.25,max=2)

  # transition probabilities: draw random probabilities in diagonal
  # -> rather large probability to stay in current state
  diag_gamma<-runif(N,0.75,0.99)
  gamma <- matrix(rep((1-diag_gamma)/(N-1),N), nrow = N, ncol = N)
  diag(gamma)<-diag_gamma

  delta <- stats::runif(N)
  delta <- delta / sum(delta)

  list(
    beta = beta.start,
    shape = shape,
    rate = rate,
    kappa = kappa,
    gamma = gamma,
    delta = delta
  )
}
