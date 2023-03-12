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

  beta.start <- matrix(rep(sample(seq(-1, 1, len = N))), ncol = N, nrow = length(covs)) # make more random

  suggest_start <- function(x, N, tolerance = 0.2, range = c(-1, 1)) {
    x + (x * tolerance) * seq(range[1], range[2], len = N)
  }

  shape <- suggest_start(amt::sl_distr(data)$params$shape, N)
  rate <- suggest_start(1 / amt::sl_distr(data)$params$scale, N)
  kappa <- suggest_start(amt::ta_distr(data)$params$kappa, N)

  gamma <- matrix(stats::runif(N^2), nrow = N, ncol = N)
  gamma <- sweep(gamma, 1, rowSums(gamma), FUN = "/")

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
