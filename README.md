
<!-- README.md is generated from README.Rmd. Please edit that file -->

# msissf

<!-- badges: start -->
<!-- badges: end -->

The goal of `msissf` is to fit Markov-Switching integrated
Step-Selection Functions.

## Installation

You can install the development version of `msissf` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jmsigner/msissf")
```

## Example

``` r
library(amt)
#> 
#> Attaching package: 'amt'
#> The following object is masked from 'package:stats':
#> 
#>     filter
library(msissf)
set.seed(123)

# Load example data set
data(deer) 
forest <- get_sh_forest()
```

Next, prepare the data set using `amt`.

``` r
rs <- deer |> 
  filter_min_n_burst(min_n = 4) |> 
  steps_by_burst() |> 
  random_steps() |> extract_covariates(forest)
```

Specify the msiSSF:

``` r
N <- 2
f <- list(
  habitat.selection = case_ ~ forest + strata(step_id_),
  step.length = ~ I(sl_ * -1) + log(sl_), 
  turn.angle = ~ cos(ta_)
)
res1 <- fit_msissf(
  formula = f,
  data = rs,
  start.values = generate_starting_values(f, N = N, rs),
  N = N,
  burst_ = "burst_", decode.states = TRUE, print.level = 0
)
```

And inspect results

``` r
res1$beta
#> [1] 0.5520501 1.4571744
res1$CI
#> $CI_up
#> [1] 2.143323 1.904994
#> 
#> $CI_low
#> [1] -1.039223  1.009354
#> 
#> $p.value
#> [1] 4.965313e-01 1.799034e-10

res1$shape
#> [1] 0.791280 1.004012
res1$rate
#> [1] 0.0002112867 0.0025676928
head(res1$decoded.states, 20)
#>  [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1

library(tidyverse)
#> ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
#> ✔ ggplot2 3.4.1     ✔ purrr   1.0.1
#> ✔ tibble  3.1.8     ✔ dplyr   1.1.0
#> ✔ tidyr   1.3.0     ✔ stringr 1.5.0
#> ✔ readr   2.1.4     ✔ forcats 1.0.0
#> ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
#> ✖ dplyr::filter() masks amt::filter(), stats::filter()
#> ✖ dplyr::lag()    masks stats::lag()
x <- 1:2400
bind_rows(
  tibble(x = x, 
         y = dgamma(x, rate = res1$rate[1], shape = res1$shape[1]), state = "1"), 
  tibble(x = x, 
         y = dgamma(x, rate = res1$rate[2], shape = res1$shape[2]), state = "2")
) |> 
  ggplot(aes(x, y, col = state)) + geom_line() + theme_light()
```

<img src="man/figures/README-example3-1.png" width="100%" />
