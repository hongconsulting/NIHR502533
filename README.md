
<!-- README.md is generated from README.Rmd. Please edit that file -->

# NIHR502533

<!-- badges: start -->
<!-- badges: end -->

Sample size simulator for NIHR502533. This version calls the original
‘BuyseTest’ package via the ‘Rcpp’ interface. The surrounding ‘WH’ code
is designed to maximize speed for use in simulation studies. As such,
input validation is minimal and the user is responsible for ensuring
valid use.

Note that formal sample size calculations were performed using 100000
replicates.

## Installation

You can install the development version of NIHR502533 from
[GitHub](https://github.com/hongconsulting/NIHR502533) with:

``` r
remotes::install_github("hongconsulting/NIHR502533")
```

## Example: Cohort A (single-threaded)

``` r
library(NIHR502533)

print(WH_NIHR502533a(n_target = 204,
                    p_arm0surv = 0.73, p_arm1surv = 0.70,
                    p_arm0tox1 = 0.25, p_arm1tox1 = 0.1,
                    p_arm0tox2 = 0.45, p_arm1tox2 = 0.05, p_LTFU = 0.1, 
                    recrate = 7.5, reps = 100))
#> NIHR502533 cohort A: scenario 1
#>  - Target N = 204
#>  - RFS at 3 years = 73% vs 70%
#>  - Tox 1 = 25% vs 10%
#>  - Tox 2 = 45% vs 5%
#>  - Average recruitment = 7.5/month
#>  - Loss to follow-up = 10% every 5 years
#>  - Recruitment duration = up to 5 years
#>  - Max duration = 3 years after the last recruitment, up to 6 years
#> Loading required package: Rcpp
#> BuyseTest version 3.2.0
#>  - Average N = 193.38
#>  - Power (NTB) = 79.00%
#>  - Average duration = 5.27 years
#> [1] 193.38000   0.79000   5.27186
```

## Example: Cohort A (chunk parallelized)

``` r
library(NIHR502533)
library(parallel)

reps <- 1000
n_cores <- parallel::detectCores()
chunk_reps <- rep(ceiling(reps / n_cores), n_cores)
chunk_reps[length(chunk_reps)] <- reps - sum(head(chunk_reps, -1))
seeds <- seq_along(chunk_reps)

sims <- do.call(rbind, parallel::mclapply(seeds, function(i) {
  WH_NIHR502533a(n_target = 204,
                p_arm0surv = 0.73, p_arm1surv = 0.70,
                p_arm0tox1 = 0.25, p_arm1tox1 = 0.1,
                p_arm0tox2 = 0.45, p_arm1tox2 = 0.05,
                p_LTFU = 0.1, recrate = 7.5, reps = chunk_reps[i],
                seed = seeds[i], echo = FALSE)
}, mc.cores = n_cores))
print(colSums(sims * rep(chunk_reps, times = 3))/reps)
#> [1] 193.328000   0.830000   5.271666
```

## Example: Cohort B

``` r
library(NIHR502533)
print(WH_NIHR502533b(n_target = 19,
                     p_nullclear = 0.1, p_arm0clear = 0.34,
                     p_LTFU = 0.15, recrate = 7.5, reps = 1e5))
#> NIHR502533 cohort B: scenario 1
#>  - Target N = 19
#>  - Clearance at 1 years = 34%
#>  - Average recruitment = 7.5/month
#>  - Loss to follow-up = 15% every 5 years
#>  - Recruitment duration = up to 5 years
#>  - Max duration = 3 years after the last recruitment, up to 6 years
#>  - Average N = 18.50
#>  - Power = 81.62%
#>  - Average duration = 3.21 years
#> [1] 18.500900  0.816180  3.212489
```
