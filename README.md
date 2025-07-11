
<!-- README.md is generated from README.Rmd. Please edit that file -->

# harmBounds

<!-- badges: start -->

[![R-CMD-full](https://github.com/CTU-Bern/harmBounds/workflows/R-CMD-full/badge.svg)](https://github.com/CTU-Bern/harmBounds/actions)
[![R-CMD-standard](https://github.com/CTU-Bern/harmBounds/workflows/R-CMD-standard/badge.svg)](https://github.com/CTU-Bern/harmBounds/actions)
<!-- badges: end -->

The harmBounds package calculates boundaries, simulates data and
generates plots for safety monitoring using an event based approach.

The idea is to do simple one sample binomial exact tests on the
proportion of events in the experimental arm.

A safety problem is claimed if there is evidence that this is higher
than what would be expected from the number of people under observation
in the two arms (e.g. 0.5 in a 1:1 randomized trial).

This can be done continuously or at pre-specified total number of
events.The nominal test-wise alpha would be calibrated to obtain desired
properties, such as a certain proportion of stopping under null and
alternative scenarios or an overall type I error control (which we do
not necessarily recommend for safety testing).

## Installation

It can be installed from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("CTU-Bern/harmBounds")
```

``` r
devtools::load_all()
#> ℹ Loading harmBounds
```

## Example

### Boundaries

To get boundaries, e.g. for a trial with 3 safety interim analyses after
10, 50 and 100 events (in total over both groups) using a nominal
test-wise alpha of 0.025.

``` r
hb<-getHarmBound(nevents=c(10,50,100), alpha_test = 0.025, pH0 = 0.5)
hb
#>       n n_treat n_control       RR   stopProb cumStopProb alphaVal
#> 10   10       9         1 9.000000 0.01074219  0.01074219    0.025
#> 50   50      33        17 1.941176 0.01490030  0.02564249    0.025
#> 100 100      61        39 1.564103 0.01194192  0.03758441    0.025
```

‘n_treat’ indicates the minimal number of events in the treatment group
that would lead to a stopping of the trial. The overall type I error is
3.8%.

If we would like to control the type I error at a specific level, We can
get the test-wise alpha using:

``` r
alphaPerTest <- getAlphaPerTest(nevents = c(10,50,100), pH0= 0.5, totalAlpha = 0.05)

alphaPerTest
#> [1] 0.03245429

getHarmBound(nevents=c(10,50,100), alpha_test = alphaPerTest, pH0 = 0.5)
#>       n n_treat n_control       RR   stopProb cumStopProb   alphaVal
#> 10   10       9         1 9.000000 0.01074219  0.01074219 0.03245429
#> 50   50      33        17 1.941176 0.01490030  0.02564249 0.03245429
#> 100 100      60        40 1.500000 0.02085634  0.04649882 0.03245429
```

The overall type I error is 3.8%, i.e. as close to 5% as possible (given
the discrete nature of the test).

### Simulation

To simulate data under the null hypothesis (H0) and an alternative with
80% events in the experimental arm (H1).

``` r
set.seed(123)

sim_safety_stop(nevents = c(10,50,100), pH1 = 0.5, alpha_test = 0.025)
#> $tests
#>   nevents n1 n0 ulim   out
#> 1      10  6  4    9 FALSE
#> 2      50 25 25   33 FALSE
#> 3     100 47 53   61 FALSE
#> 
#> $nstop
#> [1] 0
#> 
#> $tstop
#> [1] NA

sim_safety_stop(nevents = c(10,50,100), pH1 = 0.6, alpha_test = 0.025)
#> $tests
#>   nevents n1 n0 ulim   out
#> 1      10  6  4    9 FALSE
#> 2      50 29 21   33 FALSE
#> 3     100 62 38   61  TRUE
#> 
#> $nstop
#> [1] 1
#> 
#> $tstop
#> [1] 3
```

And repeat that to get the proportion of trials with a stop (more
simulations would be necessary to get precise results)

``` r

set.seed(123)

ssp0<-vapply(1:100,
    function(x) sim_safety_stop(nevents = c(10,50,100), pH0 = 0.5, pH1 = 0.5,
        alpha_test=0.025)$nstop,
    numeric(1))

mean(ssp0>0)
#> [1] 0.02

ssp1<-vapply(1:100,
    function(x) sim_safety_stop(nevents = c(10,50,100), pH0 = 0.5, pH1 = 0.6,
        alpha_test=0.025)$nstop,
    numeric(1))

mean(ssp1>0)
#> [1] 0.53
```

We would stop in 2% of the trials under the null (no safety problem) and
in 53% under the alternative (safety problem).

We could now use different test-wise alphas to define an overall design
with desired properties, i.e. enough stopping under the alternative
while not stopping too often under the null.
