source("harmBound_functions.R")

## 1:1 randomization (Vacc:Plac) -> pH0 = 1/2
## 2:1 randomization (Vacc:Plac) -> pH0 = 2/3
pH0 <- 1/2 

# total "alpha" to be spent across all tests
totalAlpha <- 0.05

## total number of events at which interim analysis should be done
nevents <- c(10,50,100)

## 'getAlphaPerTest' determines the "per-test-alpha", the nominal level at which each
## each test is performed, such that the total type-I error across the tests equals
## (as near as possible) the values provided by 'totalAlpha'
alphaPerTest <- getAlphaPerTest(nevents = nevents, pH0= p0,
	totalAlpha = totalAlpha)

## Use the value 'alphaPerTest' from previous step to generate the set of stopping bounds
harmBounds <- getHarmBound(nevents = nevents,
	alpha_test = alphaPerTest,
	pH0 = p0)

print(harmBounds)
