source("harmBound_functions.R")

## 1:1 randomization (Vacc:Plac) -> null.p = 1/2
## 2:1 randomization (Vacc:Plac) -> null.p = 2/3
null.p <- 1/2 

# total "alpha" to be spent across all tests
harmMonitorAlpha <- 0.05

## range of infection totals to perform testing over
harmMonitorRange <- c(1, 84)
#harmMonitorRange <- c(9, 105)


## 'getAlphaPerTest' determines the "per-test-alpha", the nominal level at which each
## each test is performed, such that the total type-I error across the tests equals
## (as near as possible) the values provided by arg. 'totalAlpha'
alphaPerTest <- getAlphaPerTest(harmMonitorRange, null.p,
                                totalAlpha = harmMonitorAlpha)

## Use the value 'alphaPerTest' from previous step to generate the set of stopping bounds
harmBounds <- getHarmBound(
                  N = harmMonitorRange[2],
                  per.test = alphaPerTest,
                  harmBoundRange = harmMonitorRange,
                  null.p = null.p)

print( harmBounds )

outFile <- paste0(
    "harmBounds_", paste(harmMonitorRange, collapse="_"), 
    "_p=", round(null.p, dig=2), "_alpha=", harmMonitorAlpha, ".csv")

write.csv(harmBounds, file=outFile, row.names=FALSE) 

q(save = "no")
