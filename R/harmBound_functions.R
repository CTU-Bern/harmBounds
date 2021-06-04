# Functions for continuous monitoring harm boundary
# --------------------------------------------------
# Functions intended to be used directly:
#   getHarmBounds
#   getAlphaPerTest
#
# Utility functions used by above functions, not intended for direct use:
#   pNS
#   semiConstSpending
#

#' Alpha per test
#'
#' @param harmMonitorRange harmMonitorRange
#' @param null.p null.p
#' @param totalAlpha totalAlpha
#' @param alpha.interval alpha.interval
#'
#' @return
#' @export
#'
#' @examples
getAlphaPerTest <- function(harmMonitorRange,
                            null.p,
                            totalAlpha=0.05,
                            alpha.interval=c(0.000001, 0.05) ){
    getCumAlpha <- function( alphaPerTest, harmMonitorRange, null.p) {
        harmBounds <- getHarmBound(
                          N = harmMonitorRange[2],
                          per.test = alphaPerTest,
                          harmBoundRange = harmMonitorRange,
                          null.p = null.p)
        return( harmBounds$cumStopProb[ nrow(harmBounds) ] - totalAlpha )
    }

    return( uniroot(getCumAlpha, interval = alpha.interval, tol=1e-7,
                    harmMonitorRange = harmMonitorRange, null.p=null.p)$root )
}



#' Harm boundaries for safety testing
#' We consider the total number of infected participants to be 'N' and assume apriori
#' and equal likelihood of infection for vaccinees as for placebos.  So the number of
#' infected vaccinees should follow a binomial distribution with size N and probability
#' p=\code{null.p} (where \code{null.p} is the proportion of vaccinees in the trial).  We wish to
#' have a total probability of a "type I error" (stopping the trial for 'harm' when
#' there is no real difference between vaccine and placebo) of e.g .05.  Our approach
#' will be to test for harm after each new infection, and we wish to find a fixed
#' 'alpha' value to use for each test, so that the overall prob. of a false positive
#' over the course of the trial is .05.
#'
#' To accomplish this, we iteratively choose values alpha, generate a set of
#' 'stopping bounds' corresponding to that alpha, and then estimate the overall type I
#' error rate for those bounds.  We then adjust our alpha value (up or down), as needed,
#' and repeat until the type I error is as desired. Obviously, due to the discreteness
#' of the distribution, the desired type I error is unlikely to be attainable exactly.
#' But if the \code{harmBoundRange} is not too small, a value reasonably close to that desired
#' should be possible.
#'
#' @param N the total number of infections to create bounds for
#' @param per.test the nominal alpha level to use for each test
#' @param harmBoundRange a vector of length 2, giving the range of infection counts over
#'      which the testing will be done.  In standard usage, the upper end of this range
#'      should be the same value as given in argument \code{N}
#' @param null.p the probability that an infection occurs in a vaccinee, under the null
#'      hypothesis that infection is equally likely in vaccinees and placebo recipients.
#'      Hence \code{null.p} equals the fraction of the population randomized to vaccine.
#'      So under a 1:1 randomization \code{null.p}=1/2 and for a 2:1 randomization (Vacc:Plac),
#'      \code{null.p}= 2/3.
#' @param dataDir the directory in which to store the output CSV file of bounds
#'               If NULL, then not file is output.
#' @param outFile the name of the CSV file to contain the bounds. If NULL then a default
#'               filename is constructed from inputs.  This argument is only used if
#'               \code{dataDir} is non-NULL
#' @param verbose logical value controlling printing of message giving name/location
#'               that output CSV has been saved.  Only relevant if \code{dataDir} is non-NULL
#'
#' @return
#' @export
#'
#' @examples
getHarmBound <- function(N,
                         per.test,
                         harmBoundRange,
                         null.p,
                         dataDir = NULL,
                         outFile = NULL,
                         verbose = TRUE){


  ## create data frame to store results in
  bounds <- data.frame(totInfec=1:N, vaccInfecBound=NA, alphaLevelBound=NA,
                       cutoff=NA )

  ## initialize 'bound' variable, which is iteratively updated in the 'j' loop below
  bound <- NULL

  for (j in 1:nrow(bounds))
  {
    totInfec <- bounds$totInfec[j]

    alphaVal <- semiConstSpending( totInfec, alphaVals=c(0, per.test),
                                   startHarmMonitor = harmBoundRange)

    ## we don't need to do the next few steps unless alphaVal is > 0
    if (alphaVal <= 0) next

    ## choose the lowerBound for searching for the next cutoff value.
    ## Under our framework it will always be the same as, or higher than
    ## the cutoff from the previous (smaller) value of 'totInfec'
    if ( is.null(bound) ) {
      lowerBnd <- ceiling( null.p * totInfec )
    } else lowerBnd <- bound

    # startVal <- min(totInfec, bound)
    # lowerBnd <- startVal - 2
    valSeq <- totInfec:lowerBnd

    upperTailProbs <- cumsum( dbinom(valSeq, totInfec, null.p) )
    signif <- ( upperTailProbs <= alphaVal )

    ## if we have at least one significant value then do...
    if ( isTRUE(signif[1]) )
    {
      ## get "largest" (last) index for which signif == TRUE
      largest.index <- max( which( signif ) )

      ## define 'bound' to be the infection count corresponding to
      ## the 'largest.index' (i.e. the smallest infection count for
      ## which we have significance at per-test-level 'alphaVal'
      bound <- valSeq[ largest.index ]

      bounds$vaccInfecBound[ j ] <- bound
      bounds$alphaLevelBound[ j ] <- upperTailProbs[ largest.index ]
      bounds$cutoff[ j ] <- alphaVal
    }

  }

  # Only run pNS if bounds$vaccInfecBound is not all NAs (which can happen when this
  # function is called by 'getAlphaPerTest'.  If not run, then create an alternate 'out'
  # object to use, containing necessary components.
  if ( !all(is.na(bounds$vaccInfecBound) ) ) {
      out <- pNS(Bound=bounds$vaccInfecBound, p=null.p)
  } else {
    n <- nrow(bounds)
    out <- list( Bounds= data.frame(
                           n=1:n,
                           StoppingBound = bounds$vaccInfecBound ),
                 Stop = rep(0, n),
                 cutoff = bounds$cutoff )
  }


  names(out$Bounds)[ names(out$Bounds)=="StoppingBound" ] <- "Nvacc"
  boundOut <- within(out$Bounds, {
                       Nplac <- n-Nvacc
                       RR <- round( Nvacc*(1-null.p)/(Nplac*null.p), digits=2)
                     })
  boundOut <- cbind( boundOut,
                     stopProb=round(out$Stop,4),
                     cumStopProb=round(cumsum(out$Stop),4),
                     alphaVal = bounds$cutoff )



  ## -- Remainder of code just does some renaming/reordering and then outputs -- ##

  harmBounds <- boundOut

  # Rename columns from 'old.names' to 'new.names'
  old.names <- c("n","Nvacc","Nplac")
  new.names <- c("N","V","P")
  names(harmBounds)[ match(old.names, names(harmBounds)) ] <- new.names

  ## reorder columns
  ord.columns <- c("N","V","P","RR","stopProb","cumStopProb","alphaVal")
  harmBounds <- harmBounds[, ord.columns]

  if (!is.null(dataDir)) {

      ## if outFile wasn't specified by the user, create a value for it
      if (is.null(outFile) ) {
          outFile <- sprintf("harmBounds_N=%d_alphaPerTest=%6.4f_pVacc=%4.2f.csv",
                             N, round(per.test, 4), round(null.p, 2) )
      }

      write.csv(harmBounds, file.path(dataDir, outFile), row.names=FALSE)

      if (verbose) {
          cat("Potential-harm stopping boundaries saved in:\n",
              file.path(dataDir, outFile), "\n\n")
      }
  }
  return(harmBounds)
}


#' Title
#' Function pNS creates an object containing the values of the function P(n,s)
#' defined by Breslow (1970, JASA) as, for 0 <= s <= n <= N,  the probability
#' of the binomial random walk S_n reaching S_n = s without "absorption" into
#' the rejection region (i.e. without hitting any of the stopping bounds).
#'
#'  Notation:
#'    N   = max number of samples
#'    S_n = sum{ x_i, i=1,..,n }, and the {x_i} are iid bernoulli(p)
#'
#'  To calculate P(n,s), we need to specify the probability 'p' the value 'N',
#'  and a vector 'B' of length N that specifies the stopping boundaries
#'  associated with the values of 'n' from 1:N.  Hence, B[1] will be the stopping
#'  value for n=1, B[2] the stopping value for n=2, ..., and B[N] the stopping
#'  value for n=N.
#'
#'  The vector of boundaries should contain NAs for the initial components
#'  up until you want to allow stopping (e.g. for n=1,2,3,...,8(?)
#'  how far you go will depend on the value of 'p' used.
#'
#'  The object is a list with components:
#'    Bounds - data.frame with columns 'n' = 1:N and
#'              'StoppingBound'= contents of input argument 'Bound'
#'    p - value of 'p' passed to function
#'    N - value of 'N' passed to function
#'    totalStopProb - the cumulative probability of stopping by the N-th infection
#'   'Stop' - numeric vector of length N, with Stop[i] giving the prob.
#'            of hitting the stopping bound (for first time) at n=i
#'
#'   If argument 'returnPns' is TRUE, then 'Pns' is included in the output obj.
#'       Pns is not normally needed, which is why it's been made optional.
#'   Pns =  list of length N, each sublist containing a numeric vector.
#'            Pns[[ i ]] is a numeric vector of length i+1 containing values
#'            of P(n,s), for n=i and s in 0,..,i+1
#'
#'          Note that values of s >= Bound[i] are zero as they are part of
#'          the rejection region (hence there is zero probability of reaching
#'          them without hitting the boundary).
#'
#' @param Bounds Vector of stopping bounds which correspond to the infection totals
#'           1,2,3,... For infection totals where stopping is not permitted (for
#'           small infection totals, e.g. 1-7) the 'Bound' should be set to NA.
#' @param p The randomization probability for the active treatment, when only
#'       considering it and the placebo group (in cases when other groups exist).
#'       If there are 3 active groups and 1 placebo group, and the randomization
#'       prob.s are: 2/9, 2/9, 2/9, 3/9, (the last  being for the placebo group),
#'       then the null.p value for generating bounds for one active vs. placebo
#'       would be:   2/9 /(2/9 + 3/9) = 2/5
#' @param returnPns
#' @noRd
#'
#' @return
#'
#' @examples
pNS <- function(Bounds, p=.5, returnPns=FALSE){
  N <- length(Bounds)

  if ( all(is.na(Bounds)) )
    stop("The vector provided for argument 'Bounds' contains only NAs.\n",
         "Exiting...\n\n")

  if ( any(Bounds > (1:N), na.rm=TRUE) )
      stop("The bounds provided do not appear to correspond to the infection",
           " totals (1,2,...,N).\n", "One of more of the bounds are larger",
           "than their corresponding infection total.  Exiting...\n\n")

  ## create 'Pns' and 'Stop'. Note that Stop has initial values of 0.
  Stop <- numeric(N)
  Pns  <- vector("list", length=N)

  ## 'first.bound' is the point at which stopping is first allowed
  ## We checked above to make sure that not all values of Bounds are NAs
  ## so we shouldn't have any problem with this...
  first.bound <- which( !is.na( Bounds ) )[1]

  if (first.bound == 1)
    stop("Cannot stop at the first test. Why would you want to?\n\n")


  ## Initialize the base level (i==1). 'i' indexes the total number of infections
  Pns[[ 1 ]] <- c("0" = (1-p), "1" = p)

  ## loop over remaining infection totals, building up the prob.s as we go
  for (i in 2:N) {

    ## initialize vector to store prob.s of the i+1 possible outcome of the
    ## binomial(i, p) variable (the values are: 0,...,i)
    Pns[[ i ]] <- structure( numeric(i+1), names= as.character(0:i) )

    ## Only need to compute for up to the smaller of i and Bounds[i]-1
    ## And since Bounds[i] has to be <= i, Bounds[i] is always < i, so long
    ## as Bounds[i] exists (i.e. as long as it's not NA)
    max.S <- ifelse( is.na(Bounds[i]), i, Bounds[i]-1)

    Pns[[i]]["0"] <- (1-p)*Pns[[i-1]]["0"]
    for (s in 1:max.S )
    {
      ## 'S' is character version of 's'. Used for indexing.
      S         <- as.character( s )
      S.minus.1 <- as.character( s-1 )

      if (s < i )
        Pns[[i]][ S ] <- p*Pns[[i-1]][S.minus.1] + (1-p)*Pns[[i-1]][S]
      else
        ## done when s == i (only happens if Bounds[i] is NA)
        Pns[[i]][ S ] <- p*Pns[[i-1]][S.minus.1]
    }

    ## 'Stop' is the prob. that we encounter a stopping bound for the first
    ## time at infection count 'i'.  Only computed if stopping is possible,
    ## else remains at initialized value of 0.
    if ( !is.na(Bounds[i]) )

      ## If this is the first infection count at which we allow stopping, then
      ## the stop value must be computed allowing for people using non-standard
      ## bounds (i.e. ones that that *don't* begin with bound[n] = n).
      if (i == first.bound && Bounds[i] < i ) {

        ## stop prob. is same as in the 'normal' case (below) *plus* the prob.
        ## that we were already at/above the number of "success" needed to stop
        ## at the previous infection total (i.e. i-1) (for which we should have
        ## been given a stopping bound, but weren't).
        Stop[i] <- p*Pns[[ i-1 ]][ as.character(Bounds[i]-1) ] +
                     sum( Pns[[ i-1 ]][as.character(Bounds[i]:(i-1))] )
      } else {
        ## if 'i' is not the the first total at which stopping is allowed, or
        ## if it is, but we have a proper bound (i.e. Bounds[i] = i) then we do
        ## we compute as this:
        Stop[ i ] <- p*Pns[[ i-1 ]][ as.character(Bounds[i]-1) ]
      }
  }
  outObj <- list(
              totalStopProb = sum(Stop),
              Stop = Stop,
              Bounds  = data.frame(n=1:N, StoppingBound=Bounds),
              N = N,
              p = p )

  ## If 'returnPns' was set to TRUE then add Pns to the output object
  if ( isTRUE(returnPns) )
      outObj$Pns <- Pns

  return( outObj )
}
####################### end of function 'pNS'#############################


### THIS USES 'nonConstBounds' framework to do constant-bounds.
####################### end of function 'pNS'#############################

#' Title
#'
#' @param x the total number of infections (vacc + placebo)
#' @param alphaVals  is a vector of nominal (un-adjusted) p-value thresholds
#'    to use in establishing cutoffs for harm-monitoring.  This vector
#'    must have the same length as \code{startHarmMonitor}.  The i-th value
#'    of \code{alphaVals} applies to the i-th interval defined by \code{startHarmMonitor}
#' @param startHarmMonitor gives the endpoints of all intervals.
#'    The starting point of the first interval is 1, and of the the i-th interval
#'    (for i>1) is 1 + \code{startHarmMonitor}[i-1]
#'
#' @return
#'
#' @examples
semiConstSpending <- function(x, alphaVals, startHarmMonitor )
{
  which.interval <- findInterval(x, c(1, startHarmMonitor),
                                 rightmost.closed=TRUE)
  alphaVals[ which.interval ]
}


