# Functions for monitoring harm boundary
# --------------------------------------------------

#' getAlphaPerTest
#'
#' Test-wise alpha necessary to control the overall type I error at a specified level (0.05 by default)
#'
#' @param nevents vector with number of events at which an interim analysis is done
#' @param totalAlpha Overall type I error, 0.05 by default
#' @param pH0 proportion of events in the experimental arm under the null hypothesis,
#'	typically based on randomization ratio (e.g. 0.5 for a 1:1 randomization)
#' @param alpha.interval Range for test-wise alpha, c(10^(-10),0.05) by default
#'
#' @return Test-wide alpha
#'
#' @export
#'
#' @importFrom stats uniroot
#'
#' @examples
#'	getAlphaPerTest(nevents = c(10,50,100), totalAlpha = 0.05, pH0 = 0.5)
#'
#'
getAlphaPerTest <- function(nevents,
	totalAlpha=0.05,
	pH0 = 0.5,
	alpha.interval=c(10^(-10), 0.05)){
	
    getCumAlpha <- function(alphaPerTest, nevents, pH0) {
        harmBounds <- getHarmBound(
                          nevents = nevents,
                          alpha_test = alphaPerTest,
                          pH0 = pH0)
        return(harmBounds$cumStopProb[nrow(harmBounds)] - totalAlpha)
    }

	ur<-uniroot(getCumAlpha, interval = alpha.interval, tol=1e-7,
		nevents = nevents, pH0=pH0)
			
	return(ur$root)

}



#' Harm boundaries for safety testing
#'
#' Calculates the boundaries at each interim analysis, i.e. the number of events in the experimental group
#' that would lead to a stopping of the trial based on a binomial exact test,
#' assuming that the events should be equally distributed amont both groups.
#' The indicated scenario (and all more extreme) would lead to a rejection of H0 (equal distribution) and a stopping for safety.
#'
#' @param nevents vector with number of events at which an interim analysis is done
#' @param alpha_test the nominal alpha level to use for each test
#' @param pH0 proportion of events in the experimental arm under the null hypothesis,
#'	typically based on randomization ratio (e.g. 0.5 for a 1:1 randomization)
#' @return data frame with a row for each interim analysis and with columns for
#'	number of events (n),
#'	number of events in control and experimental group that would lead to a stop (n_treat, n_control),
#'	risk ratio (RR) at that level,
#'	sopping probability (stopProb),
#'	cumulative stopping probability (cumStopProb).
#'	and the nominal alpha for each test (alphaVal)
#'
#' @export
#'
#' @importFrom stats dbinom 
#'
#' @examples
#' getHarmBound(nevents=c(10,50,100), alpha_test=0.025, pH0=0.5)

getHarmBound <- function(nevents,alpha_test,pH0){
		
	# create data frame to store results in
	bounds <- data.frame(totevents=1:max(nevents), treatBound=NA,
		alphaLevelBound=NA, cutoff=NA)
	
	bound <- NULL
	for (j in 1:nrow(bounds)) {
		
		totevents <- bounds$totevents[j]
		
		if (!(totevents %in% nevents)) {
			alphaVal<-0
		} else {
			alphaVal<-alpha_test
		}
		
		## we don't need to do the next few steps unless alphaVal is > 0
		if (alphaVal <= 0) next
		
		## choose the lowerBound for searching for the next cutoff value.
		if (is.null(bound)) {
			lowerBnd <- ceiling(pH0 * totevents)
		} else {
			lowerBnd <- bound
		}
		
		valSeq <- totevents:lowerBnd
		
		upperTailProbs <- cumsum(dbinom(valSeq, totevents, pH0))
		signif <- (upperTailProbs <= alphaVal)
		
		## if we have at least one significant value then do...
		if (isTRUE(signif[1]))	{
			## get "largest" (last) index for which signif == TRUE
			largest.index <- max(which(signif))
		
			## define 'bound' to be the infection count corresponding to
			## the 'largest.index' (i.e. the smallest infection count for
			## which we have significance at per-test-level 'alphaVal'
			bound <- valSeq[largest.index]
		
			bounds$treatBound[j] <- bound
			bounds$alphaLevelBound[j] <- upperTailProbs[largest.index]
			bounds$cutoff[j] <- alphaVal
		}
		
	}
	
	if (!all(is.na(bounds$treatBound))) {
	out <- pNS(Bounds=bounds$treatBound, pH0=pH0)
	} else {
	n <- nrow(bounds)
	out <- list( Bounds= data.frame(
							n=1:n,
							StoppingBound = bounds$treatBound ),
				Stop = rep(0, n),
				cutoff = bounds$cutoff )
	}
	
	names(out$Bounds)[ names(out$Bounds)=="StoppingBound" ] <- "n_treat"
	boundOut <- within(out$Bounds, {
						n_control <- n-n_treat
						RR <- n_treat*(1-pH0)/(n_control*pH0)
					})
	boundOut <- cbind(boundOut,
					stopProb = out$Stop,
					cumStopProb=cumsum(out$Stop),
					alphaVal = bounds$cutoff)
	
	## reorder columns
	ord.columns <- c("n","n_treat","n_control","RR","stopProb","cumStopProb","alphaVal")
	boundOut <- boundOut[, ord.columns]
	
	#only show those done:
	boundOutna<-boundOut[boundOut$n %in% nevents, ]
	
	return(boundOutna)
}


#' Title
#' Helper function which creates an object containing the values of the function P(n,s)
#' defined by Breslow (1970, JASA) as, for 0 <= s <= n <= N,  the probability
#' of the binomial random walk S_n reaching S_n = s without "absorption" into
#' the rejection region (i.e. without hitting any of the stopping bounds).
#'
#' @param Bounds Vector of stopping bounds after each event. 
#'	For event totals where stopping is not permitted the 'Bound' should be set to NA.
#' @param pH0 proportion of events in the experimental arm under the null hypothesis,
#'	typically based on randomization ratio (e.g. 0.5 for a 1:1 randomization)
#' @param returnPns Whether to inlcude the binomial random walk in the output
#'
#' @noRd
#'
#' @return list with stopping probabilities and boundaries
#'
#' @examples 
#'
#'	pNS(Bounds=c(rep(NA,5),3), pH0=0.5, returnPns=TRUE)
#'
pNS <- function(Bounds, pH0=0.5, returnPns=FALSE){
	
	N <- length(Bounds)
	
	if (all(is.na(Bounds))) {
		stop("The vector provided for argument 'Bounds' contains only NAs.\n",
			"Exiting...\n\n")
	}	 
	
	if (any(Bounds > (1:N), na.rm=TRUE)) {
		stop("The bounds provided do not appear to correspond to the",
			" total number of events (1,2,...,N).\n", "One of more of the bounds are larger",
			"than their corresponding total.  Exiting...\n\n")
	}
	
	## create 'Pns' and 'Stop'. Note that Stop has initial values of 0.
	Stop <- numeric(N)
	Pns  <- vector("list", length=N)
	
	## 'first.bound' is the point at which stopping is first allowed
	first.bound <- which(!is.na(Bounds))[1]
	
	if (first.bound == 1)
	stop("Cannot stop at the first test. Why would you want to?\n\n")
	
	
	## Initialize the base level (i==1).
	Pns[[1]] <- c("0" = (1-pH0), "1" = pH0)
	
	## loop over remaining infection totals, building up the prob.s as we go
	for (i in 2:N) {
	
	Pns[[i]] <- structure(numeric(i+1), names= as.character(0:i))
	
	max.S <- ifelse(is.na(Bounds[i]), i, Bounds[i]-1)
	
	Pns[[i]]["0"] <- (1-pH0)*Pns[[i-1]]["0"]
	
	for (s in 1:max.S ) {
		S         <- as.character(s)
		S.minus.1 <- as.character(s-1)
	
		if (s < i) {
			Pns[[i]][S] <- pH0*Pns[[i-1]][S.minus.1] + (1-pH0)*Pns[[i-1]][S]
		} else {
			Pns[[i]][S] <- pH0*Pns[[i-1]][S.minus.1]
		}	
	}
	#all probs outside the boundary are 0->trials with are stopped 
	#sum(Pns[[]]) < 1 after first stop 
	
	## 'Stop' is the prob. that we encounter a stopping bound for the first
	## time at infection count 'i'.  Only computed if stopping is possible,
	## else remains at initialized value of 0.
	if (!is.na(Bounds[i])) {
	
		## If this is the first infection count at which we allow stopping, then
		## the stop value must be computed allowing for people using non-standard
		## bounds (i.e. ones that that *don't* begin with bound[n] = n).
		if ((i == first.bound | is.na(Bounds[i-1])) && Bounds[i] < i ) {
			## stop prob. = 
			#	having one event less in the previous round * prob of hving a further event plus
			#	no has already nee
			## that we were already at/above the number of "success" needed to stop
			## at the previous infection total (i.e. i-1) (for which we should have
			## been given a stopping bound, but weren't).
			Stop[i] <- pH0*Pns[[i-1]][as.character(Bounds[i]-1)] +
						sum(Pns[[i-1]][as.character(Bounds[i]:(i-1))])
			} else {
			## if 'i' is not the the first total at which stopping is allowed, or
			## if it is, but we have a proper bound (i.e. Bounds[i] = i) then we do
			## we compute as this:
			Stop[i] <- pH0*Pns[[i-1]][as.character(Bounds[i]-1)]
			}
		}
	}
	
	outObj <- list(
		totalStopProb = sum(Stop),
		Stop = Stop,
		Bounds  = data.frame(n=1:N, StoppingBound=Bounds),
		N = N,
		pH0 = pH0 )
	
	if (isTRUE(returnPns)) {
		outObj$Pns <- Pns
	}	
	
	return(outObj)
}

