# Functions for monitoring harm boundary
# --------------------------------------------------

#' getAlphaPerTest
#'
#' Test-wise alpha necessary to control the overall type I error at a specified level (0.05 by default)
#'
#' @param nevents vector with number of events at which an interim analysis is done
#' @param totalAlpha Overall type I error, 0.05 by default
#' @param pH0 proportion of events in the intervention arm under the null hypothesis,
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
#'	apt<-getAlphaPerTest(nevents = c(10,50,100), totalAlpha = 0.05, pH0 = 0.5)
#'	apt
#'	getHarmBound(nevents = c(10,50,100),alpha_test = apt, pH0 = 0.5)
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
        return(harmBounds$opchar[1,"cum_stop_prob"] - totalAlpha)
    }

	ur<-uniroot(getCumAlpha, interval = alpha.interval, tol=1e-7,
		nevents = nevents, pH0=pH0)

	return(ur$root)

}



#' Harm boundaries for safety testing
#'
#' Calculates the boundaries at each interim analysis, i.e. the number of events in the intervention group
#' that would lead to a stopping of the trial based on a binomial exact test,
#' assuming that the events should be equally distributed amont both groups.
#' The indicated scenario (and all more extreme) would lead to a rejection of H0 (equal distribution) and a stopping for safety.
#'
#' @param nevents vector with number of events (over both arms) at which an interim analysis is done
#' @param alpha_test the nominal alpha level to use for each test
#' @param pH0 proportion of events in the intervention arm under the null hypothesis,
#'	typically based on randomization ratio (e.g. 0.5 for a 1:1 randomization)
#' @param maxevents optional maximum number of events expected for the trial (over both arms), used to calculate the expected number of events
#' @param pH1 optional alternative, numeric vector, proportion of events in the intervention arm
#' @param rrH1 alternative specification of alternative as risk ratio (intervention / control)
#' @param orH1 alternative specification of alternative as risk ratio (intervention / control). Requires the control proportion (r0).
#' @param rdH1 alternative specification of alternative as risk difference (intervention - control). Requires the control proportion (r0) and the number of participants (n).
#' @param r0 risk in the control group. Required if the alternative is given as risk difference or odds ratio.
#' @param n total number of participants. Required if the alternative is given as risk difference.
#' @return a list with 3 data.frames: bounds, stopprob and opchar.
#' bounds has a row for each interim analysis and columns for
#'	number of events (events),
#'	number of events in control and intervention group that would lead to a stop
#'	(events_intervention, events_control), and the nominal alpha for each test (alpha_test).
#'  stopprob has a row for each interim analysis and columns for
#'	number of events (events),
#'	the hypothesis (pH),
#'	the stopping probability (stop_prob), and
#'	the cumulative stopping probability (cum_stop_prob)
#'	opchar has a row for each hypothesis (null plus each alternative) and columns 
#'	for the assumed proportion of events in the intervention group (p),
#'	the cumulative stopping probabilities (cum_stop_prob) and 
#'	the expected total number of events (expected_events)
#'	for the null and each alternative.
#'
#'
#' @export
#'
#' @importFrom stats dbinom
#'
#' @examples
#' getHarmBound(nevents=c(10,50,100), alpha_test=0.025, pH0=0.5)
#' #adding an alternative
#' getHarmBound(nevents=c(10,50,100), alpha_test=0.025, pH0=0.5, pH1=0.6)
#' #assume that a total of 150 events might occur
#' getHarmBound(nevents=c(10,50,100), alpha_test=0.025, pH0=0.5, pH1=0.6, maxevents=150)
#' #or several alternatives
#' getHarmBound(nevents=c(10,50,100), alpha_test=0.025, pH0=0.5,
#'	pH1 = seq(0.6,0.8,by=0.05), maxevents=150)
#' #or as risk ratio
#' getHarmBound(nevents=c(10,50,100), alpha_test=0.025, pH0=0.5, rrH1=1.5, maxevents=150)


getHarmBound <- function(nevents,alpha_test,pH0,
	maxevents=NULL,
	pH1=NULL, 
	rrH1=NULL, orH1=NULL,rdH1=NULL,
	r0=NULL, n=NULL){

	#check alternative
	nn<-sum(!is.null(pH1) | is.null(rdH1) | is.null(rrH1) | is.null(orH1))
	if (nn>1) {
		stop("Only one of pH1, rdH1, rrH1, or rrH1 should be entered.")
	}
	
	#check maxevents 
	if (!is.null(maxevents)) {
		if (maxevents<max(nevents)) {
			stop("maxevents has to be larger or equal than the maximum of nevents")
		}
	} else {
		maxevents<-max(nevents)
	}
	
	if (!is.null(rdH1) | !is.null(rrH1) | !is.null(orH1)) {

		if (!is.null(rrH1)) {
			#for risk ratio, pH1 does not depend on r0 or n:
			pH1<-rrH1*pH0/(1-pH0)/(1+rrH1*pH0/(1-pH0))
			#cri<-convertRisks(rr=rrH1, r0=0.5, n0=(1-pH0), n1=pH0)
			#stopifnot(abs(pH1 - cri["eprop"])<10^(-10))
		}
		if (!is.null(orH1)) {
			#for odds, pH1 depends on r0
			if (is.null(r0)) {
				stop("r0 has to be given if the effect is given as odds ratio (orH1).")
			}
			if (length(orH1)!=length(r0)) {
				stop("Length of orH1 and r0 have to be the same")
			}
			rr<-orH1/(1-r0+orH1*r0)
			pH1<-rr*pH0/(1-pH0)/(1+rr*pH0/(1-pH0))
			#cri<-convertRisks(or=orH1, r0=r0, n0=(1-pH0), n1=pH0)
			#pH1 <- cri["eprop"]
		}

		if (!is.null(rdH1)) {
			#for risk difference, pH1 depends on r0 and n
			if (is.null(n)) {
				stop("n has to be given if the effect is given as risk difference (rdH1).")
			}
			if (is.null(r0)) {
				stop("r0 has to be given if the effect is given as risk difference (rdH1).")
			}
			n0 <- (1-pH0)*n
			n1 <- pH0*n
			r1<-rdH1 + r0
			pH1<-r1*n1/(r0*n0 + r1*n1)
			#cri<-convertRisks(rd=rdH1, r0=r0, n0=n0, n1=n1)
			#pH1 <- cri["eprop"]
		}
	}

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

			## define 'bound' to be the count corresponding to
			## the 'largest.index'
			## which we have significance at per-test-level 'alphaVal'
			bound <- valSeq[largest.index]

			bounds$treatBound[j] <- bound
			bounds$alphaLevelBound[j] <- upperTailProbs[largest.index]
			bounds$cutoff[j] <- alphaVal
		}

	}
	
	
	#bounds:
	assumpt<-c(pH0,pH1)
	
	stopprob<-vector(length=1+length(pH1),mode="list")
	names(stopprob)<-assumpt
	cumstop<-numeric(0)
	
	i<-1
	for (i in 1:length(assumpt)) {
		
		hyp<-ifelse(i==1,"H0","H1")
					
		if (!all(is.na(bounds$treatBound))) {	
			out <- pNS(Bounds=bounds$treatBound, pH0=assumpt[i])
			nstop<-sum(out$Bounds$n*out$Stop,(1-out$totalStopProb)*maxevents)
			outc<-data.frame(p=assumpt[i],cum_stop_prob=out$totalStopProb,
				expected_events=nstop,hyp=hyp)			
		} else {
			nr <- nrow(bounds)
			out <- list(Bounds= data.frame(n=1:nr,StoppingBound = bounds$treatBound),
				Stop = rep(0, nr),cutoff = bounds$cutoff)
			outc<-data.frame(p=pH0,cum_stop_prob=0,expected_events=maxevents,hyp=hyp)	
		}		
		
		#Boundaries
		if (i==1) {
			boundOut<-out$Bounds
			names(boundOut)[names(boundOut)=="n"] <- "events"
			names(boundOut)[names(boundOut)=="StoppingBound"] <- "events_intervention"
			boundOut$events_control <- boundOut$events - boundOut$events_intervention		
			boundOut <- cbind(boundOut,alpha_test = bounds$cutoff)
			
			stopifnot(is.na(boundOut[!boundOut$events %in% nevents,"events_intervention"]))
			boundOutna<-boundOut[boundOut$events %in% nevents, ]
			rownames(boundOutna)<-1:nrow(boundOutna)
		}
		
		#Stopping probs 
		sprobi <- data.frame(events=boundOut$events,
			pH = assumpt[i],
			hyp = hyp,
			stop_prob = out$Stop,
			cum_stop_prob=cumsum(out$Stop))
		
		stopifnot(sprobi[!sprobi$events %in% nevents,"stop_prob"]==0)
		sprobi<-sprobi[sprobi$events %in% nevents, ]
		rownames(sprobi)<-1:nrow(sprobi)
		
		stopprob[[i]]<-sprobi
		
		#summary
		if (!is.null(rrH1)) {
			outc<-cbind(outc,rr=c(1,rrH1)[i])
			outc<-outc |> dplyr::relocate(rr, .after = .data$p)
		}
		if (!is.null(orH1)) {
			outc<-cbind(outc,or=c(1,orH1)[i],r0 = .data$r0)
			outc<-outc |> 
				dplyr::relocate(.data$or, .after = .data$p) |>
				dplyr::relocate(.data$r0, .after = .data$or)
		}
		if (!is.null(rdH1)) {
			outc<-cbind(outc,rd=c(1,rdH1)[i],r0=r0,n=n)
			outc<-outc |> 
				dplyr::relocate(.data$rd, .after = .data$p) |>
				dplyr::relocate(.data$r0, .after = .data$rd) |>
				dplyr::relocate(.data$n, .after = .data$r0)
		}
			
		cumstop<-rbind(cumstop,outc)
	
	}
	
	#combine
		
	res<-list(boundOutna,stopprob,cumstop)
	names(res)<-c("bounds","stopprob","opchar")

	class(res) <- c("harmbound", class(res))

	return(res)
}



#' Title
#' Helper function which creates an object containing the values of the function P(n,s)
#' defined by Breslow (1970, JASA) as, for 0 <= s <= n <= N,  the probability
#' of the binomial random walk S_n reaching S_n = s without "absorption" into
#' the rejection region (i.e. without hitting any of the stopping bounds).
#'
#' @param Bounds Vector of stopping bounds after each event.
#'	For event totals where stopping is not permitted the 'Bound' should be set to NA.
#' @param pH0 proportion of events in the intervention arm under the null hypothesis,
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


#' Convert the proportion of events in the intervention groups to risk differences and ratios and vice versa.
#'
#' @param eprop proportion of events in intervention group
#' @param etotal total number of events
#' @param rd risk difference
#' @param rr risk ratio
#' @param or odds ratio
#' @param r0 risk in the control group
#' @param n0 number of patients in the control group
#' @param n1 number of patients in the intervention group
#'
#' @return vector with risks in control and intervention group (r0, r1),
#'	the risk difference (rd), risk ratio (rr) and odds ratio (or)
#'
#' @export
#'
#' @examples
#' convertRisks(eprop=0.5,etotal=100,n0=200)
#' convertRisks(eprop=0.6,etotal=100,n0=200)
#' convertRisks(rr=1.5,n0=200,r0=0.2)
#'
convertRisks<-function(eprop=NULL,etotal=NULL,
	rd=NULL,rr=NULL,or=NULL,
	r0=NULL,
	n0,n1=n0) {

	nn<-sum(!is.null(eprop) | !is.null(rd) | !is.null(rr) | !is.null(or))
	if (nn!=1) {
		stop("Only one of eprop, rd, rr or or should be given")
	}

	if (!is.null(eprop)) {

		if (is.null(etotal)) {
			stop("etotal has to be given for conversion to risk difference or ratios.")
		}

		r0<-(etotal-eprop*etotal)/n0
		r1<-eprop*etotal/n1
		rd<-r1-r0
		rr<-r1/r0
		or<-r1/(1-r1)/(r0/(1-r0))
	}

	if (!is.null(rd) | !is.null(rr) | !is.null(or)) {

		if (is.null(r0)) {
			stop("r0 has to be given for conversion to proportion of events in intervention group.")
		}

		if (!is.null(rd)) {
			r1 <- rd + r0
		}
		if (!is.null(rr)) {
			r1 <- rr * r0
		}
		if (!is.null(or)) {
			r1 <- (or * r0/(1-r0)) / (1 + or * r0/(1-r0))
		}

		eprop <- r1*n1 / (r0*n0 + r1*n1)

		rd<-r1-r0
		rr<-r1/r0
		or<-r1/(1-r1)/(r0/(1-r0))
		etotal<-r1*n1/eprop
	}

	res<-c(eprop,etotal,n0,n1,r0,r1,rd,rr,or)
	names(res)<-c("eprop","etotal","n0","n1","r0","r1","rd","rr","or")
	return(res)
}
