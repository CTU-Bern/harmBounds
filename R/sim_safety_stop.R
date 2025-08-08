
#' Simulate safety stopping values.
#'	The effect can be given as proportion of events in the experimental group (pH1),
#'	the risk difference (rdH1), risk ratio (rrH1) or odds ratio (orH1).
#'
#' @param nevents vector with number of events at which an interim analysis is done
#' @param pH0 proportion of events in the experimental arm under the null hypothesis,
#'	typically based on randomization ratio (e.g. 0.5 for a 1:1 randomization)
#' @param alpha_test nominal alpha level for binomial exact test
#' @param pH1 proportion of events in the experimental arm under the alternative hypothesis
#' @param rrH1 risk ratio (experimental / control).
#' @param orH1 risk ratio (experimental / control). Requires the control proportion (r0).
#' @param rdH1 risk difference (experimental - control). Requires the control proportion (r0) and the number of participants (n).
#' @param r0 risk in the control group. Required if the effect is given as risk difference or odds ratio.
#' @param n total number of participants. Required if the effect is given as risk difference.
#'
#' @return list with a dataframe with number of events in each group plus upper limit for stopping and indicator for whether stopped, plus indicators number of stops and time points at first stop
#'
#' @export
#'
#' @importFrom stats rbinom
#'
#' @examples
#'	set.seed(1)	
#'	simSafetyStop(nevents=seq(10,100,by=10),pH0 = 0.5, pH1 = 0.6,alpha_test=0.025)
#'	
#'	set.seed(1)	
#'	simSafetyStop(nevents=seq(10,100,by=10),pH0 = 0.5, rrH1 = 0.6/(1-0.6), alpha_test=0.025)
#'
#'
simSafetyStop <- function(nevents, 
	pH0=0.5,
	alpha_test=0.025,
	pH1=NULL, rrH1=NULL, orH1=NULL,rdH1=NULL,
	r0=NULL, n=NULL) {
	
	nn<-sum(!is.null(pH1) | is.null(rdH1) | is.null(rrH1) | is.null(orH1))
	if (nn!=1) {
		stop("Only one of pH1, rdH1, rrH1, or rrH1 should be entered.")
	}
	
	if (!is.null(rdH1) | !is.null(rrH1) | !is.null(orH1)) {
		
		if (!is.null(rrH1)) {	
			#for risk ratio, pH1 does not depend on r0 or n:
			pH1<-rrH1*pH0/(1-pH0)/(1+rrH1*pH0/(1-pH0))
			cri<-convertRisks(rr=rrH1, r0=0.5, n0=(1-pH0), n1=pH0)
			stopifnot(abs(pH1 - cri["eprop"])<10^(-10))	
		
		} 
		if (!is.null(orH1)) {
			#for odds, pH1 depends on r0  
			if (is.null(r0)) {
				stop("r0 has to be given if the effect is given as odds ratio (orH1).")
			}
			rr<-orH1/(1-r0+orH1*r0)
			pH1<-rr*pH0/(1-pH0)/(1+rr*pH0/(1-pH0))
			cri<-convertRisks(or=orH1, r0=r0, n0=(1-pH0), n1=pH0)
			pH1 <- cri["eprop"]	
		}	
		
		if (!is.null(rdH1)) {
			#for risk difference, pH1 depends on r0 and n 
			if (is.null(n)) {
				stop("n has to be given if the effect is given as risk difference (rdH1).")
			}
			if (is.null(r0)) {
				stop("r0 has to be given if the effect is given as risk difference (rdH1).")
			}
			n0 = pH0*n
			n1 = (1-pH0)*n
			
			cri<-convertRisks(rd=rdH1, rr=rrH1, r0=r0, n0=n0, n1=n1)
			pH1 <- cri["eprop"]	
		}	
	}
	
	group<-rbinom(max(nevents),1,pH1)
	n1<-cumsum(group)[nevents]
	n0<-nevents-n1
	
	dat<-data.frame(nevents,n1,n0)
	
	ulim<-vapply(dat$nevents,function(x) 
	findbound(x, alpha_test=alpha_test, alternative="greater", pH0=pH0),
	numeric(1))
	
	dat<-cbind(dat,ulim)
	dat$out<-with(dat,n1>=ulim)
	
	nstop<-sum(dat$out==TRUE & !is.na(dat$out))
	if (nstop>0) {
		tstop<-min(which(dat$out==TRUE & !is.na(dat$out)))
	} else {
		tstop<-NA
	}
	
	return(list(tests=dat,nstop=nstop,tstop=tstop))
}

