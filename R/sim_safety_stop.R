
#' Simulate safety stopping values
#'
#' @param nevents  vector with number of events at which an interim analysis is done
#' @param pH0 proportion of events in the experimental arm under the null hypothesis,
#'	typically based on randomization ratio (e.g. 0.5 for a 1:1 randomization)
#' @param pH1 proportion of events in the experimental arm under the alternative hypothesis
#' @param alpha_test nominal alpha level for binomial exact test
#'
#' @return list with a dataframe with number of events in each group plus upper limit for stopping and indicator for whether stopped, plus indicators number of stops and time points at first stop
#'
#' @export
#'
#' @importFrom stats rbinom
#'
#' @examples
#'	set.seed(1)	
#'	sim_safety_stop(nevents=seq(10,100,by=10),pH0 = 0.5, pH1 = 0.7,alpha_test=0.025)
#'
#'
sim_safety_stop <- function(nevents, pH0=0.5, pH1=0.5, alpha_test=0.025) {

  group<-rbinom(max(nevents),1,pH1)
  n1<-cumsum(group)[nevents]
  n0<-nevents-n1

  dat<-data.frame(nevents,n1,n0)

  ulim<-vapply(dat$nevents,function(x) 
	findbound(x,alpha_test=alpha_test,alternative="greater"),
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
