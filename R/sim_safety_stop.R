



#' Simulate safety stopping values
#'
#' @param nevents number of events at which to perform an interim analysis (vector)
#' @param pgroup proportion of events in the experimental arm
#' @param cutoff alpha level for binomial exact test
#'
#' @return list with a dataframe with number of events in each group plus upper limit for stopping and indicator for whether stopped, plus indicators number of stops and time points at first stop
#' @export
#'
#' @examples
sim_safety_stop <- function(nevents,pgroup=0.5,cutoff=0.025) {

  group<-rbinom(max(nevents),1,pgroup)

  n1<-sapply(nevents,function(x) sum(group[1:x]))
  n0<-nevents-n1

  dat<-data.frame(nevents,n1,n0)

  ulim<-unlist(lapply(dat$nevents,function(x) findbound(x,cutoff=cutoff,alternative="greater")))

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
