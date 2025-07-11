#' Find stopping boundary via binomial exact tests
#'
#' @param n total number of events
#' @param alpha_test nominal alpha for the binomial test
#' @param pH0 proportion of events in the experimental arm under the null hypothesis,
#'	typically based on randomization ratio (e.g. 0.5 for a 1:1 randomization)
#' @param alternative direction of alternative, "less" or "greater"
#'
#' @return number of events in the experimental group that would lead to a stopping
#'
#' @export
#'
#' @importFrom stats binom.test
#'
#' @examples
#'	findbound(n=20, alpha_test=0.025, pH0 = 0.5, alternative="greater")
#'	findbound(n=20, alpha_test=0.025, pH0 = 0.5, alternative="less")
#'	
#'
findbound<-function(n, alpha_test=0.025, pH0 = 0.5, alternative="greater") {
 
  pvs<-vapply(1:n,function(x) 
	binom.test(x=x, n=n, p = pH0, alternative=alternative)$p.value,
	numeric(1))
	
  xlim<-which(pvs<alpha_test)

  if (!is.null(xlim)) {
    if (alternative=="less") {
      xlim<-max(xlim)
    }
    if (alternative=="greater") {
      xlim<-min(xlim)
    }
  } else {
    xlim<-NA
  }
  return(xlim)
}
