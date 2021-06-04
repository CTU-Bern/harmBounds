#' Find boundaries via binomial exact tests
#'
#' @param n
#' @param cutoff alpha for the binomial test
#' @param alternative alternative hypothesis
#'
#' @return
#' @export
#'
#' @examples
findbound<-function(n, cutoff=0.025, alternative="less") {
  xlim<-NULL
  for (j in 1:n) {
    tr<-binom.test(x=j, n=n, p = 0.5, alternative=alternative)
    if (tr$p.value<cutoff) {
      xlim<-c(xlim,j)
    }
  }
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
