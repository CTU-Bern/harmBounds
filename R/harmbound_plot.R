#' Plot harmbounds and the observed events
#'
#' @param harmbound harmbounds objects as generated using the getHarmBound function
#' @param observed optional observed number of events,
#'	as a vector with the sequential groups in which an event occured
#' (0 for control and 1 for experimental)
#' @param colourbound vector with two colours for the bounds, in and out, default is blue and red
#' @param fill_alpha opacity of the colours for the bounds
#' @param colourobserved colour for the line with the observed events
#' @param H0line logical, whether a line to indicate the expectations should be added.
#'
#' @return plot with the bounds and optionally the observed number of events
#'
#' @export
#'
#' @import ggplot2
#' @importFrom ggplot2 ggplot geom_step geom_rect aes
#'
#' @examples
#' hb<-getHarmBound(nevents=seq(10,100,by=10),alpha_test=0.025,pH0=0.5)
#' harmboundPlot(hb)
#' set.seed(123)
#' eventgroups<-rbinom(n=100,size=1,prob=0.5)
#' harmboundPlot(hb,observed=eventgroups)
#'

harmboundPlot<-function(harmbound,
	observed=NULL,
	colourbound=c("blue","red"),
	fill_alpha=0.5,
	colourobserved="black",
	H0line=TRUE) {

	if (!is.data.frame(harmbound)) {
		harmbound<-harmbound[["bounds"]]
	}

	#limit for plot
	ymax<-max(harmbound$events_exp,na.rm=TRUE)
	xmax<-max(harmbound$events,na.rm=TRUE)

	if (!is.null(observed)) {
		obs<-data.frame(events_control=cumsum(observed==0),
			events_exp=cumsum(observed==1))
		obs<-rbind(c(0,0),obs)
		obs$events<-apply(obs,1,sum)
		ymax<-max(obs$events_exp,ymax,na.rm=TRUE)
		xmax<-max(obs$events,xmax,na.rm=TRUE)
	}

	out<-ggplot() +
		geom_rect(data=harmbound,
			mapping=aes(xmin=.data$events-0.5, xmax=.data$events+0.5,
			ymin=0, ymax=.data$events_exp-0.5),
			color=NA, fill=colourbound[1], alpha=fill_alpha) +
		geom_rect(data=harmbound,
			mapping=aes(xmin=.data$events-0.5, xmax=.data$events+0.5,
			ymin=.data$events_exp-0.5, ymax=.data$events),
			color=NA, fill=colourbound[2], alpha=fill_alpha) +
		scale_y_continuous(limits=c(0,xmax+1)) +
		scale_x_continuous(limits=c(0,xmax+1)) +
		xlab("Total number of events") +
		ylab("Number of events in experimental group")

	if (!is.null(observed)) {
		out<-out +
			geom_step(aes(x = .data$events, y = .data$events_exp),
				data=obs,
				colour = colourobserved)
	}

	if (H0line) {
		df <- data.frame(x1 = 0, x2 = max(harmbound$events),
			y1 = 0, y2 = max(harmbound$events)*harmbound$pH0[1])
		out<-out +
			geom_segment(aes(x = .data$x1, y = .data$y1, xend = .data$x2, yend = .data$y2),
				 linetype = 2, data = df)
			#geom_abline(intercept = 0, slope = harmbound$pH0[1],
			#	size = 0.5, linetype = 2)

	}

	return(out)
}


#' Plot method for harmbound objects produced by \code{getHarmBound}
#'
#' @param x harmbounds objects as generated using the getHarmBound function
#' @param ... further arguments for getHarmBound
#'
#' @export
#'
plot.harmbound <- function(x, ...){
	out<-harmboundPlot(x, ...)
	return(out)
}


