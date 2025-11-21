#' Plot method for harmbound objects produced by \code{getHarmBound}
#'
#' @param x harmbounds objects as generated using the getHarmBound function
#' @param which one of \code{"bounds"}, \code{"abs_stopping"}, \code{"cum_stopping"}, \code{"exp_n"}.
#' @param ... options passed to plot
#'
#' @export
#'
#' @examples
#' harmbound<-getHarmBound(nevents = seq(10, 100, by=10),alpha_test = 0.025,
#'	pH0 = 0.5, pH1 = seq(0.55,0.7,by=0.05), maxevents = 150)
#' plot(harmbound, which = "bounds")
#' plot(harmbound, which = "abs_stopping")
#' plot(harmbound, which = "cum_stopping")
#' plot(harmbound, which = "opchar_stop")
#' plot(harmbound, which = "opchar_n")
#'
plot.harmbound <- function(x, which = "bounds", ...){
	
	which <- match.arg(which, c("bounds", "abs_stopping", "cum_stopping", "opchar_stop", "opchar_n"))
	
	if (which == "bounds"){
		out<-harmboundPlot(x, ...)
	}
	
	if (which == "abs_stopping"){
		out<-absstopPlot(x, ...)
	}
	
	if (which == "cum_stopping"){
		out<-cumstopPlot(x, ...)
	}
	
	if (which == "opchar_stop"){
		out<-opcharStopPlot(x, ...)
	}
	
	if (which == "opchar_n"){
		out<-opcharNPlot(x, ...)
	}
	
	return(out)
}



#' Plot harmbounds and the observed events
#'
#' @param harmbound harmbounds objects as generated using the getHarmBound function
#' @param observed optional observed number of events,
#'	as a vector with the sequential groups in which an event occured
#' (0 for control and 1 for intervention)
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
#' harmbound<-getHarmBound(nevents=seq(10,100,by=10),alpha_test=0.025,pH0=0.5)
#' harmboundPlot(harmbound)
#' set.seed(123)
#' eventgroups<-rbinom(n=100,size=1,prob=0.5)
#' harmboundPlot(harmbound,observed=eventgroups)
#'

harmboundPlot<-function(harmbound,
	observed=NULL,
	colourbound=c("blue","red"),
	fill_alpha=0.5,
	colourobserved="black",
	H0line=TRUE) {

	
	hbs<-harmbound[["bounds"]]
	pH0<-harmbound$stopprob[[1]]$pH[1]
	
	#limit for plot
	ymax<-max(hbs$events_intervention,na.rm=TRUE)
	xmax<-max(hbs$events,na.rm=TRUE)

	if (!is.null(observed)) {
		obs<-data.frame(events_control=cumsum(observed==0),
			events_intervention=cumsum(observed==1))
		obs<-rbind(c(0,0),obs)
		obs$events<-apply(obs,1,sum)
		ymax<-max(obs$events_intervention,ymax,na.rm=TRUE)
		xmax<-max(obs$events,xmax,na.rm=TRUE)
	}

	out<-ggplot() +
		geom_rect(data=hbs,
			mapping=aes(xmin=.data$events-0.5, xmax=.data$events+0.5,
			ymin=0, ymax=.data$events_intervention-0.5),
			color=NA, fill=colourbound[1], alpha=fill_alpha) +
		geom_rect(data=hbs,
			mapping=aes(xmin=.data$events-0.5, xmax=.data$events+0.5,
			ymin=.data$events_intervention-0.5, ymax=.data$events),
			color=NA, fill=colourbound[2], alpha=fill_alpha) +
		scale_y_continuous(limits=c(0,xmax+1)) +
		scale_x_continuous(limits=c(0,xmax+1)) +
		xlab("Total number of events") +
		ylab("Number of events in intervention group")

	if (!is.null(observed)) {
		out<-out +
			geom_step(aes(x = .data$events, y = .data$events_intervention),
				data=obs,
				colour = colourobserved)
	}

	if (H0line) {
		df <- data.frame(x1 = 0, x2 = max(hbs$events),
			y1 = 0, y2 = max(hbs$events)*pH0)
		out<-out +
			geom_segment(aes(x = .data$x1, y = .data$y1, xend = .data$x2, yend = .data$y2),
				 linetype = 2, data = df)
			#geom_abline(intercept = 0, slope = harmbound$pH[1],
			#	size = 0.5, linetype = 2)

	}

	return(out)
}


#' Plot absolute stopping probs
#'
#' @param harmbound harmbounds objects as generated using the getHarmBound function
#'
#' @return barplot with stopping probabilities under H0 and optionally H1
#'
#' @export
#'
#' @import ggplot2
#' @importFrom ggplot2 ggplot geom_bar geom_rect aes
#'
#' @examples
#' harmbound<-getHarmBound(nevents=seq(10,100,by=10),alpha_test=0.025,pH0=0.5,pH1=c(0.6,0.7,0.8))
#' absstopPlot(harmbound)
#'
absstopPlot<-function(harmbound) {

	hbs<-harmbound[["stopprob"]]
	sprob<-do.call(rbind,hbs)
	
	out<-sprob |> 
		dplyr::mutate(na = paste0(.data$hyp,": ",.data$pH)) |>
		ggplot(aes(x=.data$events, y=.data$stop_prob)) + 
		geom_bar(stat = "identity") + 
		ylab("Stopping probability") + 
		xlab("Total number of events") + 
		facet_wrap(~na)
	
	return(out)	
}



#' Plot cumulative stopping probs
#'
#' @param harmbound harmbounds objects as generated using the getHarmBound function
#'
#' @return barplot with cumulative stopping probabilities under H0 and optionally H1
#'
#' @export
#'
#' @import ggplot2
#' @importFrom ggplot2 ggplot geom_bar geom_rect aes
#'
#' @examples
#' harmbound<-getHarmBound(nevents=seq(10,100,by=10),alpha_test=0.025,pH0=0.5,pH1=c(0.6,0.7,0.8))
#' cumstopPlot(harmbound)
#'
cumstopPlot<-function(harmbound) {

	hbs<-harmbound[["stopprob"]]
	sprob<-do.call(rbind,hbs)
	
	out<-sprob |> 
		dplyr::mutate(na = paste0(.data$hyp,": ",.data$pH)) |>
		ggplot(aes(x=.data$events, y=.data$cum_stop_prob, colour=.data$na)) + 
		geom_line() +
		geom_point() + 
		ylab("Cumulative stopping probability") + 
		xlab("Total number of events") + 
		theme(legend.title = element_blank())
		 
	
	return(out)	
}



#' Plot cumulative stopping probability by hypothesis
#'
#' @param harmbound harmbounds objects as generated using the getHarmBound function
#'
#' @return line plot with the cumulative stopping for all alternatives
#'
#' @export
#'
#' @import ggplot2
#' @importFrom ggplot2 ggplot geom_bar geom_rect aes
#'
#' @examples
#' harmbound<-getHarmBound(nevents=seq(10,100,by=10),alpha_test=0.025,
#'	pH0=0.5,pH1=seq(0.5,0.8,by=0.05),maxevents=150)
#' opcharStopPlot(harmbound)
#'
opcharStopPlot<-function(harmbound) {

	hbs<-harmbound[["opchar"]]
	
	if (nrow(hbs)==1) {
		stop("At least on alternative should be given")
	}
	
	out<-hbs |> 
		ggplot(aes(x=.data$p, y=.data$cum_stop_prob)) + 
		geom_line() +
		ylab("Cumulative stopping probability") + 
		xlab("Proportion of events in intervention group") 
	
	if (nrow(hbs)<20) {
		out<-out + geom_point() 
	}
	
	return(out)	
}



#' Plot expected number of events
#'
#' @param harmbound harmbounds objects as generated using the getHarmBound function
#'
#' @return line plot with the expected number of events for all alternatives
#'
#' @export
#'
#' @import ggplot2
#' @importFrom ggplot2 ggplot geom_bar geom_rect aes
#'
#' @examples
#' harmbound<-getHarmBound(nevents=seq(10,100,by=10),alpha_test=0.025,
#'	pH0=0.5,pH1=seq(0.5,0.8,by=0.05),maxevents=150)
#' opcharNPlot(harmbound)
#'
opcharNPlot<-function(harmbound) {

	hbs<-harmbound[["opchar"]]
	
	if (nrow(hbs)==1) {
		stop("At least on alternative should be given")
	}
	
	out<-hbs |> 
		ggplot(aes(x=.data$p, y=.data$expected_events)) + 
		geom_line() +
		ylab("Expected number of events") + 
		xlab("Proportion of events in intervention group") 
		
	if (nrow(hbs)<20) {
		out<-out + geom_point() 
	}
	
	return(out)	
}


