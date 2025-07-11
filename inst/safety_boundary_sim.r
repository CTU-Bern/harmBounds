###############################################
#Simulation of safety analysis based on binomial threshold
####################################

library(harmBounds)

#function to find boundary based on a test-specific one-sided alpha
#----------

findbound(30,alpha_test=0.025,alternative="greater")
binom.test(x=21,n=30,alternative="greater")
#gives the number of events, which would lead to a significant binomial exact test 
#	at the specified threshold
findbound(30,alpha_test=0.025,alternative="less")
binom.test(x=9,n=30,alternative="less")


#single simulation
#--------------

set.seed(123)
sim_safety_stop(nevents=seq(10,100,by=10),pH0=0.5,pH1=0.5)

set.seed(123)
sim_safety_stop(nevents=seq(10,100,by=10),pH0=0.5,pH1=0.8)


#repeated simulations
#------------

#under H0
set.seed(123)
ssp0<-vapply(1:1000,
	function(x) sim_safety_stop(nevents=seq(10,100,by=10),pH0=0.5,pH1=0.5,
		alpha_test=0.025)$nstop,
	numeric(1))
mean(ssp0>0)

#or at an alpha that would control overall type I error
apt <- getAlphaPerTest(nevents=seq(10,100,by=10),pH0=0.5,totalAlpha=0.05)
set.seed(123)
ssp0<-vapply(1:1000,
	function(x) sim_safety_stop(nevents=seq(10,100,by=10),pH0=0.5,pH1=0.5,
		alpha_test=apt)$nstop,
	numeric(1))
mean(ssp0>0)

#under H1, e.g. with pH1=0.7
set.seed(123)
ssp1<-vapply(1:1000,
	function(x) sim_safety_stop(nevents=seq(10,100,by=10),pH0=0.5,pH1=0.7,
		alpha_test=apt)$nstop,
	numeric(1))
mean(ssp1>0)



#different conditions
#---------------------

nevents <- seq(10, 100, by=10)
pH0 <- 0.5
totalAlpha <- 0.05

#one possibility it to control type I error
alphaPerTest <- getAlphaPerTest(nevents=nevents, pH0=pH0,totalAlpha = totalAlpha)
alphaPerTest

alist<-c(alphaPerTest,0.025,0.05)
nreps<-100
pgs<-c(0.1,0.2,seq(0.3,0.7,by=0.025),0.8,0.9)

grid<-expand.grid(alpha=alist,pH1=pgs)
grid$nstop<-grid$tstop<-grid$pstop<-NA
j<-1

set.seed(123)
for (j in 1:nrow(grid)) {
	
	ssp1<-vapply(1:nreps,function(x) {
		ss<-sim_safety_stop(nevents=nevents,pH0=0.5,pH1=grid$pH1[j],alpha_test=grid$alpha[j])
			return(c(ss$nstop,ss$tstop))},
		numeric(2))
	rj<-c(apply(ssp1,1,function(x) mean(x,na.rm=TRUE)),mean(ssp1[1,]>0))	
	grid[j,c("nstop","tstop","pstop")]<-rj
}


#plot
#---------

res<-grid
res$pstop<-100*res$pstop
res$pH1<-100*res$pH1
res$alpha<-round(res$alpha,3)
res$alpha<-as.factor(res$alpha)

library(ggplot2)
ggplot(data=res, aes(x=pH1, y=pstop)) +
    geom_line(aes(colour=alpha)) +
	#geom_smooth(aes(colour=alpha)) +
    xlim(0,100) + 
    scale_y_continuous(breaks = seq(0, 100, by=20)) + 
    labs(x="Events in experimental group (%)", y = "Stopped for safety (%)", 
         color="Alpha") + 
	theme_bw()	 







