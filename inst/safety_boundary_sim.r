###############################################
#Simulation of safety analysis based on binomial threshold
####################################

#function to find boundary based on a test-specific one-sided alpha
#########
library(harmBounds)



findbound(30,cutoff=0.025)
findbound(30,cutoff=0.025,alternative="greater")
findbound(30,cutoff=0.05)
findbound(30,cutoff=0.05,alternative="greater")


#function for simulation
#########

nevents<-1:100
pgroup<-0.5

#nevents: no of events after which a test is done, max(nevents) is to maximal number of events expected to be observed


set.seed(1)
sim_safety_stop(nevents=1:100,pgroup=0.5)

set.seed(1)
sim_safety_stop(nevents=seq(10,100,by=10),pgroup=0.5)


#function for harm bound

#simulation
#########

nreps <- 1000
eval <- data.frame(neventst = rep(100,nreps),
                   neventm = 1,
                   pgroup = rep(0.5, nreps),
                   cutoff = rep(0.025, nreps))
res <- apply(eval,1,function(x) sim_safety_stop(x[1],x[2],x[3])$tests)
res <- dplyr::bind_rows(res)
mean(res[,1]>0)

# mapply(func
#0.116


null.p <- 1/2
harmMonitorAlpha <- 0.05
harmMonitorRange <- c(1, 100)
alphaPerTest <- getAlphaPerTest(harmMonitorRange, null.p,
                                totalAlpha = harmMonitorAlpha)
alphaPerTest

nreps<-1000
eval<-data.frame(nevents=rep(100,nreps),pgroup=rep(0.5,nreps),cutoff=rep(alphaPerTest,nreps))
res<-apply(eval,1,function(x) sim_safety_stop(nevents=x[1],pgroup=x[2],cutoff=x[3])$tests)
res<-dplyr::bind_rows(res)
mean(res$out>0)
#0.063


#H1

nreps<-1000
eval<-data.frame(nevents=rep(100,nreps),pgroup=rep(0.6,nreps),cutoff=rep(0.025,nreps))
res<-apply(eval,1,function(x) sim_safety_stop(nevents=x[1],pgroup=x[2],cutoff=x[3])$tests)
res<-dplyr::bind_rows(res)
mean(res$out>0)



#different conditions
##################

eval <- seq(10, 100, 10)
null.p <- 1/2
harmMonitorAlpha <- 0.05
harmMonitorRange <- c(1, 100)
alphaPerTest <- getAlphaPerTest(harmMonitorRange, null.p,
                                totalAlpha = harmMonitorAlpha)
alphaPerTest

cutoffs<-c(alphaPerTest,0.015,0.025,0.05)
nreps<-1000
pgs<-c(0.1,0.2,0.3,0.4,0.45,0.5,0.55,0.6,0.75,0.8,0.9)

res<-vector(mode = "list", length = length(cutoffs))
names(res)<-cutoffs

set.seed(123)
for (i in 1:length(cutoffs)) {

	co<-cutoffs[i]
	# resi<-data.frame(pgroup=pgs,pstop=NA,tstop_av=NA)
	resi<-data.frame(cutoff = co,
	                 prop_event = pgs,
	                 stop_I = NA,
	                 n_stop_I = NA,
	                 nevents_stop = NA,
	                 interim_stop = NA)

	for (pgroup in resi$prop_event) {
	  df <- data.frame(it = rep(NA, nreps), stop_I = NA, nevents_stop = NA, interim_stop = NA)
	  for(I in 1:nreps){
		  di <- sim_safety_stop(nevents=eval,pgroup=pgroup,cutoff=co)$tests
		  df$it[I] <- I
		  df$stop_I[I] <- any(di$out)
		  df$nevents_stop[I] <- ifelse(sum(di$out) > 0, min(di$nevents[di$out], na.rm = TRUE), NA)
		  df$interim_stop[I] <- ifelse(sum(di$out) > 0, min(which(di$out), na.rm = TRUE), NA)
	  }
		resi[resi$prop_event==pgroup, "stop_I"] <- mean(df$stop_I, na.rm=TRUE)
		resi[resi$prop_event==pgroup, "n_stop_I"] <- sum(df$stop_I, na.rm=TRUE)
		resi[resi$prop_event==pgroup, "nevents_stop"]<-mean(df$nevents_stop, na.rm=TRUE)
		resi[resi$prop_event==pgroup, "interim_stop"]<-mean(df$interim_stop, na.rm=TRUE)
	}
	res[[i]]<-resi
}


save(res,file="safety_boundary_sim.RData")



#plot
########

load(file="safety_boundary_sim.RData")

ci_region<-function(x,lci,uci,cols=rgb(.1,.1,.1,.2)) {
	polygon(c(rev(x),x), c(rev(lci), uci), col=cols, border = NA)
}
colp<-list(rgb(0,0,1,0.1),rgb(1,0,0,0.1))

#dev.off()
#dev.new(width=6.5,height=5,pointsize=12)

png("safety_boundary_sim.png",width=6.5,height=5,pointsize=12,res=600,units="in")

par(mar=c(3.5,3.5,2,1))

plot(0,type="n",axes=FALSE,xlim=c(0,100),ylim=c(0,100),
	xlab="Events in experimental group (%)",ylab="Boundary breached (%)",
	mgp=c(1.8,0,0))
axis(side=1,pos=0)
axis(side=2,las=1,pos=0)

pu<-par("usr")
pu<-c(0,100,0,100)
xs<-seq(pu[1],50,l=10)
ci_region(x=xs,lci=rep(pu[3],length(xs)),uci=rep(pu[4],length(xs)),col=colp[[1]])
xs<-seq(50,pu[2],l=10)
ci_region(x=xs,lci=rep(pu[3],length(xs)),uci=rep(pu[4],length(xs)),col=colp[[2]])

for (i in c(seq(0,40,by=10),seq(60,100,by=10))) {
	lines(y=c(0,100),x=c(i,i),col=rgb(.1,.1,.1,.1))
	#abline(v=i,lty=1,col=rgb(.1,.1,.1,.1))
}

#abline(v=50,lty=2)
lines(x=c(50,50),y=c(0,100),lty=2)
mtext(side=3,line=-0.2,at=50,"H0")
mtext(side=3,line=-0.2,at=c(25,75),c("Benefit","Harm"),col=c("blue","red")) #,adj = c(0,1)


for (i in 1:length(res)) {
	ri<-res[[i]]
	ri[,c("stop_I","prop_event")]<-ri[,c("stop_I","prop_event")]*100
	lines(stop_I~prop_event ,data=ri,col=i)
}

legend(x=0,y=100,legend=paste0("alphaTest = ",c(round(as.numeric(names(res)),4))),bty="n",
	col=1:length(res),lty=1)

dev.off()






