# Task 4: Stochastic Age-structured Population Model
# Written by Jeremy Collie on Bastille Day 2023
# Modified on 3-Aug-23 to calculate rebuilding probabilities
# Modified on 5-Aug-23 to retrospectively adjust numbers at age
# Modified on 8-Aug-23 to randomize recruitment in the first year
setwd("C:/Users/Jerem/Documents/Research/NOAA_COCA/Stock Assessments/WhiteHake/Analyses")
rm(list=ls())
library(scales)
#####################################################################
# Step 1. Gather the life-history parameters for white hake
#####################################################################
load("WHT_HAKE_NE.RData")
attach(WHT_HAKE_NE)
SR <- SR              # data frame of stock and recruitment
n_years <- n_years
year1 <- year1
ages <- ages          # note that the terminal age is a plus group
select <- as.matrix(selectivity)
M <- as.matrix(M)
FS <- fracyr_spawn
maturity <- as.matrix(maturity)
WAA_SSB <- as.matrix(WAA_SSB)
WAA_catch <- as.matrix(WAA_catch)
refpt <- refpt
#bev_par <- bev_par
FCB <- FCB            # F, Catch, SSB, plus R by calendar year
detach(WHT_HAKE_NE)

# calculate the derived parameters
avec <- SR$a          # vector of a values from time-varying model
se.a <- SR$se.a       # standard error of a values
Rb <- SR$b[1]         # note that b is negative
Y <- length(avec)     # number of stock-recruitment pairs
A <- length(ages)     # number of ages
year <- FCB$year

# take averages of four most recent years for select and weights
S <-  colMeans(select[(n_years-3):(n_years), ])  
NM <- colMeans(M[(n_years-3):(n_years), ])
MA <- colMeans(maturity[(n_years-3):(n_years), ])
WS <- colMeans(WAA_SSB[(n_years-3):(n_years), ])
WC <- colMeans(WAA_catch[(n_years-3):(n_years), ])

#####################################################################
# Step 2. Specify the initial conditions
#####################################################################
source("Read.ASAP3.rep.file.R")
whk.rep <- read.ASAP3.rep.file(rep.file="Run14.rep")
Yt <- 11                        # number of years to simulate (11)
Nt <- matrix(nrow=Yt, ncol=A)   # numbers at age matrix
# specify numbers at age in the terminal year
Nt[1, ] <- as.numeric(whk.rep$Na[n_years, ])
# retrospective adjustment
retro.adjust <- c(0.797791713,0.900195342,1.030630334,0.995292268,
      0.912350489,0.854540172,0.824898744,0.806042092,0.77856153)
Nt[1, ] <- Nt[1 ,]*retro.adjust
R <- whk.rep$Na[ ,1]
# SSB and Catch are in FCB
SSBt <- rep(NA,Yt)         # SSB is also retrospectively adjusted
F.adjust <- 0.108          # terminal F with retrospective adjustment
Z <- F.adjust*S + NM
SSBt[1] <- sum(Nt[1,]*exp(-Z*FS)*MA*WS)
Ct <- rep(NA,Yt)
Ct[1] <- FCB$Catch[n_years]

# a. construct the cumulative frequency distribution for 1963-2021
#jpeg('rec_dist_Fig10.jpg', width = 3500, height = 2500, res = 2500 / 5)
lnR <- log(R)
# see the test code for this method at the end of the script
lnR.sort <- sort(lnR)
long.stand <- 0:(n_years-1)/(n_years-1)
plot(long.stand,lnR.sort,col="red",xlab="Density",ylab="log(R)",bty="l")
# fit a lowess smoother to predict R from the percentile
long.smooth <- loess(lnR.sort~long.stand,span=0.2)
lines(long.smooth,col="red")

# b. construct the cumulative frequency distribution for 1995-2021
lnR <- lnR[year>1994]
# see the test code for this method at the end of the script
lnR.sort <- sort(lnR)
short.stand <- 0:(length(lnR)-1)/(length(lnR)-1)
points(short.stand,lnR.sort,col="orange")
# fit a lowess smoother to predict R from the percentile
short.smooth <- loess(lnR.sort~short.stand,span=0.25)
lines(short.smooth,col="orange")
abline(v=0.5,lty=3)
legend("topleft",c("1963-2021","1995-2021"),lty=1,
       col=c("red","orange"),bty="n")
title("Cumulative Recruitment Distributions")
#dev.off()

# c. time-invariant Ricker model
Ric.ti <- lm(Lnrs~SSB,data=SR)
RP <- coef(Ric.ti)
RV <- summary(Ric.ti)$sigma^2      # RV is same V as in Table 1
rho <- acf(resid(Ric.ti))$acf[2]   # first-order autocorrelation coef
eta <- rep(0,Yt)                   # random deviates for Ricker model
eta.0 <- resid(Ric.ti)[Y]          # start with last observed residual

# d. time-varying Ricker model
sigmaV <- sqrt(0.0344)
sigmaW <- sqrt(0.0258)
Ra <- rep(0,Yt)
Ra.0 <- avec[Y]    # Rb was specified above
# specify the parameters of the bounded random walk
# the code comes from LBRW.R
tau <- -0.15				# mean of the random walk variable
max <- 0.5 				  # maximum of the logistic function (0.5)
alpha <- 1.2    		# range of the function (0.75)
slope <- 10				  # slope of the logistic function (10)
x <- seq(-2,2,0.05)
indic <- ifelse(x>tau, -1,1)
y <- indic*max/(1+exp(-slope*(-indic*(x-tau)-alpha)))
#jpeg('bounded_RW_Fig11.jpg', width = 2500, height = 2000, res = 2500 / 5)
plot(x,y,type="l",xlab="Ricker a(t)",ylab="Penalty",bty="l")
abline(h=0,lty=4)
abline(v=tau,lty=2)
title("Logistic Bounds")
#dev.off()

#####################################################################
# Step 3. Create a function for each recruitment scenario
#####################################################################
# specify the recruitment model: elon, esht, rinv, rvar
# randomize recruitment in the first year
# the standard deviation comes from RUN14.STD line 128
random.rec <- function(Nrec=R[n_years]) {
  # Nrec is recruitment in the terminal year of the assessment
  N1 <- exp(rnorm(1,mean=log(Nrec),sd=0.24009))
  # apply the retrospective adjustment
  N1 <- N1*0.797791713
  return(N1)
}

# a. cumulative recruitment distribution for 1963-2021 
elon <- function(Ffull=0) {
Z <- Ffull*S + NM
for (i in 2:Yt) {
  Nt[i,1] <- exp(predict(long.smooth,runif(1)))
  for (j in 2:(A-1)) {Nt[i,j] <- Nt[(i-1),(j-1)]*exp(-Z[(j-1)])}
  Nt[i,A] <- Nt[(i-1),(A-1)]*exp(-Z[(A-1)]) + Nt[(i-1),A]*exp(-Z[A])
  SSBt[i] <- sum(Nt[i,]*exp(-Z*FS)*MA*WS)
  Ct[i] <- sum(Nt[i,]*(1-exp(-Z))*WC*Ffull*S/Z)
}
rbind(Nt[ ,1],SSBt,Ct)
}

# b. cumulative recruitment distribution for 1995-2021 
esht <- function(Ffull=0) {
Z <- Ffull*S + NM
for (i in 2:Yt) {
  Nt[i,1] <- exp(predict(short.smooth,runif(1)))
  for (j in 2:(A-1)) {Nt[i,j] <- Nt[(i-1),(j-1)]*exp(-Z[(j-1)])}
  Nt[i,A] <- Nt[(i-1),(A-1)]*exp(-Z[(A-1)]) + Nt[(i-1),A]*exp(-Z[A])
  SSBt[i] <- sum(Nt[i,]*exp(-Z*FS)*MA*WS)
  Ct[i] <- sum(Nt[i,]*(1-exp(-Z))*WC*Ffull*S/Z)
}
rbind(Nt[ ,1],SSBt,Ct)
}

# c. time-invariant Ricker model
rinv <- function(Ffull=0) {
Z <- Ffull*S + NM
eta[1] <- rho*eta.0+rnorm(n=1,mean=0,sd=sqrt(RV*(1-rho^2)))
for (i in 2:Yt) {
  eta[i] <- rho*eta[i-1]+rnorm(n=1,mean=0,sd=sqrt(RV*(1-rho^2)))
  Nt[i,1] <- SSBt[i-1]*exp(RP[1]+RP[2]*SSBt[i-1]+eta[i])
  for (j in 2:(A-1)) {Nt[i,j] <- Nt[(i-1),(j-1)]*exp(-Z[(j-1)])}
  Nt[i,A] <- Nt[(i-1),(A-1)]*exp(-Z[(A-1)]) + Nt[(i-1),A]*exp(-Z[A])
  SSBt[i] <- sum(Nt[i,]*exp(-Z*FS)*MA*WS)
  Ct[i] <- sum(Nt[i,]*(1-exp(-Z))*WC*Ffull*S/Z)
}
rbind(Nt[ ,1],SSBt,Ct)
}

# d. time-varying Ricker model
rvar <- function(Ffull=0) {
Z <- Ffull*S + NM
v <- rnorm(n=Yt,sd=sigmaV)
w <- rnorm(n=Yt,sd=sigmaW)
Ra[1] <- Ra.0 + w[1]
for (i in 2:Yt) {
  Ra[i] <- Ra[i-1] + w[i]
  Nt[i,1] <- SSBt[i-1]*exp(Ra[i]+Rb*SSBt[i-1]+v[i])
  for (j in 2:(A-1)) {Nt[i,j] <- Nt[(i-1),(j-1)]*exp(-Z[(j-1)])}
  Nt[i,A] <- Nt[(i-1),(A-1)]*exp(-Z[(A-1)]) + Nt[(i-1),A]*exp(-Z[A])
  SSBt[i] <- sum(Nt[i,]*exp(-Z*FS)*MA*WS)
  Ct[i] <- sum(Nt[i,]*(1-exp(-Z))*WC*Ffull*S/Z)
}
rbind(Nt[ ,1],SSBt,Ct)
}

#####################################################################
# Step 4. Run the random simulations for each scenario separately
#####################################################################
FM <- 0.213            # Fully recruited fishing mortality (0.1605)
nsim <- 10000           # number of random simulations
Btarg <- 28.19         # current rebuilding target
# a. cumulative recruitment distribution for 1963-2021 
# prime the pump
Nt[1,1] <- random.rec()
toto <- elon(Ffull=FM)      
R.tab <- toto[1,]
SSB.tab <- toto[2,]
C.tab <- toto[3,]
# now run more replicates
for (i in 1:nsim) {
  Nt[1,1] <- random.rec()
  toto <- elon(Ffull=FM)
  R.tab <- rbind(R.tab, toto[1,])
  SSB.tab <- rbind(SSB.tab, toto[2,])
  C.tab <- rbind(C.tab, toto[3,])
}
# summarize the results
roro <- apply(X=R.tab,MARGIN=2,FUN=quantile,probs=c(0.025,0.5,0.975))
elon.R <- t(roro)
soso <- apply(X=SSB.tab,MARGIN=2,FUN=quantile,probs=c(0.025,0.5,0.975))
elon.SSB <- t(soso)
coco <- apply(X=C.tab,MARGIN=2,FUN=quantile,probs=c(0.025,0.5,0.975))
elon.C <- t(coco)
# calculate the probability of meeting the rebuilding target
SSB.sort <- sort(SSB.tab[ ,Yt])
p.elon <- 1- min(which(SSB.sort>Btarg))/nsim

# c. time-invariant Ricker model
# prime the pump
Nt[1,1] <- random.rec()
toto <- rinv(Ffull=FM)
R.tab <- toto[1,]
SSB.tab <- toto[2,]
C.tab <- toto[3,]
# now run more replicates
for (i in 1:nsim) {
  Nt[1,1] <- random.rec()
  toto <- rinv(Ffull=FM)
  R.tab <- rbind(R.tab, toto[1,])
  SSB.tab <- rbind(SSB.tab, toto[2,])
  C.tab <- rbind(C.tab, toto[3,])
}
# summarize the results
roro <- apply(X=R.tab,MARGIN=2,FUN=quantile,probs=c(0.025,0.5,0.975))
rinv.R <- t(roro)
soso <- apply(X=SSB.tab,MARGIN=2,FUN=quantile,probs=c(0.025,0.5,0.975))
rinv.SSB <- t(soso)
coco <- apply(X=C.tab,MARGIN=2,FUN=quantile,probs=c(0.025,0.5,0.975))
rinv.C <- t(coco)
# calculate the probability of meeting the rebuilding target
SSB.sort <- sort(SSB.tab[ ,Yt])
p.rinv <- 1- min(which(SSB.sort>Btarg))/nsim

# b. cumulative recruitment distribution for 1995-2021
# prime the pump
Nt[1,1] <- random.rec()
toto <- esht(Ffull=FM)
R.tab <- toto[1,]
SSB.tab <- toto[2,]
C.tab <- toto[3,]
# now run more replicates
for (i in 1:nsim) {
  Nt[1,1] <- random.rec()
  toto <- esht(Ffull=FM)
  R.tab <- rbind(R.tab, toto[1,])
  SSB.tab <- rbind(SSB.tab, toto[2,])
  C.tab <- rbind(C.tab, toto[3,])
}
# summarize the results
roro <- apply(X=R.tab,MARGIN=2,FUN=quantile,probs=c(0.025,0.5,0.975))
esht.R <- t(roro)
soso <- apply(X=SSB.tab,MARGIN=2,FUN=quantile,probs=c(0.025,0.5,0.975))
esht.SSB <- t(soso)
coco <- apply(X=C.tab,MARGIN=2,FUN=quantile,probs=c(0.025,0.5,0.975))
esht.C <- t(coco)
# calculate the probability of meeting the rebuilding target
SSB.sort <- sort(SSB.tab[ ,Yt])
p.esht <- 1 - min(which(SSB.sort>Btarg))/nsim

# d. time-varying Ricker model
# prime the pump
Nt[1,1] <- random.rec()
toto <- rvar(Ffull=FM)
R.tab <- toto[1,]
SSB.tab <- toto[2,]
C.tab <- toto[3,]
# now run more replicates
for (i in 1:nsim) {
  Nt[1,1] <- random.rec()
  toto <- rvar(Ffull=FM)
  R.tab <- rbind(R.tab, toto[1,])
  SSB.tab <- rbind(SSB.tab, toto[2,])
  C.tab <- rbind(C.tab, toto[3,])
}
# summarize the results
roro <- apply(X=R.tab,MARGIN=2,FUN=quantile,probs=c(0.025,0.5,0.975))
rvar.R <- t(roro)
soso <- apply(X=SSB.tab,MARGIN=2,FUN=quantile,probs=c(0.025,0.5,0.975))
rvar.SSB <- t(soso)
coco <- apply(X=C.tab,MARGIN=2,FUN=quantile,probs=c(0.025,0.5,0.975))
rvar.C <- t(coco)
# calculate the probability of meeting the rebuilding target
SSB.sort <- sort(SSB.tab[ ,Yt])
p.rvar <- 1- min(which(SSB.sort>Btarg))/nsim

#####################################################################
# Step 5. Plot the results
#####################################################################
yr <- year[n_years] -1 + 1:Yt
x <- c(yr,rev(yr),yr[1])
#jpeg('MSY_Fig12.jpg', width = 2500, height = 3500, res = 2500 / 5)
#par(mfrow=c(3,1),mar=c(2,4,3,1)+0.1)
# Recruitment distributions
matplot(yr,elon.R,type="n",ylim=c(0,12),bty="l",
        xlab="Year",ylab="R (millions)")
y <- c(elon.R[,1],rev(elon.R[,3]),elon.R[1,1])
polygon(x,y,col=alpha("red",alpha=0.2),border=NA)
y <- c(rinv.R[,1],rev(rinv.R[,3]),rinv.R[1,1])
polygon(x,y,col=alpha("black",alpha=0.2),border=NA)
y <- c(esht.R[,1],rev(esht.R[,3]),esht.R[1,1])
polygon(x,y,col=alpha("orange",alpha=0.2),border=NA)
y <- c(rvar.R[,1],rev(rvar.R[,3]),rvar.R[1,1])
polygon(x,y,col=alpha("blue",alpha=0.2),border=NA)

# overplot the medians
lines(yr,elon.R[,2],lwd=2,col="red")
lines(yr,rinv.R[,2],lwd=2,col="black")
lines(yr,esht.R[,2],lwd=2,col="orange")
lines(yr,rvar.R[,2],lwd=2,col="blue")
title(main="a. Recruitment Projections", adj=0)
legend("topleft",c("stationary","recent","Ricker","dlm"),lwd=2,
       col=c("red","orange","black","blue"),bty="n")

# Spawning Stock Biomass distributions
#jpeg('Rebuild_Fig13.jpg', width = 3000, height = 2500, res = 2500 / 5)
matplot(yr,elon.SSB,type="n",ylim=c(0,max(rinv.SSB[,3])),bty="l",
        lab=c(10,5,7),xlab="Year",ylab="SSB (kt)")
y <- c(elon.SSB[,1],rev(elon.SSB[,3]),elon.SSB[1,1])
polygon(x,y,col=alpha("red",alpha=0.2),border=NA)
y <- c(rinv.SSB[,1],rev(rinv.SSB[,3]),rinv.SSB[1,1])
polygon(x,y,col=alpha("black",alpha=0.2),border=NA)
y <- c(esht.SSB[,1],rev(esht.SSB[,3]),esht.SSB[1,1])
polygon(x,y,col=alpha("orange",alpha=0.2),border=NA)
y <- c(rvar.SSB[,1],rev(rvar.SSB[,3]),rvar.SSB[1,1])
polygon(x,y,col=alpha("blue",alpha=0.2),border=NA)

# overplot the medians
lines(yr,elon.SSB[,2],lwd=2,col="red")
lines(yr,rinv.SSB[,2],lwd=2,col="black")
lines(yr,esht.SSB[,2],lwd=2,col="orange")
lines(yr,rvar.SSB[,2],lwd=2,col="blue")
abline(h=Btarg,lty=2,lwd=2)
legend("bottomleft",c("stationary","recent","Ricker","dlm"),lwd=2,
       col=c("red","orange","black","blue"),bty="n",ncol=2)
title(main="b. Spawning Stock Biomass Projections", adj=0)
#dev.off()

# Catch distributions
matplot(yr,elon.C,type="n",ylim=c(0,max(rinv.C[,3])),bty="l",
        xlab="Year",ylab="Catch (kt)")
y <- c(elon.C[,1],rev(elon.C[,3]),elon.C[1,1])
polygon(x,y,col=alpha("red",alpha=0.2),border=NA)
y <- c(rinv.C[,1],rev(rinv.C[,3]),rinv.C[1,1])
polygon(x,y,col=alpha("black",alpha=0.2),border=NA)
y <- c(esht.C[,1],rev(esht.C[,3]),esht.C[1,1])
polygon(x,y,col=alpha("orange",alpha=0.2),border=NA)
y <- c(rvar.C[,1],rev(rvar.C[,3]),rvar.C[1,1])
polygon(x,y,col=alpha("blue",alpha=0.2),border=NA)

# overplot the medians
lines(yr,elon.C[,2],lwd=2,col="red")
lines(yr,rinv.C[,2],lwd=2,col="black")
lines(yr,esht.C[,2],lwd=2,col="orange")
lines(yr,rvar.C[,2],lwd=2,col="blue")
legend("topleft",c("stationary","recent","Ricker","dlm"),lwd=2,
       col=c("red","orange","black","blue"),bty="n")
title(main="c. Catch Projections",adj=0)
#dev.off()

# Create a decision table from 10-year averages of R, C, and SSB
mod.col <- c("Emp.lon","Emp.sht","Ric.inv","Ric.var")
R.col <- c(mean(elon.R[-1,2]),mean(esht.R[-1,2]),mean(rinv.R[-1,2]),mean(rvar.R[-1,2]))
C.col <- c(mean(elon.C[-1,2]),mean(esht.C[-1,2]),mean(rinv.C[-1,2]),mean(rvar.C[-1,2]))
SSB.col <- c(mean(elon.SSB[-1,2]),mean(esht.SSB[-1,2]),mean(rinv.SSB[-1,2]),mean(rvar.SSB[-1,2]))
p.col <- c(p.elon,p.esht,p.rinv,p.rvar)
decision.tab <- round(cbind(R.col,C.col,SSB.col,p.col),2)
dimnames(decision.tab)[[1]] <- mod.col
dimnames(decision.tab)[[2]] <- c("R (millions)","C (kt)","SSB (kt)","p-rebuild")
decision.tab
#####################################################################
# APPENDIX. Randomization tests
#####################################################################
# a/b. sampling from the empirical recruitment distribution
lnR.samp <- sample(lnR,replace=TRUE)
# or sample from the cumulative frequency distribution
lnR <- log(R)
lnR.sort <- sort(lnR)
lnR.stand <- 0:(n_years-1)/n_years
# should the cumsum be calculated on the sorted lnR? Si
plot(lnR.stand,lnR.sort)
# fit a lowess smoother to predict R from the percentile
lnR.smooth <- loess(lnR.sort~lnR.stand,span=0.2)
lines(lnR.smooth)
predict(lnR.smooth,0.5)
abline(h=predict(lnR.smooth,0.5),lty=2)
abline(v=0.5,lty=3)

uni <- runif(1000)
points(uni,predict(lnR.smooth,uni),pch=4,col=6)
# seems to work sufficiently well

# c. simulate 1000 correlated random deviates
n <- 1000
eta <- rep(0,n)
for (i in 2:n) {
  eta[i] <- rho*eta[i-1]+rnorm(n=1,mean=0,sd=sqrt(RV*(1-rho^2)))
}
plot(eta,type="l")
abline(h=0,lty=2) 
summary(eta)                       # eta should have mean 0
var(eta)                           # this value should equal RV
acf(eta,lag.max=10)$acf[2]         # this value should equal rho


