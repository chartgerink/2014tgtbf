# Code written by CHJ Hartgerink
# Checked by: -

# Change the object mypath to where you cloned the repository
mypath <- "C:/Users/Chris/Dropbox/CJM/Masterproject/Analyzing/"

##################################################
# DO NOT EDIT PAST HERE FOR FULL REPRODUCIBILITY #
##################################################
setwd(mypath)

# Load custom functions all at once
customFunct <- list.files('a.Functions/')
for(i in 1:length(customFunct)){
  source(
    paste0('a.Functions/', customFunct[i])
    )
}

# Read- and prepare data
dat <- read.table("1.Pilot study/copilot.txt", stringsAsFactors=F)

# The following reordering of df1 to df2 is due to t and r values essentially
# being F distributions, making for more parsimonious code. dat data was
# collected by running old version of Statcheck, hence manually done here.
# REMOVE FOR FINAL
# REMOVE FOR FINAL
# REMOVE FOR FINAL
sel = dat$Statistic == "t" | dat$Statistic == "r"
dat$df2[sel] = dat$df1[dat$Statistic == "t" | dat$Statistic == "r"]
dat$df1[sel] = NA

# Adding Journal variable to pilot data, as placeholder
# REMOVE FOR FINAL
# REMOVE FOR FINAL
# REMOVE FOR FINAL
dat$Journal <- sample(1:5, 8110, replace=T)

# Eliminate impossible values, make comma's into decimal points
# and compute effect sizes
# truncate negative adjusted effect sizes to 0
dat <- prep.statcheck(dat)

# par(mfrow=c(2,1))
# par(mai=c(.6,1.2,.2,.2))

####### 1
# Descriptives full dataset
# Table
journals <- sort(unique(dat$Journal))
for(j in 1:length(journals)){
  selJournal <- dat$Journal == journals[j]
  len <- length(dat$Computed[selJournal & !is.na(dat$Computed)])
  sigRes <- sum(dat$Computed[selJournal & !is.na(dat$Computed)] < .05)
  print(
    paste0(
      journals[j], " ", len, "\t", sigRes, "\t", len-sigRes
    )
  )
}

# Effect PDF
plot(density(dat$esComp[!is.na(dat$esComp)]),
     lty=1,
     frame.plot=T, 
     main="",
     xlim=c(0,1),
     ylim=c(0,4.5),
     xaxs="i",
     yaxs="i",
     xlab="Correlation",
     ylab = "Density",
     cex.axis=.8,
     cex.lab=1,
     col = "black", las=1)


####### 2
nsig <- dat$Computed >= .05
esR <- c(.1, .25, .4)
set.seed(1234)
simNullEs <- simNullDist(dat, n.iter=length(dat$esComp)*3, alpha=.05)
simNullEs$adjESComp[simNullEs$adjESComp < 0] <- 0
# par(mfrow=c(2,1))
# Overall
plot(ecdf(na.omit(sqrt(simNullEs$esComp))),
     lty=1,
     frame.plot=T, 
     main="",
     xlim=c(0,1),
     xaxs="i",
     yaxs="i",
     xlab="Correlation",
     ylab = "Cumulative density",
     cex.axis=.8,
     cex.lab=1,
     col = "grey", las=1)
lines(ecdf(na.omit(sqrt(dat$esComp[nsig]))))
legend(x=.65,y=.8,legend=c(expression('H'[0]), 'Observed'),
       cex=.8,lty=c(1,1),
       col = c("grey","black",2),box.lwd=0 ,lwd=2, bty='n')
for(es in esR){
x <- sqrt(dat$esComp[nsig]) <= es
horiz <- sum(x[!is.na(x)]) / length(x[!is.na(x)])
text(x=.05, y=horiz+.02, labels=round(horiz, 2), cex=.8)
clip(0, es, 0, horiz)
abline(h=horiz, v=es)
clip(0, 1, 0, 1)
}

# Overall[adj]
plot(ecdf(na.omit(sqrt(simNullEs$adjESComp))),
     lty=1,
     frame.plot=T, 
     main="",
     xlim=c(0,1),
     xaxs="i",
     yaxs="i",
     xlab="Partial eta-squared",
     ylab = "Cumulative density",
     cex.axis=.8,
     cex.lab=1,
     col = "grey", las=1)
lines(ecdf(na.omit(sqrt(dat$adjESComp[nsig]))))
legend(x=.65,y=.8,legend=c(expression('H'[0]), 'Observed'),
       cex=.8,lty=c(1,1),
       col = c("grey","black",2),box.lwd=0 ,lwd=2, bty='n')
for(es in esR){
  x <- sqrt(dat$adjESComp[nsig]) <= es
  horiz <- sum(x[!is.na(x)]) / length(x[!is.na(x)])
  text(x=.05, y=horiz+.02, labels=round(horiz, 2), cex=.8)
  clip(0, es, 0, horiz)
  abline(h=horiz, v=es)
  clip(0, 1, 0, 1)
}

# par(mfrow=c(3,2))
par(ask=T)
for(i in 1:length(sort(unique(dat$Journal)))){
  sel <- dat$Journal == sort(unique(dat$Journal))[i] & nsig
  set.seed(i)
  simNullEs <- simNullDist(dat, n.iter=length(dat$esComp[sel])*3, alpha=.05)
  plot(ecdf(na.omit(sqrt(simNullEs$esComp))),
       lty=1,
       frame.plot=T, 
       main="",
       xlim=c(0,1),
       xaxs="i",
       yaxs="i",
       xlab="Correlation",
       ylab = "Cumulative density",
       cex.axis=.8,
       cex.lab=1,
       col = "grey", las=1)
  lines(ecdf(sqrt(dat$esComp[sel])),
        lwd=.5)
  for(es in esR){
    x <- sqrt(dat$esComp[sel]) <= es
    horiz <- sum(x[!is.na(x)]) / length(x[!is.na(x)])
    text(x=.05, y=horiz+.02, labels=round(horiz, 2), cex=.8)
    clip(0, es, 0, horiz)
    abline(h=horiz, v=es)
    clip(0, 1, 0, 1)
  }
}

######### 3
# Kolmogorov-Smirnov test
print("Overall")
suppressWarnings(ks.test(simNullEs$esComp, dat$esComp[nsig], alternative='greater'))
for(i in 1:length(sort(unique(dat$Journal)))){
  print(sort(unique(dat$Journal))[i])
  print(suppressWarnings(ks.test(simNullEs$esComp,
          dat$esComp[dat$Journal == sort(unique(dat$Journal))[i] & nsig],
          alternative="greater")))
}

######## 4
# Simulation study
# For actual code, see sourced file (not included for parsimony)
set.seed(35438759)
N <- c(25, # SPECIFY!
       median(dat$df2[!is.na(dat$df2)]),
       150 # SPECIFY
       )

ES <- c(.00,
        seq(.01,.99,.01))

P <- c(seq(1, 10, 1), seq(15, 50, 5))

alpha <- .05
alphaF <- 0.10
n.iter <- 10000
source('c.Simulation/simCode.R')
# Load all files back in
files <- list.files('c.Simulation/')[-5]
if(!require(stringr)){install.packages('stringr')}
names <- str_sub(files,start=1L, end=-5L)
for(i in 1:length(files)){
  assign(x=names[i],read.csv2(paste0('c.Simulation/', files[i])))
  assign(x=names[i],t(get(x=names[i])[,-1]))
}

# Data for table
t(get(x=names[1]))
t(get(x=names[2]))
t(get(x=names[3]))


#########################################
# Simulation plots
par(mfrow=c(4,1))
plot(ES, N_25[,1], type='l',main="N = 25", col=1, ylab="Power",frame.plot=T,cex.axis=.8,cex.lab=1,las=1)
lines(ES, N_25[,2],col=2)
lines(ES, N_25[,3],col=3)
lines(ES, N_25[,4],col=4)
# lines(ES, N_25[,5],col=5)
abline(v=c(.06,.14),lty=3,col="grey")
legend(x=.65,y=.8,legend=paste(P, "results", sep=" "),
       cex=.8,lty=1,
       col = 1:4,box.lwd=0 ,lwd=2, bty='n')

plot(ES, N_64[,1], type='l',main="N = 100", col=1, ylab="Power",frame.plot=T,cex.axis=.8,cex.lab=1,las=1)
lines(ES, N_64[,2],col=2)
lines(ES, N_64[,3],col=3)
lines(ES, N_64[,4],col=4)
# lines(ES, N_64[,5],col=5)
abline(v=c(.06,.14),lty=3,col="grey")
legend(x=.65,y=.8,legend=paste(P, "results", sep=" "),
       cex=.8,lty=1,
       col = 1:4,box.lwd=0 ,lwd=2, bty='n')

plot(ES, N_150[,1], type='l',main="N = 150", col=1, ylab="Power",frame.plot=T,cex.axis=.8,cex.lab=1,las=1)
lines(ES, N_150[,2],col=2)
lines(ES, N_150[,3],col=3)
lines(ES, N_150[,4],col=4)
# lines(ES, N_150[,5],col=5)
abline(v=c(.06,.14),lty=3,col="grey")
legend(x=.65,y=.8,legend=paste(P, "results", sep=" "),
       cex=.8,lty=1,
       col = 1:4,box.lwd=0 ,lwd=2, bty='n')
#######################################################

##########
# Step 5 #
# Test on paper level
##########
# Computing the Fisher Tests
fishRes <- FisherMethod(dat$Computed, dat$Source)
# Get vector to identify which papers are from which journal
sourcejour <- dat$Journal[unique(dat$Source)]
cbind(fishRes$PFish < .1, sourcejour)

esSize <- c(0.00, seq(0.01, 0.95, 0.01))
set.seed(94438)
powerRes <- powerCalc(dat,effectSize=esSize,n.iter=1000,testAlpha=.1)
write.csv2(powerRes[[1]],'powerCalcFish.csv')
write.csv2(powerRes[[2]],'powerCalcFishCompl.csv')