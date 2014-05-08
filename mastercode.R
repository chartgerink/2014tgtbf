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

# create color vector for plotting
colVec <- c("#40004b",
            "#762a83",
            "#9970ab",
            "#c2a5cf",
            "#e7d4e8",
            "#f7f7f7",
            "#d9f0d3",
            "#a6dba0",
            "#5aae61",
            "#1b7837",
            "#00441b")

##########
# Step 1 #
# Descriptive data on the dataset
# Split for overall and non-significant
##########
# Cumulative effect size distribution (overall)
par(mfrow=c(2,1))
par(mai=c(.6,1.2,.2,.2))
## Total
plot(ecdf(dat$esComp), main="",
     xlab="", ylab="Cumulative distribution", las=1,
     cex.axis=.8, col=colVec[1], yaxs='i', xlim=c(0,1), oma=c(1,1,0,0)
     )
## Field
for(i in 1:length(sort(unique(dat$Journal)))){
  lines(ecdf(dat$esComp[dat$Journal == sort(unique(dat$Journal))[i]]), col=colVec[i+1],
        lwd=.5)
}
legend(x=.75, y=.85, legend=c("Overall", sort(unique(dat$Journal))), cex=.8,
       lty=1, col = colVec[1:11], lwd=2, bty='n')
par(mai=c(1.2,1.2,.2,.2))
## Total [nonsignificant]
plot(ecdf(dat$esComp[dat$Computed > .05]), main="",
     xlab="Effect size (eta)", ylab="Cumulative distribution", las=1,
     cex.axis=.8, col=colVec[1], yaxs='i', xlim=c(0,1), oma=c(1,1,0,0)
)
## Field [nonsignificant]
for(i in 1:length(sort(unique(dat$Journal[dat$Computed > .05])))){
  lines(ecdf(dat$esComp[dat$Computed > .05][dat$Journal[dat$Computed > .05] == sort(unique(dat$Journal[dat$Computed > .05]))[i]]), col=colVec[i+1],
        lwd=.5)
}
legend(x=.75, y=.85, legend=c("Overall", sort(unique(dat$Journal[dat$Computed > .05]))), cex=.8,
       lty=1, col = colVec[1:11], lwd=2, bty='n')

# P-value distribution
## Total
par(mai=c(.6,1.2,.2,.2))
plot(density(na.omit(dat$Computed)), main="",
     xlab="", ylab="Density", las=1,
     cex.axis=.8, yaxs='i', xlim=c(0,1), col=colVec[1], xaxs='i'
)
## Field
for(i in 1:length(sort(unique(dat$Journal)))){
  lines(density(na.omit(dat$Computed[dat$Journal == sort(unique(dat$Journal))[i]])), col=colVec[i+1],
        lwd=.5)
}
legend(x=.75, y=20, legend=c("Overall", sort(unique(dat$Journal))), cex=.8,
       lty=c(1, 1), col = colVec, box.lwd=0, bty='n', lwd=2)
par(mai=c(1.2,1.2,.2,.2))
## Total [nonsignificant]
plot(density(na.omit(dat$Computed[dat$Computed > .05])), main="",
     xlab="P-value", ylab="Density", las=1,
     cex.axis=.8, yaxs='i', xlim=c(0,1), col=colVec[1], xaxs='i'
)
## Field [nonsignificant]
for(i in 1:length(sort(unique(dat$Journal[dat$Computed > .05])))){
  lines(density(na.omit(dat$Computed[dat$Computed > .05][dat$Journal[dat$Computed > .05] == sort(unique(dat$Journal[dat$Computed > .05]))[i]])), col=colVec[i+1],
        lwd=.5)
}
legend(x=.75, y=2, legend=c("Overall", sort(unique(dat$Journal[dat$Computed > .05]))), cex=.8,
       lty=c(1, 1), col = colVec, box.lwd=0, bty='n', lwd=2)

par(mfrow=c(1,1))
# Effect size distribution (sig vs n-sig)
sig <- dat$Computed < .05
nsig <- dat$Computed >= .05
plot(ecdf(dat$esComp[sig]), main="",
     xlab="Effect size (eta)", ylab="Cumulative distribution", las=1, cex.axis=.8,
     col=colVec[1], xlim=c(0,1), yaxs='i')
lines(ecdf(na.omit(dat$esComp[nsig])), col=colVec[8])
legend(x=.75,y=.85,legend=c("Significant", "Non-significant"), cex=.8,
       lty=c(1, 1), col = c(colVec[1], colVec[8]), box.lwd=0, lwd=2, bty='n')

# Proportion significant
## Information readily pastable for table
testVal <- c("t", "F", "r")
journals <- sort(unique(dat$Journal))
for(t in 1:length(testVal)){
  selTest <- dat$Statistic == testVal[t]
  for(j in 1:length(journals)){
    selJournal <- dat$Journal == journals[j]
    len <- length(dat$Computed[selTest & selJournal])
    sigRes <- sum(dat$Computed[selTest & selJournal] < .05)
    print(
      paste0(
        journals[j], " ", testVal[t], " ", len, "\t", sigRes, "\t", len-sigRes
        )
      )
  }
}


##########
# Step 3 #
# Observed vs null effect
##########
set.seed(1234)
simNullEs <- simNullDist(dat, n.iter=length(dat$esComp)*3, alpha=.05)
# Figure showing the difference between the nil effects and the observed effect distributions
# Under non-significant test
plot(ecdf(na.omit(simNullEs$esComp)),
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
lines(ecdf(na.omit(dat$esComp)))
legend(x=.65,y=.8,legend=c(expression('H'[0]), 'Effects'),
       cex=.8,lty=c(1,1),
       col = c("grey","black"),box.lwd=0 ,lwd=2, bty='n')

# Kolmogorov-Smirnov test
ks.test(simNullEs$esComp, dat$esComp, alternative='greater')

########
# Step 4
# Simulation study results
# For actual code, see folder c.Simulations
########
# For running the simulation.
# Commented out to prevent redundant rerunning.
set.seed(35438759)
source('c.Simulation/simCode.R')
# Load all files back in
files <- list.files('c.Simulation/')[-5]
if(!require(stringr)){install.packages('stringr')}
names <- str_sub(files,start=1L, end=-5L)
for(i in 1:length(files)){
  assign(x=names[i],read.csv2(paste0('c.Simulation/', files[i])))
  assign(x=names[i],t(get(x=names[i])[,-1]))
}

# Simulation plots
par(mfrow=c(2,2))
plot(ES, N_25[,1], type='l',main="N = 25", col=1, ylab="Power",frame.plot=T,cex.axis=.8,cex.lab=1,las=1)
lines(ES, N_25[,2],col=2)
lines(ES, N_25[,3],col=3)
lines(ES, N_25[,4],col=4)
# lines(ES, N_25[,5],col=5)
abline(v=c(.06,.14),lty=3,col="grey")

plot(ES, N_50[,1], type='l',main="N = 50", col=1, ylab="Power",frame.plot=T,cex.axis=.8,cex.lab=1,las=1)
lines(ES, N_50[,2],col=2)
lines(ES, N_50[,3],col=3)
lines(ES, N_50[,4],col=4)
# lines(ES, N_50[,5],col=5)
abline(v=c(.06,.14),lty=3,col="grey")

plot(ES, N_100[,1], type='l',main="N = 100", col=1, ylab="Power",frame.plot=T,cex.axis=.8,cex.lab=1,las=1)
lines(ES, N_100[,2],col=2)
lines(ES, N_100[,3],col=3)
lines(ES, N_100[,4],col=4)
# lines(ES, N_100[,5],col=5)
abline(v=c(.06,.14),lty=3,col="grey")

plot(ES, N_150[,1], type='l',main="N = 150", col=1, ylab="Power",frame.plot=T,cex.axis=.8,cex.lab=1,las=1)
lines(ES, N_150[,2],col=2)
lines(ES, N_150[,3],col=3)
lines(ES, N_150[,4],col=4)
# lines(ES, N_150[,5],col=5)
abline(v=c(.06,.14),lty=3,col="grey")

##########
# Step 5 #
# Test on paper level
##########
# Computing the Fisher Tests
resPilot <- FisherExTest(dat$Computed, dat$Source)

esSize <- c(seq(.01,.15,.02), seq(.35,.95,.2))
set.seed(94438)
powerRes <- powerCalc(dat,effectSize=esSize,n.iter=1000,testAlpha=.1)
write.csv2(powerRes[[1]],'powerCalcFish.csv')
write.csv2(powerRes[[2]],'powerCalcFishCompl.csv')