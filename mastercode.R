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
dat$Journal <- sample(1:10, 8110, replace=T)

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
##########

# Cumulative effect size distribution (overall)
par(mai=c(1.2,1.2,.2,.2))
## Total
plot(ecdf(dat$esComp), main="",
     xlab="Effect size (eta)", ylab="Cumulative distribution", las=1,
     cex.axis=.8, col=colVec[1], yaxs='i', xlim=c(0,1), oma=c(1,1,0,0)
     )
## Field
for(i in 1:length(sort(unique(dat$Journal)))){
  lines(ecdf(dat$esComp[dat$Journal == sort(unique(dat$Journal))[i]]), col=colVec[i+1],
        lwd=.5)
}
legend(x=.75, y=.85, legend=c("Overall", sort(unique(dat$Journal))), cex=.8,
       lty=1, col = colVec[1:11], lwd=2, bty='n')

# Effect size distribution (sig vs n-sig)
sig <- dat$Computed < .05
nsig <- dat$Computed >= .05
plot(ecdf(dat$esComp[sig]), main="",
     xlab="P-value", ylab="Cumulative distribution", las=1, cex.axis=.8,
     col=colVec[1], xlim=c(0,1), yaxs='i')
lines(ecdf(na.omit(dat$esComp[nsig])), col=colVec[8])
legend(x=.7,y=.85,legend=c("Significant", "Non-significant"), cex=.8,
       lty=c(1, 1), col = c(colVec[1], colVec[8]), box.lwd=0, lwd=2, bty='n')

# P-value distribution
## Total
plot(density(na.omit(dat$Computed)), main="",
     xlab="P-value", ylab="Density", las=1,
     cex.axis=.8, yaxs='i', xlim=c(0,1), col=colVec[1], xaxs='i'
)
## Field
for(i in 1:length(sort(unique(dat$Journal)))){
  lines(density(na.omit(dat$Computed[dat$Journal == sort(unique(dat$Journal))[i]])), col=colVec[i+1],
        lwd=.5)
}
legend(x=.75, y=15, legend=c("Overall", sort(unique(dat$Journal))), cex=.8,
       lty=c(1, 1), col = colVec, box.lwd=0, bty='n')

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
# Step 2 #
# Descriptive data on the dataset
##########


##########
# Step 3 #
# Observed vs null effect
##########
set.seed(1234)
simNullEs <- simNullDist(dat, n.iter=length(dat$esComp), alpha=.05)
# Figure showing the difference between the nil effects and the observed effect distributions
# Under non-significant test
plot(ecdf(na.omit(simNullEs$esComp)),
     lty=1,
     frame.plot=F, 
     main="Effect size distributions",
     xlim=c(0,1),
     xaxs="i",
     yaxs="i",
     xlab="Partial eta-squared",
     ylab = "Cumulative density",
     cex.axis=.6,
     cex.lab=.7,
     col = "grey")
lines(ecdf(na.omit(dat$esComp)))
lines(ecdf(na.omit(dat$adjESComp)),col="blue")
legend(x=.7,y=.8,legend=c(expression('H'[0]), 'Observed effects'),cex=.7,lty=c(1,1), col = c("grey","black"),box.lwd=0)

# Kolmogorov-Smirnov test
ks.test(simNullEs$esComp, dat$esComp, alternative='greater')




##########
# Step 5 #
# Test on paper level
##########

# Computing the Fisher Tests
resPilot <- FisherExTest(dat$Computed, dat$Source)

esSize <- c(seq(.01,.15,.02), seq(.2,.9,.1))
set.seed(94438)
powerRes <- powerCalc(dat,effectSize=esSize,n.iter=1000,testAlpha=.1)
write.csv2(powerRes[[1]],'powerCalcFish.csv')
write.csv2(powerRes[[2]],'powerCalcFishCompl.csv')