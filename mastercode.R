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


# Introduction
# Assumptions calculations below
# 1. No QRPs
# 2. No violations of assumptions
# 3. Equal likelihood of null and alternative being true
Rrate <- .5
signifFindings <- c(.95, .97)
powEstimate <- c(.35, .50)
alpha <- .05
# False positive rate
Rrate*alpha + (1-Rrate)*powEstimate
# False negative rate
Rrate*(1-alpha) + (1-Rrate)*powEstimate

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
# Add effect size proportions S-M-L?
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
abline(v=c(.1,.25,.4), lty=2, col="grey")


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
  abline(h=horiz, v=es, lty=2, col="grey")
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
  abline(h=horiz, v=es, col="grey", lty=2)
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
    abline(h=horiz, v=es, col="grey", lty=2)
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

##########
# Step 5 #
# Test on paper level
##########
# Computing the Fisher Tests
fishRes <- FisherMethod(dat$Computed, dat$Source)
# Get vector to identify which papers are from which journal
sourcejour <- dat$Journal[unique(dat$Source)]
# Make vector of sig/nsig/na 
fishDF <- data.frame(FisherP=fishRes$PFish, journal=sourcejour, kRes=fishRes$CountNSig)
alphaF <- 0.10
# Compute amount of papers and proportion of sig/nsig/NA Fisher Method tests
final <- NULL
kLen <- c(1, 2, 3, 4, 5, 10, 20)

for(journals in sort(unique(dat$Journal))){
  sel <- fishDF$journal == journals
  
  # Amount of papers in a journal
  amount <- length(fishDF$FisherP[sel])
  # Proportion of significant fisher results
  amountSig <- sum(fishDF$FisherP[sel & !is.na(fishDF$FisherP)] < alphaF)
  
  for(k in 1:length(kLen)){
    if(k == length(kLen)){
      assign(paste0('amount', kLen[k]+1), length(fishDF$FisherP[sel & fishDF$kRes > kLen[k]]))
      assign(paste0('amountSig', kLen[k]+1),
             sum(fishDF$FisherP[sel & !is.na(fishDF$FisherP) & fishDF$kRes > kLen[k]] < alphaF))
    }
    else{
      assign(paste0('amount', kLen[k]), length(fishDF$FisherP[sel & fishDF$kRes <= kLen[k]]))
      if(kLen[k] == kLen[1]){
        assign(paste0('amountSig', kLen[k]),
             sum(fishDF$FisherP[sel & !is.na(fishDF$FisherP) & fishDF$kRes <= kLen[k]] < alphaF))}
      else{
        assign(paste0('amountSig', kLen[k]),
               sum(fishDF$FisherP[sel & !is.na(fishDF$FisherP) & fishDF$kRes <= kLen[k] & fishDF$kRes > kLen[k-1]] < alphaF))
      }
    }
  }
  
  # Amount of papers in a journal without significant results
  countNA <- sum(is.na(fishDF$FisherP[sel]))
  journalSet <- NULL  
  # Writing out the results
  for(k in 1:length(kLen)){
    if(kLen[k] == kLen[length(kLen)]){
      x <- get(paste0('amountSig', kLen[k]+1)) / get(paste0('amount', kLen[k]+1))
    }
    else{
      x <- get(paste0('amountSig', kLen[k])) / get(paste0('amount', kLen[k]))
    }
    journalSet <- cbind(journalSet, x)
  }
  temp <- cbind(journals,
                journalSet,
                amountSig / amount,
                countNA,
                amountSig,
                amount)
  # This is the result that goes into the table
  final <- rbind(final, temp)
}
final <- as.data.frame(final)
names(final) <- c('journals', paste0('k ', kLen), 'overall', 'countNA', 'amountSig', 'nrpapers')
round(final, 3)

# Ad hoc effect estimation
esSize <- seq(.00, .99, .01)
set.seed(9864)
powerRes <- powerCalc(dat, effectSize=esSize, n.iter=1000, alphaF=.1)
write.csv2(powerRes,'powerCalcFish.csv')
effectDat <- read.csv2('../testing/powerCalcFish.csv')
names(effectDat) <- c('X', 'Journal', paste0('ES ', esSize))

# Compute the expected number of significant fisher tests
# Overall
estimatedCorr <- NULL
expectedOverall <- apply(effectDat[,-c(1,2)], 2, sum)
plot(x=esSize,
     y=expectedOverall/length(unique(dat$Source)),
     ylim=c(0,1),
     type="l",
     xlab="Correlation",
     ylab="Proportion significant",
     xaxs="i",
     yaxs="i",
     cex.axis=.8,
     cex.lab=1,
     las=1)

observed <- sum(final$amountSig) / sum(final$nrpapers)
expected <- expectedOverall / length(effectDat$Journal)
minimum <- min(abs(observed-expected))
horiz <- observed
estimatedCorr <- esSize[expected-observed == minimum]
clip(0, estimatedCorr, 0, horiz )
abline(h=horiz, v=estimatedCorr, lty=1, col="grey")
clip(0,1,0,1)

# Per journal
for(i in 1:length(unique(effectDat$Journal))){
  assign(paste0('expected', unique(effectDat$Journal)[i]),
         apply(effectDat[effectDat$Journal == unique(effectDat$Journal)[i],-c(1,2)], 2, sum))
}
for(i in 1:length(unique(effectDat$Journal))){
  lines(esSize,
        get(paste0('expected', unique(effectDat$Journal)[i]))/sum(effectDat$Journal == unique(effectDat$Journal)[i]), lty=i+1)
}
for(i in 1:length(unique(effectDat$Journal))){
  observed <- final$amountSig[i] / final$nrpapers[i]
  expected <- get(paste0('expected', unique(effectDat$Journal)[i])) / sum(effectDat$Journal == unique(effectDat$Journal)[i])
  minimum <- min(abs(observed-expected))
  horiz <- observed
  estimatedCorr <- c(estimatedCorr, esSize[expected-observed == minimum])
  clip(0, tail(estimatedCorr, 1), 0, horiz )
  abline(h=horiz, v=tail(estimatedCorr, 1), lty=i+1, col="grey")
}

# Save the ad hoc estimations
estimatedCorr <- data.frame(journal=c('Overall', unique(effectDat$Journal)), estimatedCorr)

##########
# Step 6 #
# Relation k and significant Fisher tests
##########
zcv <- qt(alpha, df=Inf, lower.tail=F)
k <- seq(1, 500, 1)
plot(1-(zcv*(1/sqrt(k))), ylim=c(0,1), type="s", xlab="k", ylab="Proportion significant")
