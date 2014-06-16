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
# dat <- read.table("1.Pilot study/copilot.txt", stringsAsFactors=F)
dat <- read.csv2("fullFile.csv", stringsAsFactors=F)
# There are two test statistic indicators that are NA
# Manually correct these
dat$Statistic[is.na(dat$Statistic)] <- "F"

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
# I.e., alpha error control

# False positive rate
falsePos1 <- 1-((.25*.35) / ((.25*.35) + (.75*.05)))
falsePos2 <- 1-((.5*.5) / ((.5*.5) + (.5*.05)))
sort(c(falsePos1, falsePos2))

# False negative rate
falseNeg1 <- ((.25*.5) / ((.25*.5) + (.75*.95)))
falseNeg2 <- ((.5*.65) / ((.5*.65) + (.5*.95)))
sort(c(falseNeg1, falseNeg2))

####### 1
# Descriptives full dataset
# Table
journals <- sort(unique(dat$journals.jour.))
for(j in 1:length(journals)){
  selJournal <- dat$journals.jour. == journals[j]
  len <- length(dat$Computed[selJournal & !is.na(dat$Computed)])
  sigRes <- sum(dat$Computed[selJournal & !is.na(dat$Computed)] < .05)
  print(
    paste0(
      journals[j], "\t", len, "\t", sigRes, "\t", len-sigRes
    )
  )
}

# Selecting only the t, r and F values
dat <- dat[dat$Statistic == 't' | dat$Statistic == 'r' | dat$Statistic == 'F',]

# Effect PDF
# Add effect size proportions S-M-L?
plot(density(dat$esComp[!is.na(dat$esComp)]),
     lty=1,
     frame.plot=T, 
     main="",
     xlim=c(0,1),
#      ylim=c(0,4.5),
     xaxs="i",
     yaxs="i",
     xlab="Correlation",
     ylab = "Density",
     cex.axis=.8,
     cex.lab=1,
     col = "black", las=1)
abline(v=c(.1,.25,.4), lty=2, col="grey")
t1 <- sum(dat$esComp[!is.na(dat$esComp)] < .1) / length(dat$esComp[!is.na(dat$esComp)])
t2 <- sum(dat$esComp[!is.na(dat$esComp)] >= .1 & dat$esComp[!is.na(dat$esComp)] < .25) / length(dat$esComp[!is.na(dat$esComp)])
t3 <- sum(dat$esComp[!is.na(dat$esComp)] >= .25 & dat$esComp[!is.na(dat$esComp)] < .4) / length(dat$esComp[!is.na(dat$esComp)])
t4 <- sum(dat$esComp[!is.na(dat$esComp)] >= .4) / length(dat$esComp[!is.na(dat$esComp)])
text(x=.1/2, y=.5, labels=round(t1,2), cex=.8)
text(x=((.25-.1)/2)+.1, y=.5, labels=round(t2,2), cex=.8)
text(x=((.4-.25)/2)+.25, y=.5, labels=round(t3,2), cex=.8)
text(x=((1-.4)/2)+.4, y=.5, labels=round(t4,2), cex=.8)

####### 2
nsig <- dat$Computed >= .05
esR <- c(.1, .25, .4)
set.seed(1234)
simNullEs <- simNullDist(dat, n.iter=length(dat$esComp)*3, alpha=.05)
simNullEs$adjESComp[simNullEs$adjESComp < 0] <- 0
write.csv2(simNullEs, 'simNullEs.csv')
# simNullEs <- read.csv2('simNullEs.csv')
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
  clip(0, es, 0, 1)
  abline(h=horiz, v=es, lty=2, col="grey")
  h0horiz <- sum(sqrt(simNullEs$esComp[!is.na(simNullEs$esComp)]) < es) / length(simNullEs$esComp[!is.na(simNullEs$esComp)])
  abline(h=h0horiz, v=es, lty=2, col="grey")
  text(x=.05, y=h0horiz+.01, labels=round(h0horiz, 2), cex=.6, col='darkgrey')
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
  clip(0, es, 0, 1)
  abline(h=horiz, v=es, col="grey", lty=2)
  h0horiz <- sum(sqrt(simNullEs$adjESComp[!is.na(simNullEs$adjESComp)]) < es) / length(simNullEs$adjESComp[!is.na(simNullEs$adjESComp)])
  abline(h=h0horiz, v=es, lty=2, col="grey")
  text(x=.05, y=h0horiz+.01, labels=round(h0horiz, 2), cex=.6, col='darkgrey')
  clip(0, 1, 0, 1)
}

par(mfrow=c(2,4))
par(ask=F)
for(i in 1:length(sort(unique(dat$journals.jour.)))){
  sel <- dat$journals.jour. == sort(unique(dat$journals.jour.))[i] & nsig
  set.seed(i)
  simNullEs <- simNullDist(dat, n.iter=length(dat$esComp[sel])*3, alpha=.05)
  temp <- ks.test(simNullEs$esComp,
                  dat$esComp[dat$journals.jour. == sort(unique(dat$journals.jour.))[i] & nsig],
                  alternative="greater")
  print((temp))
  plot(ecdf(na.omit(sqrt(simNullEs$esComp))),
       lty=1,
       frame.plot=T, 
       main=paste0(sort(unique(dat$journals.jour.))[i], ", D=", round(temp$statistic,3),", p=",
                   ifelse(temp$p.value>.001, round(temp$p.value,3), "<.001")),
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
    clip(0, es, 0, 1)
    abline(h=horiz, v=es, col="grey", lty=2)
    h0horiz <- sum(sqrt(simNullEs$esComp[!is.na(simNullEs$esComp)]) < es) / length(simNullEs$esComp[!is.na(simNullEs$esComp)])
    abline(h=h0horiz, v=es, lty=2, col="grey")
    text(x=.05, y=h0horiz+.01, labels=round(h0horiz, 2), cex=.6, col='darkgrey')
    clip(0, 1, 0, 1)
  }
}
# par(ask=F)

######### 3
# Kolmogorov-Smirnov test
print("Overall")
suppressWarnings(ks.test(simNullEs$esComp, dat$esComp[nsig], alternative='greater'))
for(i in 1:length(sort(unique(dat$journals.jour.)))){
  print(sort(unique(dat$journals.jour.))[i])
  print(suppressWarnings(ks.test(simNullEs$esComp,
                                 dat$esComp[dat$journals.jour. == sort(unique(dat$journals.jour.))[i] & nsig],
                                 alternative="greater")))
}

######## 4
# Simulation study
# For simulation code, see sourced file (not included for parsimony)

N <- c(as.numeric(summary(dat$df2)[2]), # 25th percentile
       as.numeric(summary(dat$df2)[3]), # 50th percentile
       as.numeric(summary(dat$df2)[5]) # 75th percentile
)

ES <- c(.00,
        seq(.01,.99,.01))

P <- c(seq(1, 10, 1), seq(15, 50, 5))

alpha <- .05
alphaF <- 0.10
n.iter <- 10000
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

# Data for table
t(get(x=names[1]))
t(get(x=names[2]))
t(get(x=names[3]))

##########
# Step 5 #
# Test on paper level
##########
# Computing the Fisher Tests
fishRes <- FisherMethod(dat$Computed, as.character(dat$Source))
# Get vector to identify which papers are from which journal
sourcejour <- NULL
for(i in 1:length(unique(dat$Source))){
  sourcejour[i] <- unique(dat$journals.jour.[dat$Source==unique(dat$Source)[i]])
}

# sourcejour <- dat$journals.jour.[unique(dat$Source)]
# Make vector of sig/nsig/na 
fishDF <- data.frame(FisherP=fishRes$PFish, journal=sourcejour, kRes=fishRes$CountNSig)
alphaF <- 0.10
# Compute amount of papers and proportion of sig/nsig/NA Fisher Method tests
final <- NULL
kLen <- c(1, 2, 3, 4, 5, 10, 20)

for(journals in sort(unique(dat$journals.jour.))){
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
names(final) <- c('journals', paste0('k', kLen), 'overall', 'countNA', 'amountSig', 'nrpapers')
# round(final, 3)

# Ad hoc effect estimation
esSize <- seq(.00, .99, .01)
# These are temporarily commented out to prevent re-runnin
powerRes <- NULL
for(i in 1:length(sort(unique(dat$Source)))){
  sel <- dat$Source == sort(unique(dat$Source))[i]
  set.seed(9864+i)
  temp <- powerCalc(dat[sel,], effectSize=esSize, n.iter=1000, alphaF=.1)
  powerRes <- rbind(powerRes, temp)
  print(paste(i, "of", length(sort(unique(dat$Source)))))
}
write.csv2(powerRes,'powerCalcFish.csv')
# effectDat <- read.csv2('../testing/powerCalcFish.csv')
effectDat <- read.csv2('powerCalcFish.csv')

names(effectDat) <- c('X', 'Journal', 'k', paste0('ES', esSize))

# Compute the expected number of significant fisher tests
# Overall
estimatedCorr <- NULL
estimatedCorrJournal <- NULL
expectedOverall <- apply(effectDat[,-c(1,2,3)], 2, sum)
plot(x=esSize,
     y=expectedOverall / length(effectDat$Journal),
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
estimatedCorrOverall <- esSize[abs(observed-expected) == minimum]
clip(0, estimatedCorrOverall, 0, horiz )
abline(h=horiz, v=estimatedCorrOverall, lty=1, col="grey")
clip(0,1,0,1)

# Per journal
for(i in 1:length(unique(effectDat$Journal))){
  assign(paste0('expected', sort(unique(effectDat$Journal))[i]),
         apply(effectDat[effectDat$Journal == sort(unique(effectDat$Journal))[i],-c(1,2,3)], 2, sum))
}
for(i in 1:length(unique(effectDat$Journal))){
  lines(esSize,
        get(paste0('expected', sort(unique(effectDat$Journal))[i]))/sum(effectDat$Journal == sort(unique(effectDat$Journal))[i]), lty=i+1)
}
for(i in 1:length(unique(effectDat$Journal))){
  observed <- final$amountSig[i] / final$nrpapers[i]
  expected <- get(paste0('expected', sort(unique(effectDat$Journal))[i])) / sum(effectDat$Journal == sort(unique(effectDat$Journal))[i])
  minimum <- min(abs(observed-expected))
  horiz <- observed
  estimatedCorrJournal[i] <- esSize[abs(observed-expected) == minimum]
  clip(0, estimatedCorrJournal[i], 0, horiz )
  abline(h=horiz, v=estimatedCorrJournal[i], lty=i+1, col="grey")
  clip(0, 1, 0, 1)
}

# Save the ad hoc estimations
estimatedCorr <- data.frame(journal=c('Overall', sort(unique(effectDat$Journal))), c(estimatedCorrOverall, estimatedCorrJournal))

##########
# Step 6 #
# Relation k and significant Fisher tests
##########


# Select out all zero-k papers
fishDF <- fishDF[!fishDF$kRes == 0,]
k <- sort(unique(fishDF$kRes))
curveES <- seq(0, 1.5, .01)
zcv <- qnorm(.9, 0, 1)

fishDF$sig <- fishDF$FisherP
fishDF$sig[fishDF$sig < alphaF] <- 1
fishDF$sig[fishDF$sig >= alphaF & !fishDF$sig == 1] <- 0

r2fit = NULL
ktemp <- seq(1:max(k))
datFit <- data.frame(yi=fishDF$sig, kRes=fishDF$kRes, jour=fishDF$journal)
xi <- list(NULL)
for(i in 1:length(curveES)){
  xi[[i]] <- (1 - pnorm(zcv/sqrt(datFit$kRes), curveES[i], 1 / sqrt(datFit$kRes)))
  datFit <- cbind(datFit, xi[[i]])
  r2fit[i] <- summary(lm(datFit$yi ~ 1 + datFit[,3+i]))$r.squared
}

# Saturated model
satur <- summary(lm(datFit$yi ~ 0 + as.factor(datFit$kRes)))


plot(x=datFit$kRes, y=datFit$yi, col= 'grey')
curve((1 - pnorm(zcv/sqrt(x), curveES[which(r2fit == max(r2fit))], 1 / sqrt(x))), from=1, to=max(datFit$kRes),
      add=T)
lines(x=sort(unique(datFit$kRes)), y=satur$coefficients[,1],col='red')

for(z in 1:length(unique(datFit$jour))){
  sel <- datFit$jour == sort(unique(datFit$jour))[z]
  kMean <- mean(datFit$kRes[sel])
  propMean <- mean(datFit$yi[sel], na.rm=T)
  points(x=kMean, y=propMean, col="blue", pch=19)
}

max(r2fit)
satur$r.squared

