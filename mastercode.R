# Code written by CHJ Hartgerink
# Checked by: -

##################################################
# DO NOT EDIT PAST HERE FOR FULL REPRODUCIBILITY #
# SIMULATIONS ARE NOT RE-RUN DUE TO RUNTIME      #
# IF YOU'D LIKE TO RE-RUN YOU WILL HAVE TO UN-   #
# COMMENT SOME PARTS                             #
##################################################
setwd(choose.dir())
if(!require(latex2exp)){install.packages('latex2exp')}
library(latex2exp)

# Load custom functions all at once
customFunct <- list.files('functions/')
for(i in 1:length(customFunct)){
  source(
    paste0('functions/', customFunct[i])
  )
}

# Read- and prepare data
dat <- read.csv2("data/statcheck_full_anonymized.csv", stringsAsFactors=F, dec = ",", sep = ";")[-1]

# There are two test statistic indicators that are NA
# Manually correct these
dat$Statistic[is.na(dat$Statistic)] <- "F"

# Computing unadjusted and adjusted effect sizes (OBSERVED)
dat <- cbind(dat, esComp.statcheck(dat))
dat$adjESComp[dat$adjESComp < 0] <- 0

# Turning df1 for t and r into 1.
dat$df1[dat$Statistic == "t" | dat$Statistic == "r"] <- 1

# Select out incorrectly exttracted r values
dat <- dat[!(dat$Statistic=="r" & dat$Value > 1),]

# Select out irrefutably wrong df reporting
dat <- dat[!dat$df1 == 0,]

# select out NA computed p-values
dat <- dat[!is.na(dat$Computed),]

# Selecting only the t, r and F values
dat <- dat[dat$Statistic == 't' | dat$Statistic == 'r' | dat$Statistic == 'F',]
nsig <- dat$Computed >= .05
esR <- c(.1, .25, .4)
alpha = .05
alphaF = .1

##########
# Method #
##########

# Simulation sample sizes
c(as.numeric(summary(dat$df2)[2]), # 25th percentile
  as.numeric(summary(dat$df2)[3]), # 50th percentile
  as.numeric(summary(dat$df2)[5])) # 75th percentile

# Number of papers
length(unique(dat$Source))

###########
# Results #
########### 
# Descriptives dataset
# Table
journals <- sort(unique(dat$journals.jour.))
for(j in 1:length(journals)){
  selJournal <- dat$journals.jour. == journals[j]
  meanK <- mean(table(dat$Source[selJournal]))
  len <- length(dat$Computed[selJournal & !is.na(dat$Computed)])
  sigRes <- sum(dat$Computed[selJournal & !is.na(dat$Computed)] < .05)
  print(
    paste0(
      journals[j], "\t", len, "\t", meanK, "\t", sigRes, "\t", len-sigRes
    )
  )
  #   print(meanK)
}

# nonsignificant proportion p/year
i <- 1
sig <- NULL
nsigtemp <- NULL
kval <- NULL

for(y in 1985:2013){
  sel <- dat$years.y. == y
  sig[i] <- sum(dat$Computed[sel] < alpha) / length(dat$Computed[sel])
  nsigtemp[i] <- sum(dat$Computed[sel] > alpha) / length(dat$Computed[sel])
  kval[i] <- median(table(dat$Source[sel])) / sum(table(dat$Source[sel]))
  i <- i + 1
}

save(nsigtemp, file = 'figures/Fig4')
rm(nsigtemp)

####### 2 + 3
## Uncomment the next four lines to re-run the simulations
set.seed(1234)
simNullEs <- simNullDist(dat, n.iter = length(dat$esComp[nsig]) * 3, alpha = .05)
simNullEs$adjESComp[simNullEs$adjESComp < 0] <- 0
write.table(simNullEs, 'archive/simNullEs.csv', sep = ";", dec = ".")

simNullEs <- read.table('archive/simNullEs.csv', sep = ";", dec = ".")
temp <- ks.test(simNullEs$esComp,
                dat$esComp[nsig],
                alternative="greater")

pdf('figures/Fig5.pdf', width=11, height=7)
par(mfrow = c(1,2), mai = c(1.2,1.2,.8,.5))
# Overall
plot(ecdf(na.omit(sqrt(simNullEs$esComp))),
     lty = 1,
     frame.plot = T, 
     main = latex2exp(
       sprintf(
         "Unadjusted, $D=%s,p<2.2\\times 10^{-16}$", round(temp$statistic,2))),
     xlim = c(0, 1),
     xaxs = "i",
     yaxs = "i",
     xlab = latex2exp("Correlation ($\\eta$)"),
     ylab = "Cumulative density",
     cex.axis=.8,
     cex.lab=1,
     cex.main=1.5,
     col = "grey", las=1)
lines(ecdf(na.omit(sqrt(dat$esComp[nsig]))))
legend(x = .6, y = .2, legend = c(latex2exp("$H_0$"), 'Observed'),
       cex = 1, lty = c(1, 1),
       col = c("grey", "black", 2),
       box.lwd = 0, lwd = 2, bty = 'n')
for(es in esR){
  h0horiz <- sum(sqrt(simNullEs$esComp[!is.na(simNullEs$esComp)]) < es) / length(simNullEs$esComp[!is.na(simNullEs$esComp)])
  clip(es, 1, 0, h0horiz)
  abline(h = h0horiz, v = es, lty = 2, col = "grey")
  clip(0, 1, 0, 1)
  text(x = .9, y = h0horiz - .02, labels = round(h0horiz, 2), cex = 1, col = 'darkgrey')
  x <- sqrt(dat$esComp[nsig]) <= es
  horiz <- sum(x[!is.na(x)]) / length(x[!is.na(x)])
  text(x = .05, y = horiz + .02, labels = round(horiz, 2), cex = 1)
  clip(0, es, 0, horiz)
  abline(h = horiz, v = es, lty = 2, col = "black")
  clip(0, 1, 0, 1)
}

temp <- ks.test(simNullEs$adjESComp,
                dat$adjESComp[nsig],
                alternative="greater")
# Overall[adj]
plot(ecdf(na.omit(sqrt(simNullEs$adjESComp))),
     lty = 1,
     frame.plot = T, 
     main = latex2exp(
       sprintf(
         "Adjusted, $D=%s,p<2.2\\times 10^{-16}$", round(temp$statistic,2))),
     xlim = c(0, 1),
     xaxs = "i",
     yaxs = "i",
     xlab = latex2exp("Correlation ($\\eta$)"),
     ylab = "Cumulative density",
     cex.axis = .8,
     cex.lab = 1,
     cex.main = 1.5,
     col = "grey", las = 1)
lines(ecdf(na.omit(sqrt(dat$adjESComp[nsig]))))
legend(x = .6,y = .2,legend = c(latex2exp("$H_0$"), 'Observed'),
       cex = 1,lty = c(1, 1),
       col = c("grey", "black", 2), box.lwd = 0, lwd = 2, bty = 'n')
for(es in esR){
  h0horiz <- sum(sqrt(simNullEs$adjESComp[!is.na(simNullEs$adjESComp)]) < es) / length(simNullEs$adjESComp[!is.na(simNullEs$adjESComp)])
  clip(es, 1, 0, h0horiz)
  abline(h = h0horiz, v = es, lty = 2, col = "grey")
  clip(0, 1, 0, 1)
  text(x = .9, y = h0horiz - .02, labels = round(h0horiz, 2), cex = 1, col = 'darkgrey')
  x <- sqrt(dat$adjESComp[nsig]) <= es
  horiz <- sum(x[!is.na(x)]) / length(x[!is.na(x)])
  text(x = .05, y = horiz + .02, labels = round(horiz, 2), cex = 1)
  clip(0, es, 0, horiz)
  abline(h = horiz, v = es, col = "black", lty = 2)
  clip(0, 1, 0, 1)
}
dev.off()

pdf('figures/Fig6.pdf', onefile = TRUE, width = 11, height = 11)
par(mfrow = c(2, 2), mai = c(1.2, 1.2, .8, .5))

for(i in 1:4){
  sel <- dat$journals.jour. == sort(unique(dat$journals.jour.))[i] & nsig
  set.seed(i)
  simNullEs <- simNullDist(dat, n.iter=length(dat$esComp[sel])*3, alpha=.05)
  temp <- ks.test(simNullEs$esComp,
                  dat$esComp[dat$journals.jour. == sort(unique(dat$journals.jour.))[i] & nsig],
                  alternative="greater")
  print((temp))
  if(i == 3){
    plot(ecdf(na.omit(sqrt(simNullEs$esComp))),
         lty=1,
         frame.plot=T, 
         main=paste0(sort(unique(dat$journals.jour.))[i], ", D=", round(temp$statistic,3),", p=7.934*10^-6"),
         xlim=c(0,1),
         xaxs="i",
         yaxs="i",
         xlab="Correlation",
         ylab = "Cumulative density",
         cex.axis=.8,
         cex.lab=1,
         cex.main=1.5,
         col = "grey", las=1)
  }
  else{plot(ecdf(na.omit(sqrt(simNullEs$esComp))),
            lty=1,
            frame.plot=T, 
            main=paste0(sort(unique(dat$journals.jour.))[i], ", D=", round(temp$statistic,3),", p<2.2*10^-16"),
            xlim=c(0,1),
            xaxs="i",
            yaxs="i",
            xlab="Correlation",
            ylab = "Cumulative density",
            cex.axis=.8,
            cex.lab=1,
            cex.main=1.5,
            col = "grey", las=1)}
  lines(ecdf(sqrt(dat$esComp[sel])),
        lwd=.5)
  legend(x=.6,y=.2,legend=c(expression('H'[0]), 'Observed'),
         cex=1.2,lty=c(1,1),
         col = c("grey","black",2),box.lwd=0 ,lwd=2, bty='n')
  for(es in esR){
    h0horiz <- sum(sqrt(simNullEs$esComp[!is.na(simNullEs$esComp)]) < es) / length(simNullEs$esComp[!is.na(simNullEs$esComp)])
    clip(es, 1, 0, h0horiz)
    abline(h=h0horiz, v=es, lty=2, col="grey")
    clip(0,1,0,1)
    x <- sqrt(dat$esComp[sel]) <= es
    horiz <- sum(x[!is.na(x)]) / length(x[!is.na(x)])
    text(x=.05, y=horiz+.02, labels=round(horiz, 2), cex=1.2)
    text(x=.9, y=h0horiz-.02, labels=round(h0horiz, 2), cex=1.2, col='darkgrey')
    clip(0, es, 0, horiz)
    abline(h=horiz, v=es, col="black", lty=2)
    clip(0, 1, 0, 1)
  }
}

for(i in 5:8){
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
       main=paste0(sort(unique(dat$journals.jour.))[i], ", D=", round(temp$statistic,3),", p<2.2*10^-16"),
       xlim=c(0,1),
       xaxs="i",
       yaxs="i",
       xlab="Correlation",
       ylab = "Cumulative density",
       cex.axis=.8,
       cex.lab=1,
       cex.main=1.5,
       col = "grey", las=1)
  lines(ecdf(sqrt(dat$esComp[sel])),
        lwd=.5)
  legend(x=.6,y=.2,legend=c(expression('H'[0]), 'Observed'),
         cex=1.2,lty=c(1,1),
         col = c("grey","black",2),box.lwd=0 ,lwd=2, bty='n')
  for(es in esR){
    h0horiz <- sum(sqrt(simNullEs$esComp[!is.na(simNullEs$esComp)]) < es) / length(simNullEs$esComp[!is.na(simNullEs$esComp)])
    clip(es, 1, 0, h0horiz)
    abline(h=h0horiz, v=es, lty=2, col="grey")
    clip(0,1,0,1)
    x <- sqrt(dat$esComp[sel]) <= es
    horiz <- sum(x[!is.na(x)]) / length(x[!is.na(x)])
    text(x=.05, y=horiz+.02, labels=round(horiz, 2), cex=1.2)
    text(x=.9, y=h0horiz-.02, labels=round(h0horiz, 2), cex=1.2, col='darkgrey')
    clip(0, es, 0, horiz)
    abline(h=horiz, v=es, col="black", lty=2)
    clip(0, 1, 0, 1)
  }
}
dev.off()


tempjour <- NULL
negjour <- NULL
kjour <- NULL
tempk <- NULL
for(i in 1:length(sort(unique(dat$journals.jour.[nsig])))){
  tempjour[i] <- sort(unique(dat$journals.jour.[nsig]))[i]
  
  sel <- nsig & dat$journals.jour.==sort(unique(dat$journals.jour.[nsig]))[i]
  
  # Number of papers containing negative results
  negjour[i] <- length(unique(dat$Source[sel]))
  
  # Minimum median k
  for(j in 1:negjour[i]){
    selJ <- nsig & dat$Source==unique(dat$Source[sel])[j]
    tempk[j] <- length(dat$Computed[selJ])
    print(j)
  }
  kjour[i] <- median(tempk)
  print(i)
}
write.csv2(cbind(tempjour, negjour, kjour), 'checks.csv')

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
# set.seed(35438759)
# source('c.Simulation/simCode.R')

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

# Agresti-Coull CI
.1-qnorm(.95,0,1)*(sqrt((1/10000)*.1*.9))
.1+qnorm(.95,0,1)*(sqrt((1/10000)*.1*.9))

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

years <- NULL
for(i in 1:length(sort(unique(dat$Source)))){
  sel <- dat$Source == sort(unique(dat$Source))[i]
  years[i] <- unique(dat$years.y.[sel])
  print(i)
}

# Make vector of sig/nsig/na 
fishDF <- data.frame(FisherP=fishRes$PFish, journal=sourcejour, kRes=fishRes$CountNSig, year=years)
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
  
  # Amount of papers in a journal without significant results
  countNA <- sum(is.na(fishDF$FisherP[sel]))
  journalSet <- NULL  
  # Writing out the results
  for(k in 1:length(kLen)){
    if(k == 7){
      x <- sum(fishDF$FisherP[sel & !is.na(fishDF$FisherP) & fishDF$kRes >= kLen[k]] < alphaF) / length(fishDF$FisherP[sel  & fishDF$kRes >= kLen[k]])
      journalSet <- cbind(journalSet, x)
    }
    else if(k == 5 | k ==6){
      x <- sum(fishDF$FisherP[sel & !is.na(fishDF$FisherP) & fishDF$kRes >= kLen[k]& fishDF$kRes < kLen[k+1]] < alphaF) / length(fishDF$FisherP[sel & fishDF$kRes >= kLen[k]& fishDF$kRes < kLen[k+1]])
      journalSet <- cbind(journalSet, x)
    } else{
      x <- sum(fishDF$FisherP[sel & !is.na(fishDF$FisherP) & fishDF$kRes == kLen[k]] < alphaF) / length(fishDF$FisherP[sel & fishDF$kRes == kLen[k]])
      journalSet <- cbind(journalSet, x)}
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
temp <- "Overall"
for(i in 1:length(kLen)){
  if(i == 7){
    temp[i+1] <- sum(fishDF$FisherP[fishDF$kRes >= kLen[i] & !is.na(fishDF$FisherP)] < alphaF )/ length(fishDF$FisherP[fishDF$kRes >= kLen[i]])
  }
  else if(i == 5|i == 6){
    temp[i+1] <- sum(fishDF$FisherP[fishDF$kRes >= kLen[i] &fishDF$kRes < kLen[i+1] & !is.na(fishDF$FisherP)] < alphaF )/ length(fishDF$FisherP[fishDF$kRes >= kLen[i] & fishDF$kRes < kLen[i+1]])
  }
  else{
    temp[i+1] <- sum(fishDF$FisherP[fishDF$kRes == kLen[i] & !is.na(fishDF$FisherP)] < alphaF )/ length(fishDF$FisherP[fishDF$kRes == kLen[i]])  
  }
  
}
temp <- c(temp,
          sum(as.numeric(as.character(final[,dim(final)[2]-1])))/sum(as.numeric(as.character(final[,dim(final)[2]]))),
          sum(as.numeric(as.character(final[,dim(final)[2]-2]))),
          sum(as.numeric(as.character(final[,dim(final)[2]-1]))),
          sum(as.numeric(as.character(final[,dim(final)[2]]))))
final <- rbind(as.character(temp), final)
final <- as.data.frame(final)
names(final) <- c('journals', paste0('k', kLen), 'overall', 'countNA', 'amountSig', 'nrpapers')
# write.csv2(final, '../Writing/Tables/table4.csv', row.names=F)

# Computing the number of significant Fisher results per year
# As proportion of all papers reporting nonsignificant results
library(plyr)

fishDF$logicalP <- ifelse(fishDF$FisherP<.1, 1, 0)
fisherYear <- ddply(fishDF, .(year), summarise, propYear=mean(logicalP, na.rm=TRUE)
)

knsYear <- ddply(fishDF, .(year), summarise, kYear=mean(kRes, na.rm=TRUE)
)

library(ggplot2)

mydf <- data.frame(x = fisherYear$year,
                   y = fisherYear$propYear,
                   count = knsYear$kYear)

ggplot(mydf, aes(x = x, y = y)) + geom_point(aes(size = count)) + ylim(0, 1) + geom_smooth(method="lm") +
  xlab("Year") + ylab("Proportion significant Fisher results")

ggsave(filename = 'fisheryears.png', plot = last_plot(), width = 21, height = 9)

# tiff('../Writing/Figures/falseneg.tiff', width=1200, height=900)
# plot(fisherYear, ylim=c(0,1), type='o', ylab="Proportion",
#      main="False negatives", xlab="Year")
# abline(lm(fisherYear$propYear~fisherYear$year))
# dev.off()
# symbols(x=fisherYear$year, y=fisherYear$propYear, circles=knsYear$kYear, inches=1/4, ann=F, bg="steelblue2", fg=NULL, ylim=c(0,1))
  # to add the frequencies of Ks inspect tabular format of k per journal
# Overall
table(fishDF$kRes)
# Per journal
table(fishDF$kRes[fishDF$journal=="DP"])
table(fishDF$kRes[fishDF$journal=="FP"])
table(fishDF$kRes[fishDF$journal=="JAP"])
table(fishDF$kRes[fishDF$journal=="JCCP"])
table(fishDF$kRes[fishDF$journal=="JEPG"])
table(fishDF$kRes[fishDF$journal=="JPSP"])
table(fishDF$kRes[fishDF$journal=="PLOS"])
table(fishDF$kRes[fishDF$journal=="PS"])

medianN <- NULL
p25 <- NULL
p75 <- NULL
i <- 1
for(y in 1985:2013){
  temp <- summary(dat$df2[dat$years.y. == y])
  
  medianN[i] <- temp[3]
  p25[i] <- temp[2]
  p75[i] <- temp[5]
  
  i <- i + 1
}

tiff('../Writing/Figures/fig7.tiff', width=1200, height=900)
plot(x=1985:2013, y=medianN, type='o', col="black",
     ylab="N", xlab="Year", ylim=c(0,80), cex.lab=1.2, las=1, lwd=1, cex.axis=1.2, xaxs='i')
# lines(x=1985:2013, y=p25, type='o', col='grey')
# lines(x=1985:2013, y=p75, type='o', col='grey')
dev.off()


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

# Null model
nullmod <- lm(datFit$yi ~ 1)

# Optimal curve model searching
# [datFit$kRes <= 20]
r2fit = NULL
ktemp <- seq(1:max(k))
datFit <- data.frame(yi=fishDF$sig, kRes=fishDF$kRes, jour=fishDF$journal)
xi <- list(NULL)
for(i in 1:length(curveES)){
  xi[[i]] <- (1 - pnorm(zcv/sqrt(datFit$kRes), curveES[i], 1 / sqrt(datFit$kRes)))
  datFit <- cbind(datFit, xi[[i]])
  r2fit[i] <- summary(lm(datFit$yi ~ 1 + datFit[,3+i]))$r.squared
}
# Optimal curve fit
xi <- (1 - pnorm(zcv/sqrt(datFit$kRes), curveES[which(r2fit == max(r2fit))], 1 / sqrt(datFit$kRes)))
curvemod <- lm(datFit$yi ~ 1 + xi)
curveES[which(r2fit == max(r2fit))]


# Saturated model
satur <- summary(lm(datFit$yi ~ 0 + as.factor(datFit$kRes)))

tiff('../Writing/Figures/fig8.tiff', width=717, height=654)
par(mfrow=c(1,1), mai=c(1.2,1.2,.2,.2))
plot(x=(datFit$kRes), y=datFit$yi, col= 'white',
     xlab="Ln(k)", ylab="Estimated nr. of significant Fisher tests", xaxs="i",
     cex.axis=1.2,
     cex.lab=1.2,
     las=1, lwd=1)
curve((1 - pnorm(zcv/sqrt(x), curveES[which(r2fit == max(r2fit))], 1 / sqrt(x))), from=1, to=max(datFit$kRes),
      add=T)
lines(x=(sort(unique(datFit$kRes))), y=satur$coefficients[,1],col='black', lty=2)

for(z in 1:length(unique(datFit$jour))){
  points(x=kMean[z], y=propMean[z], col="black", pch=z)
}

abline(h=nullmod$coefficients[1], col="grey", lty=1)
legend(x=3,y=.25,legend=c(as.character(sort(unique(datFit$jour))), "Null, R2=0",
                          paste0("Curve, R2=", round(summary(curvemod)$r.squared,3)),
                          paste0("Saturated, R2=", round(satur$r.squared,3))),
       cex=.8, pch=c(1:8, NA, NA, NA), lty=c(rep(NA, 8), 1, 1, 2), 
       col=c(rep("black",8),"grey", rep("black", 2)),
       #        col = c("black", rep("blue", 5), rep("red", 3))
       box.lwd=0 ,lwd=1, bty='n', y.intersp=1)
dev.off()

# Observed true positive
# Overall
selK <- sort(unique(fishDF$kRes))
obstruesig <- NULL
obspow <- NULL
obssig <- NULL
obstruesigSat <- NULL
i <- 1
for(k in as.numeric(selK)){
  obspow[i] <- (1 - pnorm(zcv/sqrt(k), curveES[which(r2fit == max(r2fit))], 1/sqrt(k)))
  obssig[i] <- sum(fishDF$FisherP[!is.na(fishDF$FisherP) & fishDF$kRes == k] < alphaF)
  obstruesig[i] <- (obssig[i]*obspow[i])
  obstruesigSat[i] <- obssig[i]*satur$coefficients[i]
  i <- i + 1
}
sum(obssig)*sum(obstruesig)/length(fishDF$FisherP)
sum(obssig)*sum(obstruesigSat)/length(fishDF$FisherP)
sum(obssig)

sum(obstruesig)/length(fishDF$FisherP)
sum(obstruesigSat)/length(fishDF$FisherP)

j <- 1
low <- NULL
high <- NULL
for(y in 1985:2013){
  if(sort(unique(fishDF$kRes[fishDF$year == y]))[1] == 0){
    selK <- sort(unique(fishDF$kRes[fishDF$year == y]))[-1]
  }else{
    selK <- sort(unique(fishDF$kRes[fishDF$year == y]))}
  obstruesig <- NULL
  obspow <- NULL
  obssig <- NULL
  obstruesigSat <- NULL
  i <- 1
  for(k in as.numeric(selK)){
    obspow[i] <- (1 - pnorm(zcv/sqrt(k), curveES[which(r2fit == max(r2fit))], 1/sqrt(k)))
    obssig[i] <- sum(fishDF$FisherP[!is.na(fishDF$FisherP) & fishDF$kRes == k & fishDF$year == y] < alphaF)
    obstruesig[i] <- (obssig[i]*obspow[i])
    obstruesigSat[i] <- obssig[i]*satur$coefficients[i]
    i <- i + 1
  }
  low[j] <- sum(obstruesig)/length(fishDF$FisherP[fishDF$year == y])
  high[j] <-sum(obstruesigSat)/length(fishDF$FisherP[fishDF$year == y])
  
  j <- j + 1
}
year <- 1985:2013
temp <- cbind(year, low, high)
temp <- as.data.frame(temp)
names(temp) <- c("Year", 'Low [Curve]', 'High [Saturated]')

write.csv2(temp, '../Writing/Tables/table6notused.csv', row.names=F)

tiff('../Writing/Figures/fig9.tiff', width=650, height=576)
par(mfrow=c(1,1), mai=c(1.2,1.2,.2,.2))
plot(x=1985:2013, y=low, col= 'black',
     xlab="Year", ylab="Estimated false negative rate", xaxs="i",
     cex.axis=1,
     cex.lab=1,
     las=1, lwd=1,type='o',ylim=c(0,1))
lines(x=1985:2013, y=high, ylim=c(0,1), type='o', lty=2)
legend(x=2005,y=.1,legend=c("Lowerbound", "Upperbound"),
       cex=.8, lty=c(1,2), 
       box.lwd=0 ,lwd=2, bty='n', y.intersp=1)
dev.off()

#########
# Cases #
#########
setwd("D:/files/phd/toogoodtobefalse/JPSP/selected")

library(statcheck)

cases <- checkHTMLdir("D:/files/phd/toogoodtobefalse/JPSP/selected")
cases <- cases[,1:10]
write.csv2(cases, 'cases2.csv')

FisherMethod(cases$Computed, id = cases$Source, alpha=.05)

callanN <- c(64,
             218,
             83,
             190,
             59,
             367
             #              ,85,
             #              142,
             #              139,
             #              77+103
)
sum(callanN)
mean(callanN)

jungN <- c(151,
           152,
           304,
           132,
           294,
           193,
           198,
           329,
           419,
           835,
           1065)
sum(jungN)
mean(jungN)

#################
# Gender effect #
#################
setwd("D:/files/phd/toogoodtobefalse/")

genderdat <- read.csv2('datafilegender100.csv')
# For consistency, select only t F and r values
genderdat <- genderdat[genderdat$Statistic == 't' | genderdat$Statistic == 'F' | genderdat$Statistic == 'r',]
gendersample <- genderdat[genderdat$gender == TRUE,]

set.seed(123)
sampled <- sample(x = unique(gendersample$Source), size = 100, replace = FALSE)

genderrandomsample <- gendersample[gendersample$Source %in% sampled,]

write.csv2(genderrandomsample,
           'genderrandomselect.csv')

# Discussion
require(car)
iccSS <- Anova(lm(dat$Computed[nsig] ~ dat$Source[nsig]), type="III")
# Computes the ICC
iccSS$Sum[2]/(iccSS$Sum[3]+iccSS$Sum[2])


# Table 3 power computations
ser <- 1/sqrt(c(33, 62, 119)-3)
rho <- .1
zcv <- 1.282
rcv <- (exp(2*(zcv*ser))-1)/(exp(2*(zcv*ser))+1)
zrcv <- .5*log((1+rcv)/(1-rcv))
zrho <- .5*log((1+rho)/(1-rho))
round(1-pnorm(zrcv, mean=zrho, sd=ser),4)

rho <- .25
rcv <- (exp(2*(zcv*ser))-1)/(exp(2*(zcv*ser))+1)
zrcv <- .5*log((1+rcv)/(1-rcv))
zrho <- .5*log((1+rho)/(1-rho))
round(1-pnorm(zrcv, mean=zrho, sd=ser),4)


# 
# chjh
# mva 
# jmw
# 
# set1 <- as.dataframe(chjh = chjh[1:90], jmw = jmw)
# cohen.kappa(rbind(set1$chjh, set1$jmw), w=NULL,n.obs=NULL,alpha=.05)  
# set1 <- as.dataframe(chjh = chjh[91:180], mva = mva)
# cohen.kappa(rbind(set1$chjh, set1$mva), w=NULL,n.obs=NULL,alpha=.05) 

# save <- NULL
# 
# for(journal in unique(fishDF$journal)){
#   sel <- fishDF$kRes > 0 & fishDF$journal == journal
#   
#   save <- c(save, mean(fishDF$kRes[sel]))
#   
# }
# 
# cbind(as.character(unique(fishDF$journal)), save)
