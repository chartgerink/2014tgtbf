# Code written by CHJ Hartgerink
# Checked by: -

# set to run in documents folder
setwd(normalizePath('~'))

# load packages -----------------------------------------------------------

if(!require(httr)){install.packages('httr')}
library(httr)
if(!require(latex2exp)){install.packages('latex2exp')}
library(latex2exp)
if(!require(plyr)){install.packages('plyr')}
library(plyr)
if(!require(ggplot2)){install.packages('ggplot2')}
library(ggplot2)
if(!require(stringr)){install.packages('stringr')}
library(stringr)
if(!require(car)){install.packages('car')}
library(car)

# download and load dependent files ---------------------------------------

# custom functions
# individual files available in `functions/`
GET('https://github.com/chartgerink/2014tgtbf/raw/master/functions/functions.R',
    write_disk('functions.R', overwrite = TRUE))
source("functions.R")

# simulation study code
GET('https://github.com/chartgerink/2014tgtbf/raw/master/Analyzing/simCode.R',
    write_disk('simCode.R', overwrite = TRUE))

# data
# literally gets the data from the Nuijten et al paper
GET('https://osf.io/gdr4q/?action=download',
    write_disk('statcheck_full_anonymized.csv', overwrite = TRUE))
dat <- read.csv2("statcheck_full_anonymized.csv",
                 stringsAsFactors = F,
                 dec = ",",
                 sep = ";")[-1]

# the gender related results (ALL)
GET('https://raw.githubusercontent.com/chartgerink/2014tgtbf/master/data/datafilegender500_post.csv',
    write_disk('datafilegender500_post.csv', overwrite = TRUE))
# the coded gender results
GET('https://raw.githubusercontent.com/chartgerink/2014tgtbf/master/data/gendercoded%20cleaned%20and%20discussed.csv',
    write_disk('gendercoded cleaned and discussed.csv', overwrite = TRUE))



# data cleaning -----------------------------------------------------------

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

save(dat)

# methods section ---------------------------------------------------------

# Descriptives dataset
journals <- sort(unique(dat$journals.jour.))
for(j in 1:length(journals)){
  selJournal <- dat$journals.jour. == journals[j]
  meanK <- mean(table(dat$Source[selJournal]))
  len <- length(dat$Computed[selJournal & !is.na(dat$Computed)])
  sigRes <- sum(dat$Computed[selJournal & !is.na(dat$Computed)] < .05)
  print(
    paste0(
      journals[j], "\t", len, "\t", round(meanK, 1), "\t", sigRes, "\t", len-sigRes
    )
  )
}

# Simulation sample sizes
c(as.numeric(summary(dat$df2)[2]), # 25th percentile
  as.numeric(summary(dat$df2)[3]), # 50th percentile
  as.numeric(summary(dat$df2)[5])) # 75th percentile

# 27,523 gender results
# manually extracted from datafilegender500_post.csv
# number of rows (take into account header row)

# results section ---------------------------------------------------------

# nr of results -----------------------------------------------------------

dim(dat)[1]

# nonsignificant results --------------------------------------------------

dim(dat[dat$Computed > .05, ])

# density plot effects ----------------------------------------------------

pdf('Fig3.pdf',width=7, height=8)
par(mai = c(1, 1, 0, .2))

plot(density(dat$esComp[!is.na(dat$esComp)]),
     lty = 1,
     frame.plot = T, 
     main = "",
     xlim = c(0, 1),
     xaxs = "i",
     xlab = latex2exp("Correlation ($\\eta$)"),
     ylab = "Density",
     cex.axis = .8,
     cex.lab = 1,
     col = "black", las = 1, bty = 'n')
abline(v = c(.1, .25, .4), lty = 2, col = "grey")
t1 <- sum(dat$esComp[!is.na(dat$esComp)] < .1) / length(dat$esComp[!is.na(dat$esComp)])
t2 <- sum(dat$esComp[!is.na(dat$esComp)] >= .1 & dat$esComp[!is.na(dat$esComp)] < .25) / length(dat$esComp[!is.na(dat$esComp)])
t3 <- sum(dat$esComp[!is.na(dat$esComp)] >= .25 & dat$esComp[!is.na(dat$esComp)] < .4) / length(dat$esComp[!is.na(dat$esComp)])
t4 <- sum(dat$esComp[!is.na(dat$esComp)] >= .4) / length(dat$esComp[!is.na(dat$esComp)])
text(x = .1 / 2, y = .15, labels = round(t1, 2), cex = 1)
text(x = ((.25 - .1) / 2) + .1, y = .15, labels = round(t2, 2), cex = 1)
text(x = ((.4 - .25) / 2) + .25, y = .15, labels = round(t3, 2), cex = 1)
text(x = ((1 - .4) / 2) + .4, y = .15, labels = round(t4, 2), cex = 1)
dev.off()

# proportion nonsignificant results over time -----------------------------

pdf(file = 'Fig4.pdf', width = 7, height = 5)
par(mfrow = c(1, 1), mai = c(1, 1, .2, .2))
plot(x = 1985:2013,
     nsigtemp,
     ylim = c(0, .4),
     type = 'o',
     xlab = "Year", ylab = 'Proportion nonsignificant',
     yaxs = 'i',
     cex.axis = 1, las = 1, bty = 'n')
abline(lm(nsigtemp ~ c(1985:2013)), lty = 2)
dev.off()

# Comparison of observed- and expected effect distributions of non --------

# takes some time to run these, have to simulate many values
# go get coffee
set.seed(1234)
simNullEs <- simNullDist(dat, n.iter = length(dat$esComp[nsig]) * 3, alpha = .05)
simNullEs$adjESComp[simNullEs$adjESComp < 0] <- 0

temp <- ks.test(simNullEs$esComp,
                dat$esComp[nsig],
                alternative="greater")

pdf('Fig5.pdf', width=11, height=7)
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

pdf('S1Fig.pdf', onefile = TRUE, width = 11, height = 11)
par(mfrow = c(2, 2), mai = c(1.2, 1.2, .8, .5))

for(i in 1:4){
  sel <- dat$journals.jour. == sort(unique(dat$journals.jour.))[i] & nsig
  set.seed(i)
  simNullEs <- simNullDist(dat, n.iter=length(dat$esComp[sel])*3, alpha = .05)
  temp <- ks.test(simNullEs$esComp,
                  dat$esComp[dat$journals.jour. == sort(unique(dat$journals.jour.))[i] & nsig],
                  alternative = "greater")
  print((temp))
  if(i == 3){
    plot(ecdf(na.omit(sqrt(simNullEs$esComp))),
         lty = 1,
         frame.plot = T, 
         main = latex2exp(
           sprintf("%s, D=%s, $p=7.934\\times 10^{-6}$",
                   sort(unique(dat$journals.jour.))[i],
                   round(temp$statistic,3))),
         xlim = c(0,1),
         xaxs = "i",
         yaxs = "i",
         xlab = latex2exp("Correlation ($\\eta$)"),
         ylab = "Cumulative density",
         cex.axis = .8,
         cex.lab = 1,
         cex.main = 1.5,
         col = "grey", las = 1)
  }
  else{plot(ecdf(na.omit(sqrt(simNullEs$esComp))),
            lty = 1,
            frame.plot = T, 
            main = latex2exp(
              sprintf("%s, D=%s, $p=2.2\\times 10^{-16}$",
                      sort(unique(dat$journals.jour.))[i],
                      round(temp$statistic,3))),
            xlim = c(0,1),
            xaxs = "i",
            yaxs = "i",
            xlab = latex2exp("Correlation ($\\eta$)"),
            ylab = "Cumulative density",
            cex.axis = .8,
            cex.lab = 1,
            cex.main = 1.5,
            col = "grey", las = 1)}
  lines(ecdf(sqrt(dat$esComp[sel])),
        lwd = .5)
  legend(x = .6, y = .2, legend = c(latex2exp("$H_0$"), 'Observed'),
         cex = 1.2, lty = c(1, 1),
         col = c("grey", "black", 2), box.lwd = 0, lwd = 2, bty = 'n')
  for(es in esR){
    h0horiz <- sum(sqrt(simNullEs$esComp[!is.na(simNullEs$esComp)]) < es) / length(simNullEs$esComp[!is.na(simNullEs$esComp)])
    clip(es, 1, 0, h0horiz)
    abline(h = h0horiz, v = es, lty = 2, col = "grey")
    clip(0, 1, 0, 1)
    x <- sqrt(dat$esComp[sel]) <= es
    horiz <- sum(x[!is.na(x)]) / length(x[!is.na(x)])
    text(x = .05, y = horiz + .02, labels = round(horiz, 2), cex = 1.2)
    text(x = .9, y = h0horiz - .02, labels = round(h0horiz, 2), cex = 1.2, col = 'darkgrey')
    clip(0, es, 0, horiz)
    abline(h = horiz, v = es, col = "black", lty = 2)
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
       lty = 1,
       frame.plot = T, 
       main = latex2exp(
         sprintf("%s, D=%s, $p=2.2\\times 10^{-16}$",
                 sort(unique(dat$journals.jour.))[i],
                 round(temp$statistic,3))),
       xlim = c(0, 1),
       xaxs = "i",
       yaxs = "i",
       xlab = latex2exp("Correlation ($\\eta$)"),
       ylab = "Cumulative density",
       cex.axis = .8,
       cex.lab = 1,
       cex.main = 1.5,
       col = "grey", las = 1)
  lines(ecdf(sqrt(dat$esComp[sel])),
        lwd = .5)
  legend(x = .6,y = .2,legend = c(latex2exp("$H_0$"), 'Observed'),
         cex = 1.2,lty = c(1, 1),
         col = c("grey", "black", 2), box.lwd = 0 , lwd = 2, bty = 'n')
  for(es in esR){
    h0horiz <- sum(sqrt(simNullEs$esComp[!is.na(simNullEs$esComp)]) < es) / length(simNullEs$esComp[!is.na(simNullEs$esComp)])
    clip(es, 1, 0, h0horiz)
    abline(h = h0horiz, v = es, lty = 2, col = "grey")
    clip(0, 1, 0, 1)
    x <- sqrt(dat$esComp[sel]) <= es
    horiz <- sum(x[!is.na(x)]) / length(x[!is.na(x)])
    text(x = .05, y = horiz + .02, labels = round(horiz, 2), cex = 1.2)
    text(x = .9, y = h0horiz - .02, labels = round(h0horiz, 2), cex = 1.2, col = 'darkgrey')
    clip(0, es, 0, horiz)
    abline(h = horiz, v = es, col = "black", lty = 2)
    clip(0, 1, 0, 1)
  }
}
dev.off()

# Statistical properties of the Fisher method -----------------------------

# condition setting
N <- c(as.numeric(summary(dat$df2)[2]), # 25th percentile
       as.numeric(summary(dat$df2)[3]), # 50th percentile
       as.numeric(summary(dat$df2)[5]) # 75th percentile
)

ES <- c(.00,
        seq(.01, .99, .01))

P <- c(seq(1, 10, 1), seq(15, 50, 5))

alpha <- .05
alphaF <- 0.10
n.iter <- 10000

# NOTE: runs upon source and can take a while
set.seed(35438759)
source('simCode.R')
# Load all results back in
files <- list.files()[grepl("N_", list.files())]

names <- str_sub(files,start=1L, end=-5L)
for(i in 1:length(files)){
  assign(x = names[i], read.csv(files[i]))
  assign(x = names[i], t(get(x = names[i])[ ,-1]))
}

# Data for table
# rows are k
# columns are effect size
# N = 33
t(get(x = names[2]))
# N = 62
t(get(x = names[3]))
# N = 119
t(get(x = names[1]))

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


# Agresti-Coull CI
.1 - qnorm(.95, 0, 1) * (sqrt((1/10000) * .1 * .9))
.1 + qnorm(.95, 0, 1) * (sqrt((1/10000) * .1 * .9))

# False negative nonsignificant results detected with the Fisher m --------

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
save(fishDF, file = 'figures/Fig6')

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

names(final) <- c('journals', paste0('k', kLen), 'overall', 'countNA', 'amountSig', 'nrpapers', 'perc_k1')
write.csv(final, 'table4.csv', row.names=F)

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

# Median per journal
x = as.numeric(as.character(final$amountSig)) / (as.numeric(as.character(final$nrpapers)) - as.numeric(as.character(final$countNA)))
y = c(mean(fishDF$kRes[fishDF$journal=="DP"]),
      median(fishDF$kRes[fishDF$journal=="FP"]),
      median(fishDF$kRes[fishDF$journal=="JAP"]),
      median(fishDF$kRes[fishDF$journal=="JCCP"]),
      median(fishDF$kRes[fishDF$journal=="JEPG"]),
      median(fishDF$kRes[fishDF$journal=="JPSP"]),
      median(fishDF$kRes[fishDF$journal=="PLOS"]),
      median(fishDF$kRes[fishDF$journal=="PS"]))
cor(x[-1], y)

# Computing the number of significant Fisher results per year
# As proportion of all papers reporting nonsignificant results
fishDF$logicalP <- ifelse(fishDF$FisherP < .1, 1, 0)
fisherYear <- ddply(fishDF, .(year), summarise, propYear=mean(logicalP, na.rm=TRUE))

knsYear <- ddply(fishDF, .(year), summarise, kYear=mean(kRes, na.rm=TRUE))

mydf <- data.frame(x = fisherYear$year,
                   y = fisherYear$propYear,
                   count = knsYear$kYear)

mydf

ggplot(mydf, aes(x = x, y = y)) + geom_point(aes(size = count)) + ylim(0, 1) + geom_smooth(method="lm") +
  xlab("Year") + ylab("Proportion significant Fisher results")

ggsave(filename = 'Fig6.pdf', plot = last_plot(), width = 7, height = 5)

# sample size development over time
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

pdf('Fig7.pdf', width=7, height=6)
par(mar = c(4, 4, .2, 2))
plot(x=1985:2013, y=medianN, type='o', col="black",
     ylab="N", xlab="Year", ylim=c(0,150), cex.lab=1.2, las=1, lwd=1, cex.axis=1.2, xaxs='i', bty = 'n')
lines(x=1985:2013, y=p25, type='o', col='grey')
lines(x=1985:2013, y=p75, type='o', col='grey')
text(y = 100, x = 2010, "P75", col = "grey")
text(y = 50, x = 2010, "P50", col = "black")
text(y = 22, x = 2010, "P25", col = "grey")
dev.off()

# Effect PDF
pdf('Fig8.pdf',width=7, height=8)
par(mai = c(1, 1, 0, .2))

plot(density(dat$esComp[!is.na(dat$esComp) & dat$years.y. == 1985]),
     lty = 1,
     frame.plot = T, 
     main = "",
     xlim = c(0, 1),
     xaxs = "i",
     xlab = latex2exp("Correlation ($\\eta$)"),
     ylab = "Density",
     cex.axis = .8,
     cex.lab = 1,
     col = "black", las = 1, bty = 'n')
lines(density(dat$esComp[!is.na(dat$esComp) & dat$years.y. == 2013]), col = "darkgrey")
abline(v = c(.1, .25, .4), lty = 2, col = "grey")
t11 <- sum(dat$esComp[!is.na(dat$esComp) & dat$years.y == 1985] < .1) / length(dat$esComp[!is.na(dat$esComp) & dat$years.y == 1985])
t21 <- sum(dat$esComp[!is.na(dat$esComp) & dat$years.y == 1985] >= .1 & dat$esComp[!is.na(dat$esComp) & dat$years.y == 1985] < .25) / length(dat$esComp[!is.na(dat$esComp) & dat$years.y == 1985])
t31 <- sum(dat$esComp[!is.na(dat$esComp) & dat$years.y == 1985] >= .25 & dat$esComp[!is.na(dat$esComp) & dat$years.y == 1985] < .4) / length(dat$esComp[!is.na(dat$esComp) & dat$years.y == 1985])
t41 <- sum(dat$esComp[!is.na(dat$esComp) & dat$years.y == 1985] >= .4) / length(dat$esComp[!is.na(dat$esComp) & dat$years.y == 1985])
t12 <- sum(dat$esComp[!is.na(dat$esComp) & dat$years.y == 2013] < .1) / length(dat$esComp[!is.na(dat$esComp) & dat$years.y == 2013])
t22 <- sum(dat$esComp[!is.na(dat$esComp) & dat$years.y == 2013] >= .1 & dat$esComp[!is.na(dat$esComp) & dat$years.y == 2013] < .25) / length(dat$esComp[!is.na(dat$esComp) & dat$years.y == 2013])
t32 <- sum(dat$esComp[!is.na(dat$esComp) & dat$years.y == 2013] >= .25 & dat$esComp[!is.na(dat$esComp) & dat$years.y == 2013] < .4) / length(dat$esComp[!is.na(dat$esComp) & dat$years.y == 2013])
t42 <- sum(dat$esComp[!is.na(dat$esComp) & dat$years.y == 2013] >= .4) / length(dat$esComp[!is.na(dat$esComp) & dat$years.y == 2013])
text(x = .1 / 2, y = .15, labels = round(t11, 2), cex = 1)
text(x = .1 / 2, y = .35, labels = round(t12, 2), cex = 1, col = 'darkgrey')
text(x = ((.25 - .1) / 2) + .1, y = .15, labels = round(t21, 2), cex = 1)
text(x = ((.25 - .1) / 2) + .1, y = .35, labels = round(t22, 2), cex = 1, col = 'darkgrey')
text(x = ((.4 - .25) / 2) + .25, y = .15, labels = round(t31, 2), cex = 1)
text(x = ((.4 - .25) / 2) + .25, y = .35, labels = round(t32, 2), cex = 1, col = 'darkgrey')
text(x = .45, y = .15, labels = round(t41, 2), cex = 1)
text(x = .45, y = .35, labels = round(t42, 2), cex = 1, col = 'darkgrey')
legend(x = .6, y = 1, legend = c("1985", '2013'),
       cex = 1, lty = c(1, 1),
       col = c("black", "darkgrey"),
       box.lwd = 0, lwd = 2, bty = 'n')
dev.off()


# gender effect evaluated -------------------------------------------------

gend <- read.csv2('gendercoded cleaned and discussed.csv', sep = ";", dec = ".", header = TRUE)

# 1 = null expected
# 2 = effect expected
# 3 = no expectation
table(gend$significance, gend$final_code)

options(scipen = 5, digits = 8)
# sel <- gend$significance == "significant" & gend$final_code == "no expectation"
for(sig in unique(gend$significance)){
  for(code in unique(gend$final_code[!is.na(gend$final_code)])){
    sel <- gend$significance == sig & gend$final_code == code
    
    if(sig == "significant"){
      temp <- gend$Computed[sel] / .05
      pstar <- temp[!is.na(temp)]
      x$CountNSig <- length(pstar)
      x$Fish <- -2*sum(log(pstar))
      x$PFish <- pchisq(x$Fish, df = 2 * length(pstar), lower.tail = FALSE)
    }
    if(sig == "nonsignificant"){
      x <- (FisherMethod(x = gend$Computed[sel], id = 1, alpha = 0.05))
    }
    
    cat(sprintf("For %s %s, k = %s, chi2 = %s, p = %s\n", sig, code, x$CountNSig, x$Fish, x$PFish))
  }
}


# discussion --------------------------------------------------------------

iccSS <- Anova(lm(dat$Computed[nsig] ~ dat$Source[nsig]), type="III")
# Computes the ICC
iccSS$Sum[2]/(iccSS$Sum[3]+iccSS$Sum[2])
