# Code written by CHJ Hartgerink
# Checked by: -

# set to run in documents folder
setwd(normalizePath('~'))

# load packages -----------------------------------------------------------

if(!require(httr)){install.packages('httr')}
library(httr)
if(!require(TeX)){install.packages('TeX')}
library(TeX)
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

save(dat, file = 'dat')


# theoretical framework ---------------------------------------------------

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


# Application 1 -----------------------------------------------------------

# methods app1 ---------------------------------------------------------
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

# results app1 ------------------------------------------------------------
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
     xlab = TeX("Correlation (|$\\eta$|)"),
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
     main = TeX(
       sprintf(
         "Unadjusted, $D=%s,p<2.2\\times 10^{-16}$", round(temp$statistic,2))),
     xlim = c(0, 1),
     xaxs = "i",
     yaxs = "i",
     xlab = TeX("Correlation (|$\\eta$|)"),
     ylab = "Cumulative density",
     cex.axis=.8,
     cex.lab=1,
     cex.main=1.5,
     col = "grey", las=1)
lines(ecdf(na.omit(sqrt(dat$esComp[nsig]))))
legend(x = .6, y = .2, legend = c(TeX("$H_0$"), 'Observed'),
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
     main = TeX(
       sprintf(
         "Adjusted, $D=%s,p<2.2\\times 10^{-16}$", round(temp$statistic,2))),
     xlim = c(0, 1),
     xaxs = "i",
     yaxs = "i",
     xlab = TeX("Correlation (|$\\eta$|)"),
     ylab = "Cumulative density",
     cex.axis = .8,
     cex.lab = 1,
     cex.main = 1.5,
     col = "grey", las = 1)
lines(ecdf(na.omit(sqrt(dat$adjESComp[nsig]))))
legend(x = .6,y = .2,legend = c(TeX("$H_0$"), 'Observed'),
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
         main = TeX(
           sprintf("%s, D=%s, $p=7.934\\times 10^{-6}$",
                   sort(unique(dat$journals.jour.))[i],
                   round(temp$statistic,3))),
         xlim = c(0,1),
         xaxs = "i",
         yaxs = "i",
         xlab = TeX("Correlation (|$\\eta$|)"),
         ylab = "Cumulative density",
         cex.axis = .8,
         cex.lab = 1,
         cex.main = 1.5,
         col = "grey", las = 1)
  }
  else{plot(ecdf(na.omit(sqrt(simNullEs$esComp))),
            lty = 1,
            frame.plot = T, 
            main = TeX(
              sprintf("%s, D=%s, $p=2.2\\times 10^{-16}$",
                      sort(unique(dat$journals.jour.))[i],
                      round(temp$statistic,3))),
            xlim = c(0,1),
            xaxs = "i",
            yaxs = "i",
            xlab = TeX("Correlation (|$\\eta$|)"),
            ylab = "Cumulative density",
            cex.axis = .8,
            cex.lab = 1,
            cex.main = 1.5,
            col = "grey", las = 1)}
  lines(ecdf(sqrt(dat$esComp[sel])),
        lwd = .5)
  legend(x = .6, y = .2, legend = c(TeX("$H_0$"), 'Observed'),
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
       main = TeX(
         sprintf("%s, D=%s, $p=2.2\\times 10^{-16}$",
                 sort(unique(dat$journals.jour.))[i],
                 round(temp$statistic,3))),
       xlim = c(0, 1),
       xaxs = "i",
       yaxs = "i",
       xlab = TeX("Correlation (|$\\eta$|)"),
       ylab = "Cumulative density",
       cex.axis = .8,
       cex.lab = 1,
       cex.main = 1.5,
       col = "grey", las = 1)
  lines(ecdf(sqrt(dat$esComp[sel])),
        lwd = .5)
  legend(x = .6,y = .2,legend = c(TeX("$H_0$"), 'Observed'),
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
y = c(median(fishDF$kRes[fishDF$journal=="DP"]),
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
     xlab = TeX("Correlation (|$\\eta$|)"),
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


# Application 2 -----------------------------------------------------------

# 27,523 gender results
# manually extracted from datafilegender500_post.csv
# number of rows (take into account header row)

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

# Application 3 -----------------------------------------------------------

# Read in Tilburg data
info <- GET('https://osf.io/fgjvw/?action=download', write_disk('rpp_data.csv', overwrite = TRUE)) #downloads data file from the OSF
MASTER <- read.csv("rpp_data.csv")[1:167, ]
colnames(MASTER)[1] <- "ID" # Change first column name to ID to be able to load .csv file

# number of iterations per condition
# if true probability is .025 or .975, then the SE = .00494 for n.iter = 1,000
n.iter <- 10000
# alpha of the individual test
alpha <- .05
# vector of ES
es <- c(.1, .3, .5)

# set seed
set.seed(alpha * n.iter * sum(es))

# select only the available and nonsignificant replications
# call this M
M <- MASTER[MASTER$T_sign_R == 0 
            & !is.na(MASTER$T_sign_R) 
            & (MASTER$T_Test.Statistic..R. == 't'
               | MASTER$T_Test.Statistic..R. == 'F'
               | MASTER$T_Test.Statistic..R. == 'r'
               | MASTER$T_Test.Statistic..R. == 'z'
               | MASTER$T_Test.Statistic..R. == 'Chi2'), ]

Fish_M <- FisherMethod(x = M$T_pval_USE..R.,
                       id = 1,
                       alpha = .05)

# selectors for easy selection later on  
selr <- M$T_Test.Statistic..R. == 'r'
selz <- M$T_Test.Statistic..R. == 'z'
selChi2 <- M$T_Test.Statistic..R. == 'Chi2'

# ensure that we can easily compute the effect size r
M$T_df1..R.[selr] <- 1
M$T_df2..R.[selr] <- M$T_N..R.[selr] - 2  
# ensure z is tranformed to fisher z
M$sigma[selz] <- sqrt(1 / (M$T_N..R.[selz] - 3))

# object to write count to
res <- list(NULL)
res[[1]] <- rep(0, dim(M)[1] + 1)
res[[2]] <- rep(0, dim(M)[1] + 1)
res[[3]] <- rep(0, dim(M)[1] + 1)

# compute critical value for F and t values
# we assume two-tailed tests everywhere
fCV <- qf(p = alpha, df1 = M$T_df1..R., df2 = M$T_df2..R., lower.tail = FALSE)
# compute critical value for z values
zCV <- qnorm(p = .025, mean = 0, sd = M$sigma[selz], lower.tail = FALSE)  
# compute critical value for chi2 values
Chi2CV <- qchisq(p = .05, df = M$T_df1..R.[selChi2], lower.tail = FALSE)

for(es_loop in 1:length(es)){
  
  # L will be the number of nonsignificant studies sampled from studies from M, 
  # which have true effect es > 0
  for (L in 0:dim(M)[1]){
    
    # randomly generate a p-value given the effect sizes
    # the null effects will yield uniform p-values      
    for (i in 1:n.iter){
      
      ran_L <- sample(dim(M)[1], size = L, replace = FALSE)
      
      # set population effect for the randomly sampled to es.
      pop_effect <- rep(0, dim(M)[1])
      if (!L == 0) pop_effect[ran_L] <- es[es_loop]
      
      # compute the noncentrality parameter
      # for the null effects this will be 0
      f2 <- pop_effect^2 / (1 - pop_effect^2)
      ncp <- f2 * ifelse(M$T_Test.Statistic..R. == "F",
                         M$T_df1..R. + M$T_df2..R. + 1,
                         ifelse(selr,
                                M$T_df1..R. + M$T_df2..R. + 1,
                                ifelse(selChi2,
                                       M$T_N..R.,
                                       M$T_df1..R. + M$T_df2..R.)))
      
      # determine area under H1 that corresponds to nonsignificant findings (= beta = Type II error)      
      beta <- pf(q = fCV, df1 = M$T_df1..R., df2 = M$T_df2..R., ncp = ncp, lower.tail = TRUE)
      beta[selz] <- pnorm(zCV, atanh(pop_effect[selz]), M$sigma[selz], lower.tail = TRUE) - pnorm(-zCV, atanh(pop_effect[selz]), M$sigma[selz], lower.tail = TRUE) 
      beta[selChi2] <- pchisq(Chi2CV, df = M$T_df1..R.[selChi2], ncp = ncp[selChi2], lower.tail = TRUE)
      
      # compute p-value alternative dist; sample from beta
      pA <- runif(length(beta), min = 0, max = beta)
      
      # test statistics simulated under alternative, based on beta
      fA <- qf(p = pA, df1 = M$T_df1..R., df2 = M$T_df2..R., ncp = ncp, lower.tail = TRUE)
      fA[selz] <- qnorm(pA[selz] + pnorm(-zCV, atanh(pop_effect[selz]), M$sigma[selz], lower.tail = TRUE),
                        atanh(pop_effect[selz]), M$sigma[selz], lower.tail = TRUE)
      fA[selChi2] <- qchisq(pA[selChi2], df = M$T_df1..R.[selChi2], ncp = ncp[selChi2], lower.tail = TRUE)  
      
      # pValues under null
      # convert test statistics under alternative to p-values under h0
      p0 <- pf(q = fA, df1 = M$T_df1..R., df2 = M$T_df2..R., lower.tail = FALSE)
      p0[selz] <- 1 - 2 * abs(pnorm(fA[selz], 0, M$sigma[selz], lower.tail = FALSE) - .5)
      p0[selChi2] <- pchisq(fA[selChi2], M$T_df1..R.[selChi2], lower.tail = FALSE)
      
      # Result of Fisher test for this run with L true effects > 0
      simulated <- FisherMethod(p0, 1, .05)
      
      res[[es_loop]][L + 1] <- res[[es_loop]][L + 1] + ifelse(simulated$Fish > Fish_M$Fish, 1, 0)
    }
  }
}

# pop_effect = .1
# sampled from L = 0:63
# 10000 iterations
# > res[[1]] / n.iter
# [1] 0.0386 0.0428 0.0519 0.0551 0.0576 0.0638 0.0752 0.0818 0.0902 0.0921 0.1018 0.1142
# [13] 0.1287 0.1312 0.1468 0.1595 0.1669 0.1812 0.1828 0.2010 0.2186 0.2352 0.2470 0.2629
# [25] 0.2733 0.2935 0.3037 0.3222 0.3439 0.3596 0.3723 0.3950 0.4047 0.4197 0.4432 0.4584
# [37] 0.4826 0.5048 0.5080 0.5365 0.5536 0.5576 0.5955 0.6062 0.6181 0.6363 0.6505 0.6713
# [49] 0.6738 0.6996 0.7166 0.7239 0.7394 0.7540 0.7760 0.7736 0.7920 0.8046 0.8165 0.8329
# [61] 0.8328 0.8523 0.8621 0.8719

# pop_effect = .3
# sampled from L = 0:63
# 10000 iterations
# > res[[2]] / n.iter
# [1] 0.0397 0.0590 0.0901 0.1294 0.1726 0.2231 0.2856 0.3542 0.4347 0.4947 0.5769 0.6511
# [13] 0.7153 0.7817 0.8280 0.8777 0.9038 0.9339 0.9532 0.9678 0.9790 0.9851 0.9903 0.9955
# [25] 0.9975 0.9982 0.9989 0.9995 0.9995 0.9996 0.9998 0.9999 1.0000 1.0000 1.0000 1.0000
# [37] 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000
# [49] 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000
# [61] 1.0000 1.0000 1.0000 1.0000

# pop_effect = .5
# sampled from L = 0:63
# 10000 iterations
# > res[[3]] / n.iter
# [1] 0.0386 0.0777 0.1265 0.2085 0.3026 0.4241 0.5481 0.6754 0.7932 0.8589 0.9169 0.9558
# [13] 0.9776 0.9888 0.9958 0.9985 0.9991 0.9995 0.9999 0.9999 1.0000 1.0000 1.0000 1.0000
# [25] 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000
# [37] 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000
# [49] 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000
# [61] 1.0000 1.0000 1.0000 1.0000


# discussion --------------------------------------------------------------

iccSS <- Anova(lm(dat$Computed[nsig] ~ dat$Source[nsig]), type="III")
# Computes the ICC
iccSS$Sum[2]/(iccSS$Sum[3]+iccSS$Sum[2])
