# Set to main folder of project
# NOT figures folder
setwd(choose.dir())

# Load data
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

pdf('figures/S1Fig.pdf', onefile = TRUE, width = 11, height = 11)
par(mfrow = c(2, 2), mai = c(1.2, 1.2, .8, .5))

for(i in 1:4){
  sel <- dat$journals.jour. == sort(unique(dat$journals.jour.))[i] & nsig
  set.seed(i)
  simNullEs <- simNullDist(dat, n.iter=length(dat$esComp[sel])*3, alpha = .05)
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