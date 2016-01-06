# Set to main folder of project
# NOT figures folder
setwd(choose.dir())

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
     xlab = latex2exp("Correlation ($|\\eta|$)"),
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
     xlab = latex2exp("Correlation ($|\\eta|$)"),
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
