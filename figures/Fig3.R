# Choose the directory all files are cloned in
# NOT the figures folder.
setwd(choose.dir())

load('figures/Fig3')

# Effect PDF
pdf('figures/Fig3.pdf',width=7, height=8)
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