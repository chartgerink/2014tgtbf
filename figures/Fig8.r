# Choose the directory all files are cloned in
# NOT the figures folder.
setwd(choose.dir())

load('figures/Fig3')

# Effect PDF
pdf('figures/Fig8.pdf',width=7, height=8)
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