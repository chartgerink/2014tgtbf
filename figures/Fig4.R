# Set to main folder of project
# NOT figures folder
setwd(choose.dir())
load('figures/Fig4')

pdf(file = 'figures/Fig4.pdf', width = 7, height = 5)
par(mfrow = c(1, 1), mai = c(1, 1, .2, .2))
plot(x=1985:2013,
     nsigtemp,
     ylim = c(0, .4),
     type = 'o',
     xlab = "Year", ylab = 'Proportion nonsignificant',
     yaxs = 'i',
     cex.axis = 1, las = 1, bty = 'n')
abline(lm(nsigtemp ~ c(1985:2013)), lty = 2)
dev.off()