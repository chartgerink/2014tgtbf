setwd("D:/Dropbox/projects/2014tgtbf/figures")

pdf('Fig2.pdf', width = 7, height = 5)
par(mar = c(1, 2, 3, 4))
plot(dnorm(seq(-3, 3, .001), 0, 1), type = 'l',
     bty = 'n', xaxt = 'n', yaxt = 'n')
# lines(dnorm(seq(-2, 4, .001), 1, .5))
dev.off()

