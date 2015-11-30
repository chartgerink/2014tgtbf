# Set to project folder
# NOT figures folder
setwd(choose.dir())
if(!require(ggplot2)){install.packages('ggplot2')}
if(!require(plyr)){install.packages('plyr')}
library(ggplot2)
library(plyr)

load('figures/Fig6')

# Computing the number of significant Fisher results per year
# As proportion of all papers reporting nonsignificant results
fishDF$logicalP <- ifelse(fishDF$FisherP<.1, 1, 0)
fisherYear <- ddply(fishDF, .(year), summarise, propYear=mean(logicalP, na.rm=TRUE))

knsYear <- ddply(fishDF, .(year), summarise, kYear=mean(kRes, na.rm=TRUE))

mydf <- data.frame(x = fisherYear$year,
                   y = fisherYear$propYear,
                   count = knsYear$kYear)

ggplot(mydf, aes(x = x, y = y)) + geom_point(aes(size = count)) + ylim(0, 1) + geom_smooth(method="lm") +
  xlab("Year") + ylab("Proportion significant Fisher results")

ggsave(filename = 'figures/Fig6.png', plot = last_plot(), width = 7, height = 5)