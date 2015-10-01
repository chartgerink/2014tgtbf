# Set to main folder of project
# NOT figures folder
setwd(choose.dir())

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

x <- NULL
i <- 1

for(y in 1985:2013){
  x[i] <- sum(dat$esComp[!is.na(dat$esComp) & dat$years.y. == y] < .1) / length(dat$esComp[!is.na(dat$esComp) & dat$years.y. == y])
  i <- i + 1
}

pdf('figures/Fig8.pdf', width = 7, height = 5)
par(mfrow = c(1, 1), mai = c(1, 1, .2, .2))
plot(y = x, x = 1985:2013, type = 'o', ylim = c(0, .6), xlab = "Year", ylab = "Proportion small effects (r<.1)", bty = 'n')
dev.off()
