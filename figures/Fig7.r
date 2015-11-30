# Set to project folder
# NOT figures folder
setwd(choose.dir())

# Load custom functions all at once
customFunct <- list.files('functions/')
for(i in 1:length(customFunct)){
  source(
    paste0('functions/', customFunct[i])
  )
}

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

pdf('figures/Fig7.pdf', width=7, height=6)
par(mar = c(4, 4, .2, 2))
plot(x=1985:2013, y=medianN, type='o', col="black",
     ylab="N", xlab="Year", ylim=c(0,150), cex.lab=1.2, las=1, lwd=1, cex.axis=1.2, xaxs='i', bty = 'n')
lines(x=1985:2013, y=p25, type='o', col='grey')
lines(x=1985:2013, y=p75, type='o', col='grey')
text(y = 100, x = 2010, "P75", col = "grey")
text(y = 50, x = 2010, "P50", col = "black")
text(y = 22, x = 2010, "P25", col = "grey")
dev.off()