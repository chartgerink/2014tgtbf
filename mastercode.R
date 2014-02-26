# Contributor(s): Chris H.J. Hartgerink

# Change the object mypath to where you cloned the repository
# Home Computer
mypath <- "C:/Users/Chris/Dropbox/CJM/Masterproject/Analyzing/"
setwd(mypath)
# Load custom functions
source("a.Functions/FisherExTest.R")
source("a.Functions/TerminalDigits.R")
source("a.Functions/ESComp.R")


###############
# Pilot Study #
###############
## Importing and preparing datafile
copilot <- read.table("1.Pilot study/copilot.txt",stringsAsFactors=F)
# Removing out of bounds p-values
selNA <- copilot$p_value_computed>=1
sum(selNA[!is.na(selNA)])
copilot$p_value_computed[selNA] <- NA
# Replace all comma's with decimal points and make the variable numeric.
copilot$test_statistic_value <- suppressWarnings(as.numeric(sub(",",".",copilot$test_statistic_value)))
copilot$df1 <- suppressWarnings(as.numeric(sub(",",".",copilot$df1)))
copilot$df2 <- suppressWarnings(as.numeric(sub(",",".",copilot$df2)))

## Computing fisher test statistics (inexact)
resPilot <- FisherExTest(copilot$p_value_computed, copilot$pap_id)

## Computing unadjusted and adjusted effect sizes
copilot <- cbind(copilot, ESComp(copilot))

# Descriptive statistics
# Non-significant results
# Percent
mean(resPilot$PercentNonSig)
hist(resPilot$PercentNonSig,breaks=20)


###########
# Figures #
###########
# Figure 1
# Conceptual plot of p-distributions
png(file = "b.Figures/fig1.png", width=821, height=501)
curve(exp(-x+.5), from=1, to=.05, add=F, xlim=c(1,.05), ylim=c(0,2.5)
      , xlab="P-value", ylab="Density")
curve(exp(-1.5*x+.75), from=1, to=.05, add=T)
abline(h=1, lty=2)
# text(.15,.9, expression(Under~H[0]), cex=.7)
# text(.15,1.53, expression(Under~H[A]), cex=.7)
text(.93,.85, "Upperbound", cex=.8)
text(.12,1.21, "Lowerbound", cex=.8)
arrows(1,.6
       ,1,1
       ,code=2)
arrows(.05,1.568312
       ,.05,1
       ,code=2)
legend(x = .3, y = .5
       , legend = c(expression(Under~H[0]), expression(Under~H[A]))
       , lty = c(2,1)
       ,box.lwd = 0,box.col = "white",bg = "white")
dev.off()
