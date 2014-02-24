# Masterfile source code
# Last updated: january 2014
# Contributor(s): Chris H.J. Hartgerink

# Change the object mypath to where you cloned the repository
# Home Computer
mypath <- "C:/Users/Chris/Dropbox/CJM/Masterproject/"
setwd(mypath)
# Work computer
mypath <- "D:/Chris/Dropbox/CJM/Masterproject/"
setwd(mypath)

# Load custom functions
source("Analyzing//0. Functions/FisherExTest.R")
source("Analyzing//0. Functions/TerminalDigits.R")


###############
# Pilot Study #
###############
## Importing and preparing datafile
copilot <- read.table("Analyzing//1. Pilot study/copilot.txt",stringsAsFactors=F)
# Removing out of bounds p-values
selNA <- copilot$p_value_computed>=1
sum(selNA[!is.na(selNA)])
copilot$p_value_computed[selNA] <- NA
# Replace all comma's with decimal points and make the variable numeric.
copilot$test_statistic_value <- suppressWarnings(as.numeric(sub(",",".",copilot$test_statistic_value)))
copilot$df1 <- suppressWarnings(as.numeric(sub(",",".",copilot$df1)))
copilot$df2 <- suppressWarnings(as.numeric(sub(",",".",copilot$df2)))

## Computing fisher test statistics (inexact)
ResPilot <- FisherExTest(copilot$p_value_computed, copilot$pap_id)

## Computing effect sizes for all 
# M,SD of number of non-significant p-values per study
round(mean(unlist(lapply(selPtrunc,length))),2)
round(sd(unlist(lapply(selPtrunc,length))),2)
# Fisher test *
sum(pfishTest<=.05); length(pfishTest)
# Fisher test -*
sum(pfishTestCompl<=.05); length(pfishTestCompl)
# Mean proportion of significant values
mean(na.omit(ratioSig))
# Mean proportion of nonsignificant values
mean(na.omit(ratioNSig))

###########
# Figures #
###########
# Figure 1
# Conceptual plot of p-distributions
png(file = "Writing/fig1.png", width=821, height=501)
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
