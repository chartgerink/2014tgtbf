# Contributor(s): Chris H.J. Hartgerink

# Change the object mypath to where you cloned the repository
# Home Computer
mypath <- "C:/Users/Chris/Dropbox/CJM/Masterproject/Analyzing/"
setwd(mypath)
# Load custom functions
source("a.Functions/FisherExTest.R")
# source("a.Functions/TerminalDigits.R")
source("a.Functions/esComp.R")
source("a.Functions/simNullDist.R")
source("a.Functions/simEffDist.R")
source("a.Functions/powerCalc.R")

esSize <- c(seq(.01,.05,.01)[1:4],seq(.05,.15,.02)[1:5],seq(.15,.95,.05))
source("a.Functions/powerCalc.R");x <- powerCalc(copilot,effectSize=esSize,n.iter=1000,testAlpha=.1)


###############
# Pilot Study #
###############
###############################################################################
# Preliminary stuff #
# Importing and preparing datafile
copilot <- read.table("1.Pilot study/copilot.txt",stringsAsFactors=F)
# Removing out of bounds p-values
selNA <- copilot$p_value_computed>=1
sum(selNA[!is.na(selNA)])
copilot$p_value_computed[selNA] <- NA
# Replace all comma's with decimal points and make the variable numeric.
copilot$test_statistic_value <- suppressWarnings(as.numeric(sub(",",".",copilot$test_statistic_value)))
copilot$df1 <- suppressWarnings(as.numeric(sub(",",".",copilot$df1)))
copilot$df2 <- suppressWarnings(as.numeric(sub(",",".",copilot$df2)))
# Computing unadjusted and adjusted effect sizes (OBSERVED)
copilot <- cbind(copilot, esComp.statcheck(copilot))
###############################################################################
# Step 1 - observed effect distribution versus nil effect distribution
# Computing unadjusted and adjusted effect size distributions under NO effect
simNullEs <- simNullDist(copilot, n.iter=1000, alpha=.05)
# Computing fisher test statistics (inexact)
resPilot <- FisherExTest(copilot$p_value_computed, copilot$pap_id)
plot(density(na.omit(simNullEs$esComp)),col="red", frame.plot=FALSE)
lines(density(na.omit(copilot$esComp)))
###############################################################################
# Step 2 - Power calculations fisher test for papers under different ES
powerRes <- powerCalc(copilot, effectSize=seq(.05,.95,by=.05),testAlpha=.1,n.iter=1000)














###############################################################################
######################### Random other things #################################
###############################################################################
# Fisher test #
###############
# Descriptive statistics
# Non-significant results per paper
# Percent
mean(resPilot$PercentNonSig)
hist(resPilot$PercentNonSig,breaks=20, main="Percentages of non-significant results", xlab="Percentage")
# Number of statistics reported
# Overall
mean(resPilot$CountNSig+resPilot$CountSig)
sd(resPilot$CountNSig+resPilot$CountSig)
# Significant
mean(resPilot$CountSig)
sd(resPilot$CountSig)
# Non-significant
mean(resPilot$CountNSig)
sd(resPilot$CountNSig)

# Results
signSel <- resPilot$PFish < .1
signNSel <- resPilot$PFishCompl < .1
copilot$N <- nCalc(copilot)
fishTestSign <- round(na.omit(cbind(resPilot$Fish[signSel],
                                    resPilot$PFish[signSel])),3)
fishTestComplSign <- round(na.omit(cbind(resPilot$FishCompl[signNSel],
                                         resPilot$PFishCompl[signNSel])),3)
# Number of signif results
dim(fishTestSign)[1]
# Number of signif results complement
dim(fishTestComplSign)[1]

# Computing mean of the mean sample size
meanNSig <- NULL
for(i in 1:length(signSel)){
  if(is.na(signSel[i])){
    meanNSig[i] <- NA
  } else if(signSel[i]==TRUE){
    meanNSig[i] <- mean(copilot$N[copilot$pap_id==i], na.rm=T)
  } else {
    meanNSig[i] <- NA  
  }
}
meanNComplSig <- NULL
for(i in 1:length(signNSel)){
  if(is.na(signNSel[i])){
    meanNComplSig[i] <- NA
  } else if(signNSel[i]==TRUE){
    meanNComplSig[i] <- mean(copilot$N[copilot$pap_id==i], na.rm=T)
  } else {
    meanNComplSig[i] <- NA  
  }
}
mean(meanNSig, na.rm=T)
mean(meanNComplSig, na.rm=T)



# Simulating effect sizes under the null
# Create vector of selected test stats
# Create vector of equal uniform p val under alpha
# Compute test statistic under p
# Compute effect size under test statistic

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
