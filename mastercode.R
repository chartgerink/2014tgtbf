# Contributor(s): Chris H.J. Hartgerink

# Change the object mypath to where you cloned the repository
# Home Computer
mypath <- "C:/Users/Chris/Dropbox/CJM/Masterproject/Analyzing/"
setwd(mypath)
# Load custom functions
source("a.Functions/FisherExTest.R")
# source("a.Functions/TerminalDigits.R")
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

# Computing fisher test statistics (inexact)
resPilot <- FisherExTest(copilot$p_value_computed, copilot$pap_id)
# Descriptive statistics
# Non-significant results
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


# Running meta-analysis on unadjusted effect sizes
if(!require(metafor)){install.packages('metafor')}
# Computing unadjusted and adjusted effect sizes
copilot <- cbind(copilot, ESComp(copilot))
# Converting effect sizes into r
rSqES <- sqrt(copilot$esComp)
# Fisher's r to Z transformation for M-A
rES <- r2z(rSqES)
# Calculating sampling variance
sVar <- sVarComp(copilot)
# Selecting out the negative sampling variances
# Caused by either very small samples or by df1 > df2 in F tests
# See next line of commented out code for the cases that have this
# (na.omit(copilot[((sVar)<0),3:5]))
rES <- rES[!(sVar<0)]
sVar <- sVar[!sVar<0]
# Random effects MA
# Note: this takes a long time to run (under DerSimonian-Laird approx. 5-6minutes, REML takes >1h).
# That is why I comment this out.
modelPilot <- rma(yi=rES, vi=sVar, method="DL", digits=10)
# Transforming back MA results
estPilot <- z2r(modelPilot$b[1])
sePilot <- z2r(modelPilot$se)
ciLowPilot <- z2r(modelPilot$ci.lb)
ciHighPilot <- z2r(modelPilot$ci.ub)
tau2Pilot <- z2r(modelPilot$tau2)
resModelPilot <- round(cbind(estPilot, sePilot, ciLowPilot, ciHighPilot, tau2Pilot),3)



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
