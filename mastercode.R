# Contributor(s): Chris H.J. Hartgerink

# Change the object mypath to where you cloned the repository
# Home Computer
mypath <- "C:/Users/Chris/Dropbox/CJM/Masterproject/Analyzing/"
setwd(mypath)
# Load custom functions all at once
customFunct <- list.files('a.Functions/')
for(i in 1:length(customFunct)){
  source(paste0('a.Functions/',customFunct[i]))
}


# Introduction # 
# Rejection rate comparison sample
# Commented out due to runtime
# x1 <- proc.time()[1]
# y <- 0
# for(i in 1:1000){
#   set.seed(i)
#   x1 <- rnorm(5000000,1,.1) + runif(5000000,0,.002)
#   x2 <- rnorm(5000000,1,.1)+ runif(5000000,0,.001)
#   x <- t.test(x1,x2,paired=F,var.equal=F,)
#   if(x$p.value<.05){y <- y+1}}
# x2 <- proc.time()[2]
# y


###############
# Pilot Study #
###############
###############################################################################
# Preliminary stuff #
# Importing and preparing datafile
copilot <- read.table("1.Pilot study/copilot.txt",stringsAsFactors=F)
# The following reordering of df1 to df2 is due to t and r values essentially
# being F distributions, making for more parsimonious code. Copilot data was
# collected by running old version of Statcheck, hence manually done here.
copilot$df2[copilot$Statistic == "t" | copilot$Statistic == "r"] = copilot$df1[copilot$Statistic == "t" | copilot$Statistic == "r"]
copilot$df1[copilot$Statistic == "t" | copilot$Statistic == "r"] = NA
copilot <- prep.statcheck(copilot)
###############################################################################
# Step 1 - observed effect distribution versus nil effect distribution
# Computing unadjusted and adjusted effect size distributions under NO effect
set.seed(1234)
simNullEs <- simNullDist(copilot, n.iter=length(copilot$esComp), alpha=.05)
# Figure showing the difference between the nil effects and the observed effect distributions
# Under non-significant test
plot(ecdf(na.omit(simNullEs$esComp)),
     lty=1,
     frame.plot=F, 
     main="Effect size distributions",
     xlim=c(0,1),
     xaxs="i",
     yaxs="i",
     xlab="Partial eta-squared",
     ylab = "Cumulative density",
     cex.axis=.6,
     cex.lab=.7,
     col = "grey")
lines(ecdf(na.omit(copilot$esComp)))
legend(x=.7,y=.8,legend=c(expression('H'[0]), 'Observed effects'),cex=.7,lty=c(1,1), col = c("grey","black"),box.lwd=0)

# Percentages of effects
categories <- c(.01, .06, .14)
effPerc <- effectPercent(copilot$esComp, categories = categories)
clip(x1=0, x2=categories[1], y1 = 0, y2 = effPerc[1])
abline(v = categories[1], lty = 2, col = "grey")
clip(x1=0, x2=categories[2], y1 = 0, y2 = effPerc[2])
abline(v = categories[2], lty = 3, col = "grey")
clip(x1=0, x2=categories[3], y1 = 0, y2 = effPerc[3])
abline(v = categories[3], lty = 4, col = "grey")

# Kolmogorov-Smirnov test
# Moet die over de p of over de effecten?
ks.test(copilot$esComp, simNullEs$esComp)


###############################################################################
# Step 2 - Computing (inexact) fisher test statistics
resPilot <- FisherExTest(copilot$Computed, copilot$Source)


###############################################################################
# Step 3 - Power calculations fisher test for papers under different ES
esSize <- c(seq(.01,.15,.02),seq(.20,.5,.05),seq(.6,.9,.1))
set.seed(94438)
powerRes <- powerCalc(copilot,effectSize=esSize,n.iter=1,testAlpha=.1)
write.csv2(powerRes[[1]],'powerCalcFish.csv')
write.csv2(powerRes[[2]],'powerCalcFishCompl.csv')


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
    meanNSig[i] <- mean(copilot$N[copilot$Source==i], na.rm=T)
  } else {
    meanNSig[i] <- NA  
  }
}
meanNComplSig <- NULL
for(i in 1:length(signNSel)){
  if(is.na(signNSel[i])){
    meanNComplSig[i] <- NA
  } else if(signNSel[i]==TRUE){
    meanNComplSig[i] <- mean(copilot$N[copilot$Source==i], na.rm=T)
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
