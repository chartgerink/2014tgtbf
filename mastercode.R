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

###############
# Pilot Study #
###############
## Importing and preparing datafile
copilot <- read.table("Analyzing//1. Pilot study/copilot.txt")
# Removing out of bounds p-values
selNA <- copilot$p_value_computed>=1
sum(selNA[!is.na(selNA)])
copilot$p_value_computed[selNA] <- NA
# Replace all comma's with decimal points and make the variable numeric.
copilot$test_statistic_value <- suppressWarnings(as.numeric(sub(",",".",copilot$test_statistic_value)))

## Computing fisher test statistics (inexact)
# Create placeholder objects
uniqID <- unique(copilot$pap_id)
selP <- NULL
ratioSig <- NULL
ratioNSig <- NULL
fishTest <- NULL
fishTestCompl <- NULL
pfishTest <- NULL
pfishTestCompl <- NULL
meanP <- NULL
selPtrunc <- NULL
meanPtrunc <- NULL
# Calculate all Fisher values and p-values
# Both normal and complement
for(i in 1:length(uniqID)){
  selP[[i]] <- na.omit(copilot$p_value_reported[uniqID[i]==copilot$ideed])
  selPtrunc[[i]] <- selP[[i]][selP[[i]]>.05] 
  ratioSig[i] <- sum(selP[[i]]<=.05)/length(selP[[i]])
  ratioNSig[i] <- 1-ratioSig[i]
  selPStar <- (selP[[i]][selP[[i]]>.05]-.05)/.95
  selPStarCompl <- 1-selPStar
  fishTest[i] <- -sum(log(selPStar))
  fishTestCompl[i] <- -sum(log(selPStarCompl))
  pfishTest[i] <- pgamma(fishTest[i], shape=length(selP))
  pfishTestCompl[i] <- pgamma(fishTestCompl[i], shape=length(selP[[i]]))
  meanP[i] <- mean(selP[[i]])
  meanPtrunc[i] <- mean(selPtrunc[[i]])
}


Test <- NULL
TestCompl <- NULL
kPaper <- NULL
pTest <- NULL
pTestCompl <- NULL
indTest <- NULL
indTestCompl <- NULL
for(i in 1:length(uniqID)){
  if(length(selPtrunc[[i]])==0){Test[i] <- NA
                                TestCompl[i] <- NA
                                kPaper[i] <- NA} else{
  Test[i] <- sum(-log((selPtrunc[[i]]-.05)/.95))
  indTest[[i]] <- -log((selPtrunc[[i]]-.05)/.95)
  TestCompl[i] <- sum(-log(1-((selPtrunc[[i]]-.05)/.95)))
  indTestCompl[[i]] <- -log(1-((selPtrunc[[i]]-.05)/.95))
  kPaper[i] <- length(selPtrunc[[i]])}
  pTest[i] <- pgamma(Test[i], shape=kPaper[i])
  pTestCompl[i] <- pgamma(TestCompl[i], shape=kPaper[i])
  # Kan allebei hoog zijn
}

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

# EPV
# Number with mean P-value >.5 
sum(na.omit(meanP>.5))
# Range of MPV > .5
summary(meanP[meanP>.5])

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
