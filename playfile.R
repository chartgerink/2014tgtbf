# Playfile masterproject
mypath <- "C:/Users/Chris/Dropbox/CJM/Masterproject/"
setwd(mypath)
# copilot <- read.csv("Co-piloot//error_data_clean.csv", sep=',', header=T)
copilotID <- read.table("copilotID.txt")

# Densities of reported versus computed
pSelRep <- copilotID$p_value_reported>=.05
pSelCom <- copilotID$p_value_computed>=.05
plot(density(na.omit(copilotID$p_value_reported[pSelRep])), xlim=c(1,.05), main="")
lines(density(na.omit(copilotID$p_value_computed[pSelCom])), col=2)
abline(h=1, lty=2, lwd=.5)

# Removing out of bounds p-values
selNA <- copilotID$p_value_reported>1
sum(selNA[!is.na(selNA)])
copilotID$p_value_reported[selNA] <- NA


# Create placeholder objects
selP <- NULL
ratioSig <- NULL
ratioNSig <- NULL
fishTest <- NULL
fishTestCompl <- NULL
pfishTest <- NULL
pfishTestCompl <- NULL
# Calculate all Fisher values and p-values
# Both normal and complement
for(i in 1:length(uniqID)){
  uniqID <- unique(copilotID$ideed)
  selP[[i]] <- na.omit(copilotID$p_value_reported[uniqID[i]==copilotID$ideed])
  ratioSig[i] <- sum(selP[[i]]<=.05)/length(selP[[i]])
  ratioNSig <- 1-ratioSig
  selPStar <- (selP[[i]][selP[[i]]>.05]-.05)/.95
  selPStarCompl <- 1-selPStar
  fishTest[i] <- -sum(log(selPStar))
  fishTestCompl[i] <- -sum(log(selPStarCompl))
  pfishTest[i] <- pgamma(fishTest[i], shape=length(selP))
  pfishTestCompl[i] <- pgamma(fishTestCompl[i], shape=length(selP[[i]]))
}

# Meeting Marcel 29 jan
copilotID$ideed
x <- na.omit(uniqID[pTestCompl>.9])

pValSusp <- NULL
for(i in 1:length(x)){
  pValSusp[[i]] <- copilotID$p_value_computed[copilotID$ideed==x[i]] 
}
pValSusp[[1]][pValSusp[[1]]>.05]
pValSusp[[2]][pValSusp[[2]]>.05]
pValSusp[[3]][pValSusp[[3]]>.05]
pValSusp[[4]][pValSusp[[4]]>.05]
pValSusp[[5]][pValSusp[[5]]>.05]
pValSusp[[6]][pValSusp[[6]]>.05]
pValSusp[[7]][pValSusp[[7]]>.05]
pValSusp[[8]][pValSusp[[8]]>.05]
pValSusp[[9]][pValSusp[[9]]>.05]
pValSusp[[10]][pValSusp[[10]]>.05]

# Diminishing of returns sample in power calculations
n <- c(4,9,15,22,
       29,38,49,64,68,72,77,82,88,96,105,118,140)
pow <- c(seq(0.1,.7,by=.1),seq(.8,.99,.02))
plot(n,pow)


###############################################################################
# Random stuff
###############################################################################
# # Compute a numeric ID variable for the papers
# identify <- as.character(unique(copilot$source_variable))
# for(i in 1:length(identify)){
#   for(j in 1:length(copilot$source_variable)){
#     ifelse(copilot$source_variable[j]==identify[i], copilot$ideed[j] <- i, NA) 
#   }
# }
write.table(x,"copilotID3.txt")

head(copilotID[,c(57,2,4,5,6,8,10,11,12)])
# # P-uniform under H0
# pval <- NULL
# for(i in 1:10000){
#   temp <- rnorm(100,0,1)
#   x <- t.test(temp, y=NULL, mu=0,paired=F,var.equal=T)
#   pval[i] <- x$p.value}
# hist(pval,breaks=10)