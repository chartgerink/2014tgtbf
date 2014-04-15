powerCalc2 <- function(
	x,
	### A \code{statcheck} object
	effectSize = 0,
	### Vector of effect sizes to estimate power under
	testAlpha = .1,
	### Numeric indicating the alpha used for the Fisher test
	alpha = .05,
	### The alpha used to determine rejection areas for ORIGINAL tests
	# twoTailed = FALSE,
	### Logical indicating original hypothesis
	n.iter = 1000
	### Number of iterations to simulate power under
	){
	##details<<
	##Returns dataframe with the amount of rows equalling the amount of unique papers in the original statcheck object.
	##Indicates the power per paper, per effect size (each column).

	require(gtools)

	# Create object to be filled with all data
	finalDF <- NULL
	temp <- matrix(nrow=length(unique(x$Source)),ncol=length(effectSize))
	finalDF[[1]] <- as.data.frame(temp)
	finalDF[[2]] <- as.data.frame(temp)
	names(finalDF[[1]]) <- effectSize
	names(finalDF[[2]]) <- effectSize
	names(finalDF) <- c("fishTest","fishTestCompl")
	rm(temp)

	# Create empty objects
	tempMat <- NULL
	pFishTest <- NULL
	pFishTestCompl <- NULL
	propSigFishTest <- NULL
	propSigFishTestCompl <- NULL
	sel <- NULL
	fishTest <- NULL
	fishTestCompl <- NULL
	lengthRes <- NULL
	
	# Simulate effects under the effect sizes for each test statist ic to determine power of the test
	# For each paper
	for(p in 1:length(unique(x$Source))){
		# Select the stats for that paper
		selectStats <- x[x$Source==p,]
		# Create a matrix to put results of the pVal simulations in
		tempMat[[p]] <- matrix(nrow=n.iter,ncol=dim(selectStats)[1])
		lengthRes[p] <- dim(selectStats)[1]

		# Logical for determining if power is so high further running is redundant
		cont1 <- TRUE
		cont2 <- TRUE

		for(es in 1:length(effectSize)){
###############################################################################
			if(cont1 == TRUE & cont2 == TRUE){
				for(s in 1:dim(selectStats)[1]){
					res <- NULL
					resP <- NULL
					if(selectStats$Statistic[s]=="t"){
						# Determine critical value for (1tailed) alpha
						cv <- qf(alpha,1,selectStats$df1[s],lower.tail=F)
						# Calculate NCP for effectsize
						ncp <- (effectSize[es]/(1-effectSize[es]))*(1+selectStats$df1[s]+1)
						# Calculate area under curve for alternate distribution
						prob <- pf(cv,1,selectStats$df1[s],ncp=ncp,lower.tail=T)
						# Sample between 0 and prob from previous
						sampleProb <- runif(n.iter,0,prob)
						# Calculate accompanying test statistic with probability
						sampleTest <- qf(sampleProb,1,selectStats$df1[s],ncp=ncp)
						# Calculate p-value under the null given the non-significant value drawn from alternate
						resP <- pf(sampleTest,1,selectStats$df1[s],lower.tail=F)
					} 
					else if(selectStats$Statistic[s]=="F"){
						cv <- qf(alpha,selectStats$df1[s],selectStats$df2[s],lower.tail=F)
						ncp <- (effectSize[es]/(1-effectSize[es]))*(selectStats$df1[s]+selectStats$df2[s]+1)
						prob <- pf(cv,selectStats$df1[s],selectStats$df2[s],ncp=ncp,lower.tail=T)
						sampleProb <- runif(n.iter,0,prob)	
						sampleTest <- qf(sampleProb,selectStats$df1[s],selectStats$df2[s],ncp=ncp)
						resP <- pf(sampleTest,selectStats$df1[s],selectStats$df2[s],lower.tail=F)
					}
					else if(selectStats$Statistic[s]=="r"){
						cv <- qf(alpha,1,selectStats$df1[s],lower.tail=F)
						ncp <- (effectSize[es]/(1-effectSize[es]))*(1+selectStats$df1[s]+1)
						prob <- pf(cv,1,selectStats$df1[s],ncp=ncp,lower.tail=T)
						sampleProb <- runif(n.iter,0,prob)
						sampleTest <- qf(sampleProb,1,selectStats$df1[s],ncp=ncp)
						resP <- pf(sampleTest,1,selectStats$df1[s],lower.tail=F)
					}
					else {
						resP <- NA
					}
					# Here we return the p.val to state space of [0;1]
					tempMat[[p]][,s] <- (resP-alpha)/(1-alpha)
					print(paste("Still working, no worries, paper", p, "statistic", s, sep=" "))
				}
				
				if(sum(is.na(tempMat[[p]]))==(dim(tempMat[[p]])[1]*dim(tempMat[[p]])[2])){
					finalDF[[1]][p,es] <- NA
					finalDF[[2]][p,es] <- NA
				} 
				else{
					sel[[p]] <- as.matrix(as.matrix(tempMat[[p]])[,colSums(is.na(tempMat[[p]])) != nrow(tempMat[[p]])])
					# Apply fisher test
					fishTest[[p]] <- apply(sel[[p]],1,function(x) -sum(log(x)))
					fishTestCompl[[p]] <- apply(sel[[p]],1,function(x) -sum(log(1-(x))))
					# Compute fisher p value
					pFishTest[[p]] <- pgamma(fishTest[[p]],dim(sel[[p]])[2],lower.tail=F)
					pFishTestCompl[[p]] <- pgamma(fishTestCompl[[p]],dim(sel[[p]])[2],lower.tail=T)
					# Compute the number of significant fisher tests
					propSigFishTest[p] <- sum(na.omit(pFishTest[[p]]) < testAlpha)/length(na.omit(pFishTest[[p]]))
					propSigFishTestCompl[p] <- sum(na.omit(pFishTestCompl[[p]]) < testAlpha)/length(na.omit(pFishTestCompl[[p]]))
					# Save the power
					finalDF[[1]][p,es] <- round(propSigFishTest[p], 3)
					finalDF[[2]][p,es] <- round(propSigFishTestCompl[p], 3)
					
					if(finalDF[[1]][p,es] >= .995 & finalDF[[2]][p,es] >= .995){
						cont1 <- FALSE
						cont2 <- FALSE
					}
					else if(finalDF[[1]][p,es] >= .995){
						cont1 <- FALSE
						cont2 <- TRUE
					}
					else if(finalDF[[2]][p,es] >=.995){
						cont1 <- TRUE
						cont2 <- FALSE
					}
					else{
						cont1 <- TRUE
						cont2 <- TRUE
					}
				}
			}
###############################################################################
			else if(cont1 == TRUE & cont2 == FALSE){
				for(s in 1:dim(selectStats)[1]){
					res <- NULL
					resP <- NULL
					if(selectStats$Statistic[s]=="t"){
						# Determine critical value for (1tailed) alpha
						cv <- qf(alpha,1,selectStats$df1[s],lower.tail=F)
						# Calculate NCP for effectsize
						ncp <- (effectSize[es]/(1-effectSize[es]))*(1+selectStats$df1[s]+1)
						# Calculate area under curve for alternate distribution
						prob <- pf(cv,1,selectStats$df1[s],ncp=ncp,lower.tail=T)
						# Sample between 0 and prob from previous
						sampleProb <- runif(n.iter,0,prob)
						# Calculate accompanying test statistic with probability
						sampleTest <- qf(sampleProb,1,selectStats$df1[s],ncp=ncp)
						# Calculate p-value under the null given the non-significant value drawn from alternate
						resP <- pf(sampleTest,1,selectStats$df1[s],lower.tail=F)
					} 
					else if(selectStats$Statistic[s]=="F"){
						cv <- qf(alpha,selectStats$df1[s],selectStats$df2[s],lower.tail=F)
						ncp <- (effectSize[es]/(1-effectSize[es]))*(selectStats$df1[s]+selectStats$df2[s]+1)
						prob <- pf(cv,selectStats$df1[s],selectStats$df2[s],ncp=ncp,lower.tail=T)
						sampleProb <- runif(n.iter,0,prob)	
						sampleTest <- qf(sampleProb,selectStats$df1[s],selectStats$df2[s],ncp=ncp)
						resP <- pf(sampleTest,selectStats$df1[s],selectStats$df2[s],lower.tail=F)
					}
					else if(selectStats$Statistic[s]=="r"){
						cv <- qf(alpha,1,selectStats$df1[s],lower.tail=F)
						ncp <- (effectSize[es]/(1-effectSize[es]))*(1+selectStats$df1[s]+1)
						prob <- pf(cv,1,selectStats$df1[s],ncp=ncp,lower.tail=T)
						sampleProb <- runif(n.iter,0,prob)
						sampleTest <- qf(sampleProb,1,selectStats$df1[s],ncp=ncp)
						resP <- pf(sampleTest,1,selectStats$df1[s],lower.tail=F)
					}
					else {
						resP <- NA
					}
					# Here we return the p.val to state space of [0;1]
					tempMat[[p]][,s] <- (resP-alpha)/(1-alpha)
					print(paste("Still working, no worries, paper", p, "statistic", s, sep=" "))
				}
				
				if(sum(is.na(tempMat[[p]]))==(dim(tempMat[[p]])[1]*dim(tempMat[[p]])[2])){
					finalDF[[1]][p,es] <- NA
					finalDF[[2]][p,es] <- NA
				} 
				else{
					sel[[p]] <- as.matrix(as.matrix(tempMat[[p]])[,colSums(is.na(tempMat[[p]])) != nrow(tempMat[[p]])])
					# Apply fisher test
					fishTest[[p]] <- apply(sel[[p]],1,function(x) -sum(log(x)))
					fishTestCompl[[p]] <- apply(sel[[p]],1,function(x) -sum(log(1-(x))))
					# Compute fisher p value
					pFishTest[[p]] <- pgamma(fishTest[[p]],dim(sel[[p]])[2],lower.tail=F)
					pFishTestCompl[[p]] <- pgamma(fishTestCompl[[p]],dim(sel[[p]])[2],lower.tail=T)
					# Compute the number of significant fisher tests
					propSigFishTest[p] <- sum(na.omit(pFishTest[[p]]) < testAlpha)/length(na.omit(pFishTest[[p]]))
					propSigFishTestCompl[p] <- sum(na.omit(pFishTestCompl[[p]]) < testAlpha)/length(na.omit(pFishTestCompl[[p]]))
					# Save the power
					finalDF[[1]][p,es] <- round(propSigFishTest[p], 3)
					finalDF[[2]][p,es] <- 1
					
					if(finalDF[[1]][p,es] >= .995){
						cont1 <- FALSE
					}
					else{
						cont1 <- TRUE
						cont2 <- FALSE
					}
				}
			}
###############################################################################
			else if(cont1 == FALSE & cont2 == TRUE){
				for(s in 1:dim(selectStats)[1]){
					res <- NULL
					resP <- NULL
					if(selectStats$Statistic[s]=="t"){
						# Determine critical value for (1tailed) alpha
						cv <- qf(alpha,1,selectStats$df1[s],lower.tail=F)
						# Calculate NCP for effectsize
						ncp <- (effectSize[es]/(1-effectSize[es]))*(1+selectStats$df1[s]+1)
						# Calculate area under curve for alternate distribution
						prob <- pf(cv,1,selectStats$df1[s],ncp=ncp,lower.tail=T)
						# Sample between 0 and prob from previous
						sampleProb <- runif(n.iter,0,prob)
						# Calculate accompanying test statistic with probability
						sampleTest <- qf(sampleProb,1,selectStats$df1[s],ncp=ncp)
						# Calculate p-value under the null given the non-significant value drawn from alternate
						resP <- pf(sampleTest,1,selectStats$df1[s],lower.tail=F)
					} 
					else if(selectStats$Statistic[s]=="F"){
						cv <- qf(alpha,selectStats$df1[s],selectStats$df2[s],lower.tail=F)
						ncp <- (effectSize[es]/(1-effectSize[es]))*(selectStats$df1[s]+selectStats$df2[s]+1)
						prob <- pf(cv,selectStats$df1[s],selectStats$df2[s],ncp=ncp,lower.tail=T)
						sampleProb <- runif(n.iter,0,prob)	
						sampleTest <- qf(sampleProb,selectStats$df1[s],selectStats$df2[s],ncp=ncp)
						resP <- pf(sampleTest,selectStats$df1[s],selectStats$df2[s],lower.tail=F)
					}
					else if(selectStats$Statistic[s]=="r"){
						cv <- qf(alpha,1,selectStats$df1[s],lower.tail=F)
						ncp <- (effectSize[es]/(1-effectSize[es]))*(1+selectStats$df1[s]+1)
						prob <- pf(cv,1,selectStats$df1[s],ncp=ncp,lower.tail=T)
						sampleProb <- runif(n.iter,0,prob)
						sampleTest <- qf(sampleProb,1,selectStats$df1[s],ncp=ncp)
						resP <- pf(sampleTest,1,selectStats$df1[s],lower.tail=F)
					}
					else {
						resP <- NA
					}
					# Here we return the p.val to state space of [0;1]
					tempMat[[p]][,s] <- (resP-alpha)/(1-alpha)
					print(paste("Still working, no worries, paper", p, "statistic", s, sep=" "))
				}

				if(sum(is.na(tempMat[[p]]))==(dim(tempMat[[p]])[1]*dim(tempMat[[p]])[2])){
					finalDF[[1]][p,es] <- NA
					finalDF[[2]][p,es] <- NA
				}
				else{
					sel[[p]] <- as.matrix(as.matrix(tempMat[[p]])[,colSums(is.na(tempMat[[p]])) != nrow(tempMat[[p]])])
					# Apply fisher test
					fishTest[[p]] <- apply(sel[[p]],1,function(x) -sum(log(x)))
					fishTestCompl[[p]] <- apply(sel[[p]],1,function(x) -sum(log(1-(x))))
					# Compute fisher p value
					pFishTest[[p]] <- pgamma(fishTest[[p]],dim(sel[[p]])[2],lower.tail=F)
					pFishTestCompl[[p]] <- pgamma(fishTestCompl[[p]],dim(sel[[p]])[2],lower.tail=T)
					# Compute the number of significant fisher tests
					propSigFishTest[p] <- sum(na.omit(pFishTest[[p]]) < testAlpha)/length(na.omit(pFishTest[[p]]))
					propSigFishTestCompl[p] <- sum(na.omit(pFishTestCompl[[p]]) < testAlpha)/length(na.omit(pFishTestCompl[[p]]))
					# Save the power
					finalDF[[1]][p,es] <- 1
					finalDF[[2]][p,es] <- round(propSigFishTestCompl[p], 3)
					
					if(finalDF[[2]][p,es] >= .995){
						cont2 <- FALSE
					}
					else{
						cont1 <- FALSE
						cont2 <- TRUE
					}
				}
			}
###############################################################################
			else{
				finalDF[[1]][p, es] <- 1
				finalDF[[2]][p, es] <- 1
			}
		}

	}
	# Save the power into big dataframe
	finalDF[[1]] <- cbind(finalDF[[1]], lengthRes)
	finalDF[[2]] <- cbind(finalDF[[2]], lengthRes)
	return(finalDF)
}