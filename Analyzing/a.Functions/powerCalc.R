powerCalc <- function(
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
	# Calculate effects just in case the statcheck object doesn't have them yet
	x <- cbind(x,esComp.statcheck(x))

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
	# Simulate effects under the effect sizes for each test statistic to determine power of the test
	for(p in 1:length(unique(x$Source))){
		selectStats <- x[x$Source==p,]
		tempMat[[p]] <- matrix(nrow=n.iter,ncol=dim(selectStats)[1])

		for(es in 1:length(effectSize)){
			for(s in 1:dim(selectStats)[1]){
				res <- NULL
				resP <- NULL

				if(selectStats$Statistic[s]=="t"){
					cv <- qf(alpha,1,selectStats$df1[s],lower.tail=F)
					ncp <- (effectSize[es]/(1-effectSize[es]))*(1+selectStats$df1[s]+1)
					prob <- pf(cv,1,selectStats$df1[s],ncp=ncp,lower.tail=T)
					sampleProb <- runif(n.iter,0,prob)
					sampleTest <- qf(sampleProb,1,selectStats$df1[s],ncp=ncp)
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
				tempMat[[p]][,s] <- (resP-.05)/(1-alpha)
				print(paste("Still working, no worries, paper", p, "statistic", s, sep=" "))
			}
			if(sum(is.na(tempMat[[p]]))==(dim(tempMat[[p]])[1]*dim(tempMat[[p]])[2])) {
				finalDF[[1]][p,es] <- NA
				finalDF[[2]][p,es] <- NA
				} else{
					sel[[p]] <- as.matrix(as.matrix(tempMat[[p]])[,colSums(is.na(tempMat[[p]])) != nrow(tempMat[[p]])])
					fishTest[[p]] <- apply(sel[[p]],1,function(x) -sum(log(x)))
					fishTestCompl[[p]] <- apply(sel[[p]],1,function(x) -sum(log(1-(x))))
					pFishTest[[p]] <- pgamma(fishTest[[p]],dim(sel[[p]])[2],lower.tail=F)
					pFishTestCompl[[p]] <- pgamma(fishTestCompl[[p]],dim(sel[[p]])[2],lower.tail=F)
					propSigFishTest[p] <- sum(na.omit(pFishTest[[p]]) < testAlpha)/length(na.omit(pFishTest[[p]]))
					propSigFishTestCompl[p] <- sum(na.omit(pFishTestCompl[[p]]) < testAlpha)/length(na.omit(pFishTestCompl[[p]]))
					finalDF[[1]][p,es] <- propSigFishTest[p]
					finalDF[[2]][p,es] <- propSigFishTestCompl[p]
				}
			}
		}
		return(finalDF)
	}