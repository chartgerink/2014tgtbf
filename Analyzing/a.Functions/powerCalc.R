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

	# Calculate effects just in case the statcheck object doesn't have them yet
	x <- cbind(x,esComp.statcheck(x))

	# Create object to be filled with all data
	finalDF <- NULL
	temp <- matrix(nrow=length(unique(x$pap_id)),ncol=length(effectSize))
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
	for(p in 1:length(unique(x$pap_id))){
		selectStats <- x[x$pap_id==p,]
		tempMat[[p]] <- matrix(nrow=n.iter,ncol=dim(selectStats)[1])

		for(es in 1:length(effectSize)){
			for(s in 1:dim(selectStats)[1]){
				res <- NULL
				resP <- NULL

				if(selectStats$test_statistic[s]=="t"){
					cv <- qf(alpha,1,selectStats$df1[s],lower.tail=F)
					ncp <- (effectSize[es]/(1-effectSize[es]))*(1+selectStats$df1[s]+1)
					prob <- pf(cv,1,selectStats$df1[s],ncp=ncp,lower.tail=T)
					sampleProb <- runif(n.iter,0,prob)
					sampleTest <- qf(sampleProb,1,selectStats$df1[s],ncp=ncp)
					resP <- pf(sampleTest,1,selectStats$df1[s],lower.tail=F)
				} 
				else if(selectStats$test_statistic[s]=="F"){
					cv <- qf(alpha,selectStats$df1[s],selectStats$df2[s],lower.tail=F)
					ncp <- (effectSize[es]/(1-effectSize[es]))*(selectStats$df1[s]+selectStats$df2[s]+1)
					prob <- pf(cv,selectStats$df1[s],selectStats$df2[s],ncp=ncp,lower.tail=T)
					sampleProb <- runif(n.iter,0,prob)	
					sampleTest <- qf(sampleProb,selectStats$df1[s],selectStats$df2[s],ncp=ncp)
					resP <- pf(sampleTest,selectStats$df1[s],selectStats$df2[s],lower.tail=F)
				}
				else if(selectStats$test_statistic[s]=="r"){
					cv <- qf(alpha,1,selectStats$df1[s],lower.tail=F)
					ncp <- (effectSize[es]/(1-effectSize[es]))*(1+selectStats$df1[s]+1)
					prob <- pf(cv,1,selectStats$df1[s],ncp=ncp,lower.tail=T)
					sampleProb <- runif(n.iter,0,prob)
					sampleTest <- qf(sampleProb,1,selectStats$df1[s],ncp=ncp)
					resP <- pf(sampleTest,1,selectStats$df1[s],lower.tail=F)
				}


				else {
					while(length(resP) < n.iter){
						resP <- c(resP, NA)
					}
				}
				tempMat[[p]][,s] <- resP/.95
				# minimumP <- c(minimumP,min(resP))
				print(paste("Still working, no worries, paper", p, "statistic", s, sep=" "))
				# print(resP)
			}

			require(gtools)
			if(invalid(tempMat[[p]][!is.na(colSums(tempMat[[p]])),])){

				} else{
					tempMat[[p]][,!is.na(colSums(tempMat[[p]]))]}
			# if(dim(tempMat[[p]])[2]==1){tempMat[[p]] <- as.vector(tempMat[[p]])}
			# Just making sure all non t, F, r columns are removed
			# sel[p] <- sum(colSums(is.na(tempMat[[p]])) == nrow(tempMat[[p]]))==ncol(tempMat[[p]])
			# if(sel==TRUE){
			# 	fishTest[p] <- NA
			# 	fishTestCompl[p] <- NA
			# 	} else {
			# 		sel <- tempMat[[p]][,colSums(is.na(tempMat[[p]])) != nrow(tempMat[[p]])]
			fishTest[p] <- apply(tempMat[[p]],1,function(x) -sum(log(x)))
					# fishTestCompl[p] <- apply(tempMat[[p]],1,function(x) -sum(log(1-(x))))			
			# 	}
				# pFishTest <- pgamma(fishTest,dim(sel)[2],lower.tail=F)
				# pFishTestCompl <- pgamma(fishTestCompl,dim(sel)[2],lower.tail=F)
				# propSigFishTest <- sum(na.omit(pFishTest) < testAlpha)/length(na.omit(pFishTest))
				# propSigFishTestCompl <- sum(na.omit(pFishTestCompl) < testAlpha)/length(na.omit(pFishTestCompl))
				# finalDF[[1]][p,es] <- propSigFishTest
				# finalDF[[2]][p,es] <- propSigFishTestCompl
						# if(dim(na.omit(tempMat[[p]]))[1]==0){
						# 	finalDF[[1]][p,es] <- NA
						# 	finalDF[[2]][p,es] <- NA
						# 	} else{

							# }
						}
					}
					return(tempMat)
				}

		# 			for(p in 1:length(unique(x$pap_id))){
		# selectStats <- x[x$pap_id==p,]
		# tempMat[[p]] <- matrix(nrow=n.iter,ncol=dim(selectStats)[1])

		# for(es in 1:length(effectSize)){
		# 	for(s in 1:dim(selectStats)[1]){
		# 		res <- NULL
		# 		resP <- NULL

		# 		if(selectStats$test_statistic[s]=="t"){
		# 			cv <- qf((1-alpha),1,selectStats$df1[s],lower.tail=F)
		# 			ncp <- (effectSize[es]/1-effectSize[es])*(1+selectStats$df1[s]+1)
		# 			# testBeta <- pf(cv,1,selectStats$df1[s],ncp)
		# 			while(length(resP)<n.iter){
		# 				temp <- rf(1, 1, selectStats$df1[s], ncp)
		# 				if(temp < cv){
		# 					res <- c(res, temp)
		# 					resP <- c(resP, ((1-pf(temp,1,selectStats$df1[s]))))
		# 				}
		# 			}

		# 			} else if(selectStats$test_statistic[s]=="F"){
		# 				cv <- qf((1-alpha),selectStats$df1[s],selectStats$df2[s],lower.tail=F)
		# 				ncp <- (effectSize[es]/1-effectSize[es])*(selectStats$df1[s]+selectStats$df2[s]+1)
		# 				testBeta <- pf(cv,selectStats$df1[s],selectStats$df2[s],ncp)
		# 				while(length(resP)<n.iter){
		# 					temp <- rf(1, selectStats$df1[s], selectStats$df2[s], ncp)
		# 					if(temp < cv){
		# 						res <- c(res, temp)
		# 						resP <- c(resP, ((1-pf(temp,selectStats$df1[s],selectStats$df2[s]))))
		# 					}
		# 				}

		# 				} else if(selectStats$test_statistic[s]=="r"){
		# 					cv <- qf((1-alpha),1,selectStats$df1[s],lower.tail=F)
		# 					ncp <- (effectSize[es]/1-effectSize[es])*(1+selectStats$df1[s]+1)
		# 					testBeta <- pf(cv,1,selectStats$df1[s],ncp)
		# 					while(length(resP)<n.iter){
		# 						temp <- rf(1, 1, selectStats$df1[s], ncp)
		# 						if(temp < cv){
		# 							res <- c(res, temp)
		# 							resP <- c(resP, ((1-pf(temp,1,selectStats$df1[s]))))
		# 						}
		# 					}

		# 					} else {
		# 						while(length(resP) < n.iter){
		# 							resP <- c(resP, NA)
		# 						}
		# 					}
		# 					tempMat[[p]][,s] <- resP
		# 					# minimumP <- c(minimumP,min(resP))
		# 					print(paste("Still working, no worries, paper", p, "statistic", s, sep=" "))
		# 					# print(resP)
		# 				}
		# 				if(dim(na.omit(tempMat[[p]]))[1]==0){
		# 					finalDF[[1]][p,es] <- NA
		# 					finalDF[[2]][p,es] <- NA
		# 					} else{

		# 						fishTest <- apply(tempMat[[p]],1,function(x) -sum(log(x/.95)))
		# 						fishTestCompl <- apply(tempMat[[p]],1,function(x) -sum(log(1-(x/.95))))
		# 						pFishTest <- pgamma(fishTest, dim(tempMat[[p]])[2])
		# 						pFishTestCompl <- pgamma(fishTestCompl, dim(tempMat[[p]])[2])
		# 						propSigFishTest <- sum(na.omit(pFishTest<testAlpha))/length(na.omit(pFishTest))
		# 						propSigFishTestCompl <- sum(na.omit(pFishTestCompl<testAlpha))/length(na.omit(pFishTestCompl))
		# 						finalDF[[1]][p,es] <- propSigFishTest
		# 						finalDF[[2]][p,es] <- propSigFishTestCompl
		# 					}
		# 				}
		# 			}
		# 			return(c(temp,cv))
		# 		}