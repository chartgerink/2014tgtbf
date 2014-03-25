powerCalc <- function(
	x,
	### A \code{statcheck} object
	effectSize = 0,
	### Vector of effect sizes to estimate power under
	alpha = .05,
	### The alpha used to determine rejection areas for ORIGINAL tests
	twoTailed = TRUE,
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
	temp <- matrix(nrow=length(unique(x$pap_id)),ncol=length(effectSize))
	finalDF <- as.data.frame(temp)
	names(finalDF) <- effectSize
	rm(temp)

	# Simulate effects under the effect sizes for each test statistic to determine power of the test
	for(p in 1:length(unique(x$pap_id))){
		selectStats <- x[x$pap_id==p,]
		tempMat <- matrix(nrow=n.iter,ncol=dim(selectStats)[1])
		for(s in 1:dim(selectStats)[1]){
			for(es in 1:length(effectSize)){
				res <- NULL
				resP <- NULL
				if(selectStats$test_statistic[s]=="t"){
					cv <- qf((1-alpha/2),1,selectStats$df1[s])
					ncp <- (effectSize[es]/1-effectSize[es])*(1+selectStats$df1[s]+1)
					testBeta <- pf(cv,1,selectStats$df1[s],ncp)
					while(length(resP)<n.iter){
						temp <- rf(1, 1, selectStats$df1[s], ncp)
						if(temp < testBeta){
							# res <- c(res, temp)
							resP <- c(resP, pf(temp,1,selectStats$df1[s]))
						}
					}
					# tempMat[,es] <- resP
					# tempMat[s,es] <- length(resP[resP<.1])/length(resP)
					} else if(selectStats$test_statistic[s]=="F"){
						cv <- qf((1-alpha/2),selectStats$df1[s],selectStats$df2[s])
						ncp <- (effectSize[es]/1-effectSize[es])*(selectStats$df1[s]+selectStats$df2[s]+1)
						testBeta <- pf(cv,selectStats$df1[s],selectStats$df2[s],ncp)
						while(length(resP)<n.iter){
							temp <- rf(1, selectStats$df1[s], selectStats$df2[s], ncp)
							if(temp < testBeta){
								res <- c(res, temp)
								resP <- c(resP, pf(temp,selectStats$df1[s],selectStats$df2[s]))
							}
						}
						# tempMat[,es] <- resP
						# tempMat[s,es] <- length(resP[resP<.1])/length(resP)
						} else if(selectStats$test_statistic[s]=="r"){
							cv <- qf((1-alpha/2),1,selectStats$df1[s])
							ncp <- (effectSize[es]/1-effectSize[es])*(1+selectStats$df1[s]+1)
							testBeta <- pf(cv,1,selectStats$df1[s],ncp)
							while(length(resP)<n.iter){
								temp <- rf(1, 1, selectStats$df1[s], ncp)
								if(temp < testBeta){
									res <- c(res, temp)
									resP <- c(resP, pf(temp,1,selectStats$df1[s]))
								}
							}
							# tempMat[,es] <- resP
							# tempMat[s,es] <- length(resP[resP<.1])/length(resP)

							} else {
								# tempMat[s,es] <- NA
								resP <- NA

							}
						}
					}
				}
				return(resP)
			}







# function(
# 	x,
# 	### A \code{statcheck} object with computed effect sizes.
# 	effectSize = 0,
# 	### Effect size to use (eta2 ONLY!)
# 	alpha = .05,
# 	### Alpha level that functions as cutoff to calculate Beta value for the distribution.
# 	n.iter = 1000,
# 	### Number of iterations sampling from the noncentral distribution.
# 	twoTailed = TRUE
# 	### Logical used to calculate the critical values appropriately.
# 	){
# 	##seealso<<
# 	##\link{esComp}
# 	if(twoTailed==TRUE){
# 		if(x$test_statistic=="t"){
# 			cv <- qf((1-alpha/2),1,x$df1)
# 			ncp <- (x$esComp/1-x$esComp)*(1+x$df1+1)
# 			testBeta <- pf(cv,1,x$df1,ncp)
# 			res <- NULL
# 			temp <- NULL
# 			resP <- NULL
# 			while(length(res)<=n.iter){
# 				temp <- rf(1, 1, x$df1, ncp)
# 				if(temp < testBeta){
# 					res <- c(res, temp)
# 					resP <- c(resP, pf(temp,1,x$df1))
# 				}
# 			}
# 		}
# 		else if(x$test_statistic=="F"){
# 			cv <- qf((1-alpha/2),x$df1,x$df2)
# 			ncp <- (x$esComp/1-x$esComp)*(x$df1+x$df2+1)
# 			testBeta <- pf(cv,x$df1,x$df2,ncp)
# 			res <- NULL
# 			temp <- NULL
# 			while(length(res)<=n.iter){
# 				temp <- rf(1, df1, df2, ncp)
# 				if(temp < testBeta){
# 					res <- c(res, temp)
# 					resP <- c(resP, pf(temp,x$df1,x$df2))
# 				}
# 			}
# 		}
# 		else if(x$test_statistic=="r"){
# 			cv <- qf((1-alpha/2),1,x$df1)
# 			ncp <- (x$esComp/1-x$esComp)*(1+x$df1+1)
# 			testBeta <- pf(cv,1,x$df1,ncp)
# 			res <- NULL
# 			temp <- NULL
# 			while(length(res)<=n.iter){
# 				temp <- rf(1, 1, df1, ncp)
# 				if(temp < testBeta){
# 					res <- c(res, temp)
# 					resP <- c(resP, pf(temp,1,x$df1))
# 				}
# 			}
# 		}

# 	}
# 		# else{
# 			# 
# 		# }

# 		print(res)
# 	}