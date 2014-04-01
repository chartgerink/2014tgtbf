simEffDist <- function(# Simulate effect size distributions
	### This function simulates effect size distributions under the constraints of the degrees of freedom. 
	### Defaults to simulating effect size distribution under H0; i.e., ES = 0.
	x,
	### A \code{statcheck} object.
	n.iter = 100,
	### The amount of simulations to be run
	es = 0,
	### The effect size in the population to simulate under. Specify this as an eta-squared (between 0 and 1).
	alpha = .05
	### The alpha used to determine what values are non-significant.
	){
	# Generate vector of selected test statistics
	sampledStats <- sample(1:dim(x)[1],
		n.iter,
		replace=T)
	testVal <- NULL
	pVal <- NULL
	for(i in 1:n.iter){
		temp <- x[sampledStats[i],]
		if(x$Statistic[sampledStats[i]]=="F"){
			testVal[i] <- ((x$df2[i]/x$df1[sampledStats[i]])*es)/(1-es)
			pVal[i] <- qf(testVal[i], x$df1[sampledStats[i]], x$df2[sampledStats[i]])
			} else if(x$Statistic[sampledStats[i]]=="t"){
				testVal[i] <- sqrt((x$df1[sampledStats[i]]*es)/(1-es))
				pVal[i] <- qt(testVal[i], x$df1[sampledStats[i]])
				} else if(x$Statistic[i]=="r"){
					testVal[i] <- sqrt(es)/sqrt((1-sqrt(es))(x$df1[sampledStats[i]]))
					pVal[i] <- qt(testVal[i], x$df1[sampledStats[i]])
				} else{
					testVal[i] <- NA
					pVal[i] <- NA
				}
	}
return(pVal)




	# Generate vector of p-values under no effect
	# pVal <- runif(n.iter,
	# 	min = alpha,
	# 	max = 1)
	# testValType <- x$Statistic
	# # Generate a dataframe
	# ddply()
	# for(i in 1:length(sampledStats)){
	# 	tVal <- qt(
	# 		pVal[i],
	# 		x$df1)
	# }
}


# ifelse(x$Statistic=="t",
# 	qt(
# 		pVal[i],
# 		x$df1),
# 	ifelse(
# 		x$Statistic=="F",
# 		qf(
# 			pVal[i],
# 			x$df1,x$df2),
# 		ifelse(
# 			x$Statistic=="r",
# 			qt(
# 				pVal[i],
# 				x$df1),
# 			NA
# 			)
# 		)
# 	)