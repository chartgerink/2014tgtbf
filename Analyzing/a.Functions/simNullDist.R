simNullDist <- function(# Simulate null distribution of effect sizes
	### This function generates a distribution of effect sizes under the null that there is no effect size.
	### This is done by sampling the test statistic, degrees of freedom, and selecting a uniformly distributed p-value.
	### Subsequently, the accompanying test statistic is computed and its accompanying effect.
	### This function can be used to generate these effect sizes and plot these against the observed effect sizes.
	x,
	### A statcheck object
	n.iter = 100,
	### The amount of simulations to be run. Recommend this to be a high value, minimum 1 million.
	alpha = .05
	### The alpha to be used to select non-significant values
	){
	##details<<
	##Only works for t, F, and r right now.
	##seealso<<
	##\link{statcheck}
	# Generate vector of selected test statistics
	sel <- x$Computed >= alpha
	x <- x[sel,]
	
	sampledStats <- sample(1:dim(x)[1],
		n.iter,
		replace = T)
	# Generate vector of p-values under no effect (uniform)
	pVal <- runif(n.iter,
		min = alpha,
		max = 1)
	# Compute effect sizes under the generated p-values
	esComp <- NULL
	adjESComp <- NULL
	# testStat <- NULL
	# testVal <- NULL
	for(i in 1:n.iter){
		esComp[i] <- ifelse(x$Statistic[i]=="t",
			(qt(pVal[i], x$df2[i])^2*(1 / x$df2[i]))
			/ (((qt(pVal[i], x$df2[i])^2*1) / x$df2[i]) + 1),
			ifelse(
				x$Statistic[i]=="F",
				(qf(pVal[i], x$df1[i], x$df2[i])*(x$df1[i] / x$df2[i]))
				/ (((qf(pVal[i], x$df1[i], x$df2[i])*x$df1[i]) / x$df2[i]) + 1),
				ifelse(
					x$Statistic[i]=="r",
					sqrt(qt(pVal[i], x$df2[i])^2
						/(qt(pVal[i], x$df2[i])^2+x$df2[i])^2),
					NA
					)
				)
			)
		adjESComp[i] <- ifelse(x$Statistic[i]=="t",
			(qt(pVal[i], x$df2[i])^2 * (1 / x$df2[i]) - (1 / x$df2[i])) 
			/ (((qt(pVal[i], x$df2[i])^2*1) / x$df2[i]) + 1),
			ifelse(
				x$Statistic[i]=="F",
				(qf(pVal[i], x$df1[i], x$df2[i]) * (x$df1[i] / x$df2[i]) - (x$df1[i] / x$df2[i])) 
				/ (((qf(pVal[i], x$df1[i], x$df2[i])*x$df1[i]) / x$df2[i]) + 1),
				ifelse(
					x$Statistic[i]=="r",
					sqrt(qt(pVal[i], x$df2[i])^2/(qt(pVal[i], x$df2[i])^2+x$df2[i]))^2
					-((1-sqrt(qt(pVal[i], x$df2[i])^2
						/(qt(pVal[i], x$df2[i])^2+x$df2[i]))^2)/x$df2[i]),
					NA
					)
				)
			)
		print( paste0("still working, now at iteration ",i))
	}

	# Trying to make this with ddply and not looping it for speed, not functional yet
	# res <- 	ddply(copilot,.(Statistic), .fun = function(copilot) ifelse(x$Statistic=="t",
	# 	qt(pVal, x$df1),
	# 	ifelse(x$Statistic=="F",
	# 		qf(pVal, x$df1, x$df2),
	# 		ifelse(x$Statistic=="r",
	# 			qt(pVal, x$df1),
	# 			NA))))
	# apply(x, 1, function(x) ))

	# Compute effect size for test-values

	return(data.frame(
		esComp = esComp,
		adjESComp = adjESComp))
}