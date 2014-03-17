simEffDist <- function(# Simulate effect size distributions
	### This function simulates effect size distributions under the constraints of the degrees of freedom. 
	### Defaults to simulating effect size distribution under H0; i.e., ES = 0.
	x,
	### A \code{statcheck} object.
	n.iter = 100,
	### The amount of simulations to be run
	es = 0,
	### The effect size in the population to simulate under.
	alpha = .05,
	### The alpha used to determine what values are non-significant.
	){
	# Generate vector of selected test statistics
	sampledStats <- sample(1:dim(x)[1],
		n.iter,
		replace=T)
	# Generate vector of p-values under no effect
	pVal <- runif(n.iter,
		min = alpha,
		max = 1)
	testValType <- x$test_statistic
	# Generate a dataframe
	ddply()
	for(i in 1:length(sampledStats)){
		tVal <- qt(
			pVal[i],
			x$df1)
	}
}


ifelse(x$test_statistic=="t",
	qt(
		pVal[i],
		x$df1),
	ifelse(
		x$test_statistic=="F",
		qf(
			pVal[i],
			x$df1,x$df2),
		ifelse(
			x$test_statistic=="r",
			qt(
				pVal[i],
				x$df1),
			NA
			)
		)
	)