powerCalc <- function(
	x,
	### A \code{statcheck} object
	effectSize = 0,
	### Vector of effect sizes to estimate power under
	alphaF = .1,
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

	# Simulate effects under the effect sizes for each test statist ic to determine power of the test
	# For each paper
	power <- NULL
	kPaper <- NULL
	jourPaper <- NULL

	for(p in 1:length(unique(x$Source))){
		# Select the stats for that paper and nonsignificant
		selectStats <- x[x$Source==unique(x$Source)[p] & x$Computed >= alpha & (x$Statistic == "t" | x$Statistic == "F" | x$Statistic == "r"),]
		kPaper <- rbind(kPaper, length(selectStats$Computed))
		if(dim(selectStats)[1] == 0){
			jourPaper <- rbind(jourPaper, unique(x$journal.jour.[x$Source==unique(x$Source)[p]])[1])
		} else {jourPaper <- rbind(jourPaper, unique(selectStats$journal.jour.)[1])}
		
		pow <- 0
		powe <- NULL

		# Loop through the effects
		for(es in 1:length(effectSize)){
			chiF <- NULL
      		chiP <- NULL
			if(pow < .995){
				# Step 1 - critical value
				fCV <- qf(p=alpha, df1=selectStats$df1, df2=selectStats$df2, lower.tail=F)

				# Step 2 - non-centrality parameter
				f2 <- effectSize[es]^2/(1-effectSize[es]^2)
				ncp <- f2*ifelse(selectStats$Statistic=="F", selectStats$df1+selectStats$df2+1, selectStats$df1+selectStats$df2)

				# Step 3 - beta
				beta <- pf(q=fCV, df1=selectStats$df1, df2=selectStats$df2, ncp=ncp)
				
				# Simulate results n.iter times
				for(i in 1:n.iter){
					
					# Step 4 - pValues under alternative distribution
					pA <- runif(length(beta), min=0, max=beta)

					# Step 5 - fValues simulated under alternative
					fA <- qf(p=pA, df1=selectStats$df1, df2=selectStats$df2, ncp=ncp, lower.tail=T)

					# Step 6 - pValues under null
					p0 <- pf(q=fA, df1=selectStats$df1, df2=selectStats$df2, lower.tail=F)

					# Step 7 - Transform pValues back into state space [0;1]
					p0 <- (p0-alpha)/(1-alpha)

					# Step 8 - compute fisher statistic
					chiF <- -2*sum(log(p0))
					chiP[i] <- pchisq(q=chiF, df=2*length(p0), lower.tail=F)
				}
				# Step 9 - compute power
				pow <- sum(chiP < alphaF) / n.iter
				# print(paste(p, "of", length(unique(x$Source)), effectSize[es], "1-beta", pow, sep=" "))
			}
			else{
				pow <- 1
			}
			powe <- cbind(powe, pow)
		}
		power <- rbind(power, powe)
	}

	final <- cbind(jourPaper, kPaper, power)
	return(final)
}
