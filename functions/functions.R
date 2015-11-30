effectPercent <- function(# Compute percentages of categories of effects.
	### Works only for eta-squared values.
	x,
	### Object of effects
	categories = c(.01, .06, .14)
	### Cut off points
	){
	categories <- c(categories, 1)
	percent <- matrix(ncol = length(categories))
	names(percent) <- categories
	x <- na.omit(x)
	for(i in 1:length(categories)){
		select <- (x <= categories[i])
		percent[i] <- (sum(select) / length(x))
	}
	return(percent)
}

TerminalDigits <- function(# Return terminal digit of vector.
	x
	### Vector of numeric values.
	){
  	require(stringr)
  	CharX <- as.character(x)
  	NCharX <- nchar(CharX)
  	Res <- str_sub(CharX,start=NCharX,end=-1L)
  	return(as.numeric(Res))
}

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
		esComp[i] <- ifelse(x$Statistic[sampledStats[i]]=="t",
			(qt(pVal[i], x$df2[sampledStats[i]])^2*(1 / x$df2[sampledStats[i]]))
			/ (((qt(pVal[i], x$df2[sampledStats[i]])^2*1) / x$df2[sampledStats[i]]) + 1),
			ifelse(
				x$Statistic[sampledStats[i]]=="F",
				(qf(pVal[i], x$df1[sampledStats[i]], x$df2[sampledStats[i]])*(x$df1[sampledStats[i]] / x$df2[sampledStats[i]]))
				/ (((qf(pVal[i], x$df1[sampledStats[i]], x$df2[sampledStats[i]])*x$df1[sampledStats[i]]) / x$df2[sampledStats[i]]) + 1),
				ifelse(
					x$Statistic[sampledStats[i]]=="r",
					sqrt(qt(pVal[i], x$df2[sampledStats[i]])^2
						/(qt(pVal[i], x$df2[sampledStats[i]])^2+x$df2[sampledStats[i]])^2),
					NA
					)
				)
			)
		adjESComp[i] <- ifelse(x$Statistic[sampledStats[i]]=="t",
			(qt(pVal[i], x$df2[sampledStats[i]])^2 * (1 / x$df2[sampledStats[i]]) - (1 / x$df2[sampledStats[i]])) 
			/ (((qt(pVal[i], x$df2[sampledStats[i]])^2*1) / x$df2[sampledStats[i]]) + 1),
			ifelse(
				x$Statistic[sampledStats[i]]=="F",
				(qf(pVal[i], x$df1[sampledStats[i]], x$df2[sampledStats[i]]) * (x$df1[sampledStats[i]] / x$df2[sampledStats[i]]) - (x$df1[sampledStats[i]] / x$df2[sampledStats[i]])) 
				/ (((qf(pVal[i], x$df1[sampledStats[i]], x$df2[sampledStats[i]])*x$df1[sampledStats[i]]) / x$df2[sampledStats[i]]) + 1),
				ifelse(
					x$Statistic[sampledStats[i]]=="r",
					sqrt(qt(pVal[i], x$df2[sampledStats[i]])^2/(qt(pVal[i], x$df2[sampledStats[i]])^2+x$df2[sampledStats[i]]))^2
					-((1-sqrt(qt(pVal[i], x$df2[sampledStats[i]])^2
						/(qt(pVal[i], x$df2[sampledStats[i]])^2+x$df2[sampledStats[i]]))^2)/x$df2[sampledStats[i]]),
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

r2z <- function(# Transform correlation coefficients into Fisher's Z coefficients
  ### Transforming correlation (r) coefficients into Fisher's Z coefficients for testing.
  r
  ### The (vector) of correlation coefficient(s)
  ){
  ##seealso<<
  ### \code{z2r}
  return(.5*log((1+r)/(1-r)))
}

z2r <- function(# Transform Fisher's Z coefficients into correlation coefficients
  ### Transforming Fisher's Z coefficients into correlations coefficients to get backtransformed results after testing.
  z
  ### The (vector) of Fisher Z coefficient(s)
  ){
  ##seealso<<
  ### \code{z2r}
  return((exp(2*z)-1)/(exp(2*z)+1))
}

nCalc <- function(# Compute sample sizes based on dfs.
  x
### A \code{statcheck} object.
){
  testN <- ifelse(x$Statistic=="t",
   x$df2+1,
   ifelse(
     x$Statistic=="F",
     x$df2-(x$df1+1),
     ifelse(
       x$Statistic=="r",
       x$df+2,
                       #                        ifelse(
                       #                          x$Statistic=="Chi2",
                       #                          sqrt(x$Statistic_value[x$Statistic=="Chi2"]/(x$df1[x$Statistic=="Chi2"]+1)),
                       NA
                       )
     )
   )
  #   )
return(testN)
}

sVarComp <- function(# Compute sampling variances for correlation coefficients
  ### Use this function to compute the sampling variances for a meta-analysis based on correlation coefficient effect sizes. 
  ### See the "See Also" section for the functions to use the resulting sampling variance in accordance with
  x
  ### A \code{statcheck} object
  ){
  ##seealso<<
  ##\link{ESComp}, \link{r2z}, \link{z2r}
  testN <- nCalc(x)
  sVar <- 1/(testN-3)
  return(sVar)
}

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

# Written by CHJ Hartgerink
# The Fisher method applied to test for deviation from uniformity
# In NONSIGNIFICANT P-values

FisherMethod <- function(# Compute Fisher's exact test for non-significant p-values.
	### This function computes paper level Fisher test statistics, testing whether the distribution of non-significant p-values is uniform. Significant values indicate deviation from uniformity. 
	### Returns both the normal Fisher test, as well as the complement test.
	### Computations are done for p*=log(p), where p is all non-significant p-values for each identifier.
	x,
	### Vector of p-values.
    id,
    ### Vector giving paper identifiers.
    alpha = .05
    ### Indicate what alpha level is being maintained for the study results, which serves as a cut-off for selecting the non-significant p-values.
    ){
	Res <- NULL
	for(i in 1:length(unique(id)))
	{
			selP <- x[id==unique(id)[i]]
			nSigP <- (na.omit(selP[selP>alpha])-alpha)/(1-alpha)
			SigP <- na.omit(selP[selP<=alpha])
			if(!length(nSigP)==0){
				# Compute the Fisher test statistic
				FMeth <- -2*sum(log(nSigP))
				# Compute p-values analytically
				pFMeth <- pchisq(q=FMeth, df=2*length(nSigP), lower.tail=F)
				} else {
					FMeth <- NA
					pFMeth <- NA
				}
			Res <- rbind(Res, data.frame(
				Fish = FMeth,
				PFish = pFMeth,
				CountNSig = length(nSigP),
				CountSig = length(SigP),
				PercentNonSig = length(nSigP)/length(selP)))
	}
	return(Res)
}

esComp <- function(
  x,
  df1,
  df2,
  esType,
  adjusted=TRUE){
 esComp <- ifelse(esType=="t",
   (x^2*(1 / df2)) / (((x^2*1) / df2) + 1),
   ifelse(
     esType=="F",
     (x*(df1 / df2)) / (((x*df1) / df2) + 1),
     ifelse(
       esType=="r",
       x^2,
                       #                        ifelse(
                       #                          esType=="Chi2",
                       #                          sqrt(x[esType=="Chi2"]/(df1[esType=="Chi2"]+1)),
                       NA
                       )
     )
   )
  #   )
adjESComp <- ifelse(esType=="t",
  (x^2 * (1 / df2) - (1 / df2)) / (((x^2*1) / df2) + 1),
  ifelse(
    esType=="F",
    (x * (df1 / df2) - (df1 / df2)) / (((x*df1) / df2) + 1),
    ifelse(
      esType=="r",
      x^2-((1-x^2)/df2),
                          #                        ifelse(
                          #                          esType=="Chi2",
                          #                          sqrt(x[esType=="Chi2"]/(x$df1[x$esType=="Chi2"]+1)),
                          NA
                          )
    )
  )
return(cbind(esComp,adjESComp))
}

esComp.statcheck <- function(# Function to compute effect sizes for a \code{statcheck} object
  x,
  ### \code{statcheck} object
  adjusted=TRUE
  ### Logical for giving adjusted effect sizes
  ){
  esComp <- ifelse(x$Statistic=="t",
   (x$Value^2*(1 / x$df2)) / (((x$Value^2*1) / x$df2) + 1),
   ifelse(
     x$Statistic=="F",
     (x$Value*(x$df1 / x$df2)) / (((x$Value*x$df1) / x$df2) + 1),
     ifelse(
       x$Statistic=="r",
       x$Value^2,
                       #                        ifelse(
                       #                          x$Statistic=="Chi2",
                       #                          sqrt(x$Value[x$Statistic=="Chi2"]/(x$df1[x$Statistic=="Chi2"]+1)),
                       NA
                       )
     )
   )
  #   )
adjESComp <- ifelse(x$Statistic=="t",
  (x$Value^2 * (1 / x$df2) - (1 / x$df2)) / (((x$Value^2*1) / x$df2) + 1),
  ifelse(
    x$Statistic=="F",
    (x$Value * (x$df1 / x$df2) - (x$df1 / x$df2)) / (((x$Value*x$df1) / x$df2) + 1),
    ifelse(
      x$Statistic=="r",
      x$Value^2-((1-x$Value^2)/x$df2),
                          #                        ifelse(
                          #                          x$Statistic=="Chi2",
                          #                          sqrt(x$Value[x$Statistic=="Chi2"]/(x$df1[x$Statistic=="Chi2"]+1)),
                          NA
                          )
    )
  )
return(cbind(esComp,adjESComp))
}