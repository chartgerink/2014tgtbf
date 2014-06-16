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