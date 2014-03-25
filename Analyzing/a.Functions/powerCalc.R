function(
	x,
	### A \code{statcheck} object with computed effect sizes.
	effectSize = ,
	### Effect size to use (eta2 ONLY!)
	alpha = .05,
	### Alpha level that functions as cutoff to calculate Beta value for the distribution.
	n.iter = 1000,
	### Number of iterations sampling from the noncentral distribution.
	twoTailed = TRUE
	### Logical used to calculate the critical values appropriately.
	){
	##seealso<<
	##\link{esComp}
	if(twoTailed==TRUE){
		if(x$test_statistic=="t"){

		}
		else if(x$test_statistic=="F"){
		cv <- qf((1-alpha/2),x$df1,x$df2)
		testBeta <- pf(cv,x$df1,x$df2,ncp)
		ncp <- 1/(1-x$esComp)
		resF <- NULL
		temp <- NULL

		while(length(resF)<=n.iter){
			temp <- rf(1, df1, df2, ncp)
			if(temp < zBeta){
				resF <- c(resF, temp)
			}
		}
		}
		else if(x$test_statistic=="r"){

		}

		} else{

		}
	}