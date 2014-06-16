prep.statcheck <- function(# Prepare a statcheck object for analysis
	### Statcheck objects have some formatting problems, that need to be solved before analyzing.
	### E.g., comma's as decimal points whereas R reads a period as decimal point.
	### Additionally, this function immediately calls the esComp.statcheck() to compute effect sizes
	x
	### A \code{statcheck} object.
	){
	
	# Removing out of bounds p-values
	selNA <- x$Computed>=1 & is.na(x$Computed)
	sum(selNA[!is.na(selNA)])
	x$Computed[selNA] <- NA
	
	# Replace all comma's with decimal points and make the variable numeric.
	x$Value <- suppressWarnings(as.numeric(sub(",",".",x$Value)))
	x$df1 <- suppressWarnings(as.numeric(sub(",",".",x$df1)))
	x$df2 <- suppressWarnings(as.numeric(sub(",",".",x$df2)))
	x$Computed <- suppressWarnings(as.numeric(sub(",",".",x$Computed)))
	
	# Computing unadjusted and adjusted effect sizes (OBSERVED)
	x <- cbind(x, esComp.statcheck(x))
	x$adjESComp[x$adjESComp<0 ] <- 0

	# Turning df1 for t and r into 1.
	x$df1[x$Statistic == "t" | x$Statistic == "r"] <- 1

	# making sure dfs are numeric
	x$df1 <- as.numeric(x$df1)
	x$df2 <- as.numeric(x$df2)

	# correcting false correlation value reporting
	x$esComp[x$esComp > 1 & !is.na(x$esComp)] <- NA

	# Select out irrefutably wrong df reporting
	x <- x[!x$df1==0 & !is.na(x$df1==0),]

	# Select out incorrectly extracted r values
	x <- x[!(x$Statistic=="r" & x$Value > 1),]

	# select out NA computed p-values
	x <- x[!is.na(x$Computed),]

	return(x)
}