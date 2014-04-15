prep.statcheck <- function(# Prepare a statcheck object for analysis
	### Statcheck objects have some formatting problems, that need to be solved before analyzing.
	### E.g., comma's as decimal points whereas R reads a period as decimal point.
	### Additionally, this function immediately calls the esComp.statcheck() to compute effect sizes
	x
	### A \code{statcheck} object.
	){
	# Removing out of bounds p-values
	selNA <- x$Computed>=1
	sum(selNA[!is.na(selNA)])
	x$Computed[selNA] <- NA
	# Replace all comma's with decimal points and make the variable numeric.
	x$Value <- suppressWarnings(as.numeric(sub(",",".",x$Value)))
	x$df1 <- suppressWarnings(as.numeric(sub(",",".",x$df1)))
	x$df2 <- suppressWarnings(as.numeric(sub(",",".",x$df2)))
	# Computing unadjusted and adjusted effect sizes (OBSERVED)
	x <- cbind(x, esComp.statcheck(x))

	return(x)
}