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