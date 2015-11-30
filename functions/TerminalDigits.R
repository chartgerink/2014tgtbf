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