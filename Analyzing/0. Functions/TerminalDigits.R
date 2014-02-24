# Terminal Digit function
# Returns the terminal digit for a vector of numbers
#

TerminalDigits <- function(x){
  require(stringr)
  CharX <- as.character(x)
  NCharX <- nchar(CharX)
  Res <- str_sub(CharX,start=NCharX,end=-1L)
  return(as.numeric(Res))
}