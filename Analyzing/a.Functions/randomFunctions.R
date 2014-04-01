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
   x$df1+1,
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