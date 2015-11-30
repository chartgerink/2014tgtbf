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