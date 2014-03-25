esComp <- function(
  x,
  df1,
  df2,
  esType,
  adjusted=TRUE){
 esComp <- ifelse(esType=="t",
   (x^2*(1 / df1)) / (((x^2*1) / df1) + 1),
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
  (x^2 * (1 / df1) - (1 / df1)) / (((x^2*1) / df1) + 1),
  ifelse(
    esType=="F",
    (x * (df1 / df2) - (df1 / df2)) / (((x*df1) / df2) + 1),
    ifelse(
      esType=="r",
      x^2-((1-x^2)/df1),
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
  esComp <- ifelse(x$test_statistic=="t",
   (x$test_statistic_value^2*(1 / x$df1)) / (((x$test_statistic_value^2*1) / x$df1) + 1),
   ifelse(
     x$test_statistic=="F",
     (x$test_statistic_value*(x$df1 / x$df2)) / (((x$test_statistic_value*x$df1) / x$df2) + 1),
     ifelse(
       x$test_statistic=="r",
       x$test_statistic_value^2,
                       #                        ifelse(
                       #                          x$test_statistic=="Chi2",
                       #                          sqrt(x$test_statistic_value[x$test_statistic=="Chi2"]/(x$df1[x$test_statistic=="Chi2"]+1)),
                       NA
                       )
     )
   )
  #   )
adjESComp <- ifelse(x$test_statistic=="t",
  (x$test_statistic_value^2 * (1 / x$df1) - (1 / x$df1)) / (((x$test_statistic_value^2*1) / x$df1) + 1),
  ifelse(
    x$test_statistic=="F",
    (x$test_statistic_value * (x$df1 / x$df2) - (x$df1 / x$df2)) / (((x$test_statistic_value*x$df1) / x$df2) + 1),
    ifelse(
      x$test_statistic=="r",
      x$test_statistic_value^2-((1-x$test_statistic_value^2)/x$df1),
                          #                        ifelse(
                          #                          x$test_statistic=="Chi2",
                          #                          sqrt(x$test_statistic_value[x$test_statistic=="Chi2"]/(x$df1[x$test_statistic=="Chi2"]+1)),
                          NA
                          )
    )
  )
return(cbind(esComp,adjESComp))
}