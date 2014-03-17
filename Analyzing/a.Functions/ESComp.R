esComp <- function(
  x,
  df1,
  df2,
  esType,
  adjusted){
  if(esType == "t"){
(x^2*(1 / df1)) / (((x^2*1) / df1) + 1)
  } else if(esType == "F"){
(x*(df1 / df2)) / (((test_statistic_value*df1) / df2) + 1)
    } else if(esType == "r"){

    }
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