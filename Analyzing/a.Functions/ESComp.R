ESComp <- function(# Function to compute effect sizes for a \code{statcheck} object
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
                       x$test_statistic_value[x$test_statistic=="r"]^2,
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
                          x$test_statistic_value[x$test_statistic=="r"]^2,# Need to refine for adjusted effect size
                          #                        ifelse(
                          #                          x$test_statistic=="Chi2",
                          #                          sqrt(x$test_statistic_value[x$test_statistic=="Chi2"]/(x$df1[x$test_statistic=="Chi2"]+1)),
                          NA
                        )
                      )
  )
  return(cbind(esComp,adjESComp))
}