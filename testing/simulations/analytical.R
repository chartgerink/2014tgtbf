alpha <- .05
alphaF <- .10
ES <- seq(0.00, .99, 0.005)

setwd("C:/Users/Chris/Dropbox/CJM/Masterproject/testing")


fCV0 <- qf(p=alpha, df1=2-1, df2=50-2, lower.tail=F)
fCVa <- qf(p=(1-alpha)*(1-alphaF), df1=2-1, df2=50-2)


f2 <- ES^2/(1-ES^2)
ncp <- f2*50

numer <- pf(q=fCV0, df1=2-1, df2=50-2, ncp=ncp, lower.tail=T) - pf(q=fCVa, df1=2-1, df2=50-2, ncp=ncp, lower.tail=T)
denom <- pf(q=fCV0, df1=2-1, df2=50-2, ncp=ncp, lower.tail=T)


write.csv2(cbind(ES, numer/denom), 'analytical.csv')

