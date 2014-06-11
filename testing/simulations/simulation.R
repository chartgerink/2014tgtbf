N <- 25
groups <- 2
P <- 1:6
ES <- seq(0.00, .99, 0.005)
n.iter <- 10000
alpha <- .05
alphaF <- .1

setwd("C:/Users/Chris/Dropbox/CJM/Masterproject/testing")

# p-alpha
power <- NULL
for(k in 1:length(P)){
  pow <- 0
  powe <- NULL
  chiP <- NULL
  for(es in 1:length(ES)){
    if(pow < .995){
      for(i in 1:n.iter){
        # Step 1 - Critical value
        fCV <- qf(p=alpha, df1=groups-1, df2=(N*groups)-groups,lower.tail=F)
        
        # Step 2 - non-centrality parameter
        f2 <- ES[es]^2/(1-ES[es]^2)
        ncp <- f2*(N*groups)
        
        # Step 3 - beta
        beta <- pf(q=fCV, df1=groups-1, df2=(N*groups)-groups, ncp=ncp, lower.tail=T)
        
        # Step 4 - pValues under alternative distribution
        pA <- runif(n=P[k],min=0,max=beta)
        
        # Step 5 - tValues simulated under alternative
        fA <- qf(p=pA, df1=groups-1, df2=(N*groups)-groups, ncp=ncp, lower.tail=T)
        
        # Step 6 - pValues under null
        p0 <- pf(q=fA, df1=groups-1, df2=(N*groups)-groups,lower.tail=F)
        
        # Step 7 - transform p-values back into state space [0;1]
        p0 <- (p0-alpha)/(1-alpha)
        
        # Step 8 - computing fisher statistic
        chiF <- -2*sum(log(p0))
        chiP[i] <- pchisq(q=chiF, df=2*length(p0), lower.tail=F)
      }
      # Step 9 - computing the proportion significant
      pow <- sum(chiP < alphaF)/n.iter
      print(paste(N, P[k], ES[es], pow,sep=" "))
    } else{
      # For computing efficiency all subsequent effect sizes
      # are rounded to one if previous was .995 or higher
      pow <- 1
    }
#     if(pow == 0){pow <- 1}
    powe <- cbind(powe,
                  pow)
  }
  power <- rbind(power, powe)
}
filepath <- paste0('pAlpha',
                   '.csv')
write.csv2(t(power),
           file=filepath)

# 1-[p-alpha]
power <- NULL
for(k in 1:length(P)){
  pow <- 0
  powe <- NULL
  chiP <- NULL
  for(es in 1:length(ES)){
    if(pow < .995){
      for(i in 1:n.iter){
        # Step 1 - Critical value
        fCV <- qf(p=alpha, df1=groups-1, df2=(N*groups)-groups,lower.tail=F)
        
        # Step 2 - non-centrality parameter
        f2 <- ES[es]^2/(1-ES[es]^2)
        ncp <- f2*(N*groups)
        
        # Step 3 - beta
        beta <- pf(q=fCV, df1=groups-1, df2=(N*groups)-groups, ncp=ncp, lower.tail=T)
        
        # Step 4 - pValues under alternative distribution
        pA <- runif(n=P[k],min=0,max=beta)
        
        # Step 5 - tValues simulated under alternative
        fA <- qf(p=pA, df1=groups-1, df2=(N*groups)-groups, ncp=ncp, lower.tail=T)
        
        # Step 6 - pValues under null
        p0 <- pf(q=fA, df1=groups-1, df2=(N*groups)-groups,lower.tail=F)
        
        # Step 7 - transform p-values back into state space [0;1]
        p0 <- (p0-alpha)/(1-alpha)
        
        # Step 8 - computing fisher statistic
        chiF <- -2*sum(log(1-p0))
        chiP[i] <- pchisq(q=chiF, df=2*length(p0), lower.tail=T)
      }
      # Step 9 - computing the proportion significant
      pow <- sum(chiP < alphaF)/n.iter
      print(paste(N, P[k], ES[es], pow,sep=" "))
    } else{
      # For computing efficiency all subsequent effect sizes
      # are rounded to one if previous was .995 or higher
      pow <- 1
    }
#     if(pow == 0){pow <- 1}
    powe <- cbind(powe,
                  pow)
  }
  power <- rbind(power, powe)
}
filepath <- paste0('1pAlpha',
                   '.csv')
write.csv2(t(power),
           file=filepath)




# 1-p
power <- NULL
for(k in 1:length(P)){
  pow <- 0
  powe <- NULL
  for(es in 1:length(ES)){
    if(pow < .995){
      for(i in 1:n.iter){
        # Step 1 - Critical value
        fCV <- qf(p=alpha, df1=groups-1, df2=N*groups-groups,lower.tail=F)
        
        # Step 2 - non-centrality parameter
        f2 <- ES[es]^2/(1-ES[es]^2)
        ncp <- f2*(N*groups)
        
        # Step 3 - beta
        beta <- pf(q=fCV, df1=groups-1, df2=(N*groups)-groups, ncp=ncp, lower.tail=T)
        
        # Step 4 - pValues under alternative distribution
        pA <- runif(n=P[k],min=0,max=beta)
        
        # Step 5 - tValues simulated under alternative
        fA <- qf(p=pA, df1=groups-1, df2=(N*groups)-groups, ncp=ncp, lower.tail=T)
        
        # Step 6 - pValues under null
        p0 <- pf(q=fA, df1=groups-1, df2=(N*groups)-groups,lower.tail=F)
        
        # Step 7 - transform p-values back into state space [0;1]
        p0 <- (1-p0)/(1-alpha)
        
        # Step 8 - computing fisher statistic
        chiF <- -2*sum(log(p0))
        chiP[i] <- pchisq(q=chiF, df=2*length(p0), lower.tail=T)
      }
      # Step 9 - computing the proportion significant
      pow <- sum(chiP < alphaF)/n.iter
      print(paste(N, P[k], ES[es], pow,sep=" "))
    } else{
      # For computing efficiency all subsequent effect sizes
      # are rounded to one if previous was .995 or higher
      pow <- 1
    }
    #     if(pow == 0){pow <- 1}
    powe <- cbind(powe,
                  pow)
  }
  power <- rbind(power, powe)
}
filepath <- paste0('1p',
                   '.csv')
write.csv2(t(power),
           file=filepath)

# P
power <- NULL
for(k in 1:length(P)){
  pow <- 0
  powe <- NULL
  for(es in 1:length(ES)){
    if(pow < .995){
      for(i in 1:n.iter){
        # Step 1 - Critical value
        fCV <- qf(p=alpha, df1=groups-1, df2=N*groups-groups,lower.tail=F)
        
        # Step 2 - non-centrality parameter
        f2 <- ES[es]^2/(1-ES[es]^2)
        ncp <- f2*(N*groups)
        
        # Step 3 - beta
        beta <- pf(q=fCV, df1=groups-1, df2=(N*groups)-groups, ncp=ncp, lower.tail=T)
        
        # Step 4 - pValues under alternative distribution
        pA <- runif(n=P[k],min=0,max=beta)
        
        # Step 5 - tValues simulated under alternative
        fA <- qf(p=pA, df1=groups-1, df2=(N*groups)-groups, ncp=ncp, lower.tail=T)
        
        # Step 6 - pValues under null
        p0 <- pf(q=fA, df1=groups-1, df2=(N*groups)-groups,lower.tail=F)
        
        # Step 7 - transform p-values back into state space [0;1]
        p0 <- (p0)/(1-alpha)
        
        # Step 8 - computing fisher statistic
        chiF <- -2*sum(log(p0))
        chiP[i] <- pchisq(q=chiF, df=2*length(p0), lower.tail=T)
      }
      # Step 9 - computing the proportion significant
      pow <- sum(chiP < alphaF)/n.iter
      print(paste(N, P[k], ES[es], pow,sep=" "))
    } else{
      # For computing efficiency all subsequent effect sizes
      # are rounded to one if previous was .995 or higher
      pow <- 1
    }
    #     if(pow == 0){pow <- 1}
    powe <- cbind(powe,
                  pow)
  }
  power <- rbind(power, powe)
}
filepath <- paste0('p',
                   '.csv')
write.csv2(t(power),
           file=filepath)
