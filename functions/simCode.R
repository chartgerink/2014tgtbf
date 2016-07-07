#################
# Do not adjust #
#################
powe <- NULL
power <- NULL
chiP <- NULL
# Per N
for(n in 1:length(N)){
  power <- NULL
  for(k in 1:length(P)){
    pow <- 0
    powe <- NULL
    for(es in 1:length(ES)){
      if(pow < .995){
        for(i in 1:n.iter){
          # Step 1 - Critical value
          tCV <- qt(p=alpha, df=N[n]-1,lower.tail=F)
          
          # Step 2 - non-centrality parameter
          f2 <- ES[es]^2/(1-ES[es]^2)
          ncp <- f2*N[n]
          
          # Step 3 - beta
          beta <- pt(q=tCV,df=N[n]-1,ncp=ncp, lower.tail=T)
          
          # Step 4 - pValues under alternative distribution
          pA <- runif(n=P[k],min=0,max=beta)
          
          # Step 5 - tValues simulated under alternative
          tA <- qt(p=pA, df=N[n]-1, ncp=ncp, lower.tail=T)
          
          # Step 6 - pValues under null
          p0 <- pt(q=tA, df=N[n]-1,lower.tail=F)
          
          # Step 7 - transform p-values back into state space [0;1]
          p0 <- (p0-alpha)/(1-alpha)
          
          # Step 8 - computing fisher statistic
          chiF <- -2*sum(log(p0))
          chiP[i] <- pchisq(q=chiF, df=2*length(p0), lower.tail=F)
        }
        # Step 9 - computing the proportion significant
        pow <- sum(chiP < alphaF)/n.iter
        print(paste(N[n], P[k], ES[es], pow,sep=" "))
      } else{
        # For computing efficiency all subsequent effect sizes
        # are rounded to one if previous was .995 or higher
        pow <- 1
      }
      if(pow == 0){pow <- 1}
      powe <- cbind(powe,
                    pow)
    }
    power <- rbind(power, powe)
  }
  filepath <- sprintf('%sdata/N_%s.csv', 
                      run,
                      N[n])
  power <- as.data.frame(power)
  names(power) <- ES
  write.csv(cbind(P, power), row.names = FALSE,
            file=filepath)
  ### General comments
  ### Comment 1
  ### Noncentral distribution estimation yields the possibility
  ### of null beta values, which yields a power of null as all
  ### p0 values are 1. This is a limit due to the precision of
  ### numeric estimates in R/computers in general. This becomes a
  ### problem when the previous element in the row (i.e., ES|nr results)
  ### was still smaller than .995. Especially for 1 or 2 results 
  ### this yielded power values of zero. After inspecting several trial
  ### runs of the simulation, I decided to put 0 values (as they only
  ### followed high power values) to 1, to counteract this computational
  ### problem
}

