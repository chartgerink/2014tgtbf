setwd('D:/Dropbox/projects/bump')

# # Select out all unneccessary columns and rows
# dat <- read.csv2('datafilegender500.csv')
# dat <- dat[dat$gender == TRUE & dat$DecisionError == FALSE, c(2:11, 18:29)]
# # Add a random ordering
# set.seed(12345)
# dat$id <- sample(1:dim(dat)[1], dim(dat)[1], replace=FALSE)
# # Actually order the file by that id
# dat <- dat[with(dat, order(dat$id)),]
# # Remove the lines that have NA's
# dat <- dat[!(is.na(dat$Reported.Comparison) & is.na(dat$Reported.P.Value)),]
# # IMPORTANT HERE
# # Make the nonsignificant readable for in the program, otherwise it will lead to errors in counting
# # Arbitrary nonsignificant p-value selected.
# dat$Reported.P.Value[dat$Reported.Comparison == "ns"] <- .1
# # save the changes
# write.csv2(dat, 'datafilegender500_pre.csv', row.names=FALSE)

# Read in new datafile, including truncated columns and random id order.
dat <- read.csv2('datafilegender500_pre.csv')
# 
# # Create a start object for later only when nothing has been coded yet.
dat$code1 <- 'empty'
names(dat)
resume = 1
sigCode = 0
nsigCode = 0

# Or use this if already were coding
x <- read.csv2('chjhinterim.csv')
resume = x$resume+1
sigCode = x$sigCode
nsigCode = x$nsigCode
dat$code1 <- read.csv2('code1.csv')$x
for(i in 1:10){
  print(i)
}

cat("\014")
rm(i)
for(i in 1:max(dat$id)){
  print(i)
  if(dat$code1[i] == 'empty'){
    sub <- dat[i,]
    cat(sprintf("Is the original result\n\n %s \n\n about gender?\n\n", sub$Raw))
    input <- readline(prompt=cat(paste0(sub$sentences500)))
    if(input == 1){
      dat$code1[i] <- input
      dat$code1[dat$Source == sub$Source & !(dat$id == sub$id)] <- 'selected out, other paper'
      # Save the amount of coded sig and nsig
      sigCode <- sigCode + ifelse(sub$Reported.P.Value <= .05, 1, 0)
      nsigCode <- nsigCode + ifelse(sub$Reported.P.Value > .05, 1, 0)
    }
    if(input == 0){
      dat$code1[i] <- input
      cat("\014")
    } 
  }
  if(dat$code1[i] == 0){cat("\014 Finding next empty entry.")}
  if(sigCode == 90){
    dat$code1[dat$Reported.P.Value <= .05 & !(dat$code1 == 1 | dat$code1 == 0 | dat$code1 == 'selected out, other paper')] <- 'selected out, enough sig'
    sigCode = sigCode + 1
  }
  if(nsigCode == 90){
    dat$code1[dat$Reported.P.Value > .05 & !(dat$code1 == 1 | dat$code1 == 0 | dat$code1 == 'selected out, other paper')] <- 'selected out, enough nsig'
    nsigCode = nsigCode + 1
  }
  if(sigCode > 90 & nsigCode > 90){break}
  cat("\014")
  # Note that 91 is actually 90, but we add one to the number so the if condition does not loop at each iteration.
  cat(sprintf("The amount coded:\n\n %s significant\n\n %s nonsignificant\n\n", sigCode, nsigCode))
  write.csv2(cbind(sigCode, nsigCode, resume), 'chjhinterim.csv', row.names=FALSE)
  write.csv2(dat$code1, 'code1.csv', row.names=FALSE)
}
write.csv2(dat, 'datafilegender500_post.csv', row.names=FALSE)
