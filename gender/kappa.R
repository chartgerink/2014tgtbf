if(!require(psych)){install.packages('psych')}
library(psych)
# setwd(choose.dir())
setwd('D:/Dropbox/projects/2014tgtbf')

source('functions//FisherMethod.R')

gend <- read.csv("gender/gendercoded cleaned and discussed.csv",
                 header = TRUE, sep = ";", dec = ".")

cohen.kappa(gend[, 6:8])

# 1 = null expected
# 2 = effect expected
# 3 = no expectation
table(gend$significance, gend$final_code)
options(scipen = 5)
for(sig in unique(gend$significance)){
  for(code in unique(gend$final_code[!is.na(gend$final_code)])){
    sel <- gend$significance == sig & gend$final_code == code
  
    x <- (FisherMethod(x = gend$Computed[sel], id = 1, alpha = 0))
    
    cat(sprintf("For %s %s, k = %s, chi2 = %s, p = %s\n", sig, code, x$CountNSig, x$Fish, x$PFish))
  }
}
