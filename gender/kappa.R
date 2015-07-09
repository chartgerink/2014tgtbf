library(psych)

x <- read.csv("C:/Users/chjh/Dropbox/projects/2014toogoodtobefalse/gender/gendercoded cleaned.csv")

cohen.kappa(x = x[,4:6])
