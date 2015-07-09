setwd("D:/files/phd/toogoodtobefalse/journals2014")

folders <- list.files()[-7]
library(devtools)
# install_github('statcheckTEMP', 'chartgerink')

library(statcheck)
# 
# final <- NULL
# for(folder in folders){
#   temp <- checkHTMLdir(paste0(folder, '/'))
#   temp$journal <- folder
#   final <- rbind(final, temp)
# }
# 
# write.csv2(final, '../pilotgender.csv')

## Sampling 2 per journal unique papers
setwd("D:/files/phd/toogoodtobefalse/")
final <- read.csv2('pilotgender.csv')
# Remove decision errors
final <- final[final$DecisionError == FALSE & final$gender == TRUE,]

tests <- NULL
# Per journal
journals <- unique(final$journal)
for(journal in journals){
  uniques <- unique(final$Source[final$journal == journal])
  if(length(unique(final$Source[final$journal == journal]) == 1)){
    set.seed(123)
    papers <- sample(uniques, size=1)
  } else {papers <- sample(uniques, size=2)}
  for(paper in papers){
    # Randomly select one row
    set.seed(12)
    i <- sample(c(1, dim(final[final$Source == paper,])[1]), size=1)
    tests <- rbind(tests, final[final$Source == paper,][i,])
  }
}

write.csv2(tests, 'pilot.csv')
