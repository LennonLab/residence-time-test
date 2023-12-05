#! /bin/bash

library("spaa")

OTUs <- read.csv("/N/u/emamuell/Carbonate/GitHub/residence-time-test/data/OTUs_80.csv", header = TRUE)
rownames(OTUs) <- OTUs$X
OTUs <- OTUs[,-1]

niche.overlap <- c()

print(paste("initial: ", mean(niche.overlap(as.data.frame(OTUs), method = "pianka"))))

x <- 1
while(x < 100){
  for(col in 1:ncol(OTUs)){
    OTUs[,col] <- sample(OTUs[,col])
  }
  print(paste(x, " : ", mean(niche.overlap(as.data.frame(OTUs), method = "pianka"))))
  niche.overlap <- c(niche.overlap, mean(niche.overlap(as.data.frame(OTUs), method = "pianka")))
  x = x + 1
}

write.csv(as.data.frame(niche.overlap), "/N/u/emamuell/Carbonate/GitHub/residence-time-test/data/niche_overlap.csv", row.names = FALSE)
