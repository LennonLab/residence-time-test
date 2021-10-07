################################################################################
#                                                                              #
#	Microbial Carbon Traits: EcoPlate                                    #
#   This script produces EcoPlate data file from raw data                             #
#                                                                              #
################################################################################
#                                                                              #
#	Written by: Mario Muscarella                                                 #
#                                                                              #
#	Last update: 2015/08/06                                                      #
#                                                                              #
################################################################################

# Setup Work Environment
rm(list=ls())
library("here")
sem <- function(x){sd(na.omit(x))/sqrt(length(na.omit(x)))}

here()

source(here("code", "EcoPlate.R"))

file.path <- here("data", "RTLC", "EcoPlate")

file.names <- list.files(path = file.path, all.files=FALSE,
                         full.names=FALSE, ignore.case=FALSE, include.dirs=FALSE)

file.names <- file.names[grep("*_RT*", file.names)]

hours <- c()
strains <- c()
sets <- c()
for (name in file.names){
  file.name.info <- strsplit(name, "\\_") # split file name
  sp.id <- strsplit(file.name.info[[1]][4], "\\.t")[[1]][1]
  set <- strsplit(file.name.info[[1]][3], "\\.")[[1]][1]
  hour <- strsplit(file.name.info[[1]][1], "\\.")[[1]][1]
  strains <- c(strains, sp.id)
  sets <- c(sets, set)
  hours <- c(hours, hour)
  data <- read.ecoplate(input = paste(file.path,"/", name, sep=""), skip=23)
  assign(sp.id, data)
  }

strains <- unique(strains)
resource.names <- as.matrix(read.table(here("code","resource_matrix.txt")))
mol.groups <- as.matrix(read.delim(here("code", "moleculetype_matrix.txt"), header=F))
resources <- levels(as.factor(resource.names))
r.names <- as.factor(resource.names)[1:32]
c.grouping <- as.factor(mol.groups)[1:32]
group.res <- data.frame(r.names, c.grouping)[-1, ]
resources <- resources[resources != "Water"]

eco.data <- as.data.frame(matrix(NA, nrow = length(strains), ncol = 35))
colnames(eco.data) <- c("Tau", "Set", "Hours", resources, "NumRes")
eco.data$Tau <- strains
eco.data$Set <- sets
eco.data$Hours <- hours

for (i in strains){
  data <- get(i)
  avg.water <- mean(c(data[1,1], data[1,5], data[1,9]))
  sd.water <- sd(c(data[1,1], data[1,5], data[1,9]))
  data.corr <- round(data - (avg.water + 1.96 * sd.water), 2)
  data.corr[data.corr <= 0] <- 0
  for (j in resources){
    data.r <- which(eco.data$Tau  == i)
    data.c <- which(colnames(eco.data) == j)
    eco.data[data.r , data.c] <- round(mean(data.corr[resource.names == j]), 3)
    eco.data$NumRes[data.r] <- sum(as.numeric(eco.data[data.r, 4:34]) > 0)
    eco.data$Avg[data.r] <- mean(as.numeric(eco.data[data.r, 4:34]))
  }
}

write.table(eco.data, here("data", "RTLC", "eco.data_rt_S4_72.txt"), quote=FALSE,
            row.names=FALSE, sep="\t")




