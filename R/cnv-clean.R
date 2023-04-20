library(tidyr)
library(dplyr)

source("R/tools.R")

#! MAIN
cnvs.raw <- readRDS(file="data/RDS/cnvs-raw.rds")

cnvs <- clean.samples(cnvs.raw)

#! CLEAN GENE
bm <- read.table(file = "data/tables/biomart.csv", sep = ",", header = TRUE)
cnvs <- clean.genes(cnvs,bm)

saveRDS(cnvs, file = "data/RDS/cnvs-clean.rds")




