library(crayon)
library(TCGAbiolinks)

cat(red("0.5 Download Data"),"\n")
setwd("data/GDC")
query <- GDCquery(project = "TCGA-KIRC",
                  data.category = "Copy Number Variation",
                  data.type = "Gene Level Copy Number")
dat <- getResults(query)
query[[1]][[1]] <- dat[!duplicated(dat$cases),] # filter query
GDCdownload(query)
# #! gene annotation GRCh38.p13
data.raw <- GDCprepare(query)
setwd("..")
saveRDS(data.raw, file = "RDS/cnvs-raw.rds")