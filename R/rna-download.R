library(TCGAbiolinks)
library(crayon)

#####################################################
cat(red("0.5 Download Data"),"\n")
setwd("data/GDC")
query <- GDCquery(project = "TCGA-KIRC",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts")
GDCdownload(query)
#! gene annotation GRCh38.p13
rna.raw <- GDCprepare(query)
setwd("../")
saveRDS(rna.raw, file = "RDS/rna-raw.rds")
