require(NOISeq)
require(EDASeq)
require(crayon)
require(tidyr)

source("R/tools.R")

#! MAIN
rna.raw <- readRDS(file="data/RDS/rna-raw.rds")

rna.clean <- clean.samples(rna.raw)

bm <- read.table(file = "data/tables/biomart.csv", sep = ",", header = TRUE)
rna.clean <- clean.genes(rna.clean,bm)

saveRDS(rna.clean, file = "data/RDS/rna-clean.rds")

#! NORMALIZARTION
dout <- "data/plots/norm"
dir.create(dout, recursive = TRUE)
rna.norm <- rna.clean
pca.stage(rna.clean,fout = paste0(dout,"/PCA-rna-BeforeNorm.png"))
fac <- data.frame(tumor_stage=rna.clean$stage, 
                  row.names=colnames(rna.clean))
ln.data <- withinLaneNormalization(assay(rna.clean), 
                                   rowData(rna.clean)$geneLength, 
                                   which = "full")
gcn.data <- withinLaneNormalization(ln.data , rowData(rna.clean)$gcContent,
                                    which = "full")
norm.counts <- tmm(gcn.data, long = 1000, lc = 0, k = 0)
noiseqData <- NOISeq::readData( norm.counts, factors = fac)
mydata2corr1 = NOISeq::ARSyNseq(noiseqData, norm = "n",  logtransf = FALSE)
assay(rna.norm) <- round(exprs(mydata2corr1))
pca.stage(rna.norm,fout = paste0(dout,"/PCA-rna-AfterNorm.png"))

saveRDS(rna.norm, file = "data/RDS/rna-norm.rds")

# View(assay(rna.norm))
# View(assay(rna.raw))
