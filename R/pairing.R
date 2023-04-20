
cnvs <- readRDS(file="data/RDS/cnvs-clean.rds")
rna <- readRDS(file="data/RDS/rna-norm.rds")

# pairing samples
samples <- match(cnvs$sample,rna$sample)
samples <- samples[!is.na(samples)]
stopifnot(sum(duplicated(samples))==0)
length(samples)
rna.paired <- rna[,samples]
samples <- match(rna.paired$sample,cnvs$sample)
samples <- samples[!is.na(samples)]
stopifnot(sum(duplicated(samples))==0)
length(samples)
cnvs.paired <- cnvs[,samples]
dim(cnvs.paired)
# df <- data.frame(rna=rna.paired$sample,cnv=cnvs.paired$sample)
# saveRDS(cnvs.paired, file="data/RDS/cnv-paired.rds")

# pairing genes
rna.rd <- as.data.frame(rowData(rna))
rna.genes <- rna.rd$ensembl_gene_id
cnv.rd <- as.data.frame(rowData(cnvs))
cnv.genes <- cnv.rd$ensembl_gene_id
genes <- rna.genes %in% cnv.genes
rna.paired <- rna.paired[genes,]
#! all genes of cnv are in rna
stopifnot(sum(!(cnv.genes %in% rna.genes))==0)

saveRDS(cnvs.paired, file = "data/RDS/cnvs-paired.rds")
saveRDS(rna.paired, file = "data/RDS/rna-paired.rds")



