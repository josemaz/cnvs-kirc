library(estimate)
library(reshape2)
library(ggplot2)

rna <- readRDS(file="RDS/rna-paired.rds")
df <- as.data.frame(colData(rna))
df$muestra = substr(df$barcode,1,15)


# scores <- read.table( header = TRUE,
#   url("https://bioinformatics.mdanderson.org/estimate/tables/kidney_renal_clear_cell_carcinoma_RNAseqV2.txt")
# )
# write.table(scores, file = "tables/estimate-kirc.tsv", sep = "\t")
scores.kirc <- read.table(file = "tables/estimate-kirc.tsv", sep = "\t")
names(scores.kirc) <- c("ID","Stromal","Immune","ESTIMATE")
v <- match(df$muestra,scores.kirc$ID)
v <- v[!is.na(v)]
scores.kirc <- scores.kirc[v,]

summary(scores.kirc)

# scores <- read.table( header = TRUE,
#   url("https://bioinformatics.mdanderson.org/estimate/tables/uterine_corpus_endometrial_carcinoma_RNAseqV2.txt")
# )
# write.table(scores, file = "tables/estimate-ucec.tsv", sep = "\t")
scores.ucec <- read.table(file = "tables/estimate-ucec.tsv", sep = "\t")
names(scores.ucec) <- c("ID","Stromal","Immune","ESTIMATE")


df1 <- melt(scores.kirc, value.name = "Score", variable.name = "CellType")
df1$cancer <- "kirc"
df2 <- melt(scores.ucec, value.name = "Score", variable.name = "CellType")
df2$cancer <- "ucec"
df <-  rbind(df1,df2)

ggplot(data = df, aes(x = CellType, y = Score, fill=cancer)) +
  geom_boxplot() 

