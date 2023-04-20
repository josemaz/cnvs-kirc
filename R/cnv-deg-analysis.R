#! ANALYSIS
library(ggplot2)
library("RColorBrewer")

cnvs <- readRDS(file="data/RDS/cnv-paired.rds")
deg <- readRDS(file="data/RDS/rna-deg-NtvsPT.rds")

#! MATCHING
rd.cnvs <- as.data.frame(rowData(cnvs))
v <- match(rd.cnvs$ensembl_gene_id,deg$ensembl_gene_id)
deg <- deg[v,]

#! SAMPLE PHENOTYPES
canc <- cnvs[,cnvs$sample_type=="Primary Tumor"]
# norm <- cnvs[,cnvs$sample_type=="Solid Tissue Normal"]
dim(canc)
# dim(norm)

#! CNVS COUNT
m.canc <- assay(canc)
# m.norm <- assay(norm)
rd <- as.data.frame(rowData(cnvs))
# df<-data.frame(c.amps = rowSums(m.canc > 2), c.stabs = rowSums(m.canc == 2),
#                c.dels = rowSums(m.canc < 2), n.amps = rowSums(m.norm > 2),
#                n.stabs = rowSums(m.norm == 2), n.dels = rowSums(m.norm < 2))
df<-data.frame(c.amps = rowSums(m.canc > 2), c.stabs = rowSums(m.canc == 2),
               c.dels = rowSums(m.canc < 2))
df <- cbind(rd,df)
df$lfc <- deg$log2FoldChange
df$pval <- deg$padj

#! SORT AND NUMBERED
df1 <- df[order(nchar(df$chr), df$chr, df$start_position),]
# df2 <- df1 %>% group_by(chr) %>% mutate(id = 1:n())
df1$rowid <- 1:nrow(df1)
df1$chr <- factor(df1$chr, levels = c(1:22,"X","Y"))

#! EXPRESSION
cmp <- (df1$lfc > (-2.0)) & (df1$lfc < 2.0)
sum(!cmp)
df1[ cmp,]$lfc = 0
cmp = df1$c.amps>floor(ncol(canc)*0.60) #Hack
df2 <- df1[cmp,]

ggplot(data = df2, aes(x = rowid, y = lfc)) +
  geom_bar(stat = "identity",color = "#FF6666") + 
  labs(x = "Chromosome", y = "log2FoldChange") +
  facet_grid(~chr, scales = 'free_x', space = 'free_x', switch = 'x') + 
  theme(strip.text.x = element_text(size = rel(5)),
        text = element_text(size=rel(5)),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.margin = unit(0, "lines"))

# View(df2[df2$lfc != 0.0,])
df3 <- df2[df2$lfc != 0.0,]
df3 <- df3[df3$lfc > 2.0,]
# df3 <- df3[df3$lfc < (-2.0),]
paste(as.character(df3$gene_name),collapse="  ",sep="")


ggplot(data = df2, aes(x = rowid, y = lfc, fill = c.dels)) +
  geom_bar(stat = "identity") + 
  scale_fill_gradient2(low = "red", mid="black", high = "blue") + 
  labs(x = "Chromosome", y = "log2FoldChange") +
  facet_grid(~chr, scales = 'free_x', space = 'free_x', switch = 'x') + 
  # theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.margin = unit(0, "lines"))







#! QUANTIFICATION
# df1 <- df[df$c.amps>floor(ncol(canc)*0.60),] 
# df1[order(-df1$c.amps),][1:100,c(13,16)]
# dim(df1[order(-df1$c.amps),])
# 
# df1 <- df[df$c.dels>floor(ncol(canc)*0.60),] 
# df1[order(-df1$c.dels),][1:100,c(13,18)]
# dim(df1[order(-df1$c.dels),])
# 
# df2 <- df[df$n.amps>floor(ncol(norm)*0.60),] 
# df2[order(-df2$n.amps),][1:100,c(13,19)]
# dim(df2[order(-df2$n.amps),])
# # dim(merge(df1,df2,by="gene_name"))
# 
# df2 <- df[df$n.dels>floor(ncol(norm)*0.60),] 
# df2[order(-df2$n.dels),][1:100,c(13,21)]
# dim(df2[order(-df2$n.dels),])


cmp <- (df1$lfc > (-2.0)) & (df1$lfc < 2.0)
sum(!cmp)
df1[ cmp,]$lfc = 0
cmp = df1$c.dels>floor(ncol(canc)*0.60) #Hack
df2 <- df1[cmp,]













# m1 <- assay(norm)
# hist(rowMedians(m1))
# m2 <- assay(canc)
# hist(rowMedians(m2))
# 
# 
# df<-data.frame(max = rowMax(m1))
# rownames(df) <- rownames(m1)
# head(df %>% arrange(desc(max)))
# df<-data.frame(max = rowMax(m2))
# rownames(df) <- rownames(m2)
# head(df %>% arrange(desc(max)))
# View(bm)
# 
# df<-data.frame(sum = rowSums(m1))
# rownames(df) <- rownames(m1)
# head(df %>% arrange(desc(sum)))
# # View(df)
# hist(m1[rownames(m1)=="ENSG00000145934.14",])
# 
# df<-data.frame(sum = rowSums(m2))
# rownames(df) <- rownames(m2)
# head(df %>% arrange(desc(sum)))
# # View(df)
# hist(m2[rownames(m2)=="ENSG00000038274.15",])
# hist(df$sum)





#! DRAFTS
# ggplot(data = df2, aes(x = rowid, y = lfc, fill = c.amps)) +
#   geom_bar(stat = "identity") + 
#   scale_fill_gradient2(low = "red", mid="black", high = "blue") + 
#   labs(x = "Chromosome", y = "log2FoldChange") +
#   facet_grid(~chr, scales = 'free_x', space = 'free_x', switch = 'x') + 
#   # theme_classic() +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         panel.margin = unit(0, "lines"))