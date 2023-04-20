require(DESeq2)
require(EnhancedVolcano)
require(dplyr)

rna <- readRDS(file = "data/RDS/rna-norm.rds")
rna.pair <- readRDS(file = "data/RDS/rna-paired.rds")

# Hack add normal samples
v <- rna.pair$barcode
bool <- (rna$barcode %in% v) | (rna$sample_type == "Solid Tissue Normal")
rna <- rna[,bool]
print(dim(rna))
print(table(rna$stage))

rna$pheno = "PT"
rna[,rna$sample_type == "Solid Tissue Normal"]$pheno = "Nt"
rna$pheno = as.factor(rna$pheno)

#! DEseq2
ddsSE <- DESeqDataSet(rna, design = ~ pheno)
dds <- DESeq(ddsSE)
#! Combiantoria de los contrastes
# res <- results(dds1,contrast=c("grupo",x[2],x[1]))
res <- results(dds)
df <- cbind(as.data.frame(res),as.data.frame(rowData(rna)))
up <- df[ (df$log2FoldChange>2.0) & (df$padj<0.05),]
down <- df[ (df$log2FoldChange<(-2.0)) & (df$padj<0.05),]

EnhancedVolcano(res,
                lab=df$gene_name,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 1e-2,
                FCcutoff = 2.0,
                subtitle = "PT vs Nt",
                labSize = 3.0,
                legendPosition = 'right',
                legendLabSize = 12)
ggsave("data/plots/rna-volc-deg.png", width = 3117, height = 2205, units = "px",
       dpi = 300)
dev.off()


# histogram plot
df1 <- df[order(nchar(df$chr), df$chr, df$start_position),]
df2 <- df1 %>% group_by(chr) %>% mutate(id = 1:n())
df2$rowid <- 1:nrow(df2)
df2$chr <- factor(df2$chr, levels = c(1:22,"X","Y"))
df2$lfc = df2$log2FoldChange
cmp <- (df2$log2FoldChange > (-2.0)) & (df2$log2FoldChange < 2.0)
df2[ cmp,]$lfc = 0
library(ggplot2)
ggplot(data = df2, aes(x = rowid, y = lfc)) +
  geom_bar(stat = "identity",color = "#FF6666") + 
  labs(x = "Chromosome", y = "log2FoldChange") +
  facet_grid(~chr, scales = 'free_x', space = 'free_x', switch = 'x') + 
  theme_classic() + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.margin = unit(0, "lines"))
ggsave("data/plots/rna-log2fc-whole.png", width = 3117, height = 2205, units = "px",
       dpi = 300)
dev.off()

saveRDS(df, file = "data/RDS/rna-deg-NtvsPT.rds")


# library(circlize)
# n = 1000
# circos.par("track.height" = 0.1)
# circos.initialize(df$chr, x = df$start_position)
# circos.track(df$chr, y = df$log2FoldChange,
#              panel.fun = function(x, y) {
#                circos.text(CELL_META$xcenter, 
#                            CELL_META$cell.ylim[2] + mm_y(5), 
#                            CELL_META$sector.index)
#                circos.axis(labels.cex = 0.6)
#              })
# col = rep(c("#FF0000", "#00FF00"), 10)
# circos.trackPoints(df$chr, df$start_position, df$log2FoldChange, col = col, pch = 16, cex = 0.5)
# circos.text(-1, 0.5, "text", sector.index = "a", track.index = 1)
