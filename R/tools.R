require(SummarizedExperiment)


###########################################
clean.samples <- function(se){
  # 1.0 CLEAN SAMPLES
  print(dim(se))
  print(table(se$sample_type))
  cmp <- (se$sample_type == "Solid Tissue Normal") |
    (se$sample_type == "Primary Tumor")
  se.clean <- se[,cmp]
  se.clean <- se.clean[,!is.na(se.clean$sample_type)]
  se.clean <- se.clean[,!is.na(se.clean$ajcc_pathologic_stage)]
  se.clean <- se.clean[,!duplicated(se.clean$sample)] # duplicated samples
  print(dim(se.clean))
  # stopifnot(ncol(se.clean) == 474)
  se.clean$stage <- "Normal"
  cmp <- (se.clean$sample_type == "Primary Tumor") & 
    (se.clean$ajcc_pathologic_stage == "Stage I")
  se.clean[,cmp]$stage <- "StageI"
  cmp <- (se.clean$sample_type == "Primary Tumor") & 
    (se.clean$ajcc_pathologic_stage == "Stage II")
  se.clean[,cmp]$stage <- "StageII"
  cmp <- (se.clean$sample_type == "Primary Tumor") & 
    (se.clean$ajcc_pathologic_stage == "Stage III")
  se.clean[,cmp]$stage <- "StageIII"
  cmp <- (se.clean$sample_type == "Primary Tumor") & 
    (se.clean$ajcc_pathologic_stage == "Stage IV")
  se.clean[,cmp]$stage <- "StageIV"
  print(table(se.clean$stage))
  # View(assay(se.clean))
  print(dim(se.clean))
  return(se.clean)
}

##########################################
clean.genes <- function(dat,bm){
  #! CLEAN ROWS WITH ONE NA AT LEAST
  # dat <- rna.clean
  m <- assay(dat)
  stopifnot(sum(duplicated(rownames(m)))==0)
  print(dim(m))
  m <- m[rowSums(is.na(m)) == 0, ]
  print(dim(m))
  dat <- dat[rownames(dat) %in% rownames(m),]
  print(dim(dat))
  #! CLEAN ROWS WITHOUT ANNOTATION
  genes <- as.data.frame(rowData(dat))
  if("gene_id" %in% colnames(genes)) ## case of CNVS
    genes <- genes %>% separate(gene_id, c("ensembl_gene_id", "version"), '\\.')
  # v <- merge( bm, genes, by.x="ensembl_gene_id", by.y="ensmb")
  v <- match(genes$ensembl_gene_id, bm$ensembl_gene_id)
  d1 <- genes[!is.na(v),]
  print(dim(d1))
  # names(d1)[3] <- "gene_name"
  d1$gene_name <- NULL
  d2 <- bm[v[!is.na(v)],]
  annot <- cbind(d1,d2)
  dat <- dat[!is.na(v),]
  rowData(dat) <- annot
  print(dim(dat))
  stopifnot(nrow(dat[duplicated(rownames(dat)),])==0)
  return(dat)
}


pca.stage <- function(dat, fout = NULL){
  d2 <- data.frame(etapa = colData(dat)$stage)
  # d2$etapa <-as.factor(d2$etapa)
  rownames(d2) <- colnames(dat)
  mydat = NOISeq::readData( assay(dat) , factors = d2)
  # myPCA = dat(mydat, type = "PCA", norm=TRUE)
  myPCA = dat(mydat, type = "PCA", logtransf = F)
  if(!is.null(fout)){
    print(paste0("Writing in: ",fout))
    png(fout)
  }
  explo.plot(myPCA, factor = "etapa", plottype = "scores")
  if(!is.null(fout))dev.off()
}