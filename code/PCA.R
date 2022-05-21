require(Biobase)
require(ggplot2)


setwd("~/Documents/Source_Data_M1_expression")
load("RData/miRNA_Expr.rda")
load("RData/DEmiRNAs.rda")
load("RData/mRNA_ExprSet.rda")
load("RData/JCI_geneList.rda")

PerformPCA <- function (mat) {
  means <- colMeans(mat)
  meanCenteredMat <- apply(mat, 1, function (x) x - means)
  meanCenteredMat <- t(meanCenteredMat)
  covs <- cov(meanCenteredMat)
  res <- eigen(covs)
  vectors <- res$vectors[,1:3]
  datCompress <- meanCenteredMat %*% vectors
}

PlotPCAFig <- function (dat, types, title) {

  datCompress <- PerformPCA(dat)

  pca_df <- data.frame(PC1 = datCompress[,1], PC2 = datCompress[,2])
  #pca_data_perc <- round(100*pca$sdev^2/sum(pca$sdev^2),1)

  g <- ggplot(pca_df, aes(datCompress[,1], datCompress[,2], color = factor(types))) +
    geom_point(aes(color = types), size = 2) +
    labs(x=paste0("PC1 (",pca_data_perc[1],")"),
         y=paste0("PC2 (",pca_data_perc[2],")"), color = "Sample Type",
         title = title)
    labs(x = "PC1", y = "PC2", title = title)

  return(list(g, datCompress))
}

allGeneNames <- mRNA_ExprSet@featureData@data$SYMBOL

mRNA_Expr <- t(exprs(mRNA_ExprSet))
colnames(mRNA_Expr) <- allGeneNames

JCI_genes <- JCI_geneList$SYMBOL
mRNA_Expr <- mRNA_Expr[, colnames(mRNA_Expr) %in% JCI_genes]

miRNA_Expr <- miRNA_Expr[, rownames(DEmiRNAs)]

hasPree <- as.logical(mRNA_ExprSet@phenoData@data$pefull)
hasPree <- ifelse(hasPree, "preeclampsia", "normal")

vitDLevel <- mRNA_ExprSet@phenoData@data$basevitdng
vitDLevel <- ifelse(vitDLevel > 30, "high", "low")

caseTypes <- NULL
for (i in 1:157) {
  typ <- paste(hasPree[i], vitDLevel[i], sep = "_")
  caseTypes <- c(caseTypes, typ)
}
caseTypes <- factor(caseTypes)

#genes1 <- read.table("miRNA_targets/JCISetDiffNorm.csv", sep = "\n")$V1
#genes2 <- read.table("miRNA_targets/JCISetDiffPree.csv", sep = "\n")$V1

#genes <- setdiff(union(genes1, genes2), intersect(genes1, genes2))

#mRNA_Expr <- mRNA_Expr[, colnames(mRNA_Expr) %in% genes]

mat <- cbind(mRNA_Expr, miRNA_Expr)
res <- PlotPCAFig(mat, caseTypes, "Samples Plot")

res[[1]]

