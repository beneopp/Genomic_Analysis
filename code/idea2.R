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
  vectors <- res$vectors[,1:2]
  datCompress <- meanCenteredMat %*% vectors
}

PlotPCAFig <- function (dat, types, title) {

  datCompress <- PerformPCA(dat)

  pca_df <- data.frame(PC1 = datCompress[,1], PC2 = datCompress[,2])
  #pca_data_perc <- round(100*pca$sdev^2/sum(pca$sdev^2),1)

  g <- ggplot(pca_df, aes(datCompress[,1], datCompress[,2], color = factor(types))) +
    geom_point(aes(color = types), size = 2) +
    #labs(x=paste0("PC1 (",pca_data_perc[1],")"),
         #y=paste0("PC2 (",pca_data_perc[2],")"), color = "Sample Type",
         #title = title)
    labs(x = "PC1", y = "PC2", title = title)

  return(list(g, datCompress))
}

DivideSamples <- function () {
  hasPree <- as.logical(mRNA_ExprSet@phenoData@data$pefull)

  mRNA_ExprPree <- mRNA_Expr[hasPree,]
  miRNA_ExprPree <- miRNA_Expr[hasPree,]

  mRNA_ExprNorm <- mRNA_Expr[!hasPree,]
  miRNA_ExprNorm <- miRNA_Expr[!hasPree,]

  expressions <- list("Pree" = list("mRNA_Expr" = mRNA_ExprPree, "miRNA_Expr" = miRNA_ExprPree),
                      "Norm" = list("mRNA_Expr" = mRNA_ExprNorm, "miRNA_Expr" = miRNA_ExprNorm))
}

allGeneNames <- mRNA_ExprSet@featureData@data$SYMBOL

mRNA_Expr <- t(exprs(mRNA_ExprSet))
colnames(mRNA_Expr) <- allGeneNames

JCI_genes <- JCI_geneList$SYMBOL
mRNA_Expr <- mRNA_Expr[, colnames(mRNA_Expr) %in% JCI_genes]

miRNA_Expr <- miRNA_Expr[, rownames(DEmiRNAs)]

expressions <- DivideSamples()

means <- list()

for (exprType in names(expressions)) {
  avg <- colMeans(expressions[[exprType]]$mRNA_Expr)
  means[[exprType]] <- avg
}

ggplot() +
  geom_point(aes(x = means$Pree, y = means$Norm)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  labs(x = "Means in Preeclampsia", y = "Means in Normal", title = "mRNA Means")

medians <- list()

for (exprType in names(expressions)) {
  avg <- apply(expressions[[exprType]]$mRNA_Expr, 2, median)
  medians[[exprType]] <- avg
}

ggplot() +
  geom_point(aes(x = means$Pree, y = means$Norm)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  labs(x = "Medians in Preeclampsia", y = "Medians in Normal", title = "mRNA Medians")


means <- list()

for (exprType in names(expressions)) {
  avg <- colMeans(expressions[[exprType]]$miRNA_Expr)
  means[[exprType]] <- avg
}

ggplot() +
  geom_point(aes(x = means$Pree, y = means$Norm)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  labs(x = "Means in Preeclampsia", y = "Means in Normal", title = "miRNA Means")

medians <- list()

for (exprType in names(expressions)) {
  avg <- apply(expressions[[exprType]]$miRNA_Expr, 2, median)
  medians[[exprType]] <- avg
}

ggplot() +
  geom_point(aes(x = means$Pree, y = means$Norm)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  labs(x = "Medians in Preeclampsia", y = "Medians in Normal", title = "miRNA Medians")

vars <- list()

for (exprType in names(expressions)) {
  avg <- apply(expressions[[exprType]]$mRNA_Expr, 2, var)
  vars[[exprType]] <- avg
}

ggplot() +
  geom_point(aes(x = vars$Pree, y = vars$Norm)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  labs(x = "Medians in Preeclampsia", y = "Medians in Normal", title = "mRNA Variance")
