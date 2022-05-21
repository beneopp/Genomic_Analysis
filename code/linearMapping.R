library(Biobase)
library(ggplot2)

setwd("~/Documents/Source_Data_M1_expression")
load("RData/miRNA_Expr.rda")
load("RData/DEmiRNAs.rda")
load("RData/mRNA_ExprSet.rda")
load("RData/JCI_geneList.rda")

DivideSamples <- function () {
  hasPree <- as.logical(mRNA_ExprSet@phenoData@data$pefull)

  mRNA_ExprPree <- mRNA_Expr[hasPree,]

  mRNA_ExprNorm <- mRNA_Expr[!hasPree,]

  expressions <- list("Pree" = mRNA_ExprPree,
                      "Norm" = mRNA_ExprNorm)
}

allGeneNames <- mRNA_ExprSet@featureData@data$SYMBOL

mRNA_Expr <- t(exprs(mRNA_ExprSet))
colnames(mRNA_Expr) <- allGeneNames

hasPree <- as.logical(mRNA_ExprSet@phenoData@data$pefull)
hasPree <- factor(ifelse(hasPree, "preclampsia", "normal"))

mRNA <- "NEIL3"
miRNA <- "hsa-miR-144-5p"

mRNA_Expr <- mRNA_Expr[,mRNA]
miRNA_Expr <- miRNA_Expr[,miRNA]

ggplot(mapping = aes(mRNA_Expr, miRNA_Expr, color = hasPree)) +
  geom_point() +
  labs(xlab = "Expression of CAMP", ylab = "Expression of hsa-miR-31-5p")
