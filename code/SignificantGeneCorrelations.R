require(ggplot2)
require(reshape2)
require(textshape)
require(Cairo)
require(corrplot)
require(Biobase)
require(psych)
require(ComplexHeatmap)

#load required files
setwd("~/Documents/Source_Data_M1_expression")
load("RData/miRNA_Expr.rda")
load("RData/DEmiRNAs.rda")
load("RData/mRNA_ExprSet.rda")
load("RData/JCI_geneList.rda")

GetNameOfGenes <- function(geneNameList) {
  # create unique list of gene names

  if (!is.character(geneNameList)) {
    geneNameList <- as.character(geneNameList)
  } else if (any(startsWith(geneNameList, "hsa-"))) {
    stop("geneNameList in GetNameOfGenes was provided with a miRNA list")
  }

  newNameList <- NULL
  numberToAssign <- NULL

  for (transcript in geneNameList) {

    if (transcript %in% names(numberToAssign)) {
      numberToAssign[transcript] <- numberToAssign[[transcript]] + 1
    } else {
      numberToAssign[[transcript]] <- 1
    }

    newName <- paste(transcript, numberToAssign[transcript], sep="_")
    newNameList <- append(newNameList, newName)
  }

  newNameList
}

SelectmRNAGenes <- function (mat, mRNAList) {
  # select mRNA genes only in mRNAList
  if (!is.character(mRNAList)) {
    stop("mRNAList in SelectmRNAGenes is not a character vector.")
  } else if (any(startsWith(mRNAList, "hsa-"))) {
    stop("mRNAList in SelectmRNAGenes was provided with a miRNA list")
  }

  mat <- mat[, colnames(mat) %in% mRNAList]
}

SelectmiRNAGenes <- function (miRNAList) {
  # select miRNA genes only in miRNAList
  miRNA_Expr <- miRNA_Expr[, miRNAList]
}

DivideSamples <- function () {
  # divide sample types by cases and controls
  hasPree <- as.logical(mRNA_ExprSet@phenoData@data$pefull)

  mRNA_ExprPree <- mRNA_Expr[hasPree,]
  miRNA_ExprPree <- miRNA_Expr[hasPree,]

  mRNA_ExprNorm <- mRNA_Expr[!hasPree,]
  miRNA_ExprNorm <- miRNA_Expr[!hasPree,]

  expressions <- list("Pree" = list("mRNA_Expr" = mRNA_ExprPree, "miRNA_Expr" = miRNA_ExprPree),
                      "Norm" = list("mRNA_Expr" = mRNA_ExprNorm, "miRNA_Expr" = miRNA_ExprNorm))
}

DrawImage <- function (results, exprType, imageType, listType, range) {

  corr <- results$corr
  pVals <- results$pVals

  save(corr, file = paste0("matrices/", listType, "/", listType, exprType, imageType, "Corr.rda"))
  save(pVals, file = paste0("matrices/", listType, "/", listType, exprType, imageType, "_pVals.rda"))

  print(paste(imageType, "genes for", exprType, listType, ":", ncol(corr)))

  # some LCC genes are not in the JCI genes
  corr <- corr[, colnames(corr) %in% JCI_genes]
  pVals <- pVals[, colnames(pVals) %in% JCI_genes]

  # create an annotation bar on the top and side to show whether the mRNA or miRNA gene is up or down-regulated

  anno_color <- c("Up-regulated" = "blue", "Down-regulated" = "yellow", "Unknown" = "green")
  top_anno_bar <- JCI_geneList$Control.Case[match(colnames(corr), JCI_geneList$SYMBOL)]
  top_anno_bar <- ifelse(top_anno_bar == "Up", "Down-regulated", "Up-regulated")
  top_anno <- HeatmapAnnotation(Regulation = top_anno_bar, col = list(Regulation = anno_color),
                                annotation_label = "Regulation (Preeclampsia/Control)")

  side_anno_bar <- DEmiRNAs$Expression[match(rownames(corr), rownames(DEmiRNAs))]
  side_anno_bar <- ifelse(side_anno_bar == "Up", "Down-regulated", "Up-regulated")
  side_anno <- rowAnnotation(Regulation = side_anno_bar, col = list(Regulation = anno_color),
                             show_annotation_name = FALSE, show_legend = FALSE)

  # col_fun is the colors for the heatmap
  col_fun <- circlize::colorRamp2(c(min(corr), 0, max(corr)), c("red", "white", "blue"))

  fileDest <- paste0("heatmaps/", listType, "/", listType, exprType, imageType, ".png")

  heatmapWidth <- as.integer((ncol(corr) * 20)/96) + 1
  imageWidth <- (ncol(corr) * 20) + (96 * 5)
  CairoPNG(fileDest, height = nrow(corr) * 20, width = imageWidth, units = "px")
  ht <- Heatmap(corr, col = col_fun, name = "Pearson Correlation Coefficient", bottom_annotation = top_anno,
                right_annotation = side_anno, width = unit(heatmapWidth, "in"))
  draw(ht)
  dev.off()

  # create another heatmap with p-values
  fileDest <- paste0("heatmaps/", listType, "/", listType, exprType, imageType, "with_pVals.png")

  CairoPNG(fileDest, height = nrow(corr) * 20, width = ncol(corr) * 20, units = "px")
  corrplot(corr, is.corr = FALSE, col.lim = range, method = "color", p.mat = pVals, pch.cex = 1)
  dev.off()
}

DrawIntersection <- function (results, listType) {

  # get intersection of normal significant gene names
  mRNA_preeGenes <- names(results$Pree$sigGenes)
  mRNA_normGenes <- names(results$Norm$sigGenes)
  intersectGenes <- intersect(mRNA_preeGenes, mRNA_normGenes)

  if (length(intersectGenes) == 0) {
    return(intersectGenes)
  }

  write(intersectGenes, file = paste0("geneLists/", listType, "_corrLists/", listType, "IntersectGenes.csv"), sep = "\n")

  dats <- list()
  rang <- c(0, 0)
  for (exprType in names(results)) {
    corr <- results[[exprType]]$corr
    pVals <- results[[exprType]]$pVals
    sigGenes <- results[[exprType]]$sigGenes

    # get unique name of the intersect genes
    uniqueIntersectGenes <- sigGenes[intersectGenes]

    corr <- corr[, uniqueIntersectGenes]
    pVals <- pVals[, uniqueIntersectGenes]

    # rename the matrices with the normal name
    colnames(corr) <- intersectGenes
    colnames(pVals) <- intersectGenes

    dats[[exprType]] <- list(corr = corr, pVals = pVals)

    if (min(corr) < rang[1]) {
      rang[1] <- min(corr)
    }

    if (max(corr) > rang[2]) {
      rang[2] <- max(corr)
    }
  }

  for (exprType in names(dats)) {
    dat <- dats[[exprType]]
    DrawImage(dat, exprType, "Intersect", listType, rang)
  }

  intersectGenes
}

DrawSetDiff <- function (results, listType, intersectGenes) {

  for (exprType in names(results)) {
    corr <- results[[exprType]]$corr
    pVals <- results[[exprType]]$pVals
    sigGenes <- results[[exprType]]$sigGenes

    setDiffGenes <- setdiff(names(sigGenes), intersectGenes)

    if (length(setDiffGenes) == 0) {
      next
    }

    # get unique name of the setDiff genes
    uniqueSetDiffGenes <- sigGenes[setDiffGenes]

    corr <- corr[, uniqueSetDiffGenes]
    pVals <- pVals[, uniqueSetDiffGenes]

    # rename the matrices with the normal name
    colnames(corr) <- setDiffGenes
    colnames(pVals) <- setDiffGenes

    dat <- list("corr" = corr, "pVals" = pVals)
    range <- c(min(corr), max(corr))
    DrawImage(dat, exprType, "SetDiff", listType, range)

    write(setDiffGenes, paste0("geneLists/", listType, "_corrLists/", listType, "SetDiff", exprType, ".csv"), sep = "\n")
  }

}

UniqueToNormalGeneNames <- function(uniqueNames) {
  normalNames <- strsplit(uniqueNames, "_")
  normalNames <- unique(sapply(normalNames, function (x) x[[1]][1]))
}

CalculatePValues <- function (corr, numOfObs) {
  # creates p-value matrix
  tVals <- (corr * sqrt(numOfObs - 2)) / sqrt(1 - corr^2)
  pVals <- 2 * pt(tVals, df = numOfObs - 2, lower.tail = FALSE)

  if (!(identical(colnames(corr), colnames(pVals))) |
  !(identical(rownames(corr), rownames(pVals)))) {
    stop("row or column names do not correspond in
    CalculatePValues function")
  }

  pVals
}

DrawUnion <- function (expressions, listType) {
  # find genes that have significant correlations with microRNA

  allGeneNames <- colnames(expressions$Pree$mRNA_Expr)
  uniqueGeneNames <- GetNameOfGenes(allGeneNames)
  unionGenes <- NULL
  results <- list()

  for (exprType in names(expressions)) {
    mRNA_Expr <- expressions[[exprType]]$mRNA_Expr
    miRNA_Expr <- expressions[[exprType]]$miRNA_Expr

    colnames(mRNA_Expr) <- uniqueGeneNames

    corr <- cor(miRNA_Expr, mRNA_Expr)
    numOfObs <- nrow(mRNA_Expr)
    pVals <- corr.p(corr, numOfObs, adjust = "none")$p

    # significant genes have at least one p-value correlation less than 0.05
    sigGenes <- colnames(pVals)[apply(pVals, 2, function (x) any(x < 0.05))]

    # convert gene names back to original names and save them as names
    normalSigGeneNames <- UniqueToNormalGeneNames(sigGenes)
    names(sigGenes) <- normalSigGeneNames

    # get only unique gene name for each gene
    sigGenes <- sigGenes[unique(normalSigGeneNames)]

    # save correlation matrix, pValue matrix, and the significant genes
    results[[exprType]]$corr <- corr
    results[[exprType]]$pVals <- pVals
    results[[exprType]]$sigGenes <- sigGenes

    unionGenes <- union(unionGenes, normalSigGeneNames)
  }

  write(unionGenes, file = paste0("geneLists/", listType, "_corrLists/", listType, "Union.csv"))

  rang <- c(0, 0)
  dats <- list()
  for (exprType in names(expressions)) {
    # retreive data
    corr <- results[[exprType]]$corr
    pVals <- results[[exprType]]$pVals
    sigGenes <- results[[exprType]]$sigGenes

    # sometimes gene1_1 has a significant correlation, but gene1_2 does not
    # consequently, sigGenesCorr gets the probe with significant correlations
    sigGenesCorr <- corr[, sigGenes]
    sigGenesPVals <- pVals[, sigGenes]

    # change the sigGenesCorr and sigGenesPVals to normal gene namess
    colnames(sigGenesCorr) <- names(sigGenes)
    colnames(sigGenesPVals) <- names(sigGenes)

    # otherGenes are genes without significant correlation in the expression type
    otherGenes <- setdiff(unionGenes, names(sigGenes))

    # rename the gene names to normal names so otherGenesCorr and otherGenesPVals
    # can be indexed by otherGenes
    colnames(corr) <- allGeneNames
    colnames(pVals) <- allGeneNames

    otherGenesCorr <- corr[, otherGenes]
    otherGenesPVals <- pVals[, otherGenes]

    # combine the matrices
    corr <- cbind(sigGenesCorr, otherGenesCorr)
    pVals <- cbind(sigGenesPVals, otherGenesPVals)

    dats[[exprType]] <- list(corr = corr, pVals = pVals)

    # find the minimum and maximum correlation in both expression types
    if (min(corr) < rang[1]) {
      rang[1] <- min(corr)
    }

    if (max(corr) > rang[2]) {
      rang[2] <- max(corr)
    }
  }

  # save the image for each expression Type
  for (exprType in names(dats)) {
    DrawImage(dats[[exprType]], exprType, "Union", listType, rang)
  }

  results
}

DrawImages <- function (results, listType) {
  # draws intersection and setDiff images
  intersectGenes <- DrawIntersection(results, listType)
  DrawSetDiff(results, listType, intersectGenes)
}

allGeneNames <- as.character(mRNA_ExprSet@featureData@data$SYMBOL)

lccGenes <- read.csv("geneLists/LCC_new.csv", sep = "\n", header = FALSE)$V1

JCI_genes <- JCI_geneList$SYMBOL

geneLists <- list("LCC" = lccGenes, "JCI" = JCI_genes)

miRNA_Expr <- SelectmiRNAGenes()

for (listType in names(geneLists)) {
  mRNA_Expr <- t(exprs(mRNA_ExprSet))
  colnames(mRNA_Expr) <- allGeneNames

  print("Selecting DE genes")
  geneList <- geneLists[[listType]]
  mRNA_Expr <- SelectmRNAGenes(mRNA_Expr, geneList)

  miRNA_Expr <- miRNA_Expr[, rownames(DEmiRNAs)]

  print("dividing samples")
  expressions <- DivideSamples()
  print("Drawing images")
  results <- DrawUnion(expressions, listType)
  DrawImages(results, listType)
  print("\n")
}

