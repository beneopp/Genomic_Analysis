require(ggplot2)
require(reshape2)
require(textshape)
require(Cairo)
require(corrplot)
require(Biobase)
require(psych)
require(ComplexHeatmap)

GetNameOfGenes <- function(geneNameList) {

  if (!is.character(geneNameList)) {
    stop("geneNameList in GetNameOfGenes is not a character vector.")
  } else if (any(startsWith(geneNameList, "hsa-"))) {
    stop("geneNameList in GetNameOfGenes was provided with a miRNA list")
  }

  newNameList <- NULL
  numberToAssign <- replicate(length(unique(geneNameList)), 0)
  names(numberToAssign) <- unique(geneNameList)
  for (transcript in geneNameList) {
    numberToAssign[transcript] <- numberToAssign[transcript] + 1
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

  if (any(duplicated(colnames(mat)))) {
    uniqueGeneNames <- GetNameOfGenes(colnames(mat))
    colnames(mat) <- uniqueGeneNames
  }


  mat
}

DivideSamples <- function () {
  hasPree <- as.logical(mRNA_ExprSet@phenoData@data$pefull)

  mRNA_ExprPree <- mRNA_Expr[hasPree,]

  mRNA_ExprNorm <- mRNA_Expr[!hasPree,]

  expressions <- list("Pree" = mRNA_ExprPree,
                      "Norm" = mRNA_ExprNorm)
}

CalculatePValues <- function (corr, numOfObs) {
  tVals <- (abs(corr) * sqrt(numOfObs - 2)) / sqrt(1 - corr^2)
  pVals <- pt(tVals, numOfObs-2, lower.tail = FALSE) * 2
}

GetCorr <- function (exprType) {
  mRNA_Expr <- expressions[[exprType]]

  corr <- cor(mRNA_Expr)
  pVals <- CalculatePValues(corr, nrow(mRNA_Expr))

  comparisons <- (ncol(corr)^2 - ncol(corr))/2
  pValAdj <- 0.05/comparisons

  selector <- pVals < pValAdj & pVals > 0
  mRNAGenes <- colnames(corr)[apply(selector, 2, any)]

  corrStuff <- list("corr" = corr, "pVals" = pVals, "mRNAGenes" = mRNAGenes)
}

SelectSignificantCorrs <- function (expressions) {

  results <- list()

  for (exprType in names(expressions)) {
    corrStuff <- GetCorr(exprType)
    results[[exprType]] <- corrStuff
  }

  if (!identical(names(expressions), names(results))) {
    stop("something is wrong with SelectSignificantCorrs function")
  }

  results
}

CreateGeneAnnotationBar <- function (geneList, exprType) {
  anno_bar <- NULL

  nameToName <- c("Up" = "Up-regulated", "Down" = "Down-regulated")


  for (gene in geneList) {
    if (gene %in% JCI_geneList$SYMBOL) {
      reg <- JCI_geneList$Conrol.Case[JCI_geneList$SYMBOL == gene]
      anno_bar <- c(anno_bar, nameToName[reg])
    } else {
      anno_bar <- c(anno_bar, "Unknown")
    }
  }

  if (any(is.na(anno_bar))) {
    stop("A gene in JCI_geneList was not recognized")
  }

  anno_bar
}

DrawImage <- function (results, exprType, imageType, listType, range) {

  if (exprType %in% names(results)) {
    corr <- results[[exprType]]$corr
    pVals <- results[[exprType]]$pVals
  } else {
    corr <- results$corr
    pVals <- results$pVals
  }

  print(paste(imageType, "genes for", exprType, listType, ":", ncol(corr)))

  anno_color <- c("Up-regulated" = "blue", "Down-regulated" = "yellow", "Unknown" = "green")

  top_anno_bar <- CreateGeneAnnotationBar(colnames(corr), exprType)
  top_anno <- HeatmapAnnotation(Regulation = top_anno_bar, col = list(Regulation = anno_color))

  side_anno_bar <- DEmiRNAs$Expression[match(rownames(corr), rownames(DEmiRNAs))]
  if (exprType == "Pree") {
    side_anno_bar <- ifelse(side_anno_bar == "Up", "Down-regulated", "Up-regulated")
  } else {
    side_anno_bar <- ifelse(side_anno_bar == "Up", "Up-regulated", "Down-regulated")
  }
  side_anno <- rowAnnotation(Regulation = side_anno_bar, col = list(Regulation = anno_color),
                             show_annotation_name = FALSE, show_legend = FALSE)

  col_fun <- circlize::colorRamp2(c(min(corr), 0, max(corr)), c("red", "white", "blue"))

  fileDest <- paste0("heatmaps/OnlyGenes_", listType, exprType, imageType, ".png")

  CairoPNG(fileDest, height = nrow(corr) * 20, width = ncol(corr) * 20, units = "px")
  ht <- Heatmap(corr, col = col_fun, name = "Pearson Correlation Coefficient", bottom_annotation = top_anno,
                right_annotation = side_anno)
  draw(ht)
  dev.off()

  fileDest <- paste0("heatmaps/OnlyGenes_", listType, exprType, imageType, "with_pVals.png")

  CairoPNG(fileDest, height = nrow(corr) * 20, width = ncol(corr) * 20, units = "px")
  corrplot(corr, is.corr = FALSE, col.lim = range, method = "color", p.mat = pVals, pch.cex = 1)
  dev.off()
}

MakeUniqueMatrix <- function(oldGeneList) {

  indices <- NULL

  for (i in seq(length(oldGeneList))) {
    newName <- strsplit(oldGeneList[i], "_")[[1]][1]
    if (!(newName %in% names(indices))) {
      indices[newName] <- i
    }
  }

  indices
}

DrawIntersection <- function (results, listType) {

  # get intersecting set
  mRNApreeGenes <- results$Pree$mRNAGenes
  mRNAnormGenes <- results$Norm$mRNAGenes
  intersectGenes <- intersect(mRNApreeGenes, mRNAnormGenes)

  # If no intersecting set, stop function
  if (length(intersectGenes) == 0) {
    return(NULL)
  }

  geneNames <- sapply(intersectGenes, function (x) strsplit(x, "_")[[1]][[1]])
    if (!identical(names(geneNames), intersectGenes)) {
    stop("did not count intersect genes correctly")
  }

  write(unique(geneNames), sep = "/n", file = paste("csv_files/", listType, "IntersectGenes.csv", sep = ""))

  # make correlation and p-value matrices with intersecting gene set
  intersectResults <- list()
  for (exprType in names(results)) {
    corr <- results[[exprType]]$corr
    corr <- corr[intersectGenes, intersectGenes]

    pVals <- results[[exprType]]$pVals
    pVals <- pVals[rownames(corr), colnames(corr)]

    newIndices <- MakeUniqueMatrix(colnames(corr))

    corr <- corr[newIndices, newIndices]
    pVals <- pVals[newIndices, newIndices]

    colnames(corr) <- names(newIndices)
    rownames(corr) <- names(newIndices)

    colnames(pVals) <- names(newIndices)
    rownames(pVals) <- names(newIndices)

    corr <- cluster_matrix(corr)
    pVals <- pVals[rownames(corr), colnames(corr)]

    intersectResults[[exprType]] <- list("corr" = corr, "pVals" = pVals)
  }

  range <- c(min(intersectResults$Pree$corr, intersectResults$Norm$corr),
             max(intersectResults$Pree$corr, intersectResults$Norm$corr))

  for (exprType in names(results)) {
    DrawImage(intersectResults, exprType, "Intersect", listType, range)
  }


  intersectGenes
}

DrawSetDiff <- function (results, listType, intersectGenes) {

  for (exprType in names(results)) {
    corr <- results[[exprType]]$corr

    sigGenes <- results[[exprType]]$mRNAGenes
    sigGenes <- setdiff(sigGenes, intersectGenes)

    if (length(sigGenes) == 0) {
      next
    }

    corr <- corr[sigGenes,sigGenes]

    pVals <- results[[exprType]]$pVals
    pVals <- pVals[sigGenes, sigGenes]

    newIndices <- MakeUniqueMatrix(colnames(corr))

    corr <- corr[newIndices, newIndices]
    pVals <- pVals[newIndices, newIndices]

    colnames(corr) <- names(newIndices)
    rownames(corr) <- names(newIndices)

    colnames(pVals) <- names(newIndices)
    rownames(pVals) <- names(newIndices)

    corr <- cluster_matrix(corr)
    pVals <- pVals[rownames(corr), colnames(corr)]

    dat <- list("corr" = corr, "pVals" = pVals)

    range <- c(min(corr), max(corr))

    DrawImage(dat, exprType, "SetDiff", listType, range)

    write(names(newIndices), paste0("csv_files/", listType, "SetDiff", exprType, ".csv"))
  }

}

ChangeNames <- function (indices, mat) {
  mat <- mat[indices, indices]
  colnames(mat) <- names(indices)
  rownames(mat) <- names(indices)
  mat
}

DrawUnion <- function (expressions, listType) {

  unionGenes <- NULL
  results <- list()
  rang <- c(0, 0)
  for (exprType in names(expressions)) {
    mRNA_Expr <- expressions[[exprType]]

    corr <- cor(mRNA_Expr)
    pVals <- CalculatePValues(corr, nrow(mRNA_Expr))

    newIndices <- MakeUniqueMatrix(colnames(corr))

    corr <- ChangeNames(newIndices, corr)
    pVals <- ChangeNames(newIndices, pVals)

    results[[exprType]]$corr <- corr
    results[[exprType]]$p <- pVals

    sigGenes <- colnames(pVals)[apply(pVals, 2, function (x) any(x < 0.05 & x > 0))]
    unionGenes <- union(unionGenes, sigGenes)

    if (min(corr) < rang[1]) {
      rang[1] <- min(corr)
    }

    if (max(corr) > rang[2]) {
      rang[2] <- max(corr)
    }
  }

  unionGenesNames <- names(MakeUniqueMatrix(unionGenes))
  write(sort(unionGenesNames), file = paste0("csv_files/", listType, "Union.csv"))

  for (exprType in names(results)) {
    corr <- results[[exprType]]$corr
    pVals <- results[[exprType]]$p
    dat <- list("corr" = corr, "pVals" = pVals)
    DrawImage(dat, exprType, "Union", listType, rang)
  }
}

DrawImages <- function (results, listType) {
  intersectResults <- DrawIntersection(results, listType)
  DrawSetDiff(results, listType, intersectResults)
}


setwd("~/Documents/Source_Data_M1_expression")
load("RData/miRNA_Expr.rda")
load("RData/DEmiRNAs.rda")
load("RData/mRNA_ExprSet.rda")
load("RData/JCI_geneList.rda")

allGeneNames <- mRNA_ExprSet@featureData@data$SYMBOL

lccGenes <- read.csv("geneLists/LCC_new.csv", sep = "\n")

JCI_genes <- JCI_geneList$SYMBOL
geneLists <- list("LCC" = lccGenes, "JCI" = JCI_genes)

for (listType in names(geneLists)) {
  mRNA_Expr <- t(exprs(mRNA_ExprSet))
  colnames(mRNA_Expr) <- allGeneNames

  print("Selecting DE genes")
  geneList <- geneLists[[listType]]
  mRNA_Expr <- SelectmRNAGenes(mRNA_Expr, geneList)

  print("dividing samples")
  expressions <- DivideSamples()
  DrawUnion(expressions, listType)
  print("Selecting significant correlations")
  results <- SelectSignificantCorrs(expressions)
  print("Drawing images")
  DrawImages(results, listType)
  print("\n")
}

