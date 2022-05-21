library(corrplot)
library(Cairo)
library(ComplexHeatmap)
library(Biobase)

# load data
dir <- "/Users/benoppenheimer/Documents/Source_Data_M1_expression"
load(paste(dir,"Rdata", "DEmiRNAs.rda", sep = "/"))
load(paste(dir, "Rdata", "miRNA_Expr.rda", sep = "/"))
load(paste(dir, "Rdata", "mRNA_ExprSet.rda", sep = "/"))
phenoDat <- phenoData(mRNA_ExprSet)@data

# subset miRNAs by DE miRNAs
DEmiRNA_Expr <- miRNA_Expr[, rownames(DEmiRNAs)]

#divide subjects by cases and controls
hasPree <- as.logical(phenoDat$pefull)
DEmiRNA_ExprPree <- DEmiRNA_Expr[hasPree, ]

#get correlation of cases data
pearsonCorrPree <- cor(DEmiRNA_ExprPree)

anno_col <- ifelse(DEmiRNAs$Expression == "Up", "Down-regulated", "Up-regulated")
anno_color <- c("Up-regulated" = "blue", "Down-regulated" = "yellow", "Unknown" = "green")
right_anno <- rowAnnotation(Regulation = anno_col, col = list(Regulation = anno_color),
                            annotation_label = "Regulation (Preeclampsia/Control)")
col_fun <- circlize::colorRamp2(c(min(pearsonCorrPree), 0, max(pearsonCorrPree)), c("red", "white", "blue"))

fileDest <- paste(dir, "heatmaps",  "miRNA_only", "mirVsMirCasesPearson.png", sep = "/")
title <- "Pearson Correlation of miRNAs for Preeclampsia Patients"

CairoPNG(fileDest, width = 700)
h <- Heatmap(pearsonCorrPree, col = col_fun, column_title = title, name = "Pearson Correlation", right_annotation = right_anno)
draw(h)
dev.off()

testRes <- cor.mtest(DEmiRNA_ExprPree, conf.level = 0.95, method = "pearson")

fileDest <- paste(dir, "heatmaps",  "miRNA_only", "mirVsMirCasesPearsonPvals.png", sep = "/")

CairoPNG(fileDest)
corrplot(pearsonCorrPree, sig.level = 0.05, order = "hclust", method = "color",
         p.mat = testRes$p, pch.cex = 2)
dev.off()

#controls heatmap
DEmiRNA_ExprNorm <- DEmiRNA_Expr[!hasPree, ]
pearsonCorrNorm <- cor(DEmiRNA_ExprNorm)

col_fun <- circlize::colorRamp2(c(min(pearsonCorrNorm), 0, max(pearsonCorrNorm)), c("red", "white", "blue"))

fileDest <- paste(dir, "heatmaps",  "miRNA_only", "mirVsMirControlPearson.png", sep = "/")
title <- "Pearson Correlation of miRNAs for Normal Patients"
CairoPNG(fileDest, width = 700)
h <- Heatmap(pearsonCorrNorm, col = col_fun, column_title = title, name = "Pearson Correlation", right_annotation = right_anno)
draw(h)
dev.off()

testRes <- cor.mtest(DEmiRNA_ExprNorm, conf.level = 0.95, method = "pearson")

fileDest <- paste(dir, "heatmaps", "miRNA_only", "mirVsMirControlPearsonPvals.png", sep = "/")

CairoPNG(fileDest)
corrplot(pearsonCorrNorm, sig.level = 0.05, order = "hclust", method = "color",
         p.mat = testRes$p, pch.cex = 2)
dev.off()


