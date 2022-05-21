print("getting other data")
setwd("~/Documents/Source_Data_M1_expression")
lccGenes <- read.csv("geneLists/LCC_new.csv", header = FALSE)$V1

load("RData/targetPredictions.rda")

GetSubsetOf_df <- function(df, geneList) {
  selector <- df$TargetGene %in% geneList
}

upSelector <- GetSubsetOf_df(targetsOfUp_miRNAs_df, lccGenes)
lccTargetsOfUp_miRNAs_df <- targetsOfUp_miRNAs_df[upSelector,]
write.csv(lccTargetsOfUp_miRNAs_df, "miRNA_targets/lccTargetsOfUp_miRNAs.csv", quote = FALSE)

downSelector <- GetSubsetOf_df(targetsOfDown_miRNAs_df, lccGenes)
lccTargetsOfDown_miRNAs_df <- targetsOfDown_miRNAs_df[downSelector,]
write.csv(lccTargetsOfDown_miRNAs_df, "miRNA_targets/lccTargetsOfDown_miRNAs.csv", quote = FALSE)

allSelector <- GetSubsetOf_df(all_miRNAs_df, lccGenes)
lccAll_miRNAs_df <- all_miRNAs_df[allSelector,]
write.csv(lccAll_miRNAs_df, "miRNA_targets/lccAll_miRNAs.csv", quote = FALSE)

mirPrintSelector <- rownames(mirPrint) %in% lccGenes
lccMirPrint <- mirPrint[mirPrintSelector,]
write.csv(lccMirPrint, "miRNA_targets/lccMirPrint.csv", quote = FALSE)