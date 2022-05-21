require(miRNAtap)
require(miRNAtap.db)
require(biomaRt)

# get human gene database
print("retrieving human gene database")
geneDatabase <- useMart("ENSEMBL_MART_ENSEMBL")
geneDatabase <- useDataset("hsapiens_gene_ensembl", geneDatabase)

print("creating data frame")
setwd("~/Documents/Source_Data_M1_expression")
load("RData/DEmiRNAs.rda")
load("RData/JCI_geneList.rda")

Predict_miRNATargets <- function (miRNA) {

  targetPredictions <- getPredictedTargets(miRNA,
                                           sources=c("pictar","diana","targetscan","mirdb"),
                                           species="hsa", method="geom", min_src=2)

  if ( is.null(targetPredictions) ) {
    return()
  }

  targetAnnots <- getBM(mart = geneDatabase, attributes=c("entrezgene_id", "hgnc_symbol"),
                        filter="entrezgene_id", values=rownames(targetPredictions),
                        uniqueRows=TRUE)

  targetAnnots <- targetAnnots[targetAnnots[,2] != "", ]

  # rename rows of targetPredictions using gene name
  rownames(targetPredictions) <- targetAnnots[match(rownames(targetPredictions),targetAnnots[,1]), 2]

  # remove NAs and empty strings
  selector <- (!is.na(rownames(targetPredictions))) & (rownames(targetPredictions) != "")

  targetPredictions <- targetPredictions[selector, ]

  if ( is.null(targetPredictions) ) {
    return()
  }

  targetPredictionsInJCI <- rownames(targetPredictions)[rownames(targetPredictions) %in% JCI_genes]
}

JCI_genes <- JCI_geneList$SYMBOL

binaryTemplate <- replicate(length(JCI_genes), 0)
stringTemplate <- replicate(length(JCI_genes), "")

target_df <- data.frame(row.names = JCI_genes, mapped = binaryTemplate, LCC = binaryTemplate,
                        "Target of miRNA" = binaryTemplate, "Name of miRNA" = stringTemplate, check.names = FALSE)

lccGenes <- read.csv("geneLists/LCC_new.csv", sep = "\n", header = FALSE)$V1

target_df$LCC[rownames(target_df) %in% lccGenes] <- 1

print("finding targets")

for (miRNA in rownames(DEmiRNAs)) {
  targets <- Predict_miRNATargets(miRNA)
  if (is.null(targets)) {
    next
  }

  for (target in targets) {
    target_df[target, "Target of miRNA"] <- 1
    if (target_df[target, "Name of miRNA"] == "") {
      target_df[target, "Name of miRNA"] <- miRNA
    } else {
      target_df[target, "Name of miRNA"] <- paste0(target_df[target, "Name of miRNA"], ", ", miRNA)
    }
  }

}

mappedGenes <- read.csv("geneLists/mapped_genes.csv", header = FALSE)$V1

target_df[mappedGenes, "mapped"] <- 1

write.csv(target_df, file = "miRNA_targets/target_df.csv")



