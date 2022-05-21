print("downloading packages")
require(miRNAtap)
require(miRNAtap.db)

require(biomaRt)

setwd("~/Documents/Source_Data_M1_expression")

# get human gene database
print("retrieving human gene database")
geneDatabase <- useMart("ENSEMBL_MART_ENSEMBL")
geneDatabase <- useDataset("hsapiens_gene_ensembl", geneDatabase)

print("downloading data")
# load DE miRNAs

load("RData/DEmiRNAs.rda")

# create data structures to store information on gene targets

allTargetGenes <- NULL

df_template <- data.frame("miR" = character(), "TargetGene" = character(),
                         "RankProduct" = numeric(), "RankFinal" = integer())

targetsOfUp_miRNAs_df <- df_template

targetsOfDown_miRNAs_df <- df_template

all_miRNAs_df <- df_template

print("getting miRNA targets")
# get gene targets for each miRNA

for ( miRNA in row.names(DEmiRNAs) ) {
     
  targetPredictions <- getPredictedTargets(miRNA, 
                                           sources=c("pictar","diana","targetscan","mirdb"), 
                                           species="hsa", method="geom", min_src=2)
  
  if ( is.null(targetPredictions) ) {
    next
  }
  
  print(nrow(targetPredictions))
  
  targetAnnots <- getBM(mart = geneDatabase, attributes=c("entrezgene_id", "hgnc_symbol"), 
                        filter="entrezgene_id", values=rownames(targetPredictions), 
                        uniqueRows=TRUE)
  
  targetAnnots <- targetAnnots[targetAnnots[,2] != "", ] 
  
  # rename rows of targetPredictions using gene name
  rownames(targetPredictions) <- targetAnnots[match(rownames(targetPredictions),targetAnnots[,1])
                                              , 2]
  
  # remove NAs and empty strings
  selector <- (!is.na(rownames(targetPredictions))) & (rownames(targetPredictions) != "")
  
  targetPredictions <- targetPredictions[selector, ]
  
  if ( is.null(targetPredictions) ) {
    next
  }
  
  append_df <- data.frame(miR = rep(miRNA, nrow(targetPredictions)), 
                          TargetGene = rownames(targetPredictions), 
                          RankProduct = targetPredictions[, ncol(targetPredictions)-1], 
                          RankFinal = targetPredictions[, ncol(targetPredictions)])
  
  if ( DEmiRNAs[miRNA, "Expression"] == "Up" ) {
    
    targetsOfDown_miRNAs_df <- rbind(targetsOfDown_miRNAs_df, append_df)
    
  } else if ( DEmiRNAs[miRNA, "Expression"] == "Down" ) {
    
    targetsOfUp_miRNAs_df <- rbind(targetsOfUp_miRNAs_df, append_df)
    
  } else {
    
    stop(paste("Error: the miRNA", miRNA, "is neither up or down-regulated", sep = " "))
    
  }
  
  all_miRNAs_df <- rbind(all_miRNAs_df, append_df)
  
  allTargetGenes <- append(allTargetGenes, rownames(targetPredictions))
  
}

if (nrow(targetsOfDown_miRNAs_df) + nrow(targetsOfUp_miRNAs_df) != nrow(all_miRNAs_df)) {
  stop("Error in compiling data frames. 
       The sum of the rows targetsOfDown_miRNAs_df and targetsOfUp_miRNAs_df does not
       equal the number of rows in all_miRNAs_df.")
}

write.csv(targetsOfUp_miRNAs_df, "csv_files/targetsOfUp_miRNAs.csv", quote = FALSE, row.names = FALSE)

write.csv(targetsOfDown_miRNAs_df, "csv_files/targetsOfDown_miRNAs.csv", quote = FALSE, row.names = FALSE)

write.csv(all_miRNAs_df, "csv_files/all_miRNAs.csv", quote = FALSE, row.names = FALSE)

load("RData/JCI_geneList.rda")

targetsInJCIList <- allTargetGenes[allTargetGenes %in% JCI_geneList$SYMBOL]

print(paste("number of genes:", length(targetsInJCIList)))

targetsInJCIList <- unique(targetsInJCIList)

allTargetGenes <- unique(allTargetGenes)

mirPrint <- matrix(0, length(allTargetGenes), nrow(DEmiRNAs),
                   dimnames = list(allTargetGenes, rownames(DEmiRNAs)))

for (miRNA in rownames(DEmiRNAs)) {
  
  targetPredictions <- c()
  
  if (miRNA %in% all_miRNAs_df$miR) {
    targetPredictions <- all_miRNAs_df$TargetGene[all_miRNAs_df$miR == miRNA]
  } else {
    next
  }
  
  mirPrint[targetPredictions, miRNA] <- 1
  
}

write.csv(mirPrint, "csv_files/mirPrint.csv", quote = FALSE)

save(targetsOfUp_miRNAs_df, targetsOfDown_miRNAs_df, all_miRNAs_df, 
     mirPrint, file = "RData/targetPredictions.rda", row.names)
