require(igraph)

setwd("~/Documents/Source_Data_M1_expression")

# color nodes of interest to a particular color
ChangeLabel <- function (geneList, color) {
  targetedInList <- intersect(targetedGenes, geneList)
  if (length(targetedInList) == 0) {
    return(colorTable)
  }
  if (!all(targetedInList %in% colorTable$V1)) {
    stop(paste0("not all targetedInList were found"))
  }
  colorTable$V2[colorTable$V1 %in% targetedInList] <- color
  colorTable
}

GetTargetGenes <- function () {

  theoretTargetedGenes <- target_df$X[as.logical(target_df$LCC) & as.logical(target_df$Target.of.miRNA)]
  unionGenes <- read.csv("geneLists/LCC_corrLists/LCCUnion.csv", header = FALSE)$V1
  theoretTargetedGenes <- intersect(theoretTargetedGenes, unionGenes)

  load("matrices/LCC/LCCPreeUnion_pVals.rda")
  preeUnion_pVals <- pVals

  load("matrices/LCC/LCCNormUnion_pVals.rda")
  normUnion_pVals <- pVals
  rm(pVals)

  empiricalTargetedGenes <- NULL
  for (gene in theoretTargetedGenes) {
    targets <- target_df$Name.of.miRNA[target_df$X == gene]
    targets <- unlist(strsplit(targets, ", "))
    pVal1 <- min(preeUnion_pVals[targets, gene])
    pVal2 <- min(normUnion_pVals[targets, gene])

    if (pVal1 < 0.05 & pVal2 < 0.05) {
      empiricalTargetedGenes[gene] <- "intersect"
    } else if (pVal1 < 0.05) {
      empiricalTargetedGenes[gene] <- "pree"
    } else if (pVal2 < 0.05) {
      empiricalTargetedGenes[gene] <- "norm"
    }
  }

  empiricalTargetedGenes
}

#get color table
cmd <- "python package_code/change_STRING_colors.py -s graphs/network_image.svg"
system(cmd)

file.rename("color_table.tsv", "graphs/color_table.tsv")

colorTable <- read.delim("graphs/color_table.tsv", header = FALSE)

# get sets of intere
intersectGenes <- read.csv("geneLists/LCC_corrLists/LCCIntersectGenes.csv", header = FALSE)$V1
preeSetDiffGenes <- read.csv("geneLists/LCC_corrLists/LCCSetDiffPree.csv", header = FALSE)$V1
normSetDiffGenes <- read.csv("geneLists/LCC_corrLists/LCCSetDiffNorm.csv", header = FALSE)$V1

# get miRNA target data
target_df <- read.csv("miRNA_targets/target_df.csv")

targetedGenes <- GetTargetGenes()

# color all nodes black
numOfNodes <- nrow(colorTable)
colorTable$V2 <- replicate(numOfNodes, "rgb(0,0,0)")

red <- "rgb(255,0,0)"
green <- "rgb(0,255,0)"
blue <- "rgb(0,0,255)"

#colorCodes <- list(list(intersectGenes, red), list(preeSetDiffGenes, green), list(normSetDiffGenes, blue))

#for (colorCode in colorCodes) {
  #colorTable <- ChangeLabel(colorCode[[1]], colorCode[[2]])
#}

#colorCodes <- c("intersect" = red, "pree" = green, "norm" = blue)

#for (type in names(colorCodes)) {
  #genesOfInterest <- names(targetedGenes[targetedGenes == type])
  #colorTable$V2[colorTable$V1 %in% genesOfInterest] <- colorCodes[type]
#}

theoretTargetedGenes <- target_df$X[as.logical(target_df$LCC) & as.logical(target_df$Target.of.miRNA)]
write(theoretTargetedGenes, file = "geneLists/theoretTargetedGenes.csv")

write.table(colorTable, "graphs/color_table.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

cmd <- "python package_code/change_STRING_colors.py -s graphs/network_image.svg -c graphs/color_table.tsv"

system(cmd)
file.rename("graphs/network_image.new_colors.svg", "../graphs/marksTarget.svg")

load("RData/JCI_geneList.rda")

yellow <- "rgb(255,255,0)"
colorCodes <- c("Up" = blue, "Down" = yellow)

JCI_geneList$Control.Case <- ifelse(JCI_geneList$Control.Case == "Up", "Down", "Up")

#genesChanged <- NULL
#for (type in names(colorCodes)) {
  #selectedGeneNames <- JCI_geneList$SYMBOL[JCI_geneList$SYMBOL %in% theoretTargetedGenes & JCI_geneList$Control.Case == type]
  #genesChanged <- append(genesChanged, selectedGeneNames)
  #if (!all(selectedGeneNames %in% colorTable$V1)) {
    #stop("some genes were not found in colorTable")
  #}
  #colorTable$V2[colorTable$V1 %in% selectedGeneNames] <- colorCodes[type]
#}

lccOfLcc <- read.csv("geneLists/LCCOfLCC.csv", header = FALSE)$V1
colorTable$V2[colorTable$V1 %in% lccOfLcc] <- orange

if (!identical(sort(genesChanged), sort(theoretTargetedGenes))) {
  stop("some genes were not changed")
}

#orange <- "rgb(255,0,255)"
#colorCodes <- c("Up" = orange, "Down" = green)
#for (type in names(colorCodes)) {
  #lccGenes <- target_df$X[target_df$LCC == 1]
  #selectedGeneNames <- JCI_geneList$SYMBOL[JCI_geneList$SYMBOL %in% lccGenes &
                                             #!(JCI_geneList$SYMBOL %in% theoretTargetedGenes)
                                             #& JCI_geneList$Control.Case == type]
  #selectedGeneNames[selectedGeneNames == "KIR2DS5"] <- "KIR2DL3"
  #if (!all(selectedGeneNames %in% colorTable$V1)) {
    #stop("some genes were not found in colorTable")
  #}
  #colorTable$V2[colorTable$V1 %in% selectedGeneNames] <- colorCodes[type]
#}

write.table(colorTable, "graphs/color_table.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

file.copy("graphs/network_image.svg", "graphs/network_image2.svg")
cmd <- "python package_code/change_STRING_colors.py -s graphs/network_image2.svg -c graphs/color_table.tsv"

system(cmd)
file.rename("../graphs/onlyTarget.svg", "../graphs/targetAndNontarget'.svg")