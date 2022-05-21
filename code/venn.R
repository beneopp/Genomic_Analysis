require(ggVennDiagram)
require(ggplot2)

AreDisjoint <- function (sets) {
  for (i in 1:(length(sets)-1)) {
    for (j in (i+1):length(sets)) {
      intersection <- intersect(sets[[i]], sets[[j]])
      if (length(intersection) != 0) {
        return(FALSE)
      }
    }
  }
  return(TRUE)
}

MakeVenn <- function (listType) {

  intersectGenes <- read.csv(paste0("geneLists/", listType, "IntersectGenes.csv"), header = FALSE)$V1
  setDiffNormGenes <- read.csv(paste0("geneLists/", listType, "SetDiffNorm.csv"), header = FALSE)$V1
  setDiffPreeGenes <- read.csv(paste0("geneLists/", listType, "SetDiffPree.csv"), header = FALSE)$V1

  if (!AreDisjoint(list(intersectGenes, setDiffPreeGenes, setDiffNormGenes))) {
    stop("Wrong Overlapping of Genes")
  }

  normGenes <- union(intersectGenes, setDiffNormGenes)
  preeGenes <- union(intersectGenes, setDiffPreeGenes)
  unionGenes <- read.csv(paste0("geneLists/", listType, "Union.csv"), header = FALSE)$V1

  if (length(setdiff(unionGenes,union(normGenes, preeGenes))) != 0 |
    length(setdiff(union(normGenes, preeGenes), unionGenes)) != 0) {
    stop("unionGenes and union(normGenes, preeGenes) are not equal")
  }

  g <- ggVennDiagram(list('Controls' = normGenes, 'Cases' = preeGenes))
}

listTypes <- c("JCI", "LCC")

for (listType in listTypes) {
  g <- MakeVenn(listType) +
    labs(title = paste0("Venn Diagram for ", listType, " Genes"))
  ggsave(g, filename = paste0("VennDiagrams/", listType, "Venn.png"))
}

