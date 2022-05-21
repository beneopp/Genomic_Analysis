#################################
#Perform micro-RNA target analysis and generate 'mirPrint'
#################################

require(miRNAtap)
require(miRNAtap.db)

require(biomaRt)
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)

#gsrep=read.csv("data/348 repgenes_updatedpeexclusions_EntrezID.csv")
#head(gsrep)

load("/Users/rehom/Desktop/MiRNA_PE_04132018_newanalysis_VDAART/Kevin_miRNA/To Kevin/DEmiRNAS.rda")

repgenes=read.csv("/Users/rehom/Desktop/MiRNA_PE_04132018_newanalysis_VDAART/Kevin_miRNA/To Kevin/JCI_rep_genes.csv")

load("/Users/rehom/Desktop/MiRNA_PE_04132018_newanalysis_VDAART/Kevin_miRNA/To Kevin/mapped_LCC.rda")

repgenes_LCC=merge(mapped_LCC, repgenes, by="SYMBOL", all.x=TRUE)

save(repgenes_LCC, file="repgenes.LCC.rda")

diffexprs=DEmiRNAs

#Filter out genes that pass FDR
#sigGeneListUp <- rownames(subset(diffexprs, Fold.change.Patient...Control>=2 & Benjamini.Hochberg.FDR<=0.05))
#topMatrixUp <- norm[which(rownames(norm) %in% sigGeneListUp),]

sigGeneListUp <- subset(diffexprs, Exxpression=="Down")$mirpecandidates
topMatrixUp <- subset(diffexprs, Exxpression=="Down")


#sigGeneListDown <- rownames(subset(diffexprs, Fold.change.Patient...Control<=-2 & Benjamini.Hochberg.FDR<=0.05))
#topMatrixDown <- norm[which(rownames(norm) %in% sigGeneListDown),]

sigGeneListDown <- subset(diffexprs, Exxpression=="Up" | Exxpression=="UP")$mirpecandidates
topMatrixDown <- subset(diffexprs, Exxpression=="Up" | Exxpression=="UP")


#sigGeneList <- rownames(subset(diffexprs, abs(Fold.change.Patient...Control)>=2 & Benjamini.Hochberg.FDR<=0.05))
#topMatrix <- norm[which(rownames(norm) %in% sigGeneList),]

sigGeneList <- diffexprs$mirpecandidates
topMatrix <- diffexprs

rownames(topMatrix)<-diffexprs$mirpecandidates

#rownames(topMatrix)
#c("hsa-let-7d-5p","hsa-let-7f-5p","hsa-miR-103a-3p","hsa-miR-107","hsa-miR-1260a","hsa-miR-126-3p","hsa-miR-127-3p","hsa-miR-140-3p","hsa-miR-142-3p","hsa-miR-144-5p","hsa-miR-148a-3p","hsa-miR-151a-3p","hsa-miR-151a-5p","hsa-miR-181a-5p","hsa-miR-191-5p","hsa-miR-193a-5p","hsa-miR-194-5p","hsa-miR-1972","hsa-miR-199a-3p","hsa-miR-199a-5p","hsa-miR-19a-3p","hsa-miR-210-3p","hsa-miR-215-5p","hsa-miR-21-5p","hsa-miR-219a-5p","hsa-miR-221-3p","hsa-miR-222-3p","hsa-miR-22-3p","hsa-miR-22-5p","hsa-miR-26a-5p","hsa-miR-26b-3p","hsa-miR-26b-5p","hsa-miR-28-3p","hsa-miR-28-5p","hsa-miR-29a-3p","hsa-miR-29c-3p","hsa-miR-301a-3p","hsa-miR-30b-5p","hsa-miR-30c-5p","hsa-miR-320b","hsa-miR-328-3p","hsa-miR-331-3p","hsa-miR-339-5p","hsa-miR-34a-5p","hsa-miR-374a-5p","hsa-miR-374b-5p","hsa-miR-375","hsa-miR-382-5p","hsa-miR-409-3p","hsa-miR-423-5p","hsa-miR-424-5p","hsa-miR-485-3p","hsa-miR-501-3p","hsa-miR-502-3p","hsa-miR-543","hsa-miR-629-5p","hsa-miR-660-5p","hsa-miR-744-5p","hsa-miR-885-5p","hsa-miR-92a-3p","hsa-miR-99a-5p","mmu-miR-378a-3p")

#Up-regulated
targetGenesMaster <- c()
write.table("miR\tTargetGene\tRankProduct\tRankFinal", "MirTargetGenes.tsv", col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE, append=FALSE)

write.table("miR\tTargetGene\tRankProduct\tRankFinal", "MirTargetGenes.Up.tsv", col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE, append=FALSE)

write.table("miR\tTargetGene\tRankProduct\tRankFinal", "MirTargetGenes.Down.tsv", col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE, append=FALSE)

rownames(topMatrixUp)<-topMatrixUp$mirpecandidates

for (i in 1:length(rownames(topMatrixUp)))
{
  searchTerm <- rownames(topMatrixUp)[i]
  
  #Difficult mirs
  #if (searchTerm=="miR-1260a") searchTerm="miR-1260"
  #if (searchTerm=="miR-210-3p") searchTerm="miR-210"
  #if (searchTerm=="miR-215-5p") searchTerm="miR-215"
  #if (searchTerm=="miR-629-5p") searchTerm="miR-629"
  #if (searchTerm=="miR-219a-5p") searchTerm="miR-219"
  #if (searchTerm=="miR-328-3p") searchTerm="mir-328"
  
  mirTargetPredictions <- getPredictedTargets(searchTerm, sources=c("pictar","diana","targetscan","mirdb"), species="hsa", method="geom", min_src=2)
  
  #Change entrez gene target names to HUGO
  #require(biomaRt)
  #mart <- useMart("ENSEMBL_MART_ENSEMBL")
  #mart <- useDataset("hsapiens_gene_ensembl", mart)
  annots <- getBM(mart=mart, attributes=c("entrezgene", "hgnc_symbol"), filter="entrezgene", values=rownames(mirTargetPredictions), uniqueRows=TRUE)
  rownames(mirTargetPredictions) <- annots[match(rownames(mirTargetPredictions), annots[,1]),2]
  
  #Remove NAs
  mirTargetPredictions <- mirTargetPredictions[!is.na(rownames(mirTargetPredictions)),]
  
  wObject <- data.frame(rep(rownames(topMatrixUp)[i], nrow(mirTargetPredictions)), rownames(mirTargetPredictions), mirTargetPredictions[,(ncol(mirTargetPredictions)-1):ncol(mirTargetPredictions)], row.names=NULL)
  colnames(wObject) <- c("miR","TargetGene","RankProduct","RankFinal")
  write.table(wObject, "MirTargetGenes.Up.tsv", col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE, append=TRUE)
  
  targetGenesMaster <- c(targetGenesMaster, rownames(mirTargetPredictions))
}


#"miranda" removed not working

head(mirTargetPredictions)


#Down_regulated
rownames(topMatrixDown)<-topMatrixDown$mirpecandidates

for (i in 1:length(rownames(topMatrixDown)))
{
  searchTerm <- rownames(topMatrixDown)[i]
  
  #Difficult mirs
  #if (searchTerm=="miR-1260a") searchTerm="miR-1260"
  #if (searchTerm=="miR-210-3p") searchTerm="miR-210"
  #if (searchTerm=="miR-215-5p") searchTerm="miR-215"
  #if (searchTerm=="miR-629-5p") searchTerm="miR-629"
  #if (searchTerm=="miR-219a-5p") searchTerm="miR-219"
  #if (searchTerm=="miR-328-3p") searchTerm="mir-328"
  
  mirTargetPredictions <- getPredictedTargets(searchTerm, sources=c("pictar","diana","targetscan","mirdb"), species="hsa", method="geom", min_src=2)
  
  #Change entrez gene target names to HUGO
  #require(biomaRt)
  #mart <- useMart("ENSEMBL_MART_ENSEMBL")
  #mart <- useDataset("hsapiens_gene_ensembl", mart)
  annots <- getBM(mart=mart, attributes=c("entrezgene", "hgnc_symbol"), filter="entrezgene", values=rownames(mirTargetPredictions), uniqueRows=TRUE)
  rownames(mirTargetPredictions) <- annots[match(rownames(mirTargetPredictions), annots[,1]),2]
  
  #Remove NAs
  mirTargetPredictions <- mirTargetPredictions[!is.na(rownames(mirTargetPredictions)),]
  
  wObject <- data.frame(rep(rownames(topMatrixDown)[i], nrow(mirTargetPredictions)), rownames(mirTargetPredictions), mirTargetPredictions[,(ncol(mirTargetPredictions)-1):ncol(mirTargetPredictions)], row.names=NULL)
  colnames(wObject) <- c("miR","TargetGene","RankProduct","RankFinal")
  write.table(wObject, "MirTargetGenes.Down.tsv", col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE, append=TRUE)
  
  targetGenesMaster <- c(targetGenesMaster, rownames(mirTargetPredictions))
}


# Error in getBM(mart = mart, attributes = c("entrezgene", "hgnc_symbol"),  :               
#Values argument contains no data.
#In addition: Warning message:
#In getPredictedTargets(searchTerm, sources = c("pictar", "diana",  :
#no targets found for mirna miR-31-5

#removing miR-31-5P from topMatrixDown

topMatrixDown_1=topMatrixDown[!rownames((topMatrixDown))=="hsa-miR-31-5", ]




for (i in 1:length(rownames(topMatrixDown_1)))
{
  searchTerm <- rownames(topMatrixDown_1)[i]
  
  #Difficult mirs
  #if (searchTerm=="miR-1260a") searchTerm="miR-1260"
  #if (searchTerm=="miR-210-3p") searchTerm="miR-210"
  #if (searchTerm=="miR-215-5p") searchTerm="miR-215"
  #if (searchTerm=="miR-629-5p") searchTerm="miR-629"
  #if (searchTerm=="miR-219a-5p") searchTerm="miR-219"
  #if (searchTerm=="miR-328-3p") searchTerm="mir-328"
  
  mirTargetPredictions <- getPredictedTargets(searchTerm, sources=c("pictar","diana","targetscan","mirdb"), species="hsa", method="geom", min_src=2)
  
  #Change entrez gene target names to HUGO
  #require(biomaRt)
  #mart <- useMart("ENSEMBL_MART_ENSEMBL")
  #mart <- useDataset("hsapiens_gene_ensembl", mart)
  annots <- getBM(mart=mart, attributes=c("entrezgene", "hgnc_symbol"), filter="entrezgene", values=rownames(mirTargetPredictions), uniqueRows=TRUE)
  rownames(mirTargetPredictions) <- annots[match(rownames(mirTargetPredictions), annots[,1]),2]
  
  #Remove NAs
  mirTargetPredictions <- mirTargetPredictions[!is.na(rownames(mirTargetPredictions)),]
  
  wObject <- data.frame(rep(rownames(topMatrixDown_1)[i], nrow(mirTargetPredictions)), rownames(mirTargetPredictions), mirTargetPredictions[,(ncol(mirTargetPredictions)-1):ncol(mirTargetPredictions)], row.names=NULL)
  colnames(wObject) <- c("miR","TargetGene","RankProduct","RankFinal")
  write.table(wObject, "MirTargetGenes.Down.tsv", col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE, append=TRUE)
  
  targetGenesMaster <- c(targetGenesMaster, rownames(mirTargetPredictions))
}


#All combined

rownames(topMatrix)=topMatrix$mirpecandidates

topMatrix_1=topMatrix[!rownames(topMatrix)=="hsa-miR-31-5", ]

rownames(topMatrix_1)<-topMatrix_1$mirpecandidates

targetGenesMaster <- c()
write.table("miR\tTargetGene\tRankProduct\tRankFinal", "MirTargetGenes.tsv", col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE, append=FALSE)

for (i in 1:length(rownames(topMatrix_1)))
{
  
  searchTerm <- rownames(topMatrix_1)[i]
  
  #Difficult mirs
  #if (searchTerm=="miR-1260a") searchTerm="miR-1260"
  #if (searchTerm=="miR-210-3p") searchTerm="miR-210"
  #if (searchTerm=="miR-215-5p") searchTerm="miR-215"
  #if (searchTerm=="miR-629-5p") searchTerm="miR-629"
  #if (searchTerm=="miR-219a-5p") searchTerm="miR-219"
  #if (searchTerm=="miR-328-3p") searchTerm="mir-328"
  
  mirTargetPredictions <- getPredictedTargets(searchTerm, sources=c("pictar","diana","targetscan","mirdb"), species="hsa", method="geom", min_src=2)
  
  #Change entrez gene target names to HUGO
  #require(biomaRt)
  #mart <- useMart("ENSEMBL_MART_ENSEMBL")
  #mart <- useDataset("hsapiens_gene_ensembl", mart)
  annots <- getBM(mart=mart, attributes=c("entrezgene", "hgnc_symbol"), filter="entrezgene", values=rownames(mirTargetPredictions), uniqueRows=TRUE)
  rownames(mirTargetPredictions) <- annots[match(rownames(mirTargetPredictions), annots[,1]),2]
  
  #Remove NAs
  mirTargetPredictions <- mirTargetPredictions[!is.na(rownames(mirTargetPredictions)),]
  
  wObject <- data.frame(rep(rownames(topMatrix_1)[i], nrow(mirTargetPredictions)), rownames(mirTargetPredictions), mirTargetPredictions[,(ncol(mirTargetPredictions)-1):ncol(mirTargetPredictions)], row.names=NULL)
  colnames(wObject) <- c("miR","TargetGene","RankProduct","RankFinal")
  write.table(wObject, "MirTargetGenes.tsv", col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE, append=TRUE)
  
  targetGenesMaster <- c(targetGenesMaster, rownames(mirTargetPredictions))
}

#Create a unique gene list
targetGenesMaster <- unique(targetGenesMaster)
targetGenesMaster <- targetGenesMaster[!is.na(targetGenesMaster)]

#Create a data-frame that will be used to genearte the mirPrint
mirPrint <- data.frame(row.names=targetGenesMaster)

for (i in 1:length(rownames(topMatrix_1)))
{
  #searchTerm <- gsub("mmu-", "", gsub("hsa-", "", rownames(topMatrix)[i]))
  searchTerm <- rownames(topMatrix_1)[i]
  
  #Difficult mirs
  #if (searchTerm=="miR-1260a") searchTerm="miR-1260"
  #if (searchTerm=="miR-210-3p") searchTerm="miR-210"
  #if (searchTerm=="miR-215-5p") searchTerm="miR-215"
  #if (searchTerm=="miR-629-5p") searchTerm="miR-629"
  #if (searchTerm=="miR-219a-5p") searchTerm="miR-219"
  #if (searchTerm=="miR-328-3p") searchTerm="mir-328"
  
  mirTargetPredictions <- getPredictedTargets(searchTerm, sources=c("diana","pictar","targetscan","mirdb"), species="hsa", method="geom", min_src=2)
  
  #Change entrez gene target names to HUGO
  #require(biomaRt)
  #mart <- useMart("ENSEMBL_MART_ENSEMBL")
  #mart <- useDataset("hsapiens_gene_ensembl", mart)
  annots <- getBM(mart=mart, attributes=c("entrezgene", "hgnc_symbol"), filter="entrezgene", values=rownames(mirTargetPredictions), uniqueRows=TRUE)
  rownames(mirTargetPredictions) <- annots[match(rownames(mirTargetPredictions), annots[,1]),2]
  
  #Populate the mirPrint dataframe that has '1' for when mir targets a gene, and '0' when not
  for (j in 1:length(rownames(mirPrint)))
  {
    searchTerm <- paste("^", rownames(mirPrint)[j], "$", sep="")
    if (any(grepl(searchTerm, rownames(mirTargetPredictions))))
    {
      mirPrint[j,i] <- 1
    }
    else
    {
      mirPrint[j,i] <- 0
    }
  }
  colnames(mirPrint)[i] <- rownames(topMatrix_1)[i]
}


require(ComplexHeatmap)
require(circlize)
require(cluster)
require(RColorBrewer)


topMatrix_1$pfp[1:4]<-0.001

#Colour bar for -log (base 10) FDR Q value for DE mirs, and fold changes
dfMinusLog10FDR <- data.frame(-log10(topMatrix_1[match(colnames(mirPrint), rownames(topMatrix_1)),"pfp"]))
dfFoldChange <- data.frame(topMatrix_1[match(colnames(mirPrint), rownames(topMatrix_1)),"FC..class1.class2."])
dfGeneAnno <- data.frame(dfMinusLog10FDR, dfFoldChange)
colnames(dfGeneAnno) <- c("Mir significance score", "Regulation")
dfGeneAnno[,2] <- ifelse(dfGeneAnno[,2]<1, "Up-regulated", "Down-regulated")
table(dfGeneAnno[,2])
colours <- list("Regulation"=c("Up-regulated"="red3", "Down-regulated"="forestgreen"))
haGenes <- rowAnnotation(df=dfGeneAnno, col=colours, width=unit(1,"cm"))



mirPrintUp <- mirPrint[,which(dfGeneAnno[,2]=="Up-regulated")]
mirPrintDown <- mirPrint[,which(dfGeneAnno[,2]=="Down-regulated")]


table(rownames(mirPrintUp) %in% mapped_LCC$SYMBOL)

table(rownames(mirPrintDown) %in% mapped_LCC$SYMBOL)

mirPrintUp_LCC<-mirPrintUp[rownames(mirPrintUp) %in% mapped_LCC$SYMBOL, ]

mirPrintDown_LCC<-mirPrintDown[rownames(mirPrintDown) %in% mapped_LCC$SYMBOL, ]



#Remove terms with no overlapping genes
mirPrintUp_LCC_ov <- mirPrintUp_LCC[,apply(mirPrintUp_LCC, 2, mean)!=0]
mirPrintDown_LCC_ov <- mirPrintDown_LCC[,apply(mirPrintDown_LCC, 2, mean)!=0]   #All have overlapps


#mirPrintUp<-mirPrintUp[!rownames(mirPrintUp)=="NA", ]

#mirPrintDown<-mirPrintDown[!rownames(mirPrintDown)=="NA", ]


#get top 50 genes
#mirPrintUp <- mirPrintUp[names(sort(apply(mirPrintUp, 1, sum), decreasing=TRUE)[1:50]),]
#mirPrintDown <- mirPrintDown[names(sort(apply(mirPrintDown, 1, sum), decreasing=TRUE)[1:50]),]

#sort by rank 
mirPrintUp_LCC_ov <- mirPrintUp_LCC_ov[names(sort(apply(mirPrintUp_LCC_ov, 1, sum), decreasing=TRUE)[1:78]),]
mirPrintDown_LCC_ov <- mirPrintDown_LCC_ov[names(sort(apply(mirPrintDown_LCC_ov, 1, sum), decreasing=TRUE)[1:78]),]




mirPrint.tmp <- apply(mirPrintUp_LCC_ov, 2, function(x) gsub(1, "Target", x))
mirPrint.tmp <- apply(mirPrint.tmp, 2, function(x) gsub("0", "", x))
rownames(mirPrint.tmp) <- rownames(mirPrintUp_LCC_ov)
colnames(mirPrint.tmp) <- colnames(mirPrintUp_LCC_ov)
mirPrintUp2 <- mirPrint.tmp

mirPrint.tmp <- apply(mirPrintDown_LCC_ov, 2, function(x) gsub(1, "Target", x))
mirPrint.tmp <- apply(mirPrint.tmp, 2, function(x) gsub("0", "", x))
rownames(mirPrint.tmp) <- rownames(mirPrintDown_LCC_ov)
colnames(mirPrint.tmp) <- colnames(mirPrintDown_LCC_ov)
mirPrintDown2 <- mirPrint.tmp


indup=mirPrintUp2[,1:9]=="Target"
table(mirPrintUp2[,1:9]=="Target")

FALSE  TRUE 
635    67 

inddown=mirPrintDown2[,1:6]=="Target"
table(inddown)
FALSE  TRUE 
425    43 



###

#Define the colours and box-sizes for each value in the oncoprint
alter_fun <- list(
  background=function(x, y, w, h)
  {
    grid.rect(x, y, w-unit(0.5, "mm"), height=unit(0.5, "mm"), gp=gpar(fill="grey95", col="grey95"))
  },
  
  Target=function(x, y, w, h)
  {
    grid.rect(x, y, w-unit(2.5, "mm"), h-unit(2.5, "mm"), gp=gpar(fill="black", col="black"))
  }
)
cols <- c("Target"="black")

#Colour bar for -log (base 10) FDR Q value for DE mirs, and fold changes

#sigGeneListUp <- subset(diffexprs, Exxpression=="Down")$mirpecandidates

topMatrix_1$pfp[1:4]<-0.001

topMatrixUp <- subset(topMatrix_1, Exxpression=="Down")

topMatrixDown <- subset(topMatrix_1, Exxpression=="Up" | Exxpression=="UP")


colnames(mirPrintUp2)==rownames(topMatrixUp)

dfMinusLog10FDR <- data.frame(-log10(topMatrixUp[match(colnames(mirPrintUp2), rownames(topMatrixUp)),"pfp"]))
dfFoldChange <- data.frame(topMatrixUp[match(colnames(mirPrintUp2), rownames(topMatrixUp)),"FC..class1.class2."])
dfGeneAnno <- data.frame(dfMinusLog10FDR, dfFoldChange)
colnames(dfGeneAnno) <- c("Mir significance score", "Regulation")
#dfGeneAnno[,2] <- ifelse(dfGeneAnno[,2]<1, "Up-regulated", "Down-regulated")
dfGeneAnno[,2] <- ifelse(dfGeneAnno[,2]<1, "Up-regulated", "Down-regulated")
table(dfGeneAnno[,2])
colours <- list("Regulation"=c("Up-regulated"="royalblue", "Down-regulated"="yellow"))
haGenesUp <- rowAnnotation(df=dfGeneAnno, col=colours, width=unit(1,"cm"))


colnames(mirPrintDown2)==rownames(topMatrixDown)

dfMinusLog10FDR <- data.frame(-log10(topMatrixDown[match(colnames(mirPrintDown2), rownames(topMatrixDown)),"pfp"]))
dfFoldChange <- data.frame(topMatrixDown[match(colnames(mirPrintDown2), rownames(topMatrixDown)),"FC..class1.class2."])
dfGeneAnno <- data.frame(dfMinusLog10FDR, dfFoldChange)
colnames(dfGeneAnno) <- c("Mir significance score", "Regulation")
dfGeneAnno[,2] <- ifelse(dfGeneAnno[,2]>1, "Down-regulated","Up-regulated")
#dfGeneAnno[,2] <- ifelse(dfGeneAnno[,2]>1, "Up-regulated", "Down-regulated")
table(dfGeneAnno[,2])
#colours <- list("Regulation"=c("Up-regulated"="royalblue", "Down-regulated"="yellow"))
colours <- list("Regulation"=c("Down-regulated"="yellow", "Up-regulated"="royalblue"))
haGenesDown <- rowAnnotation(df=dfGeneAnno, col=colours, width=unit(1,"cm"))


row.names(mirPrintDown2)==row.names(mirPrintUp2)
length(grep("^Target", mirPrintUp2[1:78, ]))
length(grep("^Target", mirPrintDown2[1:78, ]))



pdf("mirPrint5.pdf", width=12, height=10)
mirPrintPlotUp <- oncoPrint(t(mirPrintUp2), get_type=function(x) strsplit(x, ";")[[1]],
                            
                            name="mirPrintPlotUp",
                            
                            alter_fun=alter_fun,
                            
                            col=cols,
                            
                            #Order by most significant mir
                            row_order=rownames(topMatrixUp),
                            #column_order=null,
                            
                            remove_empty_columns=TRUE,
                            
                            row_title="", row_title_side="left", row_title_gp=gpar(fontsize=20, fontface="bold"), show_row_names=TRUE, row_names_gp=gpar(fontsize=10, fontface="bold"), row_names_max_width=unit(3, "cm"),
                            
                            column_title="", column_title_side="top", column_title_gp=gpar(fontsize=20, fontface="bold"), column_title_rot=0, show_column_names=TRUE, column_names_gp=gpar(fontsize=10),
                            
                            pct_gp=gpar(fontsize=12, fontface="bold"),
                            
                            axis_gp=gpar(fontsize=15, fontface="bold"),
                            
                            #bottom_annotation=annCaucasian,
                            
                            heatmap_legend_param=list(title="Mir gene targets", at=c("Target"),
                                                      labels=c("Target"), nrow=1, title_position="topcenter"))

###

mirPrintPlotDown <- oncoPrint(t(mirPrintDown2), get_type=function(x) strsplit(x, ";")[[1]],
                              
                              name="mirPrintPlotUp",
                              
                              alter_fun=alter_fun,
                              
                              col=cols,
                              
                              #Order by most significant mir
                              row_order=rownames(topMatrixDown),
                              #column_order=null,
                              
                              remove_empty_columns=TRUE,
                              
                              row_title="", row_title_side="left", row_title_gp=gpar(fontsize=20, fontface="bold"), show_row_names=TRUE, row_names_gp=gpar(fontsize=10, fontface="bold"), row_names_max_width=unit(3, "cm"),
                              
                              column_title="", column_title_side="top", column_title_gp=gpar(fontsize=20, fontface="bold"), column_title_rot=0, show_column_names=TRUE, column_names_gp=gpar(fontsize=10),
                              
                              pct_gp=gpar(fontsize=12, fontface="bold"),
                              
                              axis_gp=gpar(fontsize=15, fontface="bold"),
                              
                              #bottom_annotation=annCaucasian,
                              
                              heatmap_legend_param=list(title="Mir gene targets", at=c("Target"),
                                                        labels=c("Target"), nrow=1, title_position="topcenter"))

###

pushViewport(viewport(layout=grid.layout(nr=2, nc=1)))

pushViewport(viewport(layout.pos.row=1, layout.pos.col=1))
draw(mirPrintPlotUp + haGenesUp, heatmap_legend_side="top", annotation_legend_side="left", newpage=FALSE)
upViewport()

pushViewport(viewport(layout.pos.row=2, layout.pos.col=1))
draw(mirPrintPlotDown + haGenesDown, heatmap_legend_side="top", annotation_legend_side="left", newpage=FALSE)
upViewport()
dev.off()


######################################
#Build networks of mir-to-gene targets
######################################

require(igraph)
require(plotrix)

df <- read.table("/Users/rehom/Desktop/MiRNA_PE_04132018_newanalysis_VDAART/Kevin_miRNA/To Kevin/MirTargetGenes.Up.tsv", sep="\t", header=TRUE)[,1:2]

ind=df$TargetGene %in% mapped_LCC$SYMBOL

table(ind)

df<-df[ind, ]

dim(df)

range(table(df$TargetGene))

df2<-unique(df)

range(table(df2$TargetGene))


three <- names(table(df2$TargetGene)[table(df2$TargetGene)==3])
two <- names(table(df2$TargetGene)[table(df2$TargetGene)==2])
one <- names(table(df2$TargetGene)[table(df2$TargetGene)==1])
df2 <- df2[which(df2$TargetGene %in% c(three, two, one)),]

#Shuffle the df in order to improve final layout

setseed(123)
iRandUp <- sample(rep(1:nrow(df2)))
df2 <- df2[iRandUp,]

gUp <- graph.edgelist(as.matrix(df2), directed=FALSE)

layoutUp <- data.frame(V(gUp)$name, layout.circle(gUp))

iMir <- grep("hsa-miR", layoutUp[,1])
layoutUp[iMir,2:3] <- layoutUp[iMir,2:3] * 18
V(gUp)[iMir]$shape <- "cir"
V(gUp)[iMir]$color = "firebrick2"
#V(gUp)[iMir]$color = "brown3"

iOne <- which(layoutUp[,1] %in% one)
layoutUp[iOne,2:3] <- layoutUp[iOne,2:3] * 12
V(gUp)[iOne]$shape <- "sphere"
V(gUp)[iOne]$color = "yellow"

iTwo <- which(layoutUp[,1] %in% two)
layoutUp[iTwo,2:3] <- layoutUp[iTwo,2:3] * 6
V(gUp)[iTwo]$shape <- "sphere"
#V(gUp)[iTwo]$color = "royalblue"
V(gUp)[iTwo]$color = "cyan3"

iThree <- which(layoutUp[,1] %in% three)
layoutUp[iThree,2:3] <- layoutUp[iThree,2:3] * 1
V(gUp)[iThree]$shape <- "sphere"
#V(gUp)[iThree]$color = "forestgreen"
V(gUp)[iThree]$color = "darkolivegreen3"

#V(gUp)$name <- gsub("miR-", "miR-\n", V(gUp)$name)

V(gUp)$vertex.frame.color <- "white"

###

df <- read.table("/Users/rehom/Desktop/MiRNA_PE_04132018_newanalysis_VDAART/Kevin_miRNA/To Kevin/MirTargetGenes.Down.tsv", sep="\t", header=TRUE)[,1:2]

ind=df$TargetGene %in% mapped_LCC$SYMBOL

table(ind)

df<-df[ind, ]

dim(df)

range(table(df$TargetGene))

df2<-unique(df)

range(table(df2$TargetGene))


#five <- names(table(df$TargetGene)[table(df$TargetGene)==5])

#three <- names(table(df2$TargetGene)[table(df2$TargetGene)==3])
two <- names(table(df2$TargetGene)[table(df2$TargetGene)==2])
one <- names(table(df2$TargetGene)[table(df2$TargetGene)==1])
df2 <- df2[which(df2$TargetGene %in% c(two, one)),]

#Shuffle the df in order to improve final layout
iRandDown <- sample(rep(1:nrow(df2)))
df2 <- df2[iRandDown,]

gDown <- graph.edgelist(as.matrix(df2), directed=FALSE)

layoutDown <- data.frame(V(gDown)$name, layout.circle(gDown))

iMir <- grep("hsa-miR", layoutDown[,1])
layoutDown[iMir,2:3] <- layoutDown[iMir,2:3] * 18
V(gDown)[iMir]$shape <- "square"
#V(gDown)[iMir]$color = "red1"
V(gDown)[iMir]$color = "firebrick2"


iOne <- which(layoutDown[,1] %in% one)
layoutDown[iOne,2:3] <- layoutDown[iOne,2:3] * 16 #14
V(gDown)[iOne]$shape <- "sphere"
V(gDown)[iOne]$color = "yellow"

iTwo <- which(layoutDown[,1] %in% two)
layoutDown[iTwo,2:3] <- layoutDown[iTwo,2:3] * 8 #10
V(gDown)[iTwo]$shape <- "sphere"
#V(gDown)[iTwo]$color = "cornflowerblue"
V(gDown)[iTwo]$color = "cyan3"
#iThree <- which(layoutDown[,1] %in% three)
#layoutDown[iThree,2:3] <- layoutDown[iThree,2:3] * 6
#V(gDown)[iThree]$shape <- "sphere"
#V(gDown)[iThree]$color = "yellow"

#iFour <- which(layoutUp[,1] %in% four)
#layoutUp[iFour,2:3] <- layoutUp[iFour,2:3] * 1
#V(gUp)[iFour]$shape <- "sphere"
#V(gUp)[iFour]$color = "royalblue"



#V(gDown)$name <- gsub("miR-", "miR-\n", V(gDown)$name)

V(gDown)$vertex.frame.color <- "white"

pdf("MirToGeneNetworkDown22.pdf", width=12, height=12)
par(mfrow=c(1,1))
plot.igraph(gUp, layout=data.matrix(layoutUp[,2:3]), edge.curved=FALSE, vertex.size=10.0, vertex.label.dist=0.1, vertex.label.color="black", asp=FALSE, vertex.label.cex=0.8, vertex.label.font=2, edge.color="grey90", edge.width=0.1, edge.arrow.mode=0, main="Up-regulated miRNAs")

#Add circles to segregate the mirs based on number of genes targetting them
#draw.circle(0.0, 0.0, 0.2, lwd=2, lty=4, border="darkblue")
draw.circle(-0.1,-0.2, 0.2, lwd=2, lty=4, border="darkblue")
text(-0.1, -0.3, "3 miRNAs", cex=0.9, font=2)
draw.circle(-0.1,-0.2, 0.45, lwd=2, lty=4, border="darkblue")
text(-0.1, 0.15, "2 miRNAs", cex=0.9, font=2)
draw.circle(-0.1, -0.2, 0.85, lwd=2, lty=4, border="darkblue")
text(-0.1, 0.5, "1 miRNAs", cex=0.9, font=2)

###

plot.igraph(gDown, layout=data.matrix(layoutDown[,2:3]), edge.curved=FALSE, vertex.size=10.0, vertex.label.dist=0.1, vertex.label.color="black", asp=FALSE, vertex.label.cex=0.8,vertex.label.font=2, edge.color="grey90", edge.width=0.1, edge.arrow.mode=0, main="Down-regulated miRNAs")

#Add circles to segregate the mirs based on number of genes targetting them

#draw.circle(-0.1,-0.2, 0.2, lwd=2, lty=4, border="darkblue")
#text(-0.1, -0.3, "3 miRNAs", cex=0.9, font=2)
draw.circle(0.0,-0.18, 0.45, lwd=2, lty=4, border="darkblue")
text(0.0, 0.15, "2 miRNAs", cex=0.9, font=2)
draw.circle(0.0, -0.18, 0.85, lwd=2, lty=4, border="darkblue")
text(0.0, 0.5, "1 miRNAs", cex=0.9, font=2)

dev.off()






#####
#END#
#####


df <- read.table("/Users/rehom/Desktop/MiRNA_PE_04132018_newanalysis_VDAART/Kevin_miRNA/To Kevin/MirTargetGenes.Down.tsv", sep="\t", header=TRUE)[,1:2]

ind=df$TargetGene %in% mapped_LCC$SYMBOL

table(ind)

df<-df[ind, ]

dim(df)

range(table(df$TargetGene))

df<-unique(df)

range(table(df$TargetGene))


five <- names(table(df$TargetGene)[table(df$TargetGene)==3])
four <- names(table(df$TargetGene)[table(df$TargetGene)==2])
three <- names(table(df$TargetGene)[table(df$TargetGene)==1])
df <- df[which(df$TargetGene %in% c(five, four, three)),]

#Shuffle the df in order to improve final layout
iRandUp <- sample(rep(1:nrow(df)))
df <- df[iRandUp,]

gUp <- graph.edgelist(as.matrix(df), directed=FALSE)

layoutUp <- data.frame(V(gUp)$name, layout.circle(gUp))

iMir <- grep("hsa-miR|mmu-miR", layoutUp[,1])
layoutUp[iMir,2:3] <- layoutUp[iMir,2:3] * 18
V(gUp)[iMir]$shape <- "square"
V(gUp)[iMir]$color = "firebrick1"

iThree <- which(layoutUp[,1] %in% three)
layoutUp[iThree,2:3] <- layoutUp[iThree,2:3] * 12
V(gUp)[iThree]$shape <- "sphere"
V(gUp)[iThree]$color = "yellow"

iFour <- which(layoutUp[,1] %in% four)
layoutUp[iFour,2:3] <- layoutUp[iFour,2:3] * 6
V(gUp)[iFour]$shape <- "sphere"
V(gUp)[iFour]$color = "royalblue"

iFive <- which(layoutUp[,1] %in% five)
layoutUp[iFive,2:3] <- layoutUp[iFive,2:3] * 1
V(gUp)[iFive]$shape <- "sphere"
V(gUp)[iFive]$color = "forestgreen"

#V(gUp)$name <- gsub("miR-", "miR-\n", V(gUp)$name)

V(gUp)$vertex.frame.color <- "white"

###

df <- read.table("../Results/MirTargetGenes.Down.tsv", sep="\t", header=TRUE)[,1:2]
five <- names(table(df$TargetGene)[table(df$TargetGene)==5])
four <- names(table(df$TargetGene)[table(df$TargetGene)==4])
three <- names(table(df$TargetGene)[table(df$TargetGene)==3])
df <- df[which(df$TargetGene %in% c(five, four, three)),]

#Shuffle the df in order to improve final layout
iRandDown <- sample(rep(1:nrow(df)))
df <- df[iRandDown,]

gDown <- graph.edgelist(as.matrix(df), directed=FALSE)

layoutDown <- data.frame(V(gDown)$name, layout.circle(gDown))

iMir <- grep("hsa-miR|mmu-miR", layoutDown[,1])
layoutDown[iMir,2:3] <- layoutDown[iMir,2:3] * 18
V(gDown)[iMir]$shape <- "square"
V(gDown)[iMir]$color = "red1"

iThree <- which(layoutDown[,1] %in% three)
layoutDown[iThree,2:3] <- layoutDown[iThree,2:3] * 12
V(gDown)[iThree]$shape <- "sphere"
V(gDown)[iThree]$color = "yellow"

iFour <- which(layoutDown[,1] %in% four)
layoutDown[iFour,2:3] <- layoutDown[iFour,2:3] * 6
V(gDown)[iFour]$shape <- "sphere"
V(gDown)[iFour]$color = "royalblue"

iFive <- which(layoutDown[,1] %in% five)
layoutDown[iFive,2:3] <- layoutDown[iFive,2:3] * 1
V(gDown)[iFive]$shape <- "sphere"
V(gDown)[iFive]$color = "forestgreen"

#V(gDown)$name <- gsub("miR-", "miR-\n", V(gDown)$name)

V(gDown)$vertex.frame.color <- "white"

pdf("MirToGeneNetwork8.pdf", width=12, height=6)
par(mfrow=c(1,2))
plot.igraph(gUp, layout=data.matrix(layoutUp[,2:3]), edge.curved=FALSE, vertex.size=7.0, vertex.label.dist=0.3, vertex.label.color="black", asp=FALSE, vertex.label.cex=0.6, edge.color="grey90", edge.width=0.1, edge.arrow.mode=0, main="Up-regulated miRNAs")

#Add circles to segregate the mirs based on number of genes targetting them
draw.circle(0.0, 0.0, 0.2, lwd=2, lty=4, border="darkblue")
text(0.0, 0.1, "3 miRNAs", cex=0.7, font=2)
draw.circle(0.0, 0.0, 0.5, lwd=2, lty=4, border="darkblue")
text(0.0, 0.5, "2 miRNAs", cex=0.7, font=2)
draw.circle(0.0, 0.0, 0.8, lwd=2, lty=4, border="darkblue")
text(0.0, 0.85, "1 miRNAs", cex=0.7, font=2)

###

plot.igraph(gDown, layout=data.matrix(layoutDown[,2:3]), edge.curved=FALSE, vertex.size=7.0, vertex.label.dist=0.3, vertex.label.color="black", asp=FALSE, vertex.label.cex=0.6, edge.color="grey90", edge.width=0.1, edge.arrow.mode=0, main="Down-regulated miRNAs")

#Add circles to segregate the mirs based on number of genes targetting them
draw.circle(0.0, 0.0, 0.15, lwd=2, lty=4, border="darkblue")
text(0.0, 0.0, "5 miRNAs", cex=0.7, font=2)
draw.circle(0.0, 0.0, 0.5, lwd=2, lty=4, border="darkblue")
text(0.0, 0.525, "4 miRNAs", cex=0.7, font=2)
draw.circle(0.0, 0.0, 0.8, lwd=2, lty=4, border="darkblue")
text(0.0, 0.85, "3 miRNAs", cex=0.7, font=2)
dev.off()


library(GeneOverlap)
go.obj <- newGeneOverlap(mapped_clu$SYMBOL, as.character(wObject$TargetGene), genome.size=20475)

go.obj <- testGeneOverlap(go.obj)
go.obj

print(go.obj)





