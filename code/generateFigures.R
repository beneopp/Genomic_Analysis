require(ComplexHeatmap)
require(Cairo)
require(textshape)

setwd("/Users/benoppenheimer/Documents/Source_Data_M1_expression")

load("Rdata/targetPredictions.rda")
load("Rdata/DEmiRNAs.rda")

LCCgenes <- read.csv("geneLists/LCC_new.csv", header = FALSE)$V1

alter_fun <- list(
  background = function(x, y, w, h) {
     grid.rect(x, y, w-unit(0.5, "mm"), height=unit(0.5, "mm"), gp=gpar(fill="grey95", col="grey95"))
  },
  Target=function(x, y, w, h) {
     grid.rect(x, y, w-unit(2.5, "mm"), h-unit(2.5, "mm"), gp=gpar(fill="black", col="black"))
  }
)

CreateSubMirPrint <- function(reg_type) {
  
  miRNAsOfInterest <- rownames(DEmiRNAs)[DEmiRNAs$Expression == reg_type]
  LCCmirPrint <- mirPrint[rownames(mirPrint) %in% LCCgenes, 
                          colnames(mirPrint) %in% miRNAsOfInterest]
  LCCmirPrint <- LCCmirPrint[rowSums(LCCmirPrint) != 0, 
                             colSums(LCCmirPrint) != 0]
  LCCmirPrint <- cluster_matrix(LCCmirPrint, method = "complete")
  
}

GetOncoPrint <- function(dat, name) {
  oncoPrint(t(dat), name=name, 
            
            show_pct = FALSE,
            
            alter_fun=alter_fun, alter_fun_is_vectorized = FALSE,
            
            col=c("Target"="black"),
            
            show_row_names=TRUE, row_names_gp=gpar(fontsize=10, fontface="bold"), row_names_max_width=unit(3, "cm"),
            
            show_column_names=TRUE, column_names_gp=gpar(fontsize=10),
                                            
            pct_gp=gpar(fontsize=12, fontface="bold"),
            
            heatmap_legend_param=list(title="Mir gene targets", at=c("Target"),
                                      labels=c("Target"), nrow=1, title_position="topleft")
  )
}


PlotMatrix <- function(reg_type, LCCmirPrint) {
  
  LCCmirPrint <- ifelse(LCCmirPrint, "Target", "")
  
  if (reg_type == "Up") {
    name <- "mirPrintPlotDown"
  } else {
    name <- "mirPrintPlotUp"
  }
  
  mirPrintPlot <- GetOncoPrint(LCCmirPrint, name)
}

CreateAnnotation <- function(reg_type, LCCmirPrint) {
  
  if (reg_type == "Up") {
    reg_type <- "Down-regulated"
  } else {
    reg_type <- "Up-regulated"
  }
  
  colors <- list("Regulation"=c("Down-regulated"="yellow", "Up-regulated"="royalblue"))
  regulation <- replicate(ncol(LCCmirPrint), reg_type)
  mirSigs <- -log10(DEmiRNAs$pfp[match(colnames(LCCmirPrint), rownames(DEmiRNAs))])
  arbitraryMax <- as.integer(max(mirSigs[is.finite(mirSigs)])) + 1
  mirSigs[is.infinite(mirSigs)] <- arbitraryMax
  dfGeneAnno <- data.frame("Mir significance Score" = mirSigs,
                           Regulation = regulation, check.names = FALSE)
  mirAnnoHeat <- rowAnnotation(df=dfGeneAnno, col=colors, width=unit(1,"cm"))
}

GetImages <- function(reg_type) {
  
  if (!(reg_type %in% c("Up", "Down"))) {
    stop("reg_type parameter not accepted for GetImages function.")
  }
  
  LCCmirPrint <- CreateSubMirPrint(reg_type)
  
  oncoPrintImage <- PlotMatrix(reg_type, LCCmirPrint)
  
  anno <- CreateAnnotation(reg_type, LCCmirPrint)
  
  results <- list(oncoPrintImage, anno)
}

upResults <- GetImages("Down")
plotUp <- upResults[[1]]
annoUp <- upResults[[2]]

downResults <- GetImages("Up")
plotDown <- downResults[[1]]
annoDown <- downResults[[2]]

CairoPNG(filename = "images/plot.png", width = 800, height = 800)
pushViewport(viewport(layout=grid.layout(nr=2, nc=1)))

pushViewport(viewport(layout.pos.row=1, layout.pos.col=1))
draw(plotUp + annoUp, heatmap_legend_side="top", annotation_legend_side="left", newpage=FALSE)
upViewport()

pushViewport(viewport(layout.pos.row=2, layout.pos.col=1))
draw(plotDown + annoDown, heatmap_legend_side="top", annotation_legend_side="left", newpage=FALSE)
upViewport()

dev.off()



