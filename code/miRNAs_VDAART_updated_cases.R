

load("miRNAs_for_hooman.RData")

dim(datExpr_miRNAs)

 mir=datExpr_miRNAs

 rn=gsub("ST.", "ST-", row.names(mir))

 load("final3regsvs_updatedpeexcludecases.rda")

 dim(final3regsvs_updatedpeexcludecases)

 pheno=pData(final3regsvs_updatedpeexcludecases)

 table(rn %in% pheno$SUBJECT)

 ind=which(rn %in% pheno$SUBJECT)

 mir_update_cases=mir[ind, ]

ind=which(rn %in% pheno$SUBJECT)

mir_update_cases=mir[ind, ]


na=apply(mir_update_cases, 2, function(x) sum(is.na(x)))
table(na)
range(na)

hist( apply(datExpr_miRNAs, 2, function(x) sum(is.na(x))), breaks=20)

hist( apply(datExpr_miRNAs, 2, function(x) sum(is.na(x))), breaks=20, freq=FALSE)

table(apply(datExpr_miRNAs, 2, function(x) sum(is.na(x)))<max(na)*0.2)

indnaless0.2na=which(apply(mir_update_cases, 2, function(x) sum(is.na(x)))<=max(na)*0.2)

indnaless0.2na=which(apply(mir_update_cases, 2, function(x) sum(is.na(x)))<50)

indnaless0.2na=which(apply(mir_update_cases, 2, function(x) sum(is.na(x)))<=max(na)*0.23)

mir_naless0.2na=mir_update_cases[,indnaless0.2na]

dim(mir_naless0.2na)

hist( apply(mir_naless0.2na, 2, function(x) sum(is.na(x))), breaks=20)

hist( apply(mir_naless0.2na, 2, function(x) sum(is.na(x))), breaks=20, freq=FALSE)

write.csv(mir_naless0.23na, file="mir_naless0.23na.csv")

mir_naless0.2na=read.csv("mir_naless0.2na.csv", header=T, row.names=1) #2 HK miRNAs removed
#dim(mir_naless0.2nawHK)


library("pcaMethods")

pc <- pca(mir_naless0.2na, nPcs=2, method="bpca")
pc
impmir <- completeObs(pc)
dim(impmir)

save(impmir, file="impir0.2.rda")

rn2=gsub("ST.", "ST-", row.names(impmir))

identical(pheno$SUBJECT,rn2)



mirdf=impmir
row.names(mirdf)<-pheno$sampleID
identical(colnames(t(mirdf)),rownames(pheno))
head(mirdf)

mirS4=final3regsvs_updatedpeexcludecases
exprs(mirS4)<-t(mirdf)



phenoData <- AnnotatedDataFrame(data=pheno)
eset <- ExpressionSet(cleany, phenoData=phenoData)
mirS4<- ExpressionSet(t(mirdf), phenoData=phenoData)

library(sva)
mod=model.matrix(~pefull+basevitdng, data=pheno)
mod0=model.matrix(~basevitdng, data=pheno)
res.sva=sva(dat=exprs(mirS4), mod=mod, mod0 = mod0)


#mod=model.matrix(~as.factor(pefull), data=pheno)
#mod0=model.matrix(~1, data=pheno)
#res.sva=sva(dat=exprs(mirS4), mod=mod, mod0 = mod0)



mod=model.matrix(~pefull+basevitdng, data=pheno)
mod0=model.matrix(~basevitdng, data=pheno)
res.sva=sva(dat=exprs(mirS4), mod=mod, mod0 = mod0)

cleaningY = function(y, mod, svaobj) {
X=cbind(mod,svaobj$sv)
Hat=solve(t(X)%*%X)%*%t(X)
beta=(Hat%*%t(y))
P=ncol(mod)
cleany=y-t(as.matrix(X[,-c(1:P)])%*%beta[-c(1:P),])
return(cleany)
}

cleany = cleaningY(exprs(mirS4),mod,res.sva)
save(cleany, file="regexrs.mirRNA.RData")
dim(cleany)


es5=mirS4
exprs(es5)<-cleany


head(colnames(exprs(es5)))

head(rownames(pData(es5)))

library("RankProd")
RP.outpemir <- RP(exprs(es5),pData(es5)$pefull, num.perm=1000, logged=TRUE, na.rm=FALSE,plot=FALSE, rand=123)

pemir=topGene(RP.outpemir,cutoff=0.05,method="pfp",logged=FALSE,logbase=2)

save(RP.outpemir, file="RP.outpemir0.2.rda")


pData(es5)$vitdcontrast2<-NA
pData(es5)$vitdcontrast2[which(es5$basevitdng<30)]<-1  #this will be considered 0
pData(es5)$vitdcontrast2[which(es5$basevitdng>=30)]<-2  #this will be considered 1
table(pData(es5)$vitdcontrast2)


pData(es5)$vitdcontrast2<-NA
pData(es5)$vitdcontrast2[which(es5$basevitdng<30)]<-1  #this will be considered 0
pData(es5)$vitdcontrast2[which(es5$basevitdng>=30)]<-0  #this will be considered 1
table(pData(es5)$vitdcontrast2)



RP.outvdmir <- RP(exprs(es5),pData(es5)$vitdcontrast2, num.perm=1000, logged=TRUE, na.rm=FALSE,plot=FALSE, rand=123)

vdmir=topGene(RP.outvdmir,cutoff=0.05,method="pfp",logged=FALSE,logbase=2)
save(RP.outvdmir, file="RP.outvdmir0.2.rda")


inttab1=intersect(vdmir$Table1[,1], pemir$Table1[,1])
inttab2=intersect(vdmir$Table2[,1], pemir$Table2[,1])

row.names(exprs(es5[inttab1, ])) # "hsa.miR.34a.3p"  #0.25-> "hsa-miR-182-5p" "hsa-miR-34a-3p"
row.names(exprs(es5[inttab2, ])) # "hsa.miR.145.3p"        "hsa.miR.424.5p"        "hsa.miR.378a.5p"       "hsa.miR.550a.5p"       "hsa.miR.375"           "hsa.miR.31.5"          "hsa.miR.886.3p_002194" "hsa.miR.218.5p"    
                                 #  0.25 ->  "hsa-miR-145-3p"  "hsa-miR-424-5p"  "hsa-miR-378a-5p" "hsa-miR-31-5" 
length(setdiff(inttab2, inttab1))
length(setdiff(inttab1, inttab2))
intersect(inttab2, inttab1)

mircand=union(setdiff(inttab2, inttab1),setdiff(inttab1, inttab2))
length(mircand)

length(intersect(vdmir$Table1[,1], pemir$Table1[,1]))


row.names(exprs(es5[dfvdmir$gene.index, ]))

length(intersect(vdmir$Table2[,1], pemir$Table2[,1]))
intersect(vdmir$Table2[,1], pemir$Table2[,1])




write.csv(row.names(exprs(es5)[mircand,]), file="mircan_updatedpecases.csv")


exprs(es5[mircand,])



 write.csv(row.names(exprs(es5)[mircand,]), file="mircan_updatedpecases_2.csv")


 dfpemir <- ldply (pemir, data.frame)

 dfvdmir <- ldply (vdmir, data.frame)

 mirpecand=row.names(exprs(es5[dfpemir$gene.index, ]))

 write.csv(mirpecand, file="mirpecand0.25.csv", row.names=FALSE)

 mirvdcan=row.names(exprs(es5[dfvdmir$gene.index, ]))

 write.csv(mirvdcan, file="mirvdcand0.25.csv", row.names=FALSE)


table(mirvdcan %in% mirpecand)

FALSE  TRUE 
   10    17 

length(intersect(dfpemir$gene.index, dfvdmir$gene.index))

 mircandmultidirec=intersect(dfpemir$gene.index, dfvdmir$gene.index)

 exprs(es5[mircandmultidirec, ])

write.csv(row.names(exprs(es5)[ mircandmultidirec,]), file="mircandmultidirec_updatedpecases0.25.csv")



###########


mirde=read.csv(" mirpecand_updatedpecases.csv")
indnaless150=which(apply(datExpr_miRNAs, 2, function(x) sum(is.na(x)))<40)
mir_naless150=datExpr_miRNAs[,indnaless150]

dim( mir_naless150)
mirnonde=colnames(mir_naless150)
table(mirdename %in% mirnonde)

mirdename[!mirdename %in% mirnonde]

table(mirpecand %in% mirde$miRNAs)


Housekeeping miRNAs
ath-miR159a_000338_(AVG)
RNU44_001094_(AVG)
RNU48_001006_(AVG)
U6 rRNA_001973_(AVG)


write.csv(mirpecan, file="mirpecan0.25.csv")
write.csv(mirvdcan, file="mirvdcan0.25.csv")
write.csv(dfvdmir, file="dfvdmir0.25.csv")
write.csv(dfpemir, file="dfpemir0.25.csv”)
pemir=topGene(RP.outpemir,cutoff=0.05,method="pfp",logged=TRUE,logbase=2)

