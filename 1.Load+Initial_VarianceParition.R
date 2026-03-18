rootdir="f:/OneDrive - Emory/Faundez_NULISA10/"  # SOMAplasmaMultibatch/Nets/UDS/3b.TestRegr_HNRNPA2B1_HBZ/"
setwd(rootdir)
dat<-read.csv(file="FaundezNULISA_forR.csv",header=TRUE,row.names=1)
assayNames<-rownames(dat)[14:nrow(dat)]

numericMeta<-as.data.frame(t(dat[1:13,]))
exprMat0<-dat[14:nrow(dat),]
exprMat0<-apply(exprMat0,2,as.numeric)
rownames(exprMat0)<-assayNames



colnames(numericMeta)
# [1] "Experiment"     "Foundation"     "LabSampleID"    "FamilyMember"   "Gene"           "MutationStatus" "SubjectID"      "Epilepsy"      
# [9] "APOE4"          "Sex"            "Age"            "GeneticVariant" "CNV.KdVS"

table(numericMeta$MutationStatus)
# Control Diseased 
#     148      317

numericMeta$Group<-numericMeta$MutationStatus
numericMeta$Group[numericMeta$MutationStatus=="Diseased"]<-numericMeta$Gene[numericMeta$MutationStatus=="Diseased"]

table(numericMeta$Group)
#CACNA1A   CDKL5 Control HNRNPH2  KANSL1   KCNQ2  MED13L  SLC2A1  STXBP1 SYNGAP1    TCF4 
#     18      23     148      37      41      22      20      43      56      43      14 

for (col in colnames(numericMeta)) print(paste0(col," : ",sum(table(numericMeta[,col]))))
#[1] "Experiment : 465"
#[1] "Foundation : 465"
#[1] "LabSampleID : 465"
#[1] "FamilyMember : 465"
#[1] "Gene : 465"
#[1] "MutationStatus : 465"
#[1] "SubjectID : 434"
#[1] "Epilepsy : 371"
#[1] "APOE4 : 465"
#[1] "Sex : 465"
#[1] "Age : 465"
#[1] "GeneticVariant : 372"
#[1] "CNV.KdVS : 45"
#[1] "Group : 465"

#We cannot have NA or missing values for VP.


## 3. Variance Partition pre adjustment as QC
regvars.vp<-data.frame(numericMeta)
regvars.vp$Age=as.numeric(numericMeta$Age)
regvars.vp$Sex<-factor(regvars.vp$Sex)
regvars.vp$Group<-factor(regvars.vp$Group)
#regvars.vp$Epilepsy<-factor(regvars.vp$Group)  # missing values
regvars.vp$APOE4<-factor(regvars.vp$APOE4)
regvars.vp$AB42<-as.numeric(exprMat0["AB42",])
#regvars.vp$tTau<-as.numeric(exprMat0["BD-MAPT",])
#regvars.vp$tauP217<-as.numeric(exprMat0["BD-pTau217",])
regvars.vp$Experiment<-factor(regvars.vp$Experiment)
regvars.vp$Foundation<-factor(regvars.vp$Foundation)

form <- ~Age+(1|Sex)+AB42+(1|Group)+(1|Experiment)+(1|Foundation)+(1|APOE4)

library(variancePartition)
parallelThreads=30

varPart1 <- fitExtractVarPartModel(exprMat0, form, regvars.vp,  BPPARAM=BiocParallel::SnowParam(workers = parallelThreads, type = "SOCK")) 
vp1 <- sortCols(varPart1,FUN=median,last= c("Residuals"))


pdf(file="1.VariancePartition-NULISA-130x465_unregr.pdf", width=15,height=11)
par(mfrow=c(1,1))
par(oma=c(2,2,3,1))

plotVarPart( vp1, main="NULISA 465 - Unregressed NPQ" )


DxSortOrder<-order(vp1[["Group"]],decreasing=TRUE)
#rownames(cleanDat.noNA.unreg)[DxSortOrder][1:50]
for (i in ls(vp1)) { vp1[[i]]<-vp1[[i]][DxSortOrder]; }
rownames(vp1)<-rownames(vp1)[DxSortOrder]

plotPercentBars( vp1[1:50,]) + ggtitle( "Top Group (Control/KANSL1/CDKL5..TCF4) Covariates" ) + theme(plot.margin = margin(t = 6, r = 1, b = 3, l = 1) )


DxSortOrder<-order(vp1[["Foundation"]],decreasing=TRUE)
#rownames(cleanDat.noNA.unreg)[DxSortOrder][1:50]
for (i in ls(vp1)) { vp1[[i]]<-vp1[[i]][DxSortOrder]; }
rownames(vp1)<-rownames(vp1)[DxSortOrder]

plotPercentBars( vp1[1:50,]) + ggtitle( "Top Foundation (Batch) Covariates" ) + theme(plot.margin = margin(t = 6, r = 1, b = 3, l = 1) )



DxSortOrder<-order(vp1[["Experiment"]],decreasing=TRUE)
#rownames(cleanDat.noNA.unreg)[DxSortOrder][1:50]
for (i in ls(vp1)) { vp1[[i]]<-vp1[[i]][DxSortOrder]; }
rownames(vp1)<-rownames(vp1)[DxSortOrder]

plotPercentBars( vp1[1:50,]) + ggtitle( "Top Experiment (Batch) Covariates" ) + theme(plot.margin = margin(t = 6, r = 1, b = 3, l = 1) )


DxSortOrder<-order(vp1[["APOE4"]],decreasing=TRUE)
#rownames(cleanDat.noNA.unreg)[DxSortOrder][1:50]
for (i in ls(vp1)) { vp1[[i]]<-vp1[[i]][DxSortOrder]; }
rownames(vp1)<-rownames(vp1)[DxSortOrder]

plotPercentBars( vp1[1:50,]) + ggtitle( "Top APOE4 Covariates" ) + theme(plot.margin = margin(t = 6, r = 1, b = 3, l = 1) )


DxSortOrder<-order(vp1[["AB42"]],decreasing=TRUE)
#rownames(cleanDat.noNA.unreg)[DxSortOrder][1:50]
for (i in ls(vp1)) { vp1[[i]]<-vp1[[i]][DxSortOrder]; }
rownames(vp1)<-rownames(vp1)[DxSortOrder]

plotPercentBars( vp1[1:50,]) + ggtitle( "Top Amyloid Beta 42 Covariates" ) + theme(plot.margin = margin(t = 6, r = 1, b = 3, l = 1) )


AgeSortOrder<-order(vp1[["Age"]],decreasing=TRUE)
#rownames(cleanDat.noNA.unreg)[AgeSortOrder][1:50]
for (i in ls(vp1)) { vp1[[i]]<-vp1[[i]][AgeSortOrder]; }
rownames(vp1)<-rownames(vp1)[AgeSortOrder]

plotPercentBars( vp1[1:50,]) + ggtitle( "Top Age Covariates" ) + theme(plot.margin = margin(t = 6, r = 1, b = 3, l = 1) )


SexSortOrder<-order(vp1[["Sex"]],decreasing=TRUE)
#rownames(cleanDat.noNA.unreg)[SexSortOrder][1:50]
for (i in ls(vp1)) { vp1[[i]]<-vp1[[i]][SexSortOrder]; }
rownames(vp1)<-rownames(vp1)[SexSortOrder]

plotPercentBars( vp1[1:50,]) + ggtitle( "Top Sex Covariates" ) + theme(plot.margin = margin(t = 6, r = 1, b = 3, l = 1) )


dev.off()


save.image("1.FaundezNULISA.RData")
