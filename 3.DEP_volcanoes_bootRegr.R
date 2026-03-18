rootdir3<-"f:/OneDrive - Emory/Faundez_NULISA10/"
setwd(rootdir3)
load("./2.FaundezNULISA_regr.RData")

outputfigs<-outputtabs<-rootdir<-rootdir3
setwd(rootdir)



table(numericMeta$Group)
#CACNA1A   CDKL5 Control HNRNPH2  KANSL1   KCNQ2  MED13L  SLC2A1  STXBP1 SYNGAP1    TCF4 
#     18      23     148      37      41      22      20      43      56      43      14

which(is.na(numericMeta$Group))
#integer(0)


which(!rownames(numericMeta)==colnames(cleanDat))
# all match/none don't


# One more check: Are Groups of genetic disease separate in MDS reduced dimension space?
par(mfrow=c(1,2))
#[1]  172 762
sampleMDS<-limma::plotMDS(cleanDat, col=paste0(gplots::col2hex(WGCNA::labels2colors(as.numeric(factor(numericMeta$Group)))),"88"),main=paste0("MDS of NPQ 130 NULISA assays N=",ncol(cleanDat)," samples;\ncolored by disease gene"),pch=16)

Group<-factor(numericMeta$Group)
groupColors=paste0(gplots::col2hex(WGCNA::labels2colors(as.numeric(factor(names(table(numericMeta$Group)))))),"88")

x.mean <- tapply(sampleMDS$x, Group, mean)
y.mean <- tapply(sampleMDS$y, Group, mean)
# guarantee ordering of means by Group (factor levels)
x.mean <- x.mean[levels(Group)]
y.mean <- y.mean[levels(Group)]
plot(
  x.mean, y.mean,
  xlab = "logFC Dimension 1",
  ylab = "logFC Dimension 2",
  pch = 16,
  col=groupColors,
  cex=4,
  main= "Mean (centroid) of All Samples\nwithin Each Group",
  xlim=c(-1.4,0.55)
)

legend("topleft", legend=names(table(numericMeta$Group)), bty='n', col=groupColors, pch=16, cex=1.5, pt.cex=2.5)  # , title="Sample Groups"

cleanDat.MDS<-recordPlot()
pdf(file=paste0("3.MDS_sampleDifferencesInNULISAdata(",nrow(cleanDat),"x",ncol(cleanDat),")_bootRegr.pdf"),height=8,width=16)
  print(cleanDat.MDS)
dev.off()

# Answer: - not drastically separated in 2D using variance of full data for 130 assays.





######################################
## ANOVA + Volcanoes + DEXstacked Barplots
source("./parANOVA.dex.R")

parallelThreads=31
outFilePrefix="3"
outFileSuffix="Faundez_NULISA-regrBatch"
outputCSV=FALSE


numericMeta.full<-numericMeta
cleanDat.full<-cleanDat


ANOVAoutList<-list()
diseaseGroups=c("CACNA1A", "CDKL5", "HNRNPH2", "KANSL1", "KCNQ2", "MED13L", "SLC2A1", "STXBP1", "SYNGAP1", "TCF4")
for (group in diseaseGroups) {

  numericMeta<-numericMeta.full[which(numericMeta.full$Group %in% c("Control",group)),]
  cleanDat<-cleanDat.full[,match(rownames(numericMeta),colnames(cleanDat.full))]
  
  Grouping=numericMeta$Group
  ANOVAoutList[[group]] <- parANOVA.dex()
}


cleanDat<-cleanDat.full
numericMeta<-numericMeta.full


## Assemble ANOVAout
third.cols <- sapply(ANOVAoutList, function(df) df[[3]])
fourth.cols <- sapply(ANOVAoutList, function(df) df[[4]])
colnames(third.cols)[1:2]<-paste0("Control-",colnames(third.cols)[1:2])
colnames(third.cols)[3:10]<-paste0(colnames(third.cols)[3:10],"-Control")
colnames(fourth.cols)<-paste0("diff ",colnames(third.cols))

ANOVAout.unadj <- data.frame(
  dummyCol1 = NA,
  dummyCol2 = NA,
  third.cols,
  fourth.cols,
  check.names = FALSE
)
rownames(ANOVAout.unadj)<-rownames(ANOVAoutList[[1]])


adj.pVal.cols<-t(apply(third.cols,1,function(x) p.adjust(x, method="BH")))

ANOVAout.adj <- data.frame(
  dummyCol1 = NA,
  dummyCol2 = NA,
  adj.pVal.cols,
  fourth.cols,
  check.names = FALSE
)
rownames(ANOVAout.adj)<-rownames(ANOVAoutList[[1]])


ANOVAout<-ANOVAout.adj
flip=c(3:4)
sameScale=TRUE
symbolsOnly=TRUE
highlightGeneProducts=c("SMOC1","NEFL","APOE","ABETA42","MAPT","pTau-217","pTDP43-409")
labelHighlighted=TRUE      # if true, highlighted spots get text labels with their rownames from ANOVAout
labelTop=5
useNETcolors=FALSE
plotVolc()                 # runs on ANOVAout as input (need not be specified).

#DEXpercentStacked()        # runs on prior function outputs as input; writes stacked bar plot(s) to PDF.


write.csv(ANOVAout,file="3.FDR(BH-corrected_p_values)_and_log2EffectSizes-NULISA_130assays-10gene-drivenDiseaseSampleSets_vs_148ControlSamples.csv")





######################################
## Alternative Correlation (to APOE4, epilepsy, NEFL), stats table with volcanoes

table(numericMeta$APOE4)
# No Yes 
#355 110  # 465 total, no missing
numericMeta$APOE4.01 <- 0
numericMeta$APOE4.01[numericMeta$APOE4=="Yes"] <- 1

table(numericMeta$Epilepsy)
# No Yes 
#222 149
numericMeta$Epilepsy.01<-numericMeta$Epilepsy
numericMeta$Epilepsy.01[numericMeta$Epilepsy=="Yes"]<- 1
numericMeta$Epilepsy.01[numericMeta$Epilepsy=="No"]<- 0
numericMeta$Epilepsy.01 <- as.numeric(numericMeta$Epilepsy.01)

numericMeta$'Oligo-SNCA'<-cleanDat["Oligo-SNCA",]
numericMeta$NEFL<-cleanDat["NEFL",]
numericMeta$HBA1<-cleanDat["HBA1",]


# These parameters are specific to trait correlation statistics generation; traits are provided as columns of the data frame stored in the provided example RData as the variable numericMeta. 
cor.traits=c("HBA1", "NEFL", "Oligo-SNCA", "APOE4.01","Epilepsy.01")                 # Molecular and quantitative Traits to correlate to in numericMeta columns (colnames)
filter.trait="Group"             # Trait on which to subset case samples
filter.trait.subsets=c("ALL",diseaseGroups) # Subsets of case samples will be used for correlation to the cell type proportion estimates
                                    # (7 separate cor.traits x 3 sample subsets = 21 total p and R value columns to generate)
corFn="bicor"                       #'bicor'; anything else will cause Pearson (cor) to be used


#source("./parANOVA.dex.R")
outputCSV=TRUE
outFilePrefix="3b"
outFileSuffix="Faundez_NULISA-regrBatch"
CORout <- trait.corStat()                      # runs on cleanDat and Grouping variables as required input.
# Correlation p + R table calculations complete. If you want to use the table with plotVolc(), set the variable corVolc=TRUE and use variable CORout to store the table generated.

dim(CORout)
#[1] 130 112


corVolc=TRUE        # changes the behavior of plotVolc, DEXpercentStacked, and GOparallel functions later in the pipeline, to use CORout
useNETcolors=FALSE
sameScale=TRUE
#highlightGeneProducts=rownames(CORout)[which(CORout$Filter.ALLsamples.FavoriteCorrStat.Sig)]  # Specifies which spots should be large
#labelHighlighted=TRUE

#Change values less than 1e-25 to 1e-10, so same scale volcanoes will not be stretched to -log10(p) of 200!
CORout[,3:57]<-apply(CORout[,3:57],2,function(x) { x1<-x; x1[x1<1e-20]<-1e-20; x1; })

rm(selectComps)
plotVolc()          # Plots PDFs and HTMLs

#DEXpercentStacked(CORout) 


save.image("3.FaundezNULISA-DEPandBicor.RData")




## NOT RUN BELOW - ONTOLOGY ENRICHMENT ONLY RELEVANT FOR UNBIASED PROTEOMICS.

source("./GOparallel-FET.R")  #Available from the repository file https://github.com/edammer/GOparallel/blob/main/GOparallel-FET.R
ANOVAgroups=TRUE
parallelThreads=31
outFilename="3b.CORsig.GO"
GOparallel(CORout)


# Going forward, do not use correlation statistics, use ANOVA groupwise stats.
corVolc=FALSE

## This section's outputs moved to subfolder: "3b.Correlation Stats"



##############################################################################
## GOparallel - v1.2 - with hitLists included in .csv outputs.

source("./GOparallel-FET.R")

modulesInMemory=TRUE            # uses cleanDat, net, and kMEdat from pipeline already in memory
ANOVAgroups=FALSE               # if true, modulesInMemory ignored. Volcano pipeline code should already have been run!
outFilename="3c.NULISA_regrNet.GO"
  
GOparallel()


# also skipped
modulesInMemory=FALSE            # uses cleanDat, net, and kMEdat from pipeline already in memory
ANOVAgroups=TRUE                 # if true, modulesInMemory ignored. Volcano pipeline code should already have been run!
testIndexMasterList=c(3,4,5)
outFilename="3d.NULISA_DEP_lists_byGroup.GO"

GOparallel()
