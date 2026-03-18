#####################################################################################
## WGCNA Coexpression Network Build - Seyfried Systems Biology Pipeline
##
## input cleanDat from:  plasma transposed 2 batch variable regressed NULISA NPQ data
##
## Eric Dammer -- Seyfried Lab Proteomics (2022)  run 1/24/2026 - 1/26/2026
#####################################################################################

rootdir1="F:/OneDrive - Emory/Faundez_NULISA10/4.subtyping/"
setwd(rootdir1)

load("../3.FaundezNULISA-DEPandBicor.RData")
#contains: cleanDat,numericMeta
rootdir<-rootdir1


dim(cleanDat)
#[1] 130 465





#=============================#
#  Check and Remove Outliers  #
#=============================#

if(!exists("numericMeta")) numericMeta<-traits
library(WGCNA)
enableWGCNAThreads()

sdout=3 #Z.k SD fold for outlier threshold
outliers.noOLremoval<-outliers.All<-vector()
cleanDat.noOLremoval<-cleanDat

breakoutFlag=FALSE
repeated=0
while (!breakoutFlag) {
repeated=repeated+1
normadj <- (0.5+0.5*bicor(cleanDat,use="pairwise.complete.obs")^2)

## Calculate connectivity
netsummary <- fundamentalNetworkConcepts(normadj)
ku <- netsummary$Connectivity
z.ku <- (ku-mean(ku))/sqrt(var(ku))  #corrected, moved parenthesis open to leading position
## Declare as outliers those samples which are more than sdout sd above the mean connectivity based on the chosen measure
outliers <- (z.ku < mean(z.ku)-sdout*sd(z.ku))  #previously had | z.ku > mean(z.ku)+sdout*sd(z.ku)) to remove outliers on high connectivity end...
print(paste0("There are ",sum(outliers)," outlier samples based on a bicor distance sample network connectivity standard deviation above ",sdout,".  [Round ",repeated,"]"))
targets.All=numericMeta

if (length(which(outliers))==0) breakoutFlag=TRUE
cleanDat <- cleanDat[,!outliers] 
numericMeta <- targets <- targets.All[!outliers,]
outliers.All<-c(outliers.All,outliers)
} #repeat until no more outliers

#All outliers removed
print(paste0("There are ",sum(outliers.All)," total outlier samples removed in ",repeated," iterations:"))
names(which(outliers.All))
outliersRemoved<-names(which(outliers.All))
#Note outliers as comment below, copied from R session.
#[1] "There are 5 total outlier samples removed in 3 iterations:"
#[1] "A_01_P2_2"   "A_01_P1_1"   "A_01_P4_4"   "C_11_P5_294" "D_09_P5_288"


## Enforce <50% missingness (1 less than half of cleanDat columns (or round down half if odd number of columns))
LThalfSamples<-length(colnames(cleanDat))/2
LThalfSamples<-LThalfSamples - if ((length(colnames(cleanDat)) %% 2)==1) { 0.5 } else { 1.0 }

## If operating on log2(FPKM) data, remove rows with >=50% originally 0 FPKM values (only if there are some rows to be removed)
#IndexHighMissing<-rowsRemoved<-zeroVarRows<-vector()
#temp2<-data.frame(ThrowOut=apply(cleanDat,1,function(x) length(x[x==log2(0+0.05)])>LThalfSamples))
#cleanDat<-cleanDat[!temp2$ThrowOut,]
#dim(cleanDat) #still have x genes, now for y total samples

## If working on log2(protein abundance or ratio) with NA missing values; Enforce <50% missingness (1 less than half of cleanDat columns (or round down half if odd number of columns))
#remove rows with >=50% missing values (only if there are some rows to be removed)
IndexHighMissing<-rowsRemoved<-zeroVarRows<-vector()
temp2<-as.data.frame(cleanDat[which(rowSums(as.matrix(is.na(cleanDat)))>LThalfSamples),])
#handle condition if temp2 is for one row of cleandat (a vector instead of a data frame)
if (ncol(temp2)==1) {
  temp2<-t(temp2)
  rownames(temp2)=rownames(cleanDat)[which(rowSums(as.matrix(is.na(cleanDat)))>LThalfSamples)]
}

if (nrow(temp2)>0) { IndexHighMissing=which(rowSums(as.matrix(is.na(cleanDat)))>LThalfSamples); rowsRemoved<-rownames(cleanDat)[IndexHighMissing]; cleanDat<-cleanDat[-IndexHighMissing,]; }

dim(cleanDat)
#130 x460


# We will get eigenproteins on the 130 protein x 460 sample (5 outliers removed), then go back to the 130 x 465 cleanDat, and run subtyping (transposed WGCNA) on that
# numericMeta.full, cleanDat.full for 465 samples in RData loaded




## Z-transform rows of cleanDat after transposed TAMPOR (no missing data)
#cleanDat <- t(apply(cleanDat.CSF.noNA.mpPlatform,1,function(x) (x-mean(x))/sd(x) ))
#no Z transform (Z score has no effect on network)








projectFilesOutputTag="NULISA10.plasma-eigenproteinWGCNA"
Grouping=numericMeta$Group
outputtabs<-outputfigs<-rootdir




## WGCNA blockwiseModules (for signed bicor coexpression network)
##############################################################################

##Above code may have borrowed server for parallel processing of bootstrap regression -- return parallel processing to local workstation
#library("doParallel")
#stopCluster(clusterLocal)
#parallelThreads=30 #set to # of threads on your computer
#clusterLocal <- makeCluster(c(rep("localhost",parallelThreads)),type="SOCK")
#registerDoParallel(clusterLocal)
enableWGCNAThreads() #speeds the pickSoftThreshold function

powers <- seq(4,28,by=1)  #initial power check -- try to get SFT.R.sq to go > 0.80
sft <- pickSoftThreshold(t(cleanDat),blockSize=nrow(cleanDat)+1000,   #always calculate power within a single block (blockSize > # of rows in cleanDat)
                         powerVector=powers,
                         corFnc="bicor",networkType="signed")

##Paste/replace output below and annotate based on graphical plot that follows
# NULISA10 Plasma 460 samples, protein WGCNA power selection:
   Power SFT.R.sq  slope truncated.R.sq mean.k. median.k. max.k.
1      4   0.0491 -0.238         0.1810  17.400  16.20000  30.50
2      5   0.5330 -0.793         0.4300  11.500   9.80000  24.70
3      6   0.7530 -1.130         0.6850   7.990   6.09000  20.90
4      7   0.7920 -1.130         0.7740   5.770   3.91000  18.20  << POWER=7
5      8   0.7820 -1.090         0.7900   4.340   2.56000  16.20
6      9   0.8290 -1.040         0.8720   3.390   1.65000  14.70
7     10   0.7650 -1.060         0.7730   2.730   1.07000  13.50
8     11   0.7090 -1.010         0.7120   2.260   0.69700  12.50
9     12   0.7780 -0.972         0.7950   1.920   0.45900  11.60
10    13   0.7300 -0.944         0.7410   1.660   0.30500  10.80
11    14   0.7600 -0.892         0.7820   1.450   0.21000  10.20
12    15   0.8120 -0.913         0.7940   1.290   0.14600   9.57
13    16   0.6850 -0.935         0.6690   1.160   0.10200   9.03
14    17   0.1040 -1.450        -0.1230   1.050   0.06900   8.53
15    18   0.7750 -0.895         0.7460   0.951   0.04760   8.07
16    19   0.8470 -0.876         0.8270   0.869   0.03290   7.65
17    20   0.8500 -0.843         0.8280   0.798   0.02340   7.26
18    21   0.7690 -0.851         0.7510   0.736   0.01650   6.89
19    22   0.0861 -1.160        -0.0971   0.681   0.01180   6.56
20    23   0.0872 -1.160        -0.0969   0.632   0.00866   6.24
21    24   0.8430 -0.856         0.8080   0.588   0.00632   5.94
22    25   0.8410 -0.834         0.8080   0.549   0.00449   5.66
23    26   0.1410 -1.410        -0.1040   0.513   0.00320   5.40
24    27   0.1400 -1.370        -0.1050   0.480   0.00228   5.16
25    28   0.1410 -1.350        -0.1040   0.450   0.00163   4.92


#plot initial SFT.R.sq vs. power curve
tableSFT<-sft[[2]]
plot(tableSFT[,1],tableSFT[,2],xlab="Power (Beta)",ylab="SFT R²")


#choose power at elbow of SFT R² curve approaching asymptote near or ideally above 0.80
power=7


## Run an automated network analysis (ds=4 and mergeCutHeight=0.07, more liberal)
# choose parameters deepSplit and mergeCutHeight to get respectively more modules and more stringency sending more low connectivity genes to grey (not in modules).
net <- blockwiseModules(t(cleanDat),power=power,deepSplit=4,minModuleSize=6,
                        mergeCutHeight=0.07,TOMdenom="mean", #detectCutHeight=0.9999,                        #TOMdenom="mean" may get more small modules here.
                        corType="bicor",networkType="signed",pamStage=TRUE,pamRespectsDendro=TRUE,
                        verbose=3,saveTOMs=FALSE,maxBlockSize=nrow(cleanDat)+1000,reassignThresh=0.05)       #maxBlockSize always more than the number of rows in cleanDat
#blockwiseModules can take 30 min+ for large numbers of gene products/proteins (10000s of rows); much quicker for smaller proteomic data sets

nModules<-length(table(net$colors))-1
modules<-cbind(colnames(as.matrix(table(net$colors))),table(net$colors))
orderedModules<-cbind(Mnum=paste("M",seq(1:nModules),sep=""),Color=labels2colors(c(1:nModules)))
modules<-modules[match(as.character(orderedModules[,2]),rownames(modules)),]
as.data.frame(cbind(orderedModules,Size=modules))

##copy R session output;
## NULISA10 Plasma 460 samples, protein WGCNA power selection; deepSplit=4; power=7; pamRespectsDentro=TRUE; MMS=10:
          Mnum     Color Size
turquoise   M1 turquoise   45
blue        M2      blue   40

## NULISA10 Plasma 460 samples, protein WGCNA power selection; deepSplit=4; power=6; pamRespectsDentro=FALSE; MMS=10:
          Mnum     Color Size
turquoise   M1 turquoise   44
blue        M2      blue   41

## Positive age-associated module 3 seen on dendrogram in preliminary Global Network Plots page 1, adjusted params to power=7, MMS=6:
          Mnum     Color Size
turquoise   M1 turquoise   45
blue        M2      blue   40
brown       M3     brown    7

# (we will stick with the power=7 module definitions)
net.ds4<-net


#we will explore the blockwiseModules() function-built network with parameter deepSplit=4 in the output for folder ii.Zscore_unionNetwork_ds4 (ds3 already explored)
net<-net.ds4
outputfigs<-outputtabs<-rootdir

#<SKIP> (M3 is 7 members, last module--may still grow further with iterative cleanup)
minModSize=10
# If necessary, return module members of small modules below size minSize=X to grey
if (enforceMMS) {
  removedModules<-orderedModules[which(modules<minModSize),"Color"]
  for(i in removedModules) { net$colors[net$colors==i] <- "grey" }
  for(i in removedModules) { net$MEs[,paste0("ME",i)] <- NULL }

  nModules<-length(table(net$colors))-1
  modules<-cbind(colnames(as.matrix(table(net$colors))),table(net$colors))
  orderedModules<-cbind(Mnum=paste("M",seq(1:nModules),sep=""),Color=labels2colors(c(1:nModules)))
  modules<-modules[match(as.character(orderedModules[,2]),rownames(modules)),]
  as.data.frame(cbind(orderedModules,Size=modules))
}
minModSize=10
#<END SKIP>


#calculate kME table up front, in case we need to correct color assignments
MEs<-tmpMEs<-data.frame()
MEList = moduleEigengenes(t(cleanDat), colors = net$colors)
MEs = orderMEs(MEList$eigengenes)
net$MEs <- MEs
colnames(MEs)<-gsub("ME","",colnames(MEs)) #let's be consistent in case prefix was added, remove it.
rownames(MEs)<-rownames(numericMeta)

tmpMEs <- MEs #net$MEs
colnames(tmpMEs) <- paste("ME",colnames(MEs),sep="")
MEs[,"grey"] <- NULL
tmpMEs[,"MEgrey"] <- NULL

kMEdat <- signedKME(t(cleanDat), tmpMEs, corFnc="bicor")


table(net$colors)["grey"]  #38; was 45 without M3; power=7; ds=4; MMS=6



##ITERATIVE until condition met that all module membes are at least 0.28 kMEintramodule.
#Go back and do final algorithm fix of module colors (remove kMEintramodule<0.28 members, reassign grey with kMEintramodule>0.35; max difference from kMEmax<0.10)

retry=TRUE;
kMEmaxDiff=0.1
reassignIfGT=0.30
greyIfLT=0.30
iter=1;
while (retry) {
  cat(paste0("\nkME table Cleanup, processing iteration ",iter,"..."))
  colorVecFixed<-colorVecBackup<-net$colors
  orderedModulesWithGrey=rbind(c("M0","grey"),orderedModules)
  kMEintramoduleVector<-apply( as.data.frame(cbind(net$colors,kMEdat)),1,function(x) as.numeric(x[which(colnames(kMEdat)==paste0("kME",x[1]))+1]) )  #all sig digits (no rounding), so max will be unique.
  colorVecFixed[kMEintramoduleVector<greyIfLT]<-"grey"
  kMEmaxVec<-apply( as.data.frame(kMEdat),1,function(x) max(x) )
  kMEmaxColorsVec<-apply( as.data.frame(cbind(kMEmaxVec,kMEdat)),1, function(x) gsub("kME","",colnames(kMEdat)[which(x==x[1])[2]-1]) )
  kMEintramoduleVector<-unlist(lapply(kMEintramoduleVector,function(x) if(length(x)==0) { 1 } else { x }))   #grey will be ignored in checking for kMEmaxVec-kMEintramoduleVector difference max
  kMEmaxDiffTooBig<-(kMEmaxVec-kMEintramoduleVector) >= kMEmaxDiff
  colorVecFixed[which( (colorVecFixed=="grey" & kMEmaxVec>reassignIfGT) | kMEmaxDiffTooBig )] <- kMEmaxColorsVec[which( (colorVecFixed=="grey" & kMEmaxVec>reassignIfGT) | kMEmaxDiffTooBig )]
  net$colors<-colorVecFixed

#  table(net$colors)["grey"]  #decreased to x


# Are colors still in rank order? -- put them in order by recoloring modules that changed rank
  sort(table(net$colors),decreasing=TRUE)[!names(sort(table(net$colors),decreasing=TRUE))=="grey"]

  oldcolors <- names(sort(table(net$colors),decreasing=TRUE)[!names(sort(table(net$colors),decreasing=TRUE))=="grey"])
  for (i in 1:length(oldcolors)) {
    net$colors[net$colors==oldcolors[i]]<-paste0("proxy",labels2colors(i))
  }
  for (i in 1:length(oldcolors)) {
    net$colors[net$colors==paste0("proxy",labels2colors(i))]<-labels2colors(i)
  }

# one can check that colors are in order by size now
  #sort(table(net$colors),decreasing=TRUE)[!names(sort(table(net$colors),decreasing=TRUE))=="grey"]

# recalculate kME table, since we have corrected color assignments
  MEs<-tmpMEs<-data.frame()
  MEList = moduleEigengenes(t(cleanDat), colors = net$colors, verbose=0)
  MEs = orderMEs(MEList$eigengenes)
  net$MEs <- MEs
  colnames(MEs)<-gsub("ME","",colnames(MEs)) #let's be consistent in case prefix was added, remove it.
  rownames(MEs)<-rownames(numericMeta)

  tmpMEs <- MEs #net$MEs
  colnames(tmpMEs) <- paste("ME",colnames(MEs),sep="")
  MEs[,"grey"] <- NULL
  tmpMEs[,"MEgrey"] <- NULL

  kMEdat <- signedKME(t(cleanDat), tmpMEs, corFnc="bicor")

# recheck min kMEintramodule and max diff from kMEmax
  nModules<-length(table(net$colors))-1
  modules<-cbind(colnames(as.matrix(table(net$colors))),table(net$colors))
  orderedModules<-cbind(Mnum=paste("M",seq(1:nModules),sep=""),Color=labels2colors(c(1:nModules)))
  orderedModulesWithGrey=rbind(c("M0","grey"),orderedModules)
  kMEsIntramoduleVector<-apply( as.data.frame(cbind(net$colors,kMEdat)),1,function(x) if(!x[1]=="grey") { paste0(round(as.numeric(x[which(colnames(kMEdat)==paste0("kME",x[1]))+1]),4)) } else { 1 } ) #grey proteins set to dummy value of 1 (ignore)

  kMEmaxVec<-apply( as.data.frame(kMEdat),1,function(x) max(x) )
  kMEintramoduleVector<-unlist(lapply(kMEintramoduleVector,function(x) if(length(x)==0) { 1 } else { x }))   #grey will be ignored in checking for kMEmaxVec-kMEintramoduleVector difference max
  kMEmaxDiffCalc<- kMEmaxVec-kMEintramoduleVector
  if (min(kMEsIntramoduleVector)>=greyIfLT & max(kMEmaxDiffCalc)<=kMEmaxDiff) { cat(paste0("\nkME table 'clean' in ",iter," iterations.")); retry=FALSE; }
  iter=iter+1
  if (iter>30) break; #**
}
#** breaks after iteration 30 if did not reach criteria.


nModules<-length(table(net$colors))-1
modules<-cbind(colnames(as.matrix(table(net$colors))),table(net$colors))
orderedModules<-cbind(Mnum=paste("M",seq(1:nModules),sep=""),Color=labels2colors(c(1:nModules)))
modules<-modules[match(as.character(orderedModules[,2]),rownames(modules)),]
as.data.frame(cbind(orderedModules,Size=modules))


# Final modules after 1 cleanup iteration

#NULISA10 460 sample cleanDat, no Z, pwr=7 ds=4
          Mnum     Color Size
turquoise   M1 turquoise   47
blue        M2      blue   40
brown       M3     brown    9


## saved image of R session after running and finalizing blockwiseModules() function WGCNA output (now includes net data structure)
save.image(paste0("4.saved.image.",projectFilesOutputTag,".Rdata"))  #overwrites



## Create additional numeric traits and cull meaningless ones
#  numericMeta$Group.sex<-paste0(numericMeta$Group, ".", gsub("0","Female",gsub("1","Male",numericMeta$Sex)))
#  
#  numericMeta$Age.CT<-numericMeta$Age
#  numericMeta$Age.CT[numericMeta$Group=="AD"]<-NA
#  numericMeta$Age.AD<-numericMeta$Age
#  numericMeta$Age.AD[numericMeta$Group=="Control"]<-NA
#  
#  numericMeta$Abeta.Tau.Ratio <- numericMeta$Abeta / numericMeta$Tau
#  
#  numericMeta<-numericMeta[,which(!colnames(numericMeta)=="MoCA.CollectDiscrep.Yr")]
#  numericMeta<-numericMeta[,which(!colnames(numericMeta)=="Study.ID")]


## Output GlobalNetworkPlots and kMEtable
####################################################################################################################
FileBaseName=paste0(projectFilesOutputTag,"_pwr",power,"_ds4")


library(Cairo)
CairoPDF(file=paste0(outputfigs,"4.GlobalNetworkPlots-",FileBaseName,".pdf"),width=16,height=20)

## Plot dendrogram with module colors and trait correlations
MEs<-tmpMEs<-data.frame()
MEList = moduleEigengenes(t(cleanDat), colors = net$colors)
MEs = orderMEs(MEList$eigengenes)
colnames(MEs)<-gsub("ME","",colnames(MEs)) #let's be consistent in case prefix was added, remove it.
rownames(MEs)<-rownames(numericMeta)

numericIndices<-unique(c( which(!is.na(apply(numericMeta,2,function(x) sum(as.numeric(x))))), which(!(apply(numericMeta,2,function(x) sum(as.numeric(x),na.rm=T)))==0) ))
#Warnings OK; This determines which traits are numeric and if forced to numeric values, non-NA values do not sum to 0

geneSignificance <- cor(sapply(numericMeta[,numericIndices],as.numeric),t(cleanDat),use="pairwise.complete.obs")
rownames(geneSignificance) <- colnames(numericMeta)[numericIndices]
geneSigColors <- t(numbers2colors(t(geneSignificance),,signed=TRUE,lim=c(-1,1),naColor="black"))
rownames(geneSigColors) <- colnames(numericMeta)[numericIndices]

plotDendroAndColors(dendro=net$dendrograms[[1]],
                    colors=t(rbind(net$colors,geneSigColors)),
                    cex.dendroLabels=1.2,addGuide=TRUE,
                    dendroLabels=FALSE,
                    groupLabels=c("Module Colors",colnames(numericMeta)[numericIndices]))

## Plot eigengene dendrogram/heatmap - using bicor
tmpMEs <- MEs #net$MEs
colnames(tmpMEs) <- paste("ME",colnames(MEs),sep="")
MEs[,"grey"] <- NULL
tmpMEs[,"MEgrey"] <- NULL

plotEigengeneNetworks(tmpMEs, "Eigengene Network", marHeatmap = c(3,4,2,2), marDendro = c(0,4,2,0),plotDendrograms = TRUE, xLabelsAngle = 90,heatmapColors=blueWhiteRed(50))


######################
## Find differences between Groups (as defined in Traits input); Finalize Grouping of Samples for ANOVA

#Set a vector of strings that represent each sample in order, calling out each sample as a member of named groups (used by GlobalNetworkPlot boxplots, and later, ANOVA DiffEx)
Grouping<-numericMeta$Group  #here, we will calculate the T-test equivalent P value for C9 expanded vs C9 2repeat animals, controlling age and sex; typically there is a column "Group" loaded as a column in the traits.csv file
#Grouping[numericMeta$Group==0]<-"Control"  #only necessary if Group was numerically encoded; does nothing if the Grouping vector has no numeric values
#Grouping[numericMeta$Group==1]<-"AD"


# This gets one-way ANOVA (if ranked, Kruskal-Wallis) nonparametric p-values for groupwise comparison of interest.
# look at numericMeta (traits data) and choose traits to use for linear model-determination of p value
head(numericMeta)


# Change below line to point to a factored and/or numeric trait(s), Factored trait Group will define groups for ANOVA determination of significance of overall difference (overall if >2 groups)
regvars <- data.frame(as.factor( numericMeta$Group ), as.numeric(numericMeta$Age), as.factor(numericMeta$Sex))
colnames(regvars) <- c("Group","Age","Sex") ## data frame with covaraites incase we want to try multivariate regression


# P Values for Groupwise Difference within each ME by ANOVA, type III ( AGE+SEX ADJUSTMENT )
pvec.ageSexAdj <- rep(NA,ncol(MEs))
for (i in 1:ncol(MEs)) {
  ##aov1 <- aov(data.matrix(MEs)~Group,data=regvars) ## ANOVA framework yields same results
  lm1 <- lm(MEs[,i]~Group+Age+Sex,data=regvars) #Sex effects are removed by the linear model
  #  f <- summary(lm1)[[i]]$fstatistic ## Get F statistics
  #  pvec[i] <- pf(f[1],f[2],f[3],lower.tail=F) ## Get the p-value corresponding to the whole model
  a <- car::Anova(lm1, type = 3)   # Type III tests
  pvec.ageSexAdj[i] <- a["Group", "Pr(>F)"] ## Get the p-value corresponding to Group after adjustment for other covariates (and if type=3 ANOVA, for interactions)
}
names(pvec.ageSexAdj) <- colnames(MEs)


# P Values for Groupwise Difference within each ME by ANOVA, type III ( SEX ADJUSTMENT )
pvec.sexAdj <- rep(NA,ncol(MEs))
for (i in 1:ncol(MEs)) {
  ##aov1 <- aov(data.matrix(MEs)~Group,data=regvars) ## ANOVA framework yields same results
  lm1 <- lm(MEs[,i]~Group+Sex,data=regvars) #Sex effects are removed by the linear model
  #  f <- summary(lm1)[[i]]$fstatistic ## Get F statistics
  #  pvec[i] <- pf(f[1],f[2],f[3],lower.tail=F) ## Get the p-value corresponding to the whole model
  a <- car::Anova(lm1, type = 3)   # Type III tests
  pvec.sexAdj[i] <- a["Group", "Pr(>F)"] ## Get the p-value corresponding to Group after adjustment for other covariates (and if type=3 ANOVA, for interactions)
}
names(pvec.sexAdj) <- colnames(MEs)


# Change below line to point to a factored trait, which will define groups for ANOVA ( NO ADJUSTMENT FOR COVARIATES )
regvars <- data.frame(as.factor( numericMeta$Group ), as.numeric(numericMeta$Age), as.factor(numericMeta$Sex))
colnames(regvars) <- c("Group","Age","Sex") ## data frame with covaraites incase we want to try multivariate regression

pvec <- rep(NA,ncol(MEs))
for (i in 1:ncol(MEs)) {
  ##aov1 <- aov(data.matrix(MEs)~Group,data=regvars) ## ANOVA framework yields same results
  lm1 <- lm(MEs[,i]~Group,data=regvars) #Group is only variable in the model
  #  f <- summary(lm1)[[i]]$fstatistic ## Get F statistics
  #  pvec[i] <- pf(f[1],f[2],f[3],lower.tail=F) ## Get the p-value corresponding to the whole model
  a <- car::Anova(lm1, type = 3)   # Type III tests
  pvec[i] <- a["Group", "Pr(>F)"] ## Get the p-value corresponding to Group after adjustment for other covariates (and if type=3 ANOVA, for interactions)
}
names(pvec) <- colnames(MEs)



######################
## Get sigend kME values
kMEdat <- signedKME(t(cleanDat), tmpMEs, corFnc="bicor")


######################
## Plot eigengene-trait correlations - using p value of bicor for heatmap scale
library(RColorBrewer)
MEcors <- bicorAndPvalue(MEs,numericMeta[,numericIndices])
moduleTraitCor <- MEcors$bicor
moduleTraitPvalue <- MEcors$p


textMatrix = apply(moduleTraitCor,2,function(x) signif(x, 2))
#textMatrix = paste(signif(moduleTraitCor, 2), " (",
#  signif(moduleTraitPvalue, 1), ")", sep = "");
#dim(textMatrix) = dim(moduleTraitCor)
par(mfrow=c(2,1))
#par(mar = c(6, 8.5, 3, 3));
par(mar=c(6, 8.5, 3, 3) )

## Display the correlation values within a heatmap plot
cexx <- if(nModules>75) { 0.8 } else { 1 }
xlabAngle <- if(nModules>75) { 90 } else { 45 }

colvec <- rep("white",1500)
colvec[1:500] <- colorRampPalette(rev(brewer.pal(8,"BuPu")[2:8]))(500)
colvec[501:1000]<-colorRampPalette(c("white",brewer.pal(8,"BuPu")[2]))(3)[2] #interpolated color for 0.05-0.1 p
labeledHeatmap(Matrix = t(apply(moduleTraitPvalue,2,as.numeric)),
               yLabels = colnames(numericMeta)[numericIndices],
               xLabels = paste0("ME",names(MEs)),
               xSymbols = names(MEs),
               colorLabels = FALSE,
               colors = colvec,
               textMatrix = t(textMatrix),
               setStdMargins = FALSE,
               cex.text = 0.5,
               cex.lab.x= cexx,
               xLabelsAngle = xlabAngle,
               verticalSeparator.x=c(rep(c(1:length(colnames(MEs))),as.numeric(ncol(MEs)))),
               verticalSeparator.col = 1,
               verticalSeparator.lty = 1,
               verticalSeparator.lwd = 1,
               verticalSeparator.ext = 0,
               horizontalSeparator.y=c(rep(c(1:ncol(numericMetaCustom)),ncol(numericMetaCustom))),
               horizontalSeparator.col = 1,
               horizontalSeparator.lty = 1,
               horizontalSeparator.lwd = 1,
               horizontalSeparator.ext = 0,
               zlim = c(0,0.15),
               main = paste("Module-trait relationships\n bicor r-value shown as text\nHeatmap scale: Student correlation p value"),
               cex.main=0.8)


######################
## Plot eigengene-trait heatmap custom - using bicor color scale

numericMetaCustom<-numericMeta[,numericIndices]
MEcors <- bicorAndPvalue(MEs,numericMetaCustom)
moduleTraitCor <- MEcors$bicor
moduleTraitPvalue <- MEcors$p

moduleTraitPvalue.txt<-signif(moduleTraitPvalue, 1)
moduleTraitPvalue.txt[moduleTraitPvalue.txt > as.numeric(0.05)]<-as.character("")

textMatrix = moduleTraitPvalue.txt; #paste(signif(moduleTraitCor, 2), " / (", moduleTraitPvalue, ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
                                #textMatrix = gsub("()", "", textMatrix,fixed=TRUE)

labelMat<-matrix(nrow=(length(names(MEs))), ncol=2,data=c(rep(1:(length(names(MEs)))),labels2colors(1:(length(names(MEs))))))
labelMat<-labelMat[match(names(MEs),labelMat[,2]),]
for (i in 1:(length(names(MEs)))) { labelMat[i,1]<-paste("M",labelMat[i,1],sep="") }
for (i in 1:length(names(MEs))) { labelMat[i,2]<-paste("ME",labelMat[i,2],sep="") }

#rowMin(moduleTraitPvalue) # if we want to resort rows by min P value in the row

#par(mar=c(16, 12, 3, 3) )
#par(mfrow=c(1,1))

bw<-colorRampPalette(c("#0058CC", "white"))
wr<-colorRampPalette(c("white", "#CC3300"))

colvec<-c(bw(50),wr(50))

labeledHeatmap(Matrix = t(moduleTraitCor)[,],
               yLabels = colnames(numericMetaCustom),
               xLabels = labelMat[,2],
               xSymbols = labelMat[,1],
               xColorLabels=TRUE,
               colors = colvec,
               textMatrix = t(textMatrix)[,],
               setStdMargins = FALSE,
               cex.text = 0.5,
               cex.lab.x = cexx,
               xLabelsAngle = xlabAngle,
               verticalSeparator.x=c(rep(c(1:length(colnames(MEs))),as.numeric(ncol(MEs)))),
               verticalSeparator.col = 1,
               verticalSeparator.lty = 1,
               verticalSeparator.lwd = 1,
               verticalSeparator.ext = 0,
               horizontalSeparator.y=c(rep(c(1:ncol(numericMetaCustom)),ncol(numericMetaCustom))),
               horizontalSeparator.col = 1,
               horizontalSeparator.lty = 1,
               horizontalSeparator.lwd = 1,
               horizontalSeparator.ext = 0,
               zlim = c(-1,1),
               main = "Module-trait Relationships\n Heatmap scale: signed bicor r-value", # \n (Signif. p-values shown as text)"),
               cex.main=0.8)


################
## Plot annotated heatmap - annotate all the metadata, plot the eigengenes!
# This is where we will first use the Grouping vector of string group descriptions we set above.
toplot <- MEs

colnames(toplot) <- colnames(MEs)
rownames(toplot) <- rownames(MEs)
toplot <- t(toplot)

toplot <- t(apply(toplot,1,function(x) (x-mean(x))/sd(x) ))  # scale rows so we can windsorize on Z
# Windsorize toplot
minMax=min(abs(min(toplot)),max(toplot))
toplot[toplot>minMax] <- minMax
toplot[toplot< -minMax] <- -minMax


pvec.heatmap <- pvec.ageSexAdj[match(names(pvec.ageSexAdj),rownames(toplot))]
#rownames(toplot) <- paste(rownames(toplot),"\np = ",signif(pvec,2),sep="")
rownames(toplot) <- paste(orderedModules[match(colnames(MEs),orderedModules[,2]),1]," ",rownames(toplot),"  |  AOV p=",signif(pvec.ageSexAdj,2),sep="")

# add any traits of interest you want to be in the legend
#Gender=as.numeric(numericMeta$Sex)
#Gender[Gender==0]<-"Female"
#Gender[Gender==1]<-"Male"
Gender=numericMeta$Sex
metdat=data.frame(Group=Grouping,Age=as.numeric(numericMeta$Age), Gender=Gender)

# set colors for the traits in the legend
heatmapLegendColors=list('Group'=c('goldenrod','darkturquoise',"white", "seagreen3","hotpink","purple", "darkorange", "red", "gold", "blue", "bisque4"),  #CACNA1A   CDKL5 Control HNRNPH2  KANSL1   KCNQ2  MED13L  SLC2A1  STXBP1 SYNGAP1    TCF4
                         'Age'=c("white","darkgreen"), #young to old
                         'Gender'=c("pink","dodgerblue"), #F, M
                         'Modules'=sort(colnames(MEs)))

# Count group sizes - for Groupwise ordering of samples
tab <- table(Grouping)
# Remove Control from the size ranking
non_control <- setdiff(names(tab), "Control")
# Order non-Control groups by decreasing size
ordered_groups <- c("Control", non_control[order(tab[non_control], decreasing = TRUE)] )
# Build the index of positions in the desired order
idx <- unlist(lapply(ordered_groups, function(g) which(Grouping == g)))

library(NMF)
#par(mfrow=c(2,1))
par(mar=c(3, 3, 3, 3) )
layout(matrix(c(1,1, 2,3), nrow = 2, ncol = 2, byrow=TRUE),
       heights = c(0.95,1.3), # Heights of the rows
       widths = c(0.88,0.12)) # Widths of the columns  -- the distance to squash the second plot to the left, (because we do not duplicate legends)

aheatmap(x=toplot[,idx], ## Numeric Matrix
         main="Plot of Eigengene-Trait Relationships - SAMPLES ORDERED BY GROUP/GENOTYPE",
         annCol=metdat[idx,],
         annRow=data.frame(Modules=colnames(MEs)),
         annColors=heatmapLegendColors,
         border=list(matrix = TRUE),
         scale="row",
         distfun="correlation",hclustfun="average", ## Clustering options
         cexRow=0.8, ## Character sizes
         cexCol=0.8,
         col=blueWhiteRed(100), ## Color map scheme
         treeheight=80,
         Rowv=TRUE, Colv=NA) ## Do not cluster columns - keep given order

aheatmap(x=toplot, ## Numeric Matrix
         main="Plot of Eigengene-Trait Relationships - SAMPLES CLUSTERED",
         annCol=metdat,
         annRow=data.frame(Modules=colnames(MEs)),
         annColors=heatmapLegendColors,
         annLegend=FALSE,
         legend=FALSE,
         border=list(matrix = TRUE),
         scale="row",
         distfun="correlation",hclustfun="average", ## Clustering options
         cexRow=0.8, ## Character sizes
         cexCol=0.8,
         col=blueWhiteRed(100), ## Color map scheme
         treeheight=80,
         Rowv=TRUE,Colv=TRUE) ## Cluster columns





######################################
## Change the below code in the for loop using the following session output

#These are your numerically coded traits:
colnames(numericMeta)[numericIndices] #choose traits for correlation scatterplots (verboseScatterplot functions below)

#These are your ANOVA sample groups and the number of samples in each
table(Grouping) #alphabetically ordered, you choose the order of groups in the boxplot function by typing them in

## Make changes after checking output on console for the above 2 lines
#CairoPDF(file=paste0(outputfigs,"/GlobalNetPlots(BoxPlots)_",FileBaseName,"-CAIRO.pdf"),width=18,height=11.25)
##pdf(file=paste0(outputfigs,"/GlobalNetPlots(BoxPlots)_",FileBaseName,".pdf"),width=18,height=11.25)

#par(mfrow=c(4,6))
#par(mar=c(6.5,6,4.5,1.5))

layout(mat = matrix(c(1,7,13, 2,8,14, 3,9,15, 4,10,16, 5,11,17, 6,12,18), nrow = 6, ncol = 3, byrow=TRUE),
       heights = c(1.4,0.9,0.9,0.9,0.9,0.9), # Heights of the six rows
       widths = c(0.9,0.9,0.9)) # Widths of the 3 columns
#sapply(c(1:20),layout.show)

library(beeswarm)
library(gplots)
groupNames.ordered=c("Control","STXBP1", "SLC2A1", "SYNGAP1", "KANSL1", "HNRNPH2", "CDKL5", "KCNQ2", "MED13L", "CACNA1A", "TCF4")
genotypeColors=c("white","gold","red","blue","hotpink","seagreen3", "darkturquoise","purple","darkorange","goldenrod","bisque4")
names(genotypeColors)=groupNames.ordered
numericMeta$GroupColor<-NA
for (group in groupNames.ordered) numericMeta$GroupColor[numericMeta$Group==group]<-genotypeColors[group]
transcolGeno=paste0(col2hex(genotypeColors),"99")
names(transcolGeno)=groupNames.ordered

for (i in 1:(nrow(toplot))) {  # grey already excluded, no -1
  thisModCol=colnames(MEs)[i]
  transcol=paste0(col2hex(thisModCol),"99")

  par(mar=c(6.5,6,4.5,1.5))
  titlecolor<-if(min(signif(pvec,2)[i],signif(pvec.ageSexAdj,2)[i]) <0.05) { "red" } else { "black" }
#  boxplot(toplot[i,]~factor(Grouping,c("Control","STXBP1", "SLC2A1", "SYNGAP1", "KANSL1", "HNRNPH2", "CDKL5", "KCNQ2", "MED13L", "CACNA1A", "TCF4")),col=colnames(MEs)[i],ylab="Eigenprotein Value",main=paste0(orderedModules[match(colnames(MEs)[i],orderedModules[,2]),1]," ",colnames(MEs)[i],"\nAge+Sex-disc. p = ",signif(pvec,2)[i],"\nSex-discounted p = ",signif(pvec.noAgeRegr,2)[i]),xlab=NULL,col.main=titlecolor)  #rotate x labs: ,las=2  #no outliers: ,outline=FALSE)
  boxplot(toplot[i,]~factor(Grouping,groupNames.ordered),col=genotypeColors,ylab="Eigenprotein Value",main=paste0(orderedModules[match(colnames(MEs)[i],orderedModules[,2]),1]," ",colnames(MEs)[i],"\nAge+Sex-adj. p = ",signif(pvec.ageSexAdj,2)[i],"\nUnadj p = ",signif(pvec,2)[i]),xlab=NULL,col.main=titlecolor, border = thisModCol, whiskcol = thisModCol, staplecol = thisModCol, medcol = thisModCol, boxlwd=3, las=2)  #rotate x labs: ,las=2  #no outliers: ,outline=FALSE)
  beeswarm(toplot[i,]~factor(Grouping,groupNames.ordered),method="swarm",add=TRUE,corralWidth=0.5,vertical=TRUE,pch=21,bg=transcolGeno,col="black",cex=0.8,corral="gutter") #more like prism

#  titlecolor<-if(signif(pvec.group.sex,2)[i] <0.05) { "red" } else { "black" }
#  boxplot(toplot[i,]~factor(numericMeta$Group.sex,c("Control.Male","Control.Female","AD.Male","AD.Female")),col=colnames(MEs)[i],ylab="Eigenprotein Value",main=paste0(orderedModules[match(colnames(MEs)[i],orderedModules[,2]),1]," ",colnames(MEs)[i],"\nANOVA p = ",signif(pvec.group.sex,2)[i],"\n"),xlab=NULL,las=2,col.main=titlecolor)  #no outliers: ,outline=FALSE)
#  transcol=paste0(col2hex(colnames(MEs)[i]),"99")
#  beeswarm(toplot[i,]~factor(numericMeta$Group.sex,c("Control.Male","Control.Female","AD.Male","AD.Female")),method="swarm",add=TRUE,corralWidth=0.5,vertical=TRUE,pch=21,bg=transcol,col="black",cex=0.8,corral="gutter") #more like prism

  par(mar=c(4.5,6,4.5,1.5))
  verboseScatterplot(x=numericMeta[,"APOE4.01"],y=toplot[i,],xlab="ApoE e4 Carrier Status",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=numericMeta$GroupColor,pch=16,main=paste0("bicor=",signif(moduleTraitCor[i,"APOE4.01"],2),", p=",signif(moduleTraitPvalue[i,"APOE4.01"],2),"\n"),col.main=if(moduleTraitPvalue[i,"APOE4.01"]<0.05) { "red" } else { "black" })
  verboseScatterplot(x=numericMeta[,"Epilepsy.01"],y=toplot[i,],xlab="Epilepsy (1) or None (0)",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=numericMeta$GroupColor,pch=16,main=paste0("bicor=",signif(moduleTraitCor[i,"Epilepsy.01"],2),", p=",signif(moduleTraitPvalue[i,"Epilepsy.01"],2),"\n"),col.main=if(moduleTraitPvalue[i,"Epilepsy.01"]<0.05) { "red" } else { "black" })
  verboseScatterplot(x=numericMeta[,"Oligo-SNCA"],y=toplot[i,],xlab="Oligo alpha-Synuclein (NULISA)",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=numericMeta$GroupColor,pch=16,main=paste0("bicor=",signif(moduleTraitCor[i,"Oligo-SNCA"],2),", p=",signif(moduleTraitPvalue[i,"Oligo-SNCA"],2),"\n"),col.main=if(moduleTraitPvalue[i,"Oligo-SNCA"]<0.05) { "red" } else { "black" })
  verboseScatterplot(x=numericMeta[,"NEFL"],y=toplot[i,],xlab="NEFL (NULISA)",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=numericMeta$GroupColor,pch=16,main=paste0("bicor=",signif(moduleTraitCor[i,"NEFL"],2),", p=",signif(moduleTraitPvalue[i,"NEFL"],2),"\n"),col.main=if(moduleTraitPvalue[i,"NEFL"]<0.05) { "red" } else { "black" })
  verboseScatterplot(x=numericMeta[,"HBA1"],y=toplot[i,],xlab="HBA1 (NULISA)",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=numericMeta$GroupColor,pch=16,main=paste0("bicor=",signif(moduleTraitCor[i,"HBA1"],2),", p=",signif(moduleTraitPvalue[i,"HBA1"],2),"\n"),col.main=if(moduleTraitPvalue[i,"HBA1"]<0.05) { "red" } else { "black" })

#  verboseScatterplot(x=numericMeta[,"Sex"],y=toplot[i,],xlab="Sex (1=male)",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"Sex"],2),", p=",signif(moduleTraitPvalue[i,"Sex"],2),"\n"),col.main=if(moduleTraitPvalue[i,"Sex"]<0.05) { "red" } else { "black" })
#  verboseScatterplot(x=numericMeta[,"Age.CT"],y=toplot[i,],xlab="Age of Control Indiv.",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"Age.CT"],2),", p=",signif(moduleTraitPvalue[i,"Age.CT"],2),"\n"),col.main=if(moduleTraitPvalue[i,"Age.CT"]<0.05) { "red" } else { "black" })
#  ...
}


##outputs sample-by-sample eigenprotein barplots (not useful for large number of samples)
#while(!par('page')) plot.new()
#for (i in 1:nrow(toplot)) {
# barplot(height=rev(toplot[i,]),width=5,col=colnames(MEs)[i],xlab=paste(colnames(MEs)[i]," Eigenprotein Relative Expression"),main=rownames(toplot)[i],ylab=NULL,las=2,space=0.4,horiz=TRUE) #las=2 for rotated 90° X-axis labels  main=rownames(toplot)[i]
## text(bargr,par("usr")[3] - 0.025, srt=45, adj =1, labels= c(colnames(toplot)),xpd=TRUE,font=2) # bargr <- barplot(... above; gives rotated 45° x-axis labels but overwrites on top of existing ones
#}

dev.off() #finishes and closes writing of globalNetworkPlots PDF


########################################
#Write Module Membership/kME table
orderedModulesWithGrey=rbind(c("M0","grey"),orderedModules)
kMEtableSortVector<-apply( as.data.frame(cbind(net$colors,kMEdat)),1,function(x) if(!x[1]=="grey") { paste0(paste(orderedModulesWithGrey[match(x[1],orderedModulesWithGrey[,2]),],collapse=" "),"|",round(as.numeric(x[which(colnames(kMEdat)==paste0("kME",x[1]))+1]),4)) } else { paste0("grey|AllKmeAvg:",round(mean(as.numeric(x[-1],na.rm=TRUE)),4)) } ) 
kMEtable=cbind(c(1:nrow(cleanDat)),rownames(cleanDat),net$colors,kMEdat,kMEtableSortVector)[order(kMEtableSortVector,decreasing=TRUE),]
write.table(kMEtable,file=paste0(outputtabs,"/4.ProteinAssay(130)-ModuleAssignments-",FileBaseName,".txt"),sep="\t",row.names=FALSE)
#(load above file in excel and apply green-yellow-red conditional formatting heatmap to the columns with kME values); then save as excel.


## saved image of R session
save.image(paste0("4.saved.image.",projectFilesOutputTag,".Rdata"))  #overwrites
#load("4.saved.image.NULISA10.plasma-eigenproteinWGCNA.Rdata")

#cleanDat.protein<-cleanDat
MEs.protein<-MEs
kMEtable.protein<-kMEtable
net.protein<-net


## 4b. Run WGCNA on proteins as samples, to get sample clusters (n=465)
dim(cleanDat)
# 130 460

numericMeta.5OLremoved<-numericMeta
cleanDat.5OLremoved<-cleanDat  # cleanDat.protein


numericMeta<-numericMeta.full
cleanDat<-cleanDat.full
dim(cleanDat)
# 130 465


#transpose cleanDat:
cleanDat<-t(cleanDat)



## WGCNA blockwiseModules (for signed bicor coexpression network)
##############################################################################

##Above code may have borrowed server for parallel processing of bootstrap regression -- return parallel processing to local workstation
#library("doParallel")
#stopCluster(clusterLocal)
#parallelThreads=30 #set to # of threads on your computer
#clusterLocal <- makeCluster(c(rep("localhost",parallelThreads)),type="SOCK")
#registerDoParallel(clusterLocal)
enableWGCNAThreads() #speeds the pickSoftThreshold function

powers <- seq(2,28,by=1)  #initial power check -- try to get SFT.R.sq to go > 0.80
sft <- pickSoftThreshold(t(cleanDat),blockSize=nrow(cleanDat)+1000,   #always calculate power within a single block (blockSize > # of rows in cleanDat)
                         powerVector=powers,
                         corFnc="bicor",networkType="signed")

##Paste/replace output below and annotate based on graphical plot that follows
# NULISA10 Plasma 460 samples, protein WGCNA power selection:
   Power SFT.R.sq   slope truncated.R.sq mean.k. median.k. max.k.
1      2 0.884000  9.0100          0.926  313.00    315.00  350.0
2      3 0.873000  6.1400          0.933  260.00    263.00  306.0
3      4 0.851000  4.5400          0.928  217.00    220.00  268.0  << POWER = 4 (but trend inverse to the expected)
4      5 0.807000  3.5000          0.890  183.00    185.00  235.0
5      6 0.804000  2.9200          0.914  155.00    156.00  207.0
6      7 0.762000  2.4100          0.865  132.00    133.00  183.0
7      8 0.734000  2.0600          0.847  112.00    114.00  162.0
8      9 0.703000  1.7900          0.845   96.50     96.80  144.0  << POWER = 9 (k, connectivity <100)
9     10 0.614000  1.4600          0.821   83.20     82.80  128.0
10    11 0.560000  1.2400          0.792   71.90     71.90  115.0
11    12 0.525000  1.0900          0.803   62.30     62.40  103.0
12    13 0.405000  0.8710          0.782   54.20     54.30   92.9
13    14 0.306000  0.7280          0.697   47.30     47.30   83.8
14    15 0.235000  0.5830          0.708   41.40     41.10   75.7
15    16 0.158000  0.4490          0.698   36.30     35.80   68.5
16    17 0.092400  0.3400          0.714   31.90     31.70   62.1
17    18 0.055300  0.2530          0.763   28.00     27.70   56.3
18    19 0.015400  0.1290          0.767   24.70     24.30   51.2
19    20 0.000641  0.0265          0.775   21.90     21.30   46.6
20    21 0.001290 -0.0373          0.828   19.40     18.70   42.4
21    22 0.017900 -0.1430          0.854   17.20     16.50   38.7
22    23 0.029700 -0.1860          0.877   15.30     14.50   35.3
23    24 0.042100 -0.2200          0.904   13.60     12.90   32.3
24    25 0.107000 -0.3710          0.904   12.10     11.50   29.5
25    26 0.140000 -0.4280          0.910   10.80     10.20   27.0
26    27 0.185000 -0.5250          0.882    9.66      9.01   24.8
27    28 0.223000 -0.5910          0.888    8.64      8.00   22.7


#plot initial SFT.R.sq vs. power curve
tableSFT<-sft[[2]]
plot(tableSFT[,1],tableSFT[,2],xlab="Power (Beta)",ylab="SFT R²")


#choose power at elbow of SFT R² curve approaching asymptote near or ideally above 0.80
power=9


## Run an automated network analysis (ds=4 and mergeCutHeight=0.07, more liberal)
# choose parameters deepSplit and mergeCutHeight to get respectively more modules and more stringency sending more low connectivity genes to grey (not in modules).
net <- blockwiseModules(t(cleanDat),power=power,deepSplit=4,minModuleSize=6,
                        mergeCutHeight=0.007,TOMdenom="mean", #detectCutHeight=0.9999,                        #TOMdenom="mean" may get more small modules here.
                        corType="bicor",networkType="signed",pamStage=TRUE,pamRespectsDendro=TRUE,
                        verbose=3,saveTOMs=FALSE,maxBlockSize=nrow(cleanDat)+1000,reassignThresh=0.05)       #maxBlockSize always more than the number of rows in cleanDat
#blockwiseModules can take 30 min+ for large numbers of gene products/proteins (10000s of rows); much quicker for smaller proteomic data sets

nModules<-length(table(net$colors))- if("grey" %in% names(table(net$colors))) { 1 } else { 0 }
modules<-cbind(colnames(as.matrix(table(net$colors))),table(net$colors))
orderedModules<-cbind(Mnum=paste("M",seq(1:nModules),sep=""),Color=labels2colors(c(1:nModules)))
modules<-modules[match(as.character(orderedModules[,2]),rownames(modules)),]
as.data.frame(cbind(orderedModules,Size=modules))

##copy R session output;
## NULISA10 Plasma 465 samples as genes for subtyping, protein WGCNA power selection; deepSplit=4; power=4; mergeCutheight=0.007; MMS=6:
          Mnum     Color Size
turquoise   M1 turquoise  119
blue        M2      blue  111
brown       M3     brown   84
yellow      M4    yellow   78

## NULISA10 Plasma 460 samples, protein WGCNA power selection; deepSplit=4; power=2; otherwise same as above
          Mnum     Color Size
turquoise   M1 turquoise  167
blue        M2      blue  153

## NULISA10 Plasma 460 samples, protein WGCNA power selection; deepSplit=4; power=9; otherwise same as above
          Mnum     Color Size
turquoise   M1 turquoise  118
blue        M2      blue  117
brown       M3     brown   90
yellow      M4    yellow   79
green       M5     green   61

# (we will stick with the power=9 module definitions)
net.samples<-net


#we will explore the blockwiseModules() function-built network with parameter deepSplit=4 in the output for folder ii.Zscore_unionNetwork_ds4 (ds3 already explored)
net<-net.samples
outputfigs<-outputtabs<-rootdir

#<OK TO RUN/SKIP, no effect>
enforceMMS=TRUE
minModSize=10
# If necessary, return module members of small modules below size minSize=X to grey
if (enforceMMS) {
  removedModules<-orderedModules[which(modules<minModSize),"Color"]
  for(i in removedModules) { net$colors[net$colors==i] <- "grey" }
  for(i in removedModules) { net$MEs[,paste0("ME",i)] <- NULL }

  nModules<-length(table(net$colors))- if("grey" %in% names(table(net$colors))) { 1 } else { 0 }
  modules<-cbind(colnames(as.matrix(table(net$colors))),table(net$colors))
  orderedModules<-cbind(Mnum=paste("M",seq(1:nModules),sep=""),Color=labels2colors(c(1:nModules)))
  modules<-modules[match(as.character(orderedModules[,2]),rownames(modules)),]
  as.data.frame(cbind(orderedModules,Size=modules))
}
minModSize=10
#<END SKIP>


#calculate kME table up front, in case we need to correct color assignments
MEs<-tmpMEs<-data.frame()
MEList = moduleEigengenes(t(cleanDat), colors = net$colors)
MEs = orderMEs(MEList$eigengenes)
net$MEs <- MEs
colnames(MEs)<-gsub("ME","",colnames(MEs)) #let's be consistent in case prefix was added, remove it.
rownames(MEs)<-colnames(cleanDat) #rownames(numericMeta)

tmpMEs <- MEs #net$MEs
colnames(tmpMEs) <- paste("ME",colnames(MEs),sep="")
MEs[,"grey"] <- NULL
tmpMEs[,"MEgrey"] <- NULL

kMEdat <- signedKME(t(cleanDat), tmpMEs, corFnc="bicor")


table(net$colors)["grey"]  #0 grey samples



##ITERATIVE until condition met that all module membes are at least 0.28 kMEintramodule.
#Go back and do final algorithm fix of module colors (remove kMEintramodule<0.28 members, reassign grey with kMEintramodule>0.35; max difference from kMEmax<0.10)

retry=TRUE;
kMEmaxDiff=0.1
reassignIfGT=0.30
greyIfLT=0.30
iter=1;
while (retry) {
  cat(paste0("\nkME table Cleanup, processing iteration ",iter,"..."))
  colorVecFixed<-colorVecBackup<-net$colors
  orderedModulesWithGrey=rbind(c("M0","grey"),orderedModules)
  kMEintramoduleVector<-apply( as.data.frame(cbind(net$colors,kMEdat)),1,function(x) as.numeric(x[which(colnames(kMEdat)==paste0("kME",x[1]))+1]) )  #all sig digits (no rounding), so max will be unique.
  colorVecFixed[kMEintramoduleVector<greyIfLT]<-"grey"
  kMEmaxVec<-apply( as.data.frame(kMEdat),1,function(x) max(x) )
  kMEmaxColorsVec<-apply( as.data.frame(cbind(kMEmaxVec,kMEdat)),1, function(x) gsub("kME","",colnames(kMEdat)[which(x==x[1])[2]-1]) )
  kMEintramoduleVector<-unlist(lapply(kMEintramoduleVector,function(x) if(length(x)==0) { 1 } else { x }))   #grey will be ignored in checking for kMEmaxVec-kMEintramoduleVector difference max
  kMEmaxDiffTooBig<-(kMEmaxVec-kMEintramoduleVector) >= kMEmaxDiff
  colorVecFixed[which( (colorVecFixed=="grey" & kMEmaxVec>reassignIfGT) | kMEmaxDiffTooBig )] <- kMEmaxColorsVec[which( (colorVecFixed=="grey" & kMEmaxVec>reassignIfGT) | kMEmaxDiffTooBig )]
  net$colors<-colorVecFixed

#  table(net$colors)["grey"]  #decreased to x


# Are colors still in rank order? -- put them in order by recoloring modules that changed rank
  sort(table(net$colors),decreasing=TRUE)[!names(sort(table(net$colors),decreasing=TRUE))=="grey"]

  oldcolors <- names(sort(table(net$colors),decreasing=TRUE)[!names(sort(table(net$colors),decreasing=TRUE))=="grey"])
  for (i in 1:length(oldcolors)) {
    net$colors[net$colors==oldcolors[i]]<-paste0("proxy",labels2colors(i))
  }
  for (i in 1:length(oldcolors)) {
    net$colors[net$colors==paste0("proxy",labels2colors(i))]<-labels2colors(i)
  }

# one can check that colors are in order by size now
  #sort(table(net$colors),decreasing=TRUE)[!names(sort(table(net$colors),decreasing=TRUE))=="grey"]

# recalculate kME table, since we have corrected color assignments
  MEs<-tmpMEs<-data.frame()
  MEList = moduleEigengenes(t(cleanDat), colors = net$colors, verbose=0)
  MEs = orderMEs(MEList$eigengenes)
  net$MEs <- MEs
  colnames(MEs)<-gsub("ME","",colnames(MEs)) #let's be consistent in case prefix was added, remove it.
  rownames(MEs)<-colnames(cleanDat)  #usually rownames(numericMeta)

  tmpMEs <- MEs #net$MEs
  colnames(tmpMEs) <- paste("ME",colnames(MEs),sep="")
  MEs[,"grey"] <- NULL
  tmpMEs[,"MEgrey"] <- NULL

  kMEdat <- signedKME(t(cleanDat), tmpMEs, corFnc="bicor")

# recheck min kMEintramodule and max diff from kMEmax
  nModules<-length(table(net$colors))- if("grey" %in% names(table(net$colors))) { 1 } else { 0 }
  modules<-cbind(colnames(as.matrix(table(net$colors))),table(net$colors))
  orderedModules<-cbind(Mnum=paste("M",seq(1:nModules),sep=""),Color=labels2colors(c(1:nModules)))
  orderedModulesWithGrey=rbind(c("M0","grey"),orderedModules)
  kMEsIntramoduleVector<-apply( as.data.frame(cbind(net$colors,kMEdat)),1,function(x) if(!x[1]=="grey") { paste0(round(as.numeric(x[which(colnames(kMEdat)==paste0("kME",x[1]))+1]),4)) } else { 1 } ) #grey proteins set to dummy value of 1 (ignore)

  kMEmaxVec<-apply( as.data.frame(kMEdat),1,function(x) max(x) )
  kMEintramoduleVector<-unlist(lapply(kMEintramoduleVector,function(x) if(length(x)==0) { 1 } else { x }))   #grey will be ignored in checking for kMEmaxVec-kMEintramoduleVector difference max
  kMEmaxDiffCalc<- kMEmaxVec-kMEintramoduleVector
  if (min(kMEsIntramoduleVector)>=greyIfLT & max(kMEmaxDiffCalc)<=kMEmaxDiff) { cat(paste0("\nkME table 'clean' in ",iter," iterations.")); retry=FALSE; }
  iter=iter+1
  if (iter>30) break; #**
}
#** breaks after iteration 30 if did not reach criteria.


nModules<-length(table(net$colors))- if("grey" %in% names(table(net$colors))) { 1 } else { 0 }
modules<-cbind(colnames(as.matrix(table(net$colors))),table(net$colors))
orderedModules<-cbind(Mnum=paste("M",seq(1:nModules),sep=""),Color=labels2colors(c(1:nModules)))
modules<-modules[match(as.character(orderedModules[,2]),rownames(modules)),]
as.data.frame(cbind(orderedModules,Size=modules))


# Final modules after 1 cleanup iteration

#NULISA10 460 sample cleanDat, no Z, pwr=7 ds=4
          Mnum     Color Size
turquoise   M1 turquoise  118
blue        M2      blue  117
brown       M3     brown   90
yellow      M4    yellow   79
green       M5     green   61


numericMeta$sampleSub<-net$colors
numericMeta$subtype1<-numericMeta$subtype2<-numericMeta$subtype3<-numericMeta$subtype4<-numericMeta$subtype5 <- 0
numericMeta$subtype1[numericMeta$sampleSub=="turquoise"]<- 1
numericMeta$subtype2[numericMeta$sampleSub=="blue"]     <- 1
numericMeta$subtype3[numericMeta$sampleSub=="brown"]    <- 1
numericMeta$subtype4[numericMeta$sampleSub=="yellow"]   <- 1
numericMeta$subtype5[numericMeta$sampleSub=="green"]    <- 1

numericMeta$sampleSub.sex<-paste0(numericMeta$sampleSub,".",numericMeta$Sex)


numericMeta$GroupColor<-NA
for (group in groupNames.ordered) numericMeta$GroupColor[numericMeta$Group==group]<-genotypeColors[group]

sexColor<-numericMeta$Sex
sexColor[numericMeta$Sex=="Male"] <- "dodgerblue"
sexColor[numericMeta$Sex=="Female"]<-"pink"

APOE4color<-numericMeta$APOE4
APOE4color[numericMeta$APOE4=="No"]<-"white"
APOE4color[numericMeta$APOE4=="Yes"]<-"darkslateblue"

epilepsyColor<-numericMeta$Epilepsy
epilepsyColor[numericMeta$Epilepsy=="Yes"]<-"red"
epilepsyColor[numericMeta$Epilepsy=="No"]<-"white"
epilepsyColor[is.na(epilepsyColor)]<-"grey"

numericMeta.5OLremoved$sampleSub<-numericMeta$sampleSub[match(rownames(numericMeta.5OLremoved),rownames(numericMeta))]


## saved image of R session after running and finalizing blockwiseModules() function WGCNA output (now includes net data structure)
projectFilesOutputTag="NULISA10.plasma-eigensampleWGCNA"
save.image(paste0("4b.saved.image.",projectFilesOutputTag,".Rdata"))  #overwrites


## Output GlobalNetworkPlots (sample subtypes) and kMEtable
####################################################################################################################
FileBaseName=paste0(projectFilesOutputTag,"_pwr",power,"_ds4")


library(Cairo)
CairoPDF(file=paste0(outputfigs,"4b.GlobalNetworkPlots-",FileBaseName,".pdf"),width=16,height=20)

## Plot dendrogram with module colors and trait correlations
MEs<-tmpMEs<-data.frame()
MEList = moduleEigengenes(t(cleanDat), colors = net$colors)
MEs = orderMEs(MEList$eigengenes)
colnames(MEs)<-gsub("ME","",colnames(MEs)) #let's be consistent in case prefix was added, remove it.
rownames(MEs)<-colnames(cleanDat)  # usually rownames(numericMeta)

numericIndices<-unique(c( which(!is.na(apply(numericMeta,2,function(x) sum(as.numeric(x))))), which(!(apply(numericMeta,2,function(x) sum(as.numeric(x),na.rm=T)))==0) ))
##Warnings OK; This determines which traits are numeric and if forced to numeric values, non-NA values do not sum to 0

#geneSignificance <- cor(sapply(numericMeta[,numericIndices],as.numeric),t(cleanDat),use="pairwise.complete.obs")
#rownames(geneSignificance) <- colnames(numericMeta)[numericIndices]
#geneSigColors <- t(numbers2colors(t(geneSignificance),,signed=TRUE,lim=c(-1,1),naColor="black"))
#rownames(geneSigColors) <- colnames(numericMeta)[numericIndices]

par(oma = c(4, 2.5, 1, 1))   # add a large top outer margin
par(mar = c(2, 2, 2, 2))   # normal inner margins

## viewport: the dendrogram + color bars (bottom 85%)
#par(fig = c(0, 1, 0, 0.85), new = TRUE, mar = c(2, 2, 2, 2))

plotDendroAndColors(dendro=net$dendrograms[[1]],
                    colors=t(rbind(net$colors, sexColor, APOE4color, epilepsyColor, numericMeta$GroupColor )),  # ,geneSigColors)),
                    cex.dendroLabels=1.2,addGuide=TRUE,
                    dendroLabels=FALSE,
                    groupLabels=c("Subtype Colors", "Male (blue)\nFemale (pink)", "APOE4 Carrier?\n(dark blue)", "Epilepsy\n(red; NA, grey)", "Mutation Group" ))  #,colnames(numericMeta)[numericIndices]))

# 3. Now open a small viewport ABOVE the plot
par(xpd = NA)  # allow drawing outside plot region

# Define a viewport in outer margin
par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0))
plot.new()

legend(
  x="bottom",
  y="bottom",
  legend = groupNames.ordered,
  fill   = genotypeColors,
  horiz  = TRUE,
  bty    = "n",
  cex    = 0.78,
  xpd    = NA
)

par(xpd = FALSE)



## Plot eigengene dendrogram/heatmap - using bicor
tmpMEs <- MEs #net$MEs
colnames(tmpMEs) <- paste("ME",colnames(MEs),sep="")
MEs[,"grey"] <- NULL
tmpMEs[,"MEgrey"] <- NULL

plotEigengeneNetworks(tmpMEs, "Eigensample Network", marHeatmap = c(3,4,2,2), marDendro = c(0,4,2,0),plotDendrograms = TRUE, xLabelsAngle = 90,heatmapColors=blueWhiteRed(50))



## Find differences between Disease Status Groups (e.g. AD, Control)--here, subtypes
regvars <- data.frame(as.factor(numericMeta[,"sampleSub"]),as.numeric(numericMeta[,"Age"]),as.factor(numericMeta[,"Sex"]))  #,as.numeric(numericMeta.original[,"pmi"])
colnames(regvars) <- c("Group","Age","Sex") #,"PMI"  ## data frame with covaraites in case we want to try multivariate regression
#aov1 <- aov(data.matrix(MEs)~Group,data=regvars) ## ANOVA framework yields same results
lm1 <- lm(data.matrix(cleanDat)~Group,data=regvars)

pvec <- rep(NA,ncol(cleanDat))
for (i in 1:ncol(cleanDat)) {
  f <- summary(lm1)[[i]]$fstatistic ## Get F statistics
  pvec[i] <- pf(f[1],f[2],f[3],lower.tail=F) ## Get the p-value corresponding to the whole model
}
names(pvec) <- colnames(cleanDat)

## Find differences between APOE Risk Groups
regvars <- data.frame(as.factor(numericMeta[,"APOE4"]),as.numeric(numericMeta[,"Age"]),as.factor(numericMeta[,"Sex"]))  #,as.numeric(numericMeta.original[,"pmi"])
colnames(regvars) <- c("Group","Age","Sex") #,"PMI"  ## data frame with covaraites in case we want to try multivariate regression
#aov1 <- aov(data.matrix(MEs)~Group,data=regvars) ## ANOVA framework yields same results
lm1 <- lm(data.matrix(cleanDat)~Group,data=regvars)

pvec.Apoe <- rep(NA,ncol(cleanDat))
for (i in 1:ncol(cleanDat)) {
  f <- summary(lm1)[[i]]$fstatistic ## Get F statistics
  pvec.Apoe[i] <- pf(f[1],f[2],f[3],lower.tail=F) ## Get the p-value corresponding to the whole model
}
names(pvec.Apoe) <- colnames(cleanDat)


## Find differences between New subtypes.Sex Groups
regvars <- data.frame(as.factor(numericMeta[,"sampleSub.sex"]),as.numeric(numericMeta[,"Age"]),as.factor(numericMeta[,"Sex"]))  #,as.numeric(numericMeta.original[,"pmi"])
colnames(regvars) <- c("Group","Age","Sex") #,"PMI"  ## data frame with covaraites in case we want to try multivariate regression
#aov1 <- aov(data.matrix(MEs)~Group,data=regvars) ## ANOVA framework yields same results
lm1 <- lm(data.matrix(cleanDat)~Group,data=regvars)

pvec.subGroup.Sex <- rep(NA,ncol(cleanDat))
for (i in 1:ncol(cleanDat)) {
  f <- summary(lm1)[[i]]$fstatistic ## Get F statistics
  pvec.subGroup.Sex[i] <- pf(f[1],f[2],f[3],lower.tail=F) ## Get the p-value corresponding to the whole model
}
names(pvec.subGroup.Sex) <- colnames(cleanDat)


## Plot eigenprotein-trait correlations - using bicor

## Calculate MEs for all 465 samples
MEs.proteinAllSamp<-tmpMEs.proteinAllSamp<-data.frame()
MEList.proteinAllSamp = moduleEigengenes(cleanDat, colors = net.protein$colors)  #already transposed cleanDat
MEs.proteinAllSamp = orderMEs(MEList.proteinAllSamp$eigengenes)
colnames(MEs.proteinAllSamp)<-gsub("ME","",colnames(MEs.proteinAllSamp)) #let's be consistent in case prefix was added, remove it.
rownames(MEs.proteinAllSamp)<-colnames(t(cleanDat))  # usually rownames(numericMeta)

tmpMEs.proteinAllSamp <- MEs.proteinAllSamp #net$MEs
colnames(tmpMEs.proteinAllSamp) <- paste("ME",colnames(MEs.proteinAllSamp),sep="")
MEs.proteinAllSamp[,"grey"] <- NULL
tmpMEs.proteinAllSamp[,"MEgrey"] <- NULL



library(RColorBrewer)
MEcors <- bicorAndPvalue(MEs.proteinAllSamp,numericMeta[,numericIndices])
moduleTraitCor <- MEcors$bicor
moduleTraitPvalue <- MEcors$p


textMatrix = apply(moduleTraitCor,2,function(x) signif(x, 2))
#textMatrix = paste(signif(moduleTraitCor, 2), " (",
#  signif(moduleTraitPvalue, 1), ")", sep = "");
#dim(textMatrix) = dim(moduleTraitCor)
par(mfrow=c(1,1))
par(mar = c(6, 8.5, 3, 3));

## Display the correlation values within a heatmap plot
colvec <- rep("white",1500)
colvec[1:500] <- colorRampPalette(rev(brewer.pal(8,"BuPu")[2:8]))(500)
colvec[501:1000]<-colorRampPalette(c("white",brewer.pal(8,"BuPu")[2]))(3)[2] #interpolated color for 0.05-0.1 p
labeledHeatmap(Matrix = cbind( apply(moduleTraitPvalue,2,as.numeric) ), #, apply(moduleTraitPvalue,2,as.numeric)),  #2 copies of same column for 2 rows of heatmap.
               xLabels = c(colnames(numericMeta)[numericIndices] ),  #, "duplicated row"),
               yLabels = paste0("ME",names(MEs.proteinAllSamp)),
               ySymbols = names(MEs.proteinAllSamp),
               colorLabels = FALSE,
               colors = colvec,
               textMatrix = cbind(textMatrix ), #, textMatrix),
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(0,0.15),
               main = paste("Sample Net Module-trait relationships\n bicor r-value shown as text\nHeatmap scale: Student correlation p value"),
               cex.main=0.8)


######################
## Plot eigengene-trait heatmap custom - using bicor

numericMetaCustom<-numericMeta[,numericIndices]
if(is.null(dim(numericMetaCustom)[2])) { numericMetaCustom<-data.frame(Column1=numericMeta[,numericIndices]); colnames(numericMetaCustom) <- colnames(numericMeta[,numericIndices]); }
colnames(numericMetaCustom)[which(colnames(numericMetaCustom)=="Treatment.numericMeta.original$sampleSub")]<-"Treatment (P=0 M=0.5 F=3)"
MEcors <- bicorAndPvalue(MEs.proteinAllSamp,numericMetaCustom)
moduleTraitCor <- MEcors$bicor
moduleTraitPvalue <- MEcors$p

moduleTraitPvalue<-signif(moduleTraitPvalue, 1)
moduleTraitPvalue[moduleTraitPvalue > as.numeric(0.05)]<-as.character("")

textMatrix = moduleTraitPvalue; #paste(signif(moduleTraitCor, 2), " / (", moduleTraitPvalue, ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
                                #textMatrix = gsub("()", "", textMatrix,fixed=TRUE)

labelMat<-matrix(nrow=(length(names(MEs.proteinAllSamp))), ncol=2,data=c(rep(1:(length(names(MEs.proteinAllSamp)))),labels2colors(1:(length(names(MEs.proteinAllSamp))))))
labelMat<-labelMat[match(names(MEs.proteinAllSamp),labelMat[,2]),]
for (i in 1:(length(names(MEs.proteinAllSamp)))) { labelMat[i,1]<-paste("M",labelMat[i,1],sep="") }
for (i in 1:length(names(MEs.proteinAllSamp))) { labelMat[i,2]<-paste("ME",labelMat[i,2],sep="") }

#rowMin(moduleTraitPvalue) # if we want to resort rows by min P value in the row

par(mar=c(16, 12, 3, 3) )
par(mfrow=c(1,1))

bw<-colorRampPalette(c("#0058CC", "white"))
wr<-colorRampPalette(c("white", "#CC3300"))

colvec<-c(bw(50),wr(50))

labeledHeatmap(Matrix = t(cbind(moduleTraitCor )),  #, moduleTraitCor)),
               yLabels = c(colnames(numericMetaCustom) ), #, "duplicated column"),
               xLabels = labelMat[,2],
               xSymbols = labelMat[,1],
               xColorLabels=TRUE,
               colors = colvec,
               textMatrix = t(cbind(textMatrix )), #, textMatrix))[,],
               setStdMargins = FALSE,
               cex.text = 0.5,
               verticalSeparator.x=c(rep(c(1:length(colnames(MEs.proteinAllSamp))),as.numeric(ncol(MEs.proteinAllSamp)))),
               verticalSeparator.col = 1,
               verticalSeparator.lty = 1,
               verticalSeparator.lwd = 1,
               verticalSeparator.ext = 0,
               horizontalSeparator.y=c(rep(c(1:ncol(numericMetaCustom)),ncol(numericMetaCustom))),
               horizontalSeparator.col = 1,
               horizontalSeparator.lty = 1,
               horizontalSeparator.lwd = 1,
               horizontalSeparator.ext = 0,
               zlim = c(-1,1),
               main = "Sample Net Module-trait Relationships\n Heatmap scale: signed bicor r-value", # \n (Signif. p-values shown as text)"),
               cex.main=0.8)


dev.off()

#numericIndices.protein<-unique(c( which(!is.na(apply(numericMeta,2,function(x) sum(as.numeric(x))))), which(!(apply(numericMeta,2,function(x) sum(as.numeric(x),na.rm=T)))==0) ))





## Define the heatmap function to split on the 5 sample tWGCNA clusters for all 465 samples
makeSupClusterHeatmap.corrColClust.moreTraits.noRowClust.5sampClusterSPLIT.moreAnnColors <- function(cleanDat=cleanDat,traits=numericMeta,ANOVAout=ANOVAout,FDRcutoff=1,scaleData=FALSE,forceSplit=FALSE, env=.GlobalEnv) {
  	# Check data for clustering
  range(cleanDat,na.rm=T)
  	#[1] -8.654356  7.939041
  range(cleanDat[which(ANOVAout$'minP'<FDRcutoff),],na.rm=T)
  	#[1] -2.163940  7.939041
  	#will return -Inf to Inf range if 0 rows!
  
   # Get sample traits
  #numericMeta<-as.data.frame(readxl::read_excel("EmHeparin36_TRAITS.xlsx"))
  #colnames(numericMeta)[which(colnames(numericMeta)=="cha_N_Nel")]<-"Channel"
  #rownames(numericMeta)<-paste0("b0",numericMeta$Batch,".",numericMeta$Channel)
  #numericMeta<-numericMeta[match(colnames(cleanDat),rownames(numericMeta)),]

  numericMeta$Epilepsy[is.na(numericMeta$Epilepsy)]<-"unknown"
    	
  	# Sample Trait annotation colors for heatmap are from traits
  annotation_col = data.frame(
  	    GenoGroup = factor(numericMeta$Group), 
#  	    Race = factor(numericMeta$Race),
#  	    Abeta = as.numeric(numericMeta$Abeta),
#  	    tTau = as.numeric(numericMeta$Tau),
#  	    pTau = as.numeric(numericMeta$pTau),
#  	    MoCA = as.numeric(numericMeta$MoCA),
#  	    e4.Risk = as.numeric(numericMeta$APOE.Risk),
            age = as.numeric(numericMeta$Age),
  	    sex = as.factor(numericMeta$Sex),
  	    e4.carrier = as.factor(numericMeta$APOE4),
  	    epilepsy = as.factor(numericMeta$Epilepsy),
  	    sampleSubtype = factor(numericMeta$sampleSub) #,
#  	    sampleSub.ADonly = numericMeta$sampleSub.ADonly  #note - factor after NAs replaced with dummy value.
  	#    plasma.pTau181 = as.numeric(numericMeta$'Plasma pTau181')
  	    )
  rownames(annotation_col) = rownames(numericMeta)
#  annotation_col$sampleSub.ADonly[which(is.na(annotation_col$sampleSub.ADonly))]<- "grey30"  #cannot have NA in traits here
#  annotation_col$sampleSub.ADonly<-factor(annotation_col$sampleSub.ADonly)

  #ADDED TO HANDLE NA values for Epilepsy
  annotation_col[is.na(annotation_col)]<-0
  
  DEP.status.KANSL1vs.CT = rep(c("unchanged"), nrow(ANOVAout))
  DEP.status.KANSL1vs.CT[which(ANOVAout[,'diff KANSL1-Control']<0)]<-"upregulated KANSL1 vs CT"
  DEP.status.KANSL1vs.CT[which(ANOVAout[,'diff KANSL1-Control']>0)]<-"downregulated KANSL1 vs CT"
#  DEP.status1<-DEP.status.KANSL1vs.CT[which(ANOVAout$'FDR (BH)'<FDRcutoff)]

  DEP.status.SYNGAP1vs.CT = rep(c("unchanged"), nrow(ANOVAout))
  DEP.status.SYNGAP1vs.CT[which(ANOVAout[,'diff SYNGAP1-Control']<0)]<-"upregulated SYNGAP1 vs CT"
  DEP.status.SYNGAP1vs.CT[which(ANOVAout[,'diff SYNGAP1-Control']>0)]<-"downregulated SYNGAP1 vs CT"
#  DEP.status2<-DEP.status.SYNGAP1vs.CT[which(ANOVAout$'FDR (BH)'<FDRcutoff)]

  DEP.status.SLC2A1vs.CT = rep(c("unchanged"), nrow(ANOVAout))
  DEP.status.SLC2A1vs.CT[which(ANOVAout[,'diff SLC2A1-Control']<0)]<-"upregulated SLC2A1 vs CT"
  DEP.status.SLC2A1vs.CT[which(ANOVAout[,'diff SLC2A1-Control']>0)]<-"downregulated SLC2A1 vs CT"
#  DEP.status3<-DEP.status.SLC2A1vs.CT[which(ANOVAout$'FDR (BH)'<FDRcutoff)]

  annotation_row = data.frame(
  	 KANSL1vs.CT.trend=factor(DEP.status.KANSL1vs.CT),   # factor(DEP.status1),
  	 SYNGAP1vs.CT.trend=factor(DEP.status.SYNGAP1vs.CT),   # factor(DEP.status2),
  	 SLC2A1vs.CT.trend=factor(DEP.status.SLC2A1vs.CT),  # factor(DEP.status3),
  	 Module= ANOVAout$NETcolors  # factor(gsub("^ME","",rownames(ANOVAout)))
  )
  rownames(annotation_row) = rownames(ANOVAout) #[which(ANOVAout$'FDR (BH)'<FDRcutoff)]

  #NOTE -- make sure your named vector elements giving trait colors (both to columns and rows) match the names used above!	
  ann_colors = list(
#  	 Diagnosis = c(AD="darkslateblue", Control="turquoise"),
#  	 Race = c("American Indian or Alaska Native" = "grey60", Asian="tan", "Black or African American" = "black", "Caucasian or White" = "lightyellow"),
#  	 ABeta = c("firebrick1","white"),
#  	 tTau = c("white","darkgreen"),
#  	 pTau = c("white","darkred"),
#  	 MoCA = c("red","white"),
#  	 e4.Risk = c("white","darkviolet"),
#  	 e4.Dose = c("white","darkgreen"),
#  	 sampleSub.ADandCT = c("turquoise"= "turquoise", "blue"= "blue", "brown"= "brown", "yellow"= "yellow", "green"= "green", "red"="red", "black"="black", "pink"="pink", "magenta"="magenta", "purple"="purple", "grey"= "grey"),
#  	 sampleSub.ADonly = c("turquoise"= "turquoise", "blue"= "blue", "brown"= "brown", "yellow"= "yellow", "green"= "green", "red"="red", "grey"= "grey", "grey30"= "grey30"),
        GenoGroup = c("CACNA1A"='goldenrod',"CDKL5"='darkturquoise',"Control"="white","HNRNPH2"="seagreen3","KANSL1"="hotpink","KCNQ2"="purple","MED13L"="darkorange","SLC2A1"="red","STXBP1"="gold","SYNGAP1"="blue","TCF4"="bisque4"),
        age= c("white","darkgreen"),
        sex= c("Female"="pink","Male"="dodgerblue"),
        e4.carrier= c("No"="white","Yes"="darkslateblue"),
        epilepsy= c("unknown"="grey","No"="white","Yes"="red"),
        sampleSubtype = c("turquoise"= "turquoise", "blue"= "blue", "brown"= "brown", "yellow"= "yellow", "green"= "green"),
  	#    plasma.pTau181 = c("white","darkviolet"),
  	 KANSL1vs.CT.trend = c(unchanged = "grey", "upregulated KANSL1 vs CT" = "#E7298A", "downregulated KANSL1 vs CT" = "#66A61E"),
  	 SYNGAP1vs.CT.trend = c(unchanged = "grey", "upregulated SYNGAP1 vs CT" = "#E7298A", "downregulated SYNGAP1 vs CT" = "#66A61E"),
  	 SLC2A1vs.CT.trend = c(unchanged = "grey", "upregulated SLC2A1 vs CT" = "#E7298A", "downregulated SLC2A1 vs CT" = "#66A61E"),
  	 Module = c("grey"="grey","blue"="blue","turquoise"="turquoise","brown"="brown")  #net.protein$colors  # ANOVAout$NETcolors  # gsub("^ME","",rownames(ANOVAout))  # when MEs were run through ANOVAout<-parANOVA.dex()
  )
  # Module color vector must have elements named with the exact same colors.
  # avoids ComplexHeatmap() "Error: elements in `col` should be named vectors."
#  names(ann_colors[["Module"]])<- net.protein$colors  #ANOVAout$NETcolors  # gsub("^ME","",rownames(ANOVAout))  # avoids ComplexHeatmap() "Error: elements in `col` should be named vectors."
  	
  	# Finalize heatmap data matrix
  cleanDat.input <- if (!scaleData) { cleanDat } else { data<-t(apply(cleanDat,1,function(x) x-median(x,na.rm=T))); rownames(data)=rownames(cleanDat); colnames(data)=colnames(cleanDat); data; }
  	
  supervised.cleanDat<-as.matrix(cleanDat.input[which(ANOVAout$'minP'<FDRcutoff),])
  #handle single row with FDR < FDRcutoff -- transpose matrix output from above line
  if(ncol(supervised.cleanDat)<=1) { supervised.cleanDat<-t(supervised.cleanDat); rownames(supervised.cleanDat) <- rownames(cleanDat.input)[which(ANOVAout$'minP'<FDRcutoff)]; }
  	# make data extremes symmetrical same absolute val for max and min of data  - and max z=10
  dataMinAbsExtreme=10  # min(abs(range(supervised.cleanDat,na.rm=T)))
  supervised.cleanDat.windsorized<-supervised.cleanDat
  supervised.cleanDat.windsorized[supervised.cleanDat< -dataMinAbsExtreme] <- -dataMinAbsExtreme
  supervised.cleanDat.windsorized[supervised.cleanDat> dataMinAbsExtreme] <- dataMinAbsExtreme
  	
  	# Generate Graphic and preview
  	# Manual: https://jokergoo.github.io/ComplexHeatmap-reference/book/integrate-with-other-packages.html
  heatmap <- if(!length(which(ANOVAout$'minP'<FDRcutoff))==0) { ComplexHeatmap::pheatmap(supervised.cleanDat.windsorized,scale="none",
#  	 main=paste0("Supervised Clustering of Heparin36 Plasma Samples by Proteins with FDR < ",FDRcutoff,if(scaleData) { " [Scaled Data]" }), #las=1, #angle_col=0,
  	 main=paste0("Supervised Clustering of FNIH296 CSF Samples by Coexpression Network 130 Protein Assays\n...with Select DiffEx Annotation in tWGCNA sample network 5 subtypes, and Extended Trait Annotation. ",if(scaleData) { " [Scaled Data]" }), #las=1, #angle_col=0,
  	 col = c("navy", "white", "firebrick3"), #heatmap main color scheme -- comment this row for default deep blue, yellow, deep red scale
  	 annotation_col = annotation_col, annotation_row = annotation_row, 
  	 annotation_colors = ann_colors,

  	 border_color = NA,
  	 clustering_distance_cols = "correlation",   #*** changed from default "euclidean" 
  	 
	 show_rownames = FALSE,  # default is TRUE ***changed here only
	 cluster_rows = FALSE,  # default is TRUE *** NEW change here.

  	 # SPLIT CLUSTERS BY KEY TRAITS (forcing clustering only within the split groups)	  
  	 column_split = if(forceSplit) { annotation_col$sampleSubtype } else { NULL },  #only available in pheatmap from ComplexHeatmap package, not pheatmap package.
#***  	 row_split = if(forceSplit) {annotation_row$DEP.status1 } else { NULL }#,   #only available in pheatmap from ComplexHeatmap package, not pheatmap package.
  	 row_split = NULL  # only split columns, by DxGroup (AD/Control) in this special use case. #*** commented above line
  	 #above is intended last uncommented line -- no comma at end.

  	  #column_names_gp = grid::gpar(fontsize = 8),
  	  #row_names_gp = grid::gpar(fontsize = 5, fontfamily = "sans", fontface = "bold")
  	  #                    gpar(fontsize = 5, fontfamily = "sans", fontface = "bold")
  	  #labels_row=""
  	  #cluster_rows=FALSE
  	  ) } else { c() } #when there are no rows to plot!
  	#all parameters for tweaking output:   ComplexHeatmap::ht_opt()
  	
  	## Tweak rowname font size if needed (https://www.biostars.org/p/9470610/) -- does not work for ComplexHeatmap package implementation of pheatmap fn.
  	#cols = supervised.cleanDat.windsorized[order(match(rownames(supervised.cleanDat.windsorized), heatmap$gtable$grobs[[5]]$label)), ]$colors  #Assuming row labels are in grob 5
  	#p$gtable$grobs[[5]]$gp = gpar(col = cols, fontsize = 15, fontface = "bold")
  	
  	#Show on Console plot graphics window
  dev.new()
  print(heatmap)  #heatmap.graphics<-recordPlot()
  #Check for no rows meeting FDRcutoff criteria, and return empty vector instead of heatmap if so
  if(length(which(ANOVAout$'minP'<FDRcutoff))==0) { return(list(heatmap=heatmap, data=cleanDat[which(ANOVAout$'minP'<FDRcutoff),])) } else { return(list(heatmap=heatmap, data=supervised.cleanDat.windsorized, FDRcutoff=FDRcutoff)) }
}
####################  



ANOVAout[,2]<-apply(ANOVAout[,3:12],1,min)
colnames(ANOVAout)[2]<-"minP"
ANOVAout$NETcolors<-net.protein$colors


forceSplit=TRUE
scaleData=FALSE
FDRcutoff=1


cleanDat.relatednessOrder<-t(cleanDat[,kMEtable[,1] ])

## Loops to build heatmap with all parameter combinations  [correlation-based column (sample) clustering]
#for (forceSplit in c(TRUE)) {   ## easier to see within Dx columns and within rows grouped by (UP or DOWNREGULATED status) before clustering (default=TRUE)
#  for (scaleData in c(FALSE)) {  ## rowwise center cleanDat on each proteins' median log2 abundance across all samples (default=FALSE)
#    for (FDRcutoff in c(1) ) {  ## Output cluster heatmap at different FDR thresholds for supervision of clustering  (default=0.0001)
      cat(paste0("Producing Heatmap for FDR cutoff <",FDRcutoff, " | scaling rows to median center = ",scaleData," | forcing split groups before clustering = ",forceSplit,"...\n"))

      cleanDat.relatednessOrder2<-t(as.matrix(apply( as.matrix(apply(cleanDat.relatednessOrder,2,as.numeric)) ,1,function(x) (x-mean(x,na.rm=T)) / sd(x,na.rm=T) )))
      rownames(cleanDat.relatednessOrder2)<-rownames(cleanDat.relatednessOrder)
#      numericMeta<-numericMeta.original
      this.heatmap <- makeSupClusterHeatmap.corrColClust.moreTraits.noRowClust.5sampClusterSPLIT.moreAnnColors(cleanDat=cleanDat.relatednessOrder2,traits=numericMeta,ANOVAout=ANOVAout[kMEtable[,1], ],FDRcutoff=FDRcutoff,scaleData=scaleData,forceSplit=forceSplit)

#      cleanDat<-cleanDatMarkers.relatednessOrder.top10.noM1[,which(numericMeta.original$DxGroup=="AD")]
#      numericMeta<-numericMeta.original[which(numericMeta.original$DxGroup=="AD"),]
#      this.heatmap.ADonly <- makeSupClusterHeatmap.corrColClust.moreTraits.noRowClust.5sampClusterSPLIT(cleanDat=cleanDat,traits=numericMeta,ANOVAout=ANOVAout.relatednessOrder,FDRcutoff=FDRcutoff,scaleData=scaleData,forceSplit=forceSplit)

      if(nrow(this.heatmap[["data"]])<110) { pdfWidth=26; pdfHeight=13; } else { pdfWidth=26; pdfHeight=26; }  #helps get row_label font size readable. (every other text item also scaled relative to the page size of the PDF)
pdfWidth=26
pdfHeight=17
      fileName= paste0("4b.Ztransform-All.130Assays-Supervised Clustering Heatmap--SAMPLES_ONLY_5subtypeSPLIT--[FDR lt ",FDRcutoff,"]",if(scaleData) { "scaled" },if(forceSplit) { "_forcedSplits" },"_correlSampClustering.pdf")
      pdf(file=fileName,width=pdfWidth,height=pdfHeight)
        if(nrow(this.heatmap$data)==0) { plot.new(); mtext(paste0("NO HITS IN cleanDat with FDR < ",FDRcutoff),cex=3); } else { print(this.heatmap[["heatmap"]]) }

#        if(nrow(this.heatmap.ADonly$data)==0) { plot.new(); mtext(paste0("NO HITS IN cleanDat (ADonly) with FDR < ",FDRcutoff),cex=3); } else { print(this.heatmap.ADonly[["heatmap"]]) }

      dev.off()
      #Sys.sleep(1.5)

#      cleanDat<-cleanDatMarkers.relatednessOrder.top10.noM1
#      numericMeta<-numericMeta.original
#    }
#  }
#}




## Not run below here 1/27/2026






data1<-cleanDatMarkers.relatednessOrder.top10.noM1
colorVec.samples<-numericMeta.original$sampleSub
colorVec.AD.control<-rep("hotpink",nrow(numericMeta.original))
colorVec.AD.control[numericMeta.original$Group=="Control"]<-"cyan"
colorVec.330<-ANOVAout.relatednessOrder$NETcolors

##transposed (sample tSNE) 2D, all 296 samples
rtsne_out2D.sample <- Rtsne(as.matrix(t(data1)), dims=2, pca = FALSE, verbose = TRUE)


## plot 2D t-SNE projection
library(Cairo)
pdf(file=paste0("./5f.SAMPLE.POINTS_tSNE.2D-FNIH300CSF.PDfinal33of34modules.noM1-top10hubs_byKme.pdf"),height=14,width=14)
#CairoPDF(file=paste0("tSNE.2D-MEGA419final13modules-top",topPercent,"percentile_byKme-Cairo.pdf"),height=14,width=14)
par(mfrow=c(1,1))
par(mar = c(6, 8.5, 3, 3))
plot(rtsne_out2D.sample$Y, asp = 1, pch = 21, col = colorVec.AD.control, bg = colorVec.samples, 
     cex = 2.0, cex.axis = 1.25, cex.lab = 1.25, cex.main = 1.5, 
     xlab = "t-SNE dimension 1", ylab = "t-SNE dimension 2", 
     main = "2D t-SNE projection")
legend("bottomleft",c("Control","AD"), bg="white",col=c("cyan","hotpink"), pch=21, cex=2.5)
dev.off()


## Repeat 2D SAMPLE clustering with selected path/cog corr modules' hubs only before transposing.

data2<-data1[which(colorVec.330 %in% c("royalblue","tan","yellow","darkorange","darkolivegreen")),]
dim(data2)
#[1] 50 296


# run Rtsne (Barnes-Hut-SNE algorithm)
# without PCA step (see Amir et al. 2013, Online Methods, "viSNE analysis")
rtsne_out2D <- Rtsne(as.matrix(t(data2)), dims=2, pca = FALSE, verbose = TRUE)


## plot 2D t-SNE projection
library(Cairo)
pdf(file=paste0("./5f.SAMPLE.POINTS_tSNE.2D-FNIH300CSF.PDfinal33of34modules-top10hubs_byKme[DEMEs_M20.M12,M4.M26.M33].pdf"),height=14,width=14)
#CairoPDF(file=paste0("tSNE.2D-MEGA419final13modules-top",topPercent,"percentile_byKme-Cairo.pdf"),height=14,width=14)
par(mfrow=c(1,1))
par(mar = c(6, 8.5, 3, 3))
plot(rtsne_out2D$Y, asp = 1, pch = 21, col = colorVec.AD.control, bg = colorVec.samples, 
     cex = 2.0, cex.axis = 1.25, cex.lab = 1.25, cex.main = 1.5, 
     xlab = "t-SNE dimension 1", ylab = "t-SNE dimension 2", 
     main = "2D t-SNE projection")
dev.off()



#save.image("supervisedClusterHeatmap.FINAL.03-27-2023.RData")
## we loaded this original saved image to start output of 5f, after changing trait Groups/Dx to 140/160 ct/AD
