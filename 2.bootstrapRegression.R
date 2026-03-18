rootdir2<-"f:/OneDrive - Emory/Faundez_NULISA10/"
setwd(rootdir2)
load("./1.FaundezNULISA.RData")

outputfigs<-outputtabs<-rootdir<-rootdir2
setwd(rootdir)


#######################################################
# Regress exprMat0 for Batch effects

cleanDat.unreg<-exprMat0
#numericMeta<-targets
dim(numericMeta)
#[1] 465  14


library("doParallel")
parallelThreads=31 #max is number of processes that can run on your computer at one time
#stopCluster(clusterLocal)
clusterLocal <- makeCluster(c(rep("localhost",parallelThreads)),type="SOCK")

registerDoParallel(clusterLocal)


library(boot)
boot <- TRUE
numboot <- 1000
bs <- function(formula, data, indices) {
  d <- data[indices,] # allows bootstrap function to select samples
  fit <- lm(formula, data=d)
  return(coef(fit))
}  


condition<-factor(numericMeta$Group) 
condition<-relevel(condition,"Control")

Experiment=as.numeric(as.factor(numericMeta$Experiment))-1
Foundation=as.numeric(as.factor(numericMeta$Foundation))-1

table(Experiment)
#  0   1 
# 52 413
 
table(Foundation)
#  0   1 
#434  31

# not protecting Disease/Control Groups


#regvars <- data.frame(condition=condition, bloodEP.Soma=numericMeta$bloodEigenprotein.SOMA.9genes)
regvars <- data.frame(Experiment=Experiment, Foundation=Foundation)

## Run the regression
normExpr.reg <- matrix(NA,nrow=nrow(cleanDat.unreg),ncol=ncol(cleanDat.unreg))
rownames(normExpr.reg) <- rownames(cleanDat.unreg)
colnames(normExpr.reg) <- colnames(cleanDat.unreg)
coefmat <- matrix(NA,nrow=nrow(cleanDat.unreg),ncol=ncol(regvars)+2) ## change this to ncol(regvars)+2 when condition has 2 levels if BOOT=TRUE, +1 if BOOT=FALSE

#another RNG seed set for reproducibility
set.seed(8675309);

##precheck df.reg positions for subtraction of coefficient*term; this should be the same model matrix set in the foreach below
thisexp=as.numeric(cleanDat.unreg[1,])
#df.reg<-model.matrix(~thisexp +condition +Batch, data=data.frame(thisexp,regvars))
df.reg<-model.matrix(~thisexp +Experiment+Foundation, data=data.frame(thisexp,regvars))
colnames(df.reg)
#[1] "(Intercept)"    "thisexp"        "Experiment" "Foundation"

set.seed(8675309)
cat('[bootstrap-PARALLEL] Working on ORDINARY NONPARAMETRIC BOOTSTRAP regression with ', parallelThreads, ' threads over ', nrow(cleanDat.unreg), ' iterations.\n Estimated time to complete:', round(120/parallelThreads*nrow(cleanDat.unreg)/2736,1), ' minutes.\n') #intermediate progress printouts would not be visible in parallel mode
coefmat <- foreach (i=1:nrow(cleanDat.unreg), .combine=rbind) %dopar% {
  set.seed(8675309)
  options(stringsAsFactors=FALSE)
  library(boot)
  thisexp <- as.numeric(cleanDat.unreg[i,])

#  df.reg<-model.matrix(~thisexp +condition+Batch, data=data.frame(thisexp,regvars))
  df.reg<-model.matrix(~thisexp +Experiment+Foundation, data=data.frame(thisexp,regvars))

  bs.results <- boot(data=as.data.frame(df.reg),statistic=bs,
                   R=numboot, formula=thisexp ~ .)  ## run 1000 resamplings

apply(bs.results$t,2,function(x) median(x,na.rm=TRUE))  #this is bs.stats (output to coefmat 1 row at a time)
}
coefmat.noNA<-coefmat
coefmat.noNA[is.na(coefmat)]<-0

#** NEW (w/ model.matrix) bs.stats has NA for thisexp (position 2) now, and a variable for any nonreference factor split-out; if in doubt about order and numbering/index check the df.reg column names above.
normExpr.reg <-foreach (i=1:nrow(cleanDat.unreg), .combine=rbind) %dopar% { (cleanDat.unreg[i,] - coefmat.noNA[i,3]*Experiment - coefmat.noNA[i,4]*Foundation) }


## Sanity Check -- Did regression do something unexpected to our abundance data?
quantile(cleanDat.unreg[,440],c(0,0.025,0.25,0.5,0.75,0.975,1),na.rm=TRUE)
#       0%      2.5%       25%       50%       75%     97.5%      100% 
# 5.926346  7.638857 10.541860 11.751133 12.741634 14.749210 16.049166

quantile(normExpr.reg[,440],c(0,0.025,0.25,0.5,0.75,0.975,1),na.rm=TRUE)
#       0%      2.5%       25%       50%       75%     97.5%      100% 
# 4.819811  7.544088 11.140436 12.045598 13.042994 14.504197 18.833918

##Write cleanDat (candidate) with regressed data
############################################
cleanDat<-normExpr.reg
rownames(cleanDat)<-rownames(cleanDat.unreg)

#numericMeta.Soma.CSF$RegrBloodEigenprotein.Soma9genes<-NA
#numericMeta.Soma.CSF$RegrBloodEigenprotein.Soma9genes[match(colnames(cleanDat),rownames(numericMeta.Soma.CSF))]<-bloodEP.Soma
############################################



## 2b. Variance Partition pre adjustment as QC
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

varPart2 <- fitExtractVarPartModel(cleanDat, form, regvars.vp,  BPPARAM=BiocParallel::SnowParam(workers = parallelThreads, type = "SOCK")) 
vp2 <- sortCols(varPart2,FUN=median,last= c("Residuals"))


pdf(file="2.VariancePartition-NULISA-130x465_Regressed(Experiment+Foundation).pdf", width=15,height=11)
par(mfrow=c(1,1))
par(oma=c(2,2,3,1))

plotVarPart( vp2, main="NULISA 465 - Batch Regressed NPQ (Experiment+Foundation)" )


DxSortOrder<-order(vp2[["Group"]],decreasing=TRUE)
#rownames(cleanDat.noNA.unreg)[DxSortOrder][1:50]
for (i in ls(vp2)) { vp2[[i]]<-vp2[[i]][DxSortOrder]; }
rownames(vp2)<-rownames(vp2)[DxSortOrder]

plotPercentBars( vp2[1:50,]) + ggtitle( "Top Group (Control/KANSL1/CDKL5..TCF4) Covariates" ) + theme(plot.margin = margin(t = 6, r = 1, b = 3, l = 1) )


DxSortOrder<-order(vp2[["Foundation"]],decreasing=TRUE)
#rownames(cleanDat.noNA.unreg)[DxSortOrder][1:50]
for (i in ls(vp2)) { vp2[[i]]<-vp2[[i]][DxSortOrder]; }
rownames(vp2)<-rownames(vp2)[DxSortOrder]

plotPercentBars( vp2[1:50,]) + ggtitle( "Top Foundation (Batch) Covariates" ) + theme(plot.margin = margin(t = 6, r = 1, b = 3, l = 1) )



DxSortOrder<-order(vp2[["Experiment"]],decreasing=TRUE)
#rownames(cleanDat.noNA.unreg)[DxSortOrder][1:50]
for (i in ls(vp2)) { vp2[[i]]<-vp2[[i]][DxSortOrder]; }
rownames(vp2)<-rownames(vp2)[DxSortOrder]

plotPercentBars( vp2[1:50,]) + ggtitle( "Top Experiment (Batch) Covariates" ) + theme(plot.margin = margin(t = 6, r = 1, b = 3, l = 1) )


DxSortOrder<-order(vp2[["APOE4"]],decreasing=TRUE)
#rownames(cleanDat.noNA.unreg)[DxSortOrder][1:50]
for (i in ls(vp2)) { vp2[[i]]<-vp2[[i]][DxSortOrder]; }
rownames(vp2)<-rownames(vp2)[DxSortOrder]

plotPercentBars( vp2[1:50,]) + ggtitle( "Top APOE4 Covariates" ) + theme(plot.margin = margin(t = 6, r = 1, b = 3, l = 1) )


DxSortOrder<-order(vp2[["AB42"]],decreasing=TRUE)
#rownames(cleanDat.noNA.unreg)[DxSortOrder][1:50]
for (i in ls(vp2)) { vp2[[i]]<-vp2[[i]][DxSortOrder]; }
rownames(vp2)<-rownames(vp2)[DxSortOrder]

plotPercentBars( vp2[1:50,]) + ggtitle( "Top Amyloid Beta 42 Covariates" ) + theme(plot.margin = margin(t = 6, r = 1, b = 3, l = 1) )


AgeSortOrder<-order(vp2[["Age"]],decreasing=TRUE)
#rownames(cleanDat.noNA.unreg)[AgeSortOrder][1:50]
for (i in ls(vp2)) { vp2[[i]]<-vp2[[i]][AgeSortOrder]; }
rownames(vp2)<-rownames(vp2)[AgeSortOrder]

plotPercentBars( vp2[1:50,]) + ggtitle( "Top Age Covariates" ) + theme(plot.margin = margin(t = 6, r = 1, b = 3, l = 1) )


SexSortOrder<-order(vp2[["Sex"]],decreasing=TRUE)
#rownames(cleanDat.noNA.unreg)[SexSortOrder][1:50]
for (i in ls(vp2)) { vp2[[i]]<-vp2[[i]][SexSortOrder]; }
rownames(vp2)<-rownames(vp2)[SexSortOrder]

plotPercentBars( vp2[1:50,]) + ggtitle( "Top Sex Covariates" ) + theme(plot.margin = margin(t = 6, r = 1, b = 3, l = 1) )


dev.off()



save.image("2.FaundezNULISA_regr.RData")
