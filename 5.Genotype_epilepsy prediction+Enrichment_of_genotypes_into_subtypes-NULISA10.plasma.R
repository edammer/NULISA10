#####################################################################################
## WGCNA Coexpression Network Build - Seyfried Systems Biology Pipeline
##
## input cleanDat from:  plasma transposed 2 batch variable regressed NULISA NPQ data
##
## Eric Dammer -- Seyfried Lab Proteomics (2026)  run 02/03/2026
#####################################################################################

rootdir2="F:/OneDrive - Emory/Faundez_NULISA10/5.prediction/"
setwd(rootdir2)

load("../4.subtyping/4b.saved.image.NULISA10.plasma-eigensampleWGCNA.Rdata")
#contains: cleanDat,numericMeta
rootdir<-rootdir2


dim(cleanDat)
#[1] 465 130
#^transposed (465 samples as rows)


source("binaryPredictor.R")

library(caret)
library(glmnet)
library(xgboost)
library(ranger)
library(progressr)
library(dplyr)
library(future)
library(doFuture)
library(doRNG)


genotypes<-names(table(numericMeta$Group))

fit_fns <- lapply(genotypes, function(gt)
    binaryPredictor(expr = t(na.omit(t(cleanDat))),
    #                    = t(training.cleanDat.noNA[rankedProteins.prior[[gt]][1:50,"feature"],]),     # your 15 000 × p matrix
                          APOE_gt= numericMeta$Group,   # training.gt.APOE,  # factor of length 15 000
                          target = gt,
                          ncores=8,  # parallel::detectCores() - 1,
                          target_ppv=0.95,  # HERE WE USE DEFAULT, 0.95; in stage 2, with only 50 features, we set this to 0.50 (but the floor is 0.80 in the helper)
                          seed   = 1))

names(fit_fns) <- genotypes
#fit_fns_0.95ppv<-fit_fns
## (new vars in environment after running): binaryPredictionMetrics, rankedProteins

head(rankedProteins[["Control"]],20)
       feature importance
NEFH      NEFH  6.0037542
NEFL      NEFL  3.3432740
CALB2    CALB2  1.8935397
ACHE      ACHE  1.8898445
NPTX1    NPTX1  1.8007407
IGFBP7  IGFBP7  1.4621550
IL10      IL10  1.3930647
CST3      CST3  1.3584147
HBA1      HBA1  1.3439774
TNF        TNF  1.2572841
IL4        IL4  1.1924287
GDNF      GDNF  1.0516146
CX3CL1  CX3CL1  1.0061055
CHI3L1  CHI3L1  0.9121146
CCL26    CCL26  0.9033741
POSTN    POSTN  0.8947011
FLT1      FLT1  0.8865642
CCL2      CCL2  0.7765876
IL2        IL2  0.7404990
IL7        IL7  0.6647963

predFeats<-matrix(NA,nrow=20,ncol=0)
for (gt in genotypes) if(dim(head(rankedProteins[[gt]],20))[1]==20) {
  this.set=data.frame(head(rankedProteins[[gt]],20))
  colnames(this.set)<-c(gt,paste0(gt,"_importance"))
  predFeats<-cbind(predFeats,this.set)
}

predFeats  # CDKL5 missing features.
       CACNA1A CACNA1A_importance Control Control_importance  HNRNPH2 HNRNPH2_importance     KANSL1 KANSL1_importance     KCNQ2 KCNQ2_importance   MED13L MED13L_importance  SLC2A1 SLC2A1_importance STXBP1 STXBP1_importance SYNGAP1 SYNGAP1_importance       TCF4 TCF4_importance
NEFL      NEFL          5.2915924    NEFH          6.0037542    PSEN1          4.4605323       NEFH       5.165994069    CD40LG        6.7290100     SNCB         2.5076877   VSNL1         7.0955947   IL1B         3.7424252    BDNF          3.9309473       NEFH      7.96418354
UBB        UBB          4.2816438    NEFL          3.3432740    TAFA5          3.7242592 pTDP43-409       3.447578604      ACHE        5.7133334 pTau-217         2.4247526   FABP3         3.8733060  TIMP3         3.5947967   BASP1          2.1452402     PDLIM5      4.25710509
TREM1    TREM1          3.6924752   CALB2          1.8935397    CALB2          3.2431422 Oligo-SNCA       3.201581017      GOT1        1.9442503     AB38         2.1233648  CHI3L1         2.3437335  CXCL1         3.1081088     CRH          2.1201650      TIMP3      3.91864412
ACHE      ACHE          3.1206444    ACHE          1.8898445 pTau-217          2.8079151       IL1B       1.573507813      PGK1        1.9231625     IL15         2.1191578    NEFL         2.2456070   CCL2         3.0626224   BACE1          2.0925404     IGFBP7      1.89111512
FCN2      FCN2          2.4049713   NPTX1          1.8007407      CRH          2.3815976   pTau-231       1.460107493     BACE1        1.8235871    CALB2         1.6636011   ANXA5         1.7474058 IGFBP7         2.5779721    REST          1.8829680      FOLR1      1.77453100
VEGFA    VEGFA          2.3987611  IGFBP7          1.4621550   CX3CL1          2.0413678       GDI1       1.076820593    PDLIM5        1.2707299     IL10         1.6616311 S100A12         1.6420637  BACE1         2.3319462   CXCL8          1.5417842       FLT1      1.72452528
CXCL8    CXCL8          1.6863896    IL10          1.3930647    NPTX1          2.0104698        PTN       0.750141793    IGFBP7        1.1013558    CNTN2         1.3208302     IL9         1.6238759   FCN2         2.3126460  CX3CL1          1.2283086       AGRN      1.07455144
CALB2    CALB2          1.2410919    CST3          1.3584147     REST          1.8064106       ENO2       0.494951342     YWHAZ        1.0466741    POSTN         1.2601911   BACE1         1.1423248  FOLR1         1.4214414     NGF          1.2017586      UCHL1      0.68509732
TEK        TEK          1.0503637    HBA1          1.3439774      PGF          1.7190386     CD40LG       0.460004363      CCL4        0.6643345     AB42         1.2451300   SFTPD         1.1231755   HBA1         1.3964719    NEFH          1.1867161   pTau-231      0.68212257
TREM2    TREM2          0.7114837     TNF          1.2572841    POSTN          1.5056064      ANXA5       0.407086721      KLK6        0.6065184   CX3CL1         1.2451289   BASP1         1.0248886  FABP3         1.2053695    CST3          1.1704872     CXCL10      0.44318665
SFRP1    SFRP1          0.4709559     IL4          1.1924287      DDC          1.4391228      UCHL1       0.401471600       HTT        0.5869145      NPY         1.2414272  SQSTM1         0.9627896    MME         1.1994609     PGF          1.0600391        KDR      0.26253680
IL13      IL13          0.4597654    GDNF          1.0516146    ANXA5          1.3197242   pTau-181       0.380760880       IL7        0.5634276    TREM1         1.2373973    FLT1         0.7925428   NEFL         1.0896976    IL10          0.9828391        NGF      0.16521765
IL6R      IL6R          0.4251811  CX3CL1          1.0061055     IL6R          0.9754099      PRDX6       0.327424717 pSNCA-129        0.5485958      IL6         1.1692538    APOE         0.7352735   CCL4         0.9664445    ENO2          0.9364604      CCL11      0.14568875
VEGFD    VEGFD          0.3693961  CHI3L1          0.9121146    PARK7          0.8672998       NEFL       0.261196110      MAPT        0.5220691    NPTX1         1.1397189     IL4         0.6907885 CX3CL1         0.9310355    PGK1          0.8811274     PDGFRB      0.07388892
SNCB      SNCB          0.3352320   CCL26          0.9033741     IL1B          0.7039370   pTau-217       0.255468290      IFNG        0.4559406    GDF15         1.1352665    FGF2         0.6820495 TARDBP         0.8316095   PSEN1          0.8659593 Oligo-SNCA      0.02838328
NGF        NGF          0.3039833   POSTN          0.8947011     MAPT          0.6851668      PARK7       0.226444328      IL10        0.4454872      IL9         1.1236300     UBB         0.6711440    IL7         0.7722445 S100A12          0.8440853      CCL13      0.01617336
DDC        DDC          0.2714338    FLT1          0.8865642    CCL11          0.6062340        UBB       0.179994684     IL17A        0.4362125     GFAP         1.0885736     NPY         0.6378021    VGF         0.6385978  CXCL10          0.7871419     CHI3L1     -0.01605786
SMOC1    SMOC1          0.2603274    CCL2          0.7765876     CD63          0.5260840      YWHAZ       0.174479506      MSLN        0.4028235     MDH1         1.0269985   IL17A         0.4572562   IL10         0.6208635   CCL17          0.7833760      CALB2     -0.05413069
PDLIM5  PDLIM5          0.2210531     IL2          0.7404990      IL5          0.3840412       MDH1       0.077525152     CXCL1        0.3653579      IL5         0.9893688    HBA1         0.4356635    NPY         0.5452600    CCL3          0.7212026      SFRP1     -0.07168851
UCHL1    UCHL1          0.1982996     IL7          0.6647963    NPTX2          0.3558421        CRH       0.006502523  pTau-181        0.3566668 pTau-181         0.9530451    GDI1         0.4084633  PSEN1         0.3326864    GOT1          0.6608235   pTau-217     -0.07548408


rankedProteins.prior<-rankedProteins
## ^ importance-ranked assays are input for stage 2


###########################################
## Stage 2 Reprediction on top 10 important features using 10 binary predictor ensemble, target PPV 80%

genotypes <- genotypes[which(!genotypes %in% c("CDKL5"))]

## force the same ordering as the `genotypes` vector
rankedProteins.prior <- rankedProteins.prior[genotypes]


## (1)  train the final model on *all* data for real-world deployment on unknown samples
nulisa10_full <- train_NULISA10_full(t(na.omit(t(cleanDat))),
                               numericMeta$Group,
                               rankedProteins.prior,
                               target_ppv = 0.80,    # preferred PPV (floor is 0.80 in internal helper function)
                               topX = 10,
                               ncores = 8)

#[CACNA1A]  CV Precision NaN +/- NA | Recall 0.000 +/- 0.000
#[Control]  CV Precision 0.765 +/- 0.072 | Recall 0.145 +/- 0.113
#[HNRNPH2]  CV Precision 1.000 +/- 0.000 | Recall 0.089 +/- 0.107
#[KANSL1]  CV Precision NaN +/- NA | Recall 0.000 +/- 0.000
#[KCNQ2]  CV Precision 1.000 +/- 0.000 | Recall 0.020 +/- 0.069
#[MED13L]  CV Precision NaN +/- NA | Recall 0.000 +/- 0.000
#[SLC2A1]  CV Precision 0.833 +/- 0.302 | Recall 0.139 +/- 0.157
#[STXBP1]  CV Precision 0.450 +/- 0.497 | Recall 0.018 +/- 0.037
#[SYNGAP1]  CV Precision NaN +/- NA | Recall 0.000 +/- 0.000
#[TCF4]  CV Precision NaN +/- NA | Recall 0.000 +/- 0.000

nulisa10_full.top10<-nulisa10_full


nulisa10_full.top5 <- train_NULISA10_full(t(na.omit(t(cleanDat))),
                               numericMeta$Group,
                               rankedProteins.prior,
                               target_ppv = 0.80,    # preferred PPV (floor is 0.80 in internal helper function)
                               topX = 5,
                               ncores = 25)

#[CACNA1A]  CV Precision NaN +/- NA | Recall 0.000 +/- 0.000
#[Control]  CV Precision 0.460 +/- 0.389 | Recall 0.034 +/- 0.055
#[HNRNPH2]  CV Precision 0.650 +/- 0.474 | Recall 0.044 +/- 0.078
#[KANSL1]  CV Precision 1.000 +/- NA | Recall 0.005 +/- 0.025
#[KCNQ2]  CV Precision 1.000 +/- NA | Recall 0.010 +/- 0.050
#[MED13L]  CV Precision NaN +/- NA | Recall 0.000 +/- 0.000
#[SLC2A1]  CV Precision 0.817 +/- 0.267 | Recall 0.179 +/- 0.122
#[STXBP1]  CV Precision 0.125 +/- 0.354 | Recall 0.004 +/- 0.018
#[SYNGAP1]  CV Precision NaN +/- NA | Recall 0.000 +/- 0.000
#[TCF4]  CV Precision NaN +/- NA | Recall 0.000 +/- 0.000


## (4)  predict unknown samples
# unknown sample prediction using predictor trained on full noNA data
nulisa.10pred <- nulisa10_full.top5(t(na.omit(t(cleanDat))) , ncores=25)  #[,names(gt.APOE)[which(is.na(gt.APOE))] ]), ncores=8)
table(nulisa.10pred)
#CACNA1A         Control HNRNPH2  KANSL1   KCNQ2  MED13L  SLC2A1  STXBP1 SYNGAP1    TCF4 
#     19             160      39      45      22      21      44      57      44      14
table(numericMeta$Group)
#CACNA1A   CDKL5 Control HNRNPH2  KANSL1   KCNQ2  MED13L  SLC2A1  STXBP1 SYNGAP1    TCF4 
#     18      23     148      37      41      22      20      43      56      43      14

names(nulisa.10pred)<-rownames(numericMeta)  # names(gt.APOE)[which(is.na(gt.APOE))]
length(which(nulisa.10pred==numericMeta$Group))
# 442


######################################
### Restart Stage 1 - predict epilepsy

epilepsy<-numericMeta$Epilepsy[which(!is.na(numericMeta$Epilepsy))]
names(epilepsy)<-rownames(numericMeta)[which(!is.na(numericMeta$Epilepsy))]
epilepsy[epilepsy=="Yes"]<-"Epilepsy"
table(epilepsy)
#Epilepsy       No 
#     149      222

length(epilepsy)
# 371


#known epilepsy sample data, no missing values:
#cleanDat.epilepsy=t(na.omit(t(cleanDat)))[names(epilepsy),]

epilepsy.fit_fn <- lapply(c("Epilepsy","No"), function(gt)
    binaryPredictor(expr =t(na.omit(t(cleanDat)))[names(epilepsy),] ,
    #                    = t(training.cleanDat.noNA[rankedProteins.prior[[gt]][1:50,"feature"],]),     # your 15 000 × p matrix
                          APOE_gt= epilepsy,   # training.gt.APOE,  # factor of length 15 000
                          target = gt,
                          ncores=8,  # parallel::detectCores() - 1,
                          target_ppv=0.95,  # HERE WE USE DEFAULT, 0.95; in stage 2, with only 50 features, we set this to 0.50 (but the floor is 0.80 in the helper)
                          seed   = 1))

#(previously) CV Yes - Precision 0.887 +/- 0.148 | Recall 1.000 +/- 0.000
#CV Epilepsy - Precision 0.931 +/- 0.149 | Recall 1.000 +/- 0.000
#CV No - Precision 0.876 +/- 0.141 | Recall 1.000 +/- 0.000

names(epilepsy.fit_fn) <- c("positive","negative")
## (new vars in environment after running): binaryPredictionMetrics, rankedProteins

lastColNum=length(names(rankedProteins))
names(rankedProteins)[c((lastColNum-1):lastColNum)] <- c("Epilepsy","No")

head(rankedProteins[["Epilepsy"]],20)

predFeats.epi<-matrix(NA,nrow=20,ncol=0)
for (gt in c("Epilepsy","No")) if(dim(head(rankedProteins[[gt]],20))[1]==20) {
  this.set=data.frame(head(rankedProteins[[gt]],20))
  colnames(this.set)<-c(gt,paste0(gt,"_importance"))
  predFeats.epi<-cbind(predFeats.epi,this.set)
}

predFeats.epi
#       Epilepsy Epilepsy_importance     No No_importance
#HBA1       HBA1           4.5929784   HBA1     4.1053412
#IGFBP7   IGFBP7           3.8952989   KLK6     3.2597948
#CRH         CRH           3.0479254    DDC     3.2357337
#NEFL       NEFL           3.0331622   FLT1     3.1115648
#FLT1       FLT1           2.7012729 IGFBP7     2.7293426
#ANXA5     ANXA5           2.1313775   CCL4     2.6181717
#DDC         DDC           2.0796786   NEFL     2.4537298
#IL10       IL10           1.7124297   IL10     2.0205561
#VCAM1     VCAM1           1.6837726  ANXA5     1.6069774
#CCL4       CCL4           1.6501709    CRH     1.5085497
#CRP         CRP           1.3544448    MME     1.0655302
#MAPT       MAPT           1.1239118   CST3     1.0509662
#CST3       CST3           0.9071813  NPTX1     0.8762381
#GDF15     GDF15           0.7421629   MAPT     0.7716189
#UCHL1     UCHL1           0.7088304  VCAM1     0.7699722
#MME         MME           0.6836618    TNF     0.6800827
#CHIT1     CHIT1           0.6279624   IL13     0.6631550
#TREM1     TREM1           0.5653328  CCL13     0.5302359
#IL33       IL33           0.5333263 PDLIM5     0.5065959
#ACHE       ACHE           0.4561012    IL9     0.4416678




###########################################
## Epilepsy prediction - Stage 2 prediction on top 12 important features using binary predictor ensemble, target PPV 80%

genotypes <- c("Epilepsy","No")

epilepsy[epilepsy=="Yes"]<-"Epilepsy"
table(epilepsy)
#Epilepsy       No 
#     149      222

## force the same ordering as the `genotypes` vector
rankedProteins.epi <- rankedProteins[genotypes]


source("binaryPredictor.R")
## (1)  train the final model on *all* data for real-world deployment on unknown samples
nulisa10_epi.top12 <- train_NULISA10_full(t(na.omit(t(cleanDat)))[names(epilepsy), ],  # [which(!rownames(cleanDat) %in% names(epilepsy)),] ,
                               epilepsy,
                               rankedProteins.epi,
                               target_ppv = 0.80,    # preferred PPV (floor is 0.80 in internal helper function)
                               topX = 12,
                               ncores = 25)

#[Epilepsy]  CV Precision 0.893 +/- 0.102 | Recall 0.216 +/- 0.156  # top 12
#[No]  CV Precision 0.841 +/- 0.076 | Recall 0.424 +/- 0.157

#[Epilepsy]  CV Precision 0.874 +/- 0.147 | Recall 0.193 +/- 0.132  # top 10


## (4)  predict unknown samples
# unknown sample prediction using predictor trained on full noNA data
nulisa10_epi.pred <- nulisa10_epi.top12(t(na.omit(t(cleanDat))) , ncores=25)  #[,names(gt.APOE)[which(is.na(gt.APOE))] ]), ncores=8)
names(nulisa10_epi.pred)<-rownames(cleanDat)  # sample names
table(nulisa10_epi.pred[which(names(nulisa10_epi.pred) %in% names(epilepsy))])
#Epilepsy       No 
#     149      222
table(numericMeta$Epilepsy)
# No Yes 
#222 149    # 371 total

#unknown samples
table(nulisa10_epi.pred[which(!names(nulisa10_epi.pred) %in% names(epilepsy))])
#Epilepsy       No 
#      33       61 

nulisa10_epi.pred<-as.character(nulisa10_epi.pred)
nulisa10_epi.pred[nulisa10_epi.pred=="Epilepsy"]<-"Yes"
length(which(nulisa10_epi.pred==numericMeta$Epilepsy))
# 371 / 371 samples predicted correctly.

numericMeta$Epilepsy.predicted <- nulisa10_epi.pred

numericMeta[,c("Group","Epilepsy","Epilepsy.predicted")]
# review predictions for individuals with epilepsy vs genotype

numericMeta[numericMeta$Group=="Control",c("Group","Epilepsy","Epilepsy.predicted")]
# only 1 unknown Control is predicted to have Epilepsy.  C_06_P3_6A  Some (n=3) controls are known Epileptic.

numericMeta[is.na(numericMeta$Epilepsy),c("Age","Sex","Group","Epilepsy","Epilepsy.predicted")]



#######################################################
## Enrichment of samples by genotype into subtypes


# split rownames by Group
rn_by_group <- split(rownames(numericMeta), numericMeta$Group)

# find the maximum group size to pad shorter groups
max_len <- max(lengths(rn_by_group))

# pad each group with NA so all columns have equal length
rn_by_group_padded <- lapply(rn_by_group, function(x) {
  length(x) <- max_len
  x
})

# build the final data.frame
sampGeno <- as.data.frame(rn_by_group_padded, stringsAsFactors = FALSE)

dim(sampGeno)
# 148 x 11
write.csv(sampGeno,file="samplesByGenotype_forFET.csv",row.names=FALSE)
# removed NA values in excel.


##########
# split names(net$colors) by net$colors
rn_by_group <- split(names(net$colors), net$colors)

# find the maximum group size to pad shorter groups
max_len <- max(lengths(rn_by_group))

# pad each group with NA so all columns have equal length
rn_by_group_padded <- lapply(rn_by_group, function(x) {
  length(x) <- max_len
  x
})

# build the final data.frame
sampSubtype <- as.data.frame(rn_by_group_padded, stringsAsFactors = FALSE)

dim(sampSubtype)
# 118 x 5
write.csv(sampSubtype,file="samplesBySubtype.csv",row.names=FALSE)
# removed NA values in excel.


#################### CONFIGURATION PARAMETERS FULL LIST #############################
##             WITH SAMPLE VALUES GIVEN FOR SAMPLE DATA PROVIDED                   ##

heatmapScale="minusLogFDR"                        				# Accepted options are "p.unadj" or "minusLogFDR"
heatmapTitle="NULISA 5 Subtypes Overlaps with 11 Known Genotypes of Samples in Network"	# What are your categories (or WGCNA) list of lists based on?
                                                                                        # And What gene lists are your reference lists?
paletteColors="YlGnBu"                                          # See valid palettes using RColorBrewer::display.brewer.all()
                                                                # Can be a vector if there are more than 1 refDataFiles (heatmaps to generate)

FileBaseName="5NULISA_Subtypes_FET_to_11knownGenotypes"
refDataDescription="11knownGenotypes"				# One Description of reference Data list(s) specified in PDF file name below
# File Names of Reference List(s): You will get one output PDF page per file
refDataFiles <- c(      "samplesByGenotype_forFET.csv")	# HUMAN gene symbols have been pre-converted from the original below MOUSE lists.

speciesCode=c("hsapiens") 				# species code(s) for biomaRt (one for each refDataFile)#one for each .csv in refDataFiles


# Use Modules in memory OR a .csv file with your input gene lists.
modulesInMemory=FALSE                              		    # Load modules as categories? (If TRUE, categoriesFile not used, but you need cleanDat, net[["colors"]] and numericMeta variables)
categoriesFile="samplesBySubtype.csv"	# File Name of Categories (Lists of Fly genes), only loaded if modulesInMemory=FALSE
								                            # NOTE this file format has a column for official gene symbols of each module or cluster, with the cluster name/ID as column names in row 1
categorySpeciesCode="hsapiens"				            # What species are the gene sybmols in categoriesFile?
#bkgrFileForCategories="Dmel_geneBackground.csv"             # 12/2025: Augment background of categoriesFile symbols with a single-column no-header list (in a no comma .CSV file)

# Other Options
allowDuplicates=TRUE				# Allow duplicate symbols across different row lists (refDataFiles) for overlap?
                                    # if FALSE, row lists' duplicte entries removed entirely (since removing only extra instance(s) would have to choose to keep in one row)
						            # (should be true if you have general cell type lists and e.g. disease-associated phenotype cell type lists)
strictSymmetry=FALSE                # remove duplicated symbols across columns (modules or categoriesFile, background lists)
                                    # if TRUE, and all row list symbols are present in column lists, background, then swapping row and column input will give similar or identical stats
resortListsDecreasingSize=FALSE		# resort categories/modules and reference data lists? (decreasing size order)
heatmapScale="FDR"                  # "FDR" (or still accepted previous default "minusLogFDR") will FDR-correct all p values used for plotting. Otherwise "p" or "p.unlog". This is independent of the legendScale setting.
legendScale="minusLog"              # Heatmap legend will be plotted with increasing significance and color intensity for higher -log10(p or FDR). "unlog" option has increasing significance for lower p or FDR.
barOption=FALSE					    # Draw bar charts for each list overlap instead of a heatmap.
asterisksOnly=FALSE                 # 12/2025: new (TRUE, default); if FALSE, print p (or FDR) values within heatmap cells if below maxPcolor and expand width of columns to accomodate this text
maxPcolor=0.20                      # 12/2025: p (or FDR) values below this value will get some heatmap color in the scale; above the value, cells will be WHITE. Plots with no values below this value only show color in the legend.
adjustFETforLookupEfficiency=FALSE	# adjust p FET input for cross-species lookup inefficiency/loss of list member counts?
verticalCompression=1				# DEPRECATED for values >1. Plot(s) may be squeezed into 1 row out of this many in each PDF page, compressing the heatmap tracks vertically (or the bar chart heights) for each reference list)
                                    # (PDF dimensions and margins for plots are now dynamically calculated to fit all plot elements and text optimally, one per page.)
reproduceHistoricCalc=FALSE			# should be FALSE unless trying to reproduce exact calculations of prior publications listed.
#####################################################################################


## Generate Sample Outputs

# Load Seyfried/Emory pipeline FET as function geneListFET() having all the parameters described above, many with defaults used.
source("./geneListFET.R")


## output enrichment significance of provided test network data in memory as -log10(FDR) heatmap; enrichments checked are in provided 5 brain cell type marker gene lists
geneListFET(FileBaseName=FileBaseName,
            heatmapTitle=heatmapTitle,
            maxPcolor=maxPcolor, asterisksOnly=asterisksOnly, paletteColors="plasma",  # palettes from viridisLite package valid, as well as RcolorBewer
            modulesInMemory=modulesInMemory,categoriesFile="samplesBySubtype.csv", categorySpeciesCode="hsapiens",  # use network in memory; what species code are the symbols in cleanDat rownames? In case symbol interconversion across species is needed...
            refDataFiles=refDataFiles,speciesCode=c("hsapiens"),refDataDescription=refDataDescription)  #file(s) with columns of reference gene lists to check for overlap in; what are the species code(s) for symbols in each file?




















					  memLimit=10*1024^3
					  options(future.globals.maxSize= memLimit)  #4GB Total size of all global objects that need to be exported - up from 500MB
					  Sys.setenv(R_FUTURE_GLOBALS_MAXSIZE=memLimit) #inherited by workers
					
					
					## (2)  predict training samples
					train.APOEpred<-apoe6_full(t(training.cleanDat.noNA), ncores=32)
					
					## Naive statistics; mean accuracy and confusion matrix:
					mean(train.APOEpred == training.gt.APOE)
					#0.9888196 (original run); rerun: 0.9892262 due to non-determinism
					table(train.APOEpred, training.gt.APOE)  #confusion matrix -- should (nearly) match below before rebuild of above function system
					#               training.gt.APOE
					#train.APOEpred e2/e2 e2/e3 e2/e4 e3/e3 e3/e4 e4/e4
					#         e2/e4     0     0   351     1     0     0
					#         e2/e2    50     0     0     0     0     0
					#         e4/e4     0     1     0     4     0   841
					#         e2/e3     0  1368     0     1     5     1
					#         e3/e4     0     1     4    39  4676     7
					#         e3/e3     0    17     5  7307    67    12




save.image("saved.image02-03-2026.RData")


##############
## 2 KANSL1 enriched subtypes -- here are the samples for each:

rootdir2="F:/OneDrive - Emory/Faundez_NULISA10/5.prediction/"
setwd(rootdir2)
load("../4.subtyping/4b.saved.image.NULISA10.plasma-eigensampleWGCNA.Rdata")
#contains: cleanDat,numericMeta
outputfigs=rootdir2

KANSL1.yellow=c("G_04_Sample.31_AB000457","H_04_Sample.32_AB001604","A_09_Sample.33_AB000543","B_11_Sample.36_AB001651","D_09_Sample.38_AB001686","C_11_Sample.39_AB000277","G_11_Sample.47_AB000473","H_11_Sample.48_AB000517","B_10_Sample.50_AB000509","C_12_Sample.55_AB000300","E_10_Sample.57_AB001555","G_10_Sample.61_AB000500","A_07_Sample.67_AB000329","B_07_Sample.68_AB001524")
KANSL1.blue=c("A_11_Sample.35_AB000442","D_11_Sample.40_AB000404","F_09_Sample.42_AB000293","E_11_Sample.43_AB001630","F_11_Sample.44_AB001663","G_09_Sample.45_AB001738","H_09_Sample.46_AB001591","A_10_Sample.49_AB000423","A_12_Sample.51_AB000368","B_12_Sample.52_AB000534","C_10_Sample.53_AB001544","D_12_Sample.56_AB000362","E_12_Sample.59_AB000384","H_10_Sample.62_AB000317","H_12_Sample.64_AB000525","A_05_Sample.65_AB000346","B_05_Sample.66_AB001559","C_07_Sample.71_AB001669")
KANSL1.both=c(KANSL1.yellow,KANSL1.blue)

Grouping=c(rep("KANSL1.yellow",length(KANSL1.yellow)), rep("KANSL1.blue",length(KANSL1.blue)))

#cleanDat<-t(cleanDat[c(KANSL1.yellow,KANSL1.blue),which(!colnames(cleanDat) %in% c("BD-MAPT","BD-pTau181","BD-pTau217","BD-pTau231"))])
cleanDat<-t(cleanDat[c(KANSL1.yellow,KANSL1.blue),])
#32 samples out of 41 total KANSL1 in these blue and yellow subtypes.

parallelThreads=31
outFilePrefix="5"
outFileSuffix="Faundez_NULISA-KANSL1_2subtypes-rerun"

ANOVAout<-parANOVA.dex()
ANOVAout$NETcolors=net.protein$colors #[rownames(cleanDat)]

labelTop=5
corVolc=FALSE
useNETcolors=TRUE

plotVolc()

