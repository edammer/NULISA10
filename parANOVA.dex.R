######################################################################################################################
## ANOVA / DiffEx Functions
##
## parANOVA.dex()     - output ANOVAout dataframe for all pairwise comparisons of each cleanDat gene product (row)
##                      among specified sample groups in the vector Grouping
## plotVolc(ANOVAout) - using output of parANOVA.dex(), output PDF and plotly HTML volcano plots
##
## DEXpercentStacked(ANOVAout)
##                    - using output of parANOVA.dex(), plotVolc settings, plus modules in memory,
##                      output pdf with a page for each pairwise comparison showing % of module members meeting
##                      differential expression criteria in volcano.
##
## trait.corStat()    - creates alternative output to ANOVAout, using only correlation R and p (effect size and significance)
##                      correlation will be to one or more traits available in numericMeta, across all samples (default),
##                      or subset(s) of samples defined by the values in numericMeta for a second categorical, discrete, trait.
##
#################################################
## By Eric Dammer, Duc Duong, and Qiudong Deng ##
#################################################


######################################
## Required variables as Parameters ##
######################################


parANOVA.dex <- function(dummyVar="",
## All parameters are set as variables of the global environment, or if not set, fallback to defaults in the function.
#                         cleanDat.=cleanDat,                             # cleanDat, data table of 'cleaned' relative abundance with rows (genes or proteins) and columns (samples)
#                         Grouping=numericMeta$Group,                     # Named groups (N>=2) for comparison of difference of means, in sample (column) order of cleanDat
#                         parallelThreads=30,                             # number of CPU threads to speed calculation (recommended, 2 or more)
#                         NETcolors=net$colors,                           # a vector containing module color assignments
#                                                                         # length of vector must be equal to number of rows in cleanDat.
#                                                                         # Colors are only added to table output if available; function tries to use NETcolors variable , with fallback to net$colors
#                         twoGroupCorrMethod="BH",                        # default method for full FDR correction when only 2 groups present is Benjamini-Hochberg
#                         outputCSV=TRUE,                                 # Output Full Table of statistics?  TRUE/FALSE
#                         outFilePrefix="4",                              # typically "4", or step # in pipeline being run
#                         outFileSuffix=FileBaseName,                     # A description of the project, used as a filename suffix
#                         fallbackIfSmallTukeyP=TRUE,
                         env=.GlobalEnv) {

############################################
## PARAMETER VARIABLES CHECKED / DEFAULTS ##
############################################
  if (!exists("cleanDat")) stop("cleanDat variable must exist, holding gene product (rows) X sample (columns) data in the form of log2(relative abundance).")
  if (!exists("Grouping")) { if ("Group" %in% colnames(numericMeta)) { Grouping=numericMeta[,"Group"] } else { stop("Grouping vector not found. Must specify a group for every sample column in cleanDat.\n") }}

  if (length(which(grepl("\\-",Grouping)))>0) { cat("- Grouping values have '-'; replacing with '.'...\n"); Grouping<-gsub("\\-",".",Grouping); }
  if (!exists("parallelThreads")) { cat("- parallelThreads variable not set. Running with 2 threads only.\n"); parallelThreads=2; }
  if (!exists("NETcolors")) if(exists("net")) { if ("colors" %in% names(net)) { NETcolors=net$colors } else { NETcolors=c() } } else { NETcolors=c() }
  if (!length(NETcolors)==nrow(cleanDat)) { cat("- Network color assignment vector not supplied or not of length in rows of cleanDat; will not be included in output table and data frame.\n") }
  if (!exists("twoGroupCorrMethod")) { cat("- twoGroupCorrMethod variable not set to a correction method for p.adjust when only 2 groups of samples specified in Grouping. Using Benjamini Hochberg 'BH' FDR.\n"); twoGroupCorrMethod="BH"; }
  if (!exists("outputCSV")) { outputCSV=TRUE } else { if(!is.logical(outputCSV)) { outputCSV=TRUE } }
  if (!exists("outFilePrefix")) { outFilePrefix="" } else { if (nchar(outFilePrefix)>0) outFilePrefix=paste0(outFilePrefix,".") }
  if (!exists("outFileSuffix")) { if (exists("FileBaseName")) { outFileSuffix=paste0("-",FileBaseName) } else { outFileSuffix="-unspecified_study" }} else { if (nchar(outFileSuffix)>0) outFileSuffix=paste0("-",outFileSuffix) }
  if (!exists("fallbackIfSmallTukeyP")) { cat("- fallbackIfSmallTukeyP variable not set. Using recommended Bonferroni t-test FDR for unreliable Tukey p values <10^-8.5\n"); fallbackIfSmallTukeyP=TRUE; }
  if (!is.logical(fallbackIfSmallTukeyP)) { cat("- fallbackIfSmallTukeyP variable not TRUE/FALSE. Using recommended Bonferroni t-test FDR for unreliable Tukey p values <10^-8.5.\n"); fallbackIfSmallTukeyP=TRUE; }
  if (!exists("Tunequal")) { Tunequal=FALSE } else { if(!is.logical(Tunequal)) { Tunequal=TRUE } }

## Set up parallel backend as a local cluster using n (parallelThreads) CPU threads (in the .GlobalEnv outside of this function!)
  require(doParallel, quietly=TRUE)
  require(parallel, quietly=TRUE)
  # if(exists("clusterLocal")) stopCluster(clusterLocal) #in case already started.
  clusterLocal <- makeCluster(c(rep("localhost",parallelThreads)),type="SOCK")
  registerDoParallel(clusterLocal)

  ANOVAoutList<-list()
  caseSubset="ALL"

  data = as.data.frame(cbind(as.character(colnames(cleanDat)), as.character(Grouping), t(cleanDat)))
  colnames(data)[1:2]<-c("CODE","SampleType")
  #test run gets column headers for output
  i=3
  aov<-aov(as.double(rnorm(n=nrow(data)))~SampleType, data=data)
  anovaresult<-anova(aov)
  tuk <- TukeyHSD(aov)
  tukresult1<-data.frame(tuk$SampleType) #With random (rnorm fn output), NO PAIRWISE COMPARISONS ARE MISSING--this is the template for comparison columns in ANOVAoutList[[caseSubset]]
  j=length(rownames(tukresult1))
  comparisonList<-rownames(tukresult1)

  line = c(paste("Protein", "F-Value", "Pr(>F)", sep=","))
  for (a in 1:length(comparisonList)) {
    line=c(paste(line,comparisonList[a],sep=","))
  }
  for (a in 1:length(comparisonList)) {
    line=c(paste(line,paste0("diff ",comparisonList[a]),sep=","))
  }

## Fast parallel code
  dataclipped <- data[, 3:ncol(data)] # as.matrix(as.numeric(data[,3:ncol(data)]))
  SampleType <- as.character(data$SampleType)
  if (!Tunequal) {
    parallel::clusterExport(cl=clusterLocal, list("SampleType","tukresult1","fallbackIfSmallTukeyP"), envir=environment())  # below helperfn not running in .GlobalEnv
    parts <- splitIndices(ncol(dataclipped), length(clusterLocal))
    dataclippedParts <- lapply(parts, function(i) dataclipped[,i,drop=FALSE])
#    resParts <- parallel::parApply(cl=clusterLocal, dataclippedParts, 2, function(x) {  #parApply call from within function does not parallelize.
    resParts <- parallel::clusterApply(cl=clusterLocal, dataclippedParts, helperfn)
    ANOVAoutList[[caseSubset]] <- do.call(cbind,resParts)
    ANOVAoutList[[caseSubset]] <- t(ANOVAoutList[[caseSubset]])
    countZeroFallbacks=sum(ANOVAoutList[[caseSubset]][,ncol(ANOVAoutList[[caseSubset]])])
    if (fallbackIfSmallTukeyP) { cat(paste0("\n...Tukey p<10^-8.5 Fallback calculations using Bonferroni corrected T test: ", countZeroFallbacks, " [",signif(countZeroFallbacks/(nrow(tukresult1)*nrow(ANOVAoutList[[caseSubset]]))*100,2),"%]\n")) } else { cat(paste0(countZeroFallbacks," [",signif(countZeroFallbacks/(nrow(tukresult1)*nrow(ANOVAoutList[[caseSubset]]))*100,2),"%] of all p values are below reliable calculation threshold of 10^-8.5.\nIf you see a ceiling effect in volcanoes, consider setting fallbackIfSmallTukeyP=TRUE.\nNote that volcano fallback for Tukey p=0 will be underestimated (and -log(p) displayed will be overestimated).\n")) }
    ANOVAoutList[[caseSubset]] <- ANOVAoutList[[caseSubset]][,-ncol(ANOVAoutList[[caseSubset]])]
    ANOVAcols <- as.vector(data.frame(do.call("rbind", strsplit(as.character(line), "[,]"))))
    ANOVAcols <- ANOVAcols[2:length(ANOVAcols)]
  }
  if (length(unique(SampleType))==2 & !Tunequal) { #for single pairwise comparison, essentially Pr(>F) from ANOVA is equivalent to T-test result (except 1 vs 2 non-NA measurements are allowed)
    ANOVAoutList[[caseSubset]][,3] <- p.adjust(ANOVAoutList[[caseSubset]][,2],method=twoGroupCorrMethod) #get BH FDR for ANOVA/T-test (equivalent) p values
    ANOVAoutList[[caseSubset]][,2:3]<-ANOVAoutList[[caseSubset]][,c(3,2)]
    ANOVAcols[2] <- paste0("FDR (",twoGroupCorrMethod,")")
  }
  if (Tunequal) {
    if (!length(unique(SampleType))==2) {
      stop("\nNumber of comparison groups is not 2. T test (unequal variance), option Tunequal=TRUE cannot be used!\n\n")
    } else {
      parallel::clusterExport(cl=clusterLocal, list("SampleType","tukresult1"), envir=environment())  # below helperfn not running in .GlobalEnv
      parts <- splitIndices(ncol(dataclipped), length(clusterLocal))
      dataclippedParts <- lapply(parts, function(i) dataclipped[,i,drop=FALSE])
      resParts <- parallel::clusterApply(cl=clusterLocal, dataclippedParts, TunequalHelperfn)
      ANOVAoutList[[caseSubset]] <- do.call(cbind,resParts)
      ANOVAoutList[[caseSubset]] <- t(ANOVAoutList[[caseSubset]])
      ANOVAoutList[[caseSubset]][,2] <- p.adjust(ANOVAoutList[[caseSubset]][,3],method=twoGroupCorrMethod) #get BH FDR for ANOVA/T-test (equivalent) p values
      ANOVAcols <- as.vector(data.frame(do.call("rbind", strsplit(as.character(line), "[,]"))))
      ANOVAcols <- ANOVAcols[2:length(ANOVAcols)]
      ANOVAcols[1] <- paste0("Tstatistic")
      ANOVAcols[2] <- paste0("FDR (",twoGroupCorrMethod,")")
    }
  }
  colnames(ANOVAoutList[[caseSubset]]) <- ANOVAcols
  ANOVAoutList[[caseSubset]] <- as.data.frame(ANOVAoutList[[caseSubset]])

  numComp=(ncol(ANOVAoutList[[caseSubset]])-2)/2

  ## Add network module colors, if they can be found for all output rows/gene products.
  if (length(NETcolors)==nrow(ANOVAoutList[[caseSubset]])) { ANOVAoutList[[caseSubset]]$NETcolors <- NETcolors }

  #*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
  if(outputCSV)  write.csv(ANOVAoutList[[caseSubset]],file=paste0(outFilePrefix,"ANOVA_diffEx-",caseSubset,outFileSuffix,".csv"))
  #*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+

  stopCluster(clusterLocal)
  return(ANOVAoutList[[caseSubset]])
} # end function parANOVA.dex


## Per solution to spur parallel use from within function (worked with parApply above only outside of function):  https://stackoverflow.com/questions/43269142/r-parlapply-not-parallel/45447830#45447830
helperfn <- function(xx) {
    apply(xx,2,function(x) {
    if(length(unique(SampleType[which(!is.na(x))]))<2) { #<2, handles "only one contrast level..." 
      x=data.frame(x=rep(1,length(x)), SampleType = SampleType)  #handles too many missing values
    } else {
      x <- data.frame(x = as.double(x), SampleType = SampleType) #as.double(x) instead of x corrects:   Error in lm.fit(x, y,... NA/NaN/Inf in 'y'
    }
    aov <- aov(as.double(x) ~ SampleType, data = x)
    anovaresult <- anova(aov)
    tuk <- TukeyHSD(aov)
    tukresult <- data.frame(tuk$SampleType)
    if (!length(rownames(tukresult)) == length(rownames(tukresult1))) tukresult <- tukresult[match(rownames(tukresult1), rownames(tukresult)), ]
    if (fallbackIfSmallTukeyP) {
      # Fallback to Bonferroni correction of t.test output if Tukey p returned as <10^-8.5 (complex Tukey approximation integral calculation has a ceiling effect or zero value output, beyond 10^-8.5) per https://stackoverflow.com/questions/16470404/tukeyhsd-adjusted-p-value-is-0-0000000
      zeroFallbacks=length(which(tukresult[,"p.adj"]<10^-8.5))
      if (zeroFallbacks>0) for (comp in which(tukresult[,"p.adj"]<10^-8.5)) {
        this.comp=rownames(tukresult)[comp]
        grp1=gsub("^(.*)\\-.*","\\1",this.comp)
        grp2=gsub("^.*\\-(.*)$","\\1",this.comp)
        lowTuk.p.estimate <- tryCatch( p.adjust(t.test(x[SampleType==grp1,"x"],x[SampleType==grp2,"x"],alternative="two.sided",var.equal=FALSE)$p.value, method="bonferroni",n=nrow(tukresult)), error=function(e) c(1) )
        tukresult[comp,"p.adj"] <- lowTuk.p.estimate
      }
    }

    c(anovaresult$F[1], anovaresult$Pr[1], as.vector(tukresult[, "p.adj"]), as.vector(tukresult[, "diff"]), zeroFallbacks)

    })  #close apply()

#  }, chunk.size = round(ncol(dataclipped)/parallelThreads,0)) #end parApply; this failed to parallelize
  } # end function helperfn


TunequalHelperfn  <- function(xx) {
    apply(xx,2,function(x) {
      x <- data.frame(x = as.double(x), SampleType = SampleType) #as.double(x) instead of x corrects:   Error in lm.fit(x, y,... NA/NaN/Inf in 'y'
      twoGroups=sort(unique(SampleType))
      tt <- tryCatch( t.test(as.double(x$x[x$SampleType==twoGroups[1]]),as.double(x$x[x$SampleType==twoGroups[2]]), alternative="two.sided", var.equal=FALSE), error=function(e) data.frame(statistic=0, p.value=1) )
      log2FC=mean(as.double(x$x[x$SampleType==twoGroups[2]]),na.rm=TRUE) - mean(as.double(x$x[x$SampleType==twoGroups[1]]),na.rm=TRUE)
      if (is.nan(log2FC)) { log2FC<-NA; tt$statistic[1]<-NA; }
      c(tt$statistic[1], tt$p.value[1], tt$p.value[1], log2FC)  #second alphabetically sorted group minus first for column 4
    })  #close apply()

  } # end function TunequalHelperfn

                               




################################
## Alternative to ANOVA stats ##
################################

trait.corStat <- function(dummyVar="",
## All parameters are set as variables of the global environment, or if not set, fallback to defaults in the function.
#                         cleanDat.=cleanDat,                             # cleanDat, data table of 'cleaned' relative abundance with rows (genes or proteins) and columns (samples)
#                         NETcolors=net$colors,                           # a vector containing module color assignments
#                                                                         # length of vector must be equal to number of rows in cleanDat.
#                                                                         # Colors are only added to table output if available; function tries to use NETcolors variable , with fallback to net$colors
#                         outputCSV=TRUE,                                 # Output Full Table of statistics?  TRUE/FALSE
#                         outFilePrefix="4",                              # typically "4", or step # in pipeline being run
#                         outFileSuffix=FileBaseName,                     # A description of the project, used as a filename suffix
## Specific to this function:
#                         cor.traits=c("GRGR","GAGA","GAGR")              # column name(s) of numericMeta (traits, sample metadata); columns to be used for correlation to cleanDat abundance values
#                         filter.trait="Day"                              # a single column name of numericMeta to use for filtering/subsetting samples to correlate to cor.traits
#                         filter.trait.subsets=c("ALL","Day5","Day10","Day15") # value(s) in filter.trait column to use for subsetting samples when correlating cleanDat rows to cor.trait(s)
#                         corFn="bicor"                                   #'spearman', 'kendall', or 'bicor'; anything else will cause Pearson (cor) to be used
                         env=.GlobalEnv) {

############################################
## PARAMETER VARIABLES CHECKED / DEFAULTS ##
############################################
  if (!exists("cleanDat")) stop("cleanDat variable must exist, holding gene product (rows) X sample (columns) data in the form of log2(relative abundance).")
  if (!exists("NETcolors")) if(exists("net")) { if ("colors" %in% names(net)) { NETcolors=net$colors } else { NETcolors=c() } } else { NETcolors=c() }
  if (!length(NETcolors)==nrow(cleanDat)) { cat("- Network color assignment vector not supplied or not of length in rows of cleanDat; will not be included in output table and data frame.\n") }
  if (!exists("outputCSV")) { outputCSV=TRUE } else { if(!is.logical(outputCSV)) { outputCSV=TRUE } }
  if (!exists("outFilePrefix")) { outFilePrefix="" } else { if (nchar(outFilePrefix)>0) outFilePrefix=paste0(outFilePrefix,".") }
  if (!exists("outFileSuffix")) { if (exists("FileBaseName")) { outFileSuffix=paste0("-",FileBaseName) } else { outFileSuffix="-unspecified_study" }} else { if (nchar(outFileSuffix)>0) outFileSuffix=paste0("-",outFileSuffix) }

  if (!length(cor.traits)==length(which(cor.traits %in% colnames(numericMeta)))) stop("Not all columns specified in variable cor.traits were found in traits data (numericMeta data frame).\n")
  if (!length(which(filter.trait %in% colnames(numericMeta)))==1) stop("The column of trait data specified by variable filter.trait not found, or more than one filter.trait was specified.\n")
  if (!exists("filter.trait.subsets")) filter.trait.subsets=="ALL"
  filter.trait.subsets<-gsub("^[Aa][Ll][Ll]$","ALL",filter.trait.subsets)
  if ("ALL" %in% filter.trait.subsets) { if(length(filter.trait.subsets)>1) { filter.trait.subsets.notALL=filter.trait.subsets[which(!filter.trait.subsets=="ALL")] } else { filter.trait.subsets.notALL=c() }} else { filter.trait.subsets.notALL=filter.trait.subsets }
  if (length(filter.trait.subsets.notALL)>0) if (!length(filter.trait.subsets.notALL)==length(which(filter.trait.subsets.notALL %in% names(table(numericMeta[,filter.trait])))) | min(table(numericMeta[,filter.trait])[filter.trait.subsets.notALL])<3) stop(paste0("Not all values specified in variable filter.trait.subsets found in numericMeta column ",filter.trait," or less than 3 total values (samples) found for such a value.\n"))
  if (!exists("corFn")) corFn="cor"
  if (corFn=="bicor") { this.corFn="bicor" } else if (corFn=="kendall") { cat("- Using Kendall's Tau correlation from psych::corr.test(). And enforcing n>2 pairwise complete correlation.\n"); this.corFn="kendall"; } else if (corFn=="spearman") { cat("- Using Spearman's rank correlation from psych::corr.test(). And enforcing n>2 pairwise complete correlation.\n"); this.corFn="spearman"; } else { cat("- Using Pearson correlation (corFn not set or not='bicor' or 'kendall').\n"); this.corFn="cor"; }

  CORstatList<-list()
  for(cor.trait in cor.traits) {
    for(filter.trait.value in filter.trait.subsets) {
        subset.idx= if(!filter.trait.value=="ALL") { which(numericMeta[,filter.trait]==filter.trait.value) } else { c(1:nrow(numericMeta)) }
        rawStats= if(this.corFn=="bicor") { apply(cleanDat,1,function(x) WGCNA::bicorAndPvalue(x[subset.idx],numericMeta[subset.idx,cor.trait],use="pairwise.complete.obs",alternative="two.sided")) } else if(this.corFn=="kendall" | this.corFn=="spearman") {
                                           apply(cleanDat,1,function(x) {
                                           x.pre<-x[subset.idx]
                                           y.pre<-numericMeta[subset.idx,cor.trait]
                                           keepComps=intersect(which(!is.na(x.pre)),which(!is.na(y.pre)))
                                           x.noNApairwise<-x.pre[keepComps]
                                           y.noNApairwise<-y.pre[keepComps]
                                           if(length(keepComps)>2) { psych::corr.test(x.noNApairwise, y.noNApairwise,use="pairwise.complete.obs",method=this.corFn) } else {
                                             list(r=0,p=1,se=NA,t=NA,n=0) }
                                           }) } else {
                                           apply(cleanDat,1,function(x) WGCNA::corAndPvalue(x[subset.idx],numericMeta[subset.idx,cor.trait],use="pairwise.complete.obs",alternative="two.sided")) }
      CORstatList[[paste0(gsub("\\.","_",cor.trait),".",filter.trait.value)]]=if(this.corFn=="kendall" | this.corFn=="spearman") {
        as.data.frame(matrix( unlist(lapply(rawStats,function(x) c(x["r"],x["p"],x["se"],x["t"],x["n"]))), nrow=nrow(cleanDat),ncol=5,byrow=TRUE) ) } else {
        as.data.frame(matrix( unlist(lapply(rawStats,function(x) c(x[this.corFn],x["p"],x["Z"],x["t"],x["nObs"]))), nrow=nrow(cleanDat),ncol=5,byrow=TRUE) )
      }
      this.list=paste0(gsub("\\.","_",cor.trait),".",filter.trait.value)  #cor.trait has '.' changed to '_' for plotVolc() run with corVolc=TRUE to parse column names consistently & correctly
      if (this.corFn=="kendall" | this.corFn=="spearman") {
        colnames(CORstatList[[this.list]])<-c("cor","p","SE","t","nObs")
      } else {
        colnames(CORstatList[[this.list]])<-c(this.corFn,"p","Z","t","nObs")
      }
      rownames(CORstatList[[this.list]])<-rownames(cleanDat)
      cat(paste0("Finished CORstatList[[",this.list,"]]\n"))
    }
  }

  if(!this.corFn=="bicor") { this.corlabel="cor" } else { this.corlabel="bicor" }

  CORout<-matrix(NA,nrow=nrow(cleanDat),ncol=2)
  colnames(CORout)<-c("Unused.Col","Unused.Col2")
  for(columnX in c("p",this.corlabel)) {
    this.list<-lapply(CORstatList,function(statTable) statTable[,columnX])
    #this.mat<-matrix(NA,nrow=nrow(cleanDat),ncol=0)
    #for(this.corr in names(this.list)) this.mat<-cbind(this.mat,this.list[[this.corr]])
    this.mat<-matrix(unlist(lapply(CORstatList,function(statTable) statTable[,columnX])),nrow=nrow(cleanDat),byrow=FALSE)
    #colnames(this.mat)=paste0(columnX," ",colnames(this.mat))
    colnames(this.mat)=paste0(columnX," ",names(CORstatList))
    CORout<- cbind(CORout, this.mat)
  }
  rownames(CORout)<-rownames(cleanDat)
  CORout<-as.data.frame(CORout)
  #CORout$NETcolors=net$colors

  numComp=(ncol(CORout)-2)/2

  ## Add network module colors, if they can be found for all output rows/gene products.
  if (length(NETcolors)==nrow(CORout)) { CORout$NETcolors <- NETcolors }

  #*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
  if(outputCSV)  write.csv(CORout,file=paste0(outFilePrefix,"CORstats-",outFileSuffix,".csv"))
  #*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+

  #assign("CORout",CORout, envir=.GlobalEnv)
  cat(" Correlation p + R table calculations complete. If you want to use the table with plotVolc(), set the variable corVolc=TRUE and use variable CORout to store the table generated.\n\n")
  return(CORout)
} # end function trait.corStat




##############################
## Volcano Plotter function ##
##############################
## 
## requires ggplot2, plotly packages.
##

plotVolc<- function(dummyVar="",
## All parameters are set as variables of the global environment, or if not set, fallback to defaults in the function.
#                    FCmin=0,                     # 0.25 for 25%, 0 for no threshold (vertical minimum FC threshold dashed lines)
#                    selectComps=selectComps,     # "ALL" for volcano output(s) on all pairwise comparisons in ANOVAout
#                    flip=c(),                    # p value column numbers in which to swap denominator of pair for x axis range (gene products high in denominator, will be on left)
#                    signifP=0.05,                # p value threshold for counting Differential Expression points
#                    useNETcolors=TRUE,           # use module colors saved to ANOVAout, if available; otherwise specify downColor upColor, and NCcolor
#                    downColor="royalblue",       # significant points above/beyond thresholds on the upper left are this color if useNETcolors=FALSE
#                    upColor="red",               # significant points above/beyond thresholds on the upper right are this color if useNETcolors=FALSE
#                    NCcolor="grey",              # points not significant are this color if useNETcolors=FALSE
#                    splitColors=FALSE,           # create a separate volcano plot(s) for each color in an outputfigs/splitVolcanoes subfolder (folder created if it does not exist)
#                    highlightGeneProducts=c(),   # c("APP|P05067","MAPT|P10636","APOE|P02649") ; a list of uniqueID rownames to highlight as larger gold points. If symbolsOnly=TRUE, this can be a list of symbols, like c("APP","SMOC1","MAPT")
#                    labelHighlighted=FALSE,      # if true, highlighted spots get text labels with their rownames from ANOVAout
#                    symbolsOnly=FALSE,           # for mouse-over HTML plots and the above highlight callouts, consider only displaying and using official gene symbol from first part of UniqueID rownames of ANOVAout.
#                    labelTop=0,                  # maximum p below which to label all points in the PDF output; OR an integer number of top ranked most significant points to label
#                    labelSize=4.5,               # text label font size, if any labels are found (when labelHighlighted=TRUE or labelTop>0)
#                    sameScale=FALSE,             # When multiple plots are drawn, should they all be the same scale with min and max x and y ranges?
#                    centerYaxis=FALSE,           # Creates plot(s) with x=0 (y-axis) centered between max absolute value x limits; cooperates with sameScale=TRUE to set balanced x limits to max of all plots if both TRUE
#                    HTMLout=TRUE,                # output interactive HTML copies that can be opened in browser. Requires plotly package.
#                    outFilePrefix="4",           # typically "4", or step # in pipeline being run
#                    outFileSuffix=FileBaseName,  # A description of the project, used as a filename suffix
#                    outputfigs=getwd(),          # "drive:/folder/to/location" to save output PDFs and possible HTML files

#                    corVolc=TRUE                 # Flag - correlation p values and R values in data frame CORout will be used, instead of ANOVA p values and differences of group means in ANOVAout
                    env=.GlobalEnv) {

############################################
## PARAMETER VARIABLES CHECKED / DEFAULTS ##
############################################

require(ggplot2,quietly=TRUE)

if(!exists("corVolc")) { corVolc=FALSE }
if(corVolc & exists("CORout")) { cat("- corVolc=TRUE so plotting volcanoes using correlation statistics in the table stored in variable CORout.\n"); ANOVAout<-CORout; }
if(corVolc & !exists("CORout")) if (!length(dummyVar)==1) { cat("- corVolc=TRUE, but CORout not in memory, using input provided to this function.\n"); ANOVAout=as.data.frame(dummyVar); } else { stop("Variable CORout not found or no input was provided.\nPlease run trait.corStat() function first, and save output to CORout variable or pass its output to this function.\n\n") }

if (!exists("ANOVAout")) if (!length(dummyVar)==1) { cat("- ANOVAout not in memory, using input provided to this function.\n"); ANOVAout=as.data.frame(dummyVar); } else { stop("Variable ANOVAout not found or no input was provided.\nPlease run parANOVA.dex() function first, and save output to ANOVAout variable or pass its output to this function.\n\n") }
if (!ncol(ANOVAout)>3) stop("- Input or ANOVAout variable contents are not a data (frame) with at least 4 columns. It is not valid output from the parANOVA.dex() or traits.corStat() functions.\n  Please generate input table first.\n\n")
	
numberOfNonComparisonColumns=length(colnames(ANOVAout)) - if(!corVolc) { length(which(grepl("^diff ",colnames(ANOVAout))))*2 } else { length(which(grepl("^p ",colnames(ANOVAout))))*2 }
numComp <- (length(colnames(ANOVAout)) - numberOfNonComparisonColumns) / 2 # of columns separating comparisons from matched column of log2(diffs), i.e. # of comparisons
	
if (!exists("selectComps")) { cat("- No comparison p value columns selected in selectComps. Using ALL comparisons.\n"); selectComps="ALL"; }
if (max(selectComps)>numComp+2 | min(selectComps)<3) {
  cat(" - selectComps may not reference valid integer p value column indexes of ANOVAout (or CORout).\n   Output will be for all comparisons or correlation(s).\n")
  selectComps="ALL"
}	
if (selectComps[1]=="ALL" | selectComps[1]=="all" | selectComps[1]=="All") selectComps=c(3:(numComp+2))
selectComps=as.integer(selectComps)

if(corVolc & exists("flip")) { cat("- Plotting correlation volcano(es). Variable flip will be ignored so positive correlations are to the right of x=0, consistently.\n"); flip=c(); }
if(!exists("flip")) { cat("- No comparisons selected for flipping numerator and denominator. Variable flip=c().\n"); flip=c(); }
# What columns (column numbers) of ANOVAout (T-test or Tukey p values) do we plot volcanoes for? 
testIndexMasterList <- selectComps

if(corVolc & exists("FCmin")) { cat("- Plotting correlation R and p in volcanoes. Variable FCmin will be ignored.\n"); FCmin=0; }
if (!exists("FCmin")) { cat("- No minimum fold change set in FCmin. Using 0 minimum change.\n"); FCmin=0; }

# log2(1) means NO change minimum to be counted in the volcano bookends; log2(1.25) for 25% FC min.
cutoff <- log2(1+FCmin)

# p value cutoff for Volcano significant hit counting; dashed line at -log10(sigCutoff)
if (!exists("signifP")) { cat("- No minimum significant p value threshold specified. Using signifP < 0.05.\n"); signifP=0.05; }
sigCutoff <- signifP

if (!exists("useNETcolors")) { if ("NETcolors" %in% colnames(ANOVAout)) { cat("- useNETcolors not set. NETcolors found in ANOVAout, so we will color points accordingly.\n"); useNETcolors=TRUE; } else { cat("- useNETcolors not set or not TRUE/FALSE. NETcolors not found in ANOVAout, so 3 color volcano(es) will be drawn.\n"); useNETcolors=FALSE; }}
if (!is.logical(useNETcolors)) { if ("NETcolors" %in% colnames(ANOVAout)) { cat("- useNETcolors not TRUE/FALSE. NETcolors found in ANOVAout, so we will color points accordingly.\n"); useNETcolors=TRUE; } else { cat("- useNETcolors not TRUE/FALSE. NETcolors not found in ANOVAout, so 3 color volcano(es) will be drawn.\n"); useNETcolors=FALSE; }}

if (!exists("NETcolors") & useNETcolors) {
  if ("NETcolors" %in% colnames(ANOVAout)) {
    cat("- NETcolors variable not set. But, NETcolors was found in the input ANOVA table. Using this.\n")
    NETcolors=ANOVAout$NETcolors
  } else {
    if ("colors" %in% names(net)) {
      cat("- NETcolors found as net$colors. Using this.\n")
      NETcolors=net$colors
    } else {
      stop("\nYou specified to use NETcolors, but I could not find a 'NETcolors' column in ANOVAout or the colors slot of the WGCNA::blockwiseModules() output 'net' list.\n\n")
    }
  }
}

if (!exists("downColor")) downColor="royalblue"
if (!exists("upColor")) upColor="red"
if (!exists("NCcolor")) NCcolor="grey"

if (!exists("splitColors")) splitColors=FALSE
if (!exists("symbolsOnly")) symbolsOnly=FALSE
if (!exists("HTMLout")) HTMLout=TRUE
if (!exists("outFilePrefix")) { outFilePrefix="" } else { if (nchar(outFilePrefix)>0) outFilePrefix=paste0(outFilePrefix,".") }
if (!exists("outFileSuffix")) { if (exists("FileBaseName")) { outFileSuffix=paste0("-",FileBaseName) } else { outFileSuffix="-unspecified_study" }} else { if (nchar(outFileSuffix)>0) outFileSuffix=paste0("-",outFileSuffix) }
if (!exists("outputfigs")) { cat(paste0("- Variable outputfigs not specified. Saving volcano plots to ",getwd()," .\n")); outputfigs=getwd(); }
if (!dir.exists(outputfigs)) { cat(paste0("- Directory ",outputfigs," not found. Attempting to create it.\n")); dir.create(outputfigs); }

# Gene products to highlight
if (!exists("highlightGeneProducts")) highlightGeneProducts=c()
BIGspots <- highlightGeneProducts

cutoff=log2(1+FCmin)
# shows what your cutoff for log2(FC) calculates as
if(!corVolc) { 
	cat(paste0("\n...Thresholds:\n   [x] Applying a ", FCmin*100,"% minimum fold change threshold at + and - x=", signif(cutoff,2)," .\n"))
	cat(paste0("   [y] with minimum significance cutoff p < ", sigCutoff,", equivalent to -log10(p) y=", round(-log10(sigCutoff),2)," .\n\n"))
} else {
	cat(paste0("\n...Thresholds:\n   [x] No minimum R value is needed or used with correlation statistics.\n"))
	cat(paste0("   [y] minimum Student's significance of correlation p < ", sigCutoff,", equivalent to -log10(p) y=", round(-log10(sigCutoff),2)," .\n\n"))
}

n <- nrow(ANOVAout)
dexComps <- list()
iter <- length(testIndexMasterList) + 1
comparisonIDs <- data.frame(dfVariable = rep(NA, length(testIndexMasterList)), Comparison = rep(NA, length(testIndexMasterList)))
for (i in testIndexMasterList) {
  iter <- iter - 1
  # dexRows<-which(ANOVAout[,i]<sigCutoff) #choose rows where the DEX p<sigCutoff
  if(!corVolc) {
    comparisonIDs[iter, ] <- as.vector(c(paste0("dexTargets.", gsub("-", ".", colnames(ANOVAout)[i])), paste0(as.character(gsub("-", " vs ", colnames(ANOVAout)[i])))))
  } else { 
    comparisonIDs[iter, ] <- as.vector(c(paste0("dexTargets.", gsub("^p ", "", colnames(ANOVAout)[i])), paste0(as.character(gsub("^(.*)[' '](.*)\\.(.*)$", "\\1 (\\2 in \\3 samples)", colnames(ANOVAout)[i+numComp])))))
  }
  dexComps[[comparisonIDs[iter, 1]]] <- ANOVAout
  if (!is.na(match(i, flip))) {
    dexComps[[comparisonIDs[iter, 1]]][, i + numComp] <- -1 * as.numeric(dexComps[[comparisonIDs[iter, 1]]][, i + numComp])
    comparisonIDs[iter, 2] <- gsub("(*.*) vs (*.*)", "\\2 vs \\1", comparisonIDs[iter, 2]) # flip label "vs" in comParisonIDs$Comparison[iter]
  }
}
comparisonIDs # list element names and Logical comparisons for those retrievable Dex measurements in the list elements
ls(dexComps) # list elements are dataframes with the DEX entries for that comparison


## volcano plots with (module or other) colors, SPLITTABLE on COLOR

pointColorsVectorListForPlots <- list()
volcListModColorsWeb <- list()
volcListModColors <- list()
dfListModColors <- list()
xRange<-yRange <- list()

iter <- length(testIndexMasterList) + 1
for (testIndex in testIndexMasterList) {
  iter <- iter - 1
  df <- eval(parse(text = "dexComps[[comparisonIDs$dfVariable[iter]]]"))
  cat(paste0("\rProcessing ANOVA column ", testIndex, " (", comparisonIDs$Comparison[iter], ") for volcano ...   "))

  # correct 0 Tukey pValues to ANOVA p (in column 2); it's better than taking -log10 of 0 in the next step
  if(!corVolc) {
    if (length(which(df[, testIndex] == 0))>0) print(paste0("Found ",length(which(df[, testIndex] == 0))," p=0 value(s) for this comparison! Substituted -log10(p) using one-way ANOVA overall p value (Consider fallback=TRUE in ANOVAout calculation.)"))
    df[which(df[, testIndex] == 0), testIndex] <- as.numeric(df[which(df[, testIndex] == 0), 2])
  } else {
    if (length(which(df[, testIndex] == 0))>0) print(paste0("Found ",length(which(df[, testIndex] == 0))," p=0 value(s) for this comparison! Substituted -log10(p) value of 200 for these correlations."))
    df[which(df[, testIndex] == 0), testIndex] <- 1e-200
  }    
 
  ## Check if ANOVA pVal is Significant and above FC cutoff defined above. Thresholds are used to set volcano point colors
  df[, testIndex][is.na(df[, testIndex])] <- 1 # p=0.9999 instead of NA
  df[, testIndex + numComp][is.na(df[, testIndex + numComp])] <- 0 # log2(difference)=0 instead of NA
  
  df$negLogP <- -log10(as.numeric(df[, testIndex]))
  
  df$threshold1 <- as.numeric(rep(3, n))
  ## Significance (y) and x Threshold categorization of each gene product
  for (i in 1:n) {
    if (as.numeric(df[i, testIndex + numComp]) > cutoff & as.numeric(df[i, testIndex]) < sigCutoff) {
      df$threshold1[i] <- 1
    } else {
      if (as.numeric(df[i, testIndex + numComp]) < -cutoff & as.numeric(df[i, testIndex]) < sigCutoff) {
        df$threshold1[i] <- 2
      }
    }
  }
  df$threshold1 <- as.factor(df$threshold1)
  
  df$Symbol <- rownames(df)
  if (symbolsOnly) df$Symbol <- suppressWarnings(as.data.frame(do.call("rbind", strsplit( paste0(as.data.frame(do.call("rbind", strsplit(as.character(rownames(df)), "[|]")))[, 1],";"), "[;]")))[,1])
  
  ## Color Interesting Gene Product Spots DIFFERENTLY as 4th color if doing blue/red/green (no module colors) -- (4=gold1 below)
  if(useNETcolors) { 
    if(!"NETcolors" %in% colnames(df)) df$NETcolors <- NETcolors
    df$color1<-df$NETcolors
    df$threshold2 <- as.numeric(df$threshold1)
    df$threshold2[match(intersect(df$Symbol, BIGspots), df$Symbol)] <- 4 
  } else {
    df$color1 <- as.numeric(as.character(df$threshold1))
    df$color1[match(intersect(df$Symbol, BIGspots), df$Symbol)] <- 4
    df$threshold2 <- as.numeric(as.character(df$threshold1))
    df$threshold2[match(intersect(df$Symbol, BIGspots), df$Symbol)] <- 4 
  }
  df$color1 <- as.factor(df$color1)

  if(useNETcolors) {
    df$size1 <- as.numeric(df$threshold2)
  } else {
    df$size1 <- as.numeric(df$color1)
  }
  df$size1[df$size1 < 4] <- 2.5
  df$size1[df$size1 == 4] <- 6

  #df$color2: actual spot color as.character() ; df$color1 is a factorable dummy
  if(useNETcolors) {
    df$color2<-as.character(df$color1)
    df$color2[df$threshold2 == 4] <- "gold1" #for BIGspots
  } else {
    df$color2 <- as.numeric(as.character(df$color1))
    df$color2[df$color2 == 1] <- upColor
    df$color2[df$color2 == 2] <- downColor
    df$color2[df$color2 == 3] <- NCcolor
    df$color2[df$color2 == 4] <- "gold1" #for BIGspots
  }                              #gold1 is also the 435th, last WGCNA unique color

  #df$color3 is outline color, where outlined pch symbols (21) used
  df$color3 <- df$color2
  df$color3[df$color3 == "gold1"] <- "black" #for BIGspots outline
  df$pch <- as.numeric(df$threshold2)
  df$pch[df$pch < 4] <- 16 # unfilled circles (use color2)
  df$pch[df$pch == 4] <- 21 # filled, outlined circles (border uses color3)

  #put gold1 back to module color for fill of BIGspots if (useNETcolors)
  if (useNETcolors) { df$color2[df$color2=="gold1"] <- as.character(df$color1[df$color2=="gold1"]) }

  df <- df[order(df$size1, decreasing = FALSE), ] # puts larger dots on top (at bottom of df)
  

  #splitColors TRUE/FALSE: make one volcano with all colors (FALSE), or make volcanoes for each color (TRUE)
  # SPLIT DATA FRAME FOR VOLCANO PLOT BY COLORS (if multiple eachColorSplit items)
  df.AllColors <- df
  eachColorSplit <- if (splitColors) {
    unique(df.AllColors$NETcolors)
  } else {
    c("allcolors")
  }
  for (eachColor in eachColorSplit) {
    if (splitColors) {
      df.oneColor <- df.AllColors[which(df.AllColors$NETcolors == eachColor), ]
      df.oneColor$color1 <- factor(df.oneColor$NETcolors)
    } else {
      df.oneColor <- df.AllColors
    } # end if (splitColors)
    
    names(df.oneColor)[testIndex + numComp] <- "xdata" # x=as.numeric(df.oneColor[,testIndex+numComp])
    list_element <- paste(comparisonIDs$dfVariable[iter], eachColor, sep = ".") # colnames(df)[testIndex]
    pointColorsVectorListForPlots[[list_element]] <- data.frame(color1 = factor(as.integer(df.oneColor$color1)), color2 = as.character(df.oneColor$color2), color3 = as.character(df.oneColor$color3), size = as.numeric(df.oneColor$size), pch = as.numeric(df.oneColor$pch)) #*** df.oneColor$NETcolors
    
    downregulated.text=if(!corVolc) { "Downregulated" } else { "Negative Correlations" }
    upregulated.text=if(!corVolc) { "Upregulated" } else { "Positive Correlations" }

    xlab.expr=if(!corVolc) { as.expression(bquote("Difference, log"[2] ~ .(comparisonIDs$Comparison[iter]))) } else { as.expression(comparisonIDs$Comparison[iter]) }
    xlab.text=if(!corVolc) { paste0("Difference, log2 ", comparisonIDs$Comparison[iter]) } else { comparisonIDs$Comparison[iter] }

    volcano1 <- ggplot(data = df.oneColor, aes(x = xdata, y = negLogP, color = color1, text = Symbol)) +
      # scale_colour_manual(values = unique(data.frame(col1=df.oneColor$color1,col2=df.oneColor$color2))[order(unique(data.frame(col1=df.oneColor$color1,col2=df.oneColor$color2))[,1]),2] )+ #THIS COLOR(S) IS LOOKED UP ACTIVELY DURING GRAPHICS OUTPUT IN THE VARIABLE, SO WE'VE USED A LIST ELEMENT THAT IS NEVER CHANGED
      geom_point(aes(fill = pointColorsVectorListForPlots[[list_element]][, "color3"]), alpha = 0.66, size = pointColorsVectorListForPlots[[list_element]]$size, pch = pointColorsVectorListForPlots[[list_element]][, "pch"], color = pointColorsVectorListForPlots[[list_element]][, "color3"]) +
      theme(legend.position = "none") +
      xlim(c(min(as.numeric(df.oneColor[, testIndex + numComp])), max(as.numeric(df.oneColor[, testIndex + numComp])))) + ylim(c(0, max(df.oneColor$negLogP))) +
      xlab(xlab.expr) + # colnames(df.oneColor)[testIndex]
      ylab(as.expression(bquote("-log"[10] ~ "p value"))) +
      theme(axis.title.x = element_text(size = rel(1.8), angle = 00)) +
      theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) +
      
      geom_hline(yintercept = -log10(sigCutoff), linetype = "dashed", color = "black", size = 1.2) +
       # annotate("text", label=paste0("y=",signif(-log10(sigCutoff),3),"; p=",sigCutoff), x=xRange[[list_element]][1], y=-log10(sigCutoff), vjust = -1, hjust=0) +
      geom_vline(xintercept = cutoff, linetype = "dashed", color = "black", size = 1.2) +
      geom_vline(xintercept = -cutoff, linetype = "dashed", color = "black", size = 1.2) +
      annotate("text", x = min(as.numeric(df.oneColor[, testIndex + numComp])) / 2, y = max(df.oneColor$negLogP) * .95, size = 5, label = paste0(downregulated.text,": ", bquote(.(length(which(as.numeric(as.character(df.oneColor$threshold1)) == 2)))))) +
      annotate("text", x = max(as.numeric(df.oneColor[, testIndex + numComp])) / 2, y = max(df.oneColor$negLogP) * .95, size = 5, label = paste0(upregulated.text,": ", bquote(.(length(which(as.numeric(as.character(df.oneColor$threshold1)) == 1)))))) +
      
      theme(
        # axis.text = element_text(size = 14),
        # legend.key = element_rect(fill = "navy"),
        # legend.background = element_rect(fill = "white"),
        # legend.position = c(0.14, 0.80),
        panel.grid.major = element_line(color = "darkgrey", linetype = "dashed"),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white")
      )
    
    
    # web version doesn't use as.expression! plotly fails with those, so we rebuild the volcano for the web; we end up using this version for both PDF and HTML.
    if(!exists("labelTop")) labelTop=0
    if(labelTop>0) {
      if(labelTop<1) {
        finalLabelCountUp=length(which(df.oneColor$negLogP>= -log10(labelTop) & df.oneColor$threshold1==1))
        finalLabelCountDown=length(which(df.oneColor$negLogP>= -log10(labelTop) & df.oneColor$threshold1==2))
      } else {
        finalLabelCount=as.integer(labelTop)
        finalLabelCountUp= if(length(which(df.oneColor$threshold1==1))<finalLabelCount) { length(which(df.oneColor$threshold1==1)) } else { finalLabelCount }
        finalLabelCountDown= if(length(which(df.oneColor$threshold1==2))<finalLabelCount) { length(which(df.oneColor$threshold1==2)) } else { finalLabelCount }
      }
      df.oneColor$label=rep("",length(df.oneColor$Symbol))
      if(finalLabelCountUp>0) df.oneColor$label[intersect(order(df.oneColor$negLogP,decreasing=TRUE), which(df.oneColor$threshold1==1))[1:finalLabelCountUp]] <- df.oneColor$Symbol[intersect(order(df.oneColor$negLogP,decreasing=TRUE), which(df.oneColor$threshold1==1))[1:finalLabelCountUp]]
      if(finalLabelCountDown>0) df.oneColor$label[intersect(order(df.oneColor$negLogP,decreasing=TRUE), which(df.oneColor$threshold1==2))[1:finalLabelCountDown]] <- df.oneColor$Symbol[intersect(order(df.oneColor$negLogP,decreasing=TRUE), which(df.oneColor$threshold1==2))[1:finalLabelCountDown]]

      require(ggrepel,quietly=TRUE)
    }

    if(!exists("labelHighlighted")) labelHighlighted=FALSE
    if(!exists("highlightGeneProducts")) highlightGeneProducts=c()
    if(labelHighlighted & length(highlightGeneProducts)>0) {
      if(!"label" %in% colnames(df.oneColor)) df.oneColor$label=rep("",length(df.oneColor$Symbol))
      df.oneColor$label[which(df.oneColor$Symbol %in% intersect(df.oneColor$Symbol,highlightGeneProducts))] <- df.oneColor$Symbol[which(df.oneColor$Symbol %in% intersect(df.oneColor$Symbol,highlightGeneProducts))]
    }      

    if(!exists("sameScale")) sameScale=FALSE
    if(!exists("centerYaxis")) centerYaxis=FALSE
    if(sameScale) {
      ANOVAout.all.negLogP<-t(apply(ANOVAout,1,function(x) { y=x[testIndexMasterList]; y[which(y==0)] <- x[2]; -log10(as.numeric(y)); }))
      yRange[[list_element]]=c(0,max(ANOVAout.all.negLogP,na.rm=TRUE)*1.01)
      if(!length(testIndexMasterList)==1) {
        ANOVAout.all.log2FC<-ANOVAout[,testIndexMasterList+numComp]
      } else {
        ANOVAout.all.log2FC=data.frame(col1=ANOVAout[,testIndexMasterList+numComp])
        rownames(ANOVAout.all.log2FC)<-rownames(ANOVAout)
        colnames(ANOVAout.all.log2FC)<-colnames(ANOVAout)[testIndexMasterList+numComp]
      }
      if(length(flip)>0) for (column in flip) if(!length(testIndexMasterList)==1) { ANOVAout.all.log2FC[,column-2] = ANOVAout.all.log2FC[,column-2]*(-1) } else { ANOVAout.all.log2FC=ANOVAout.all.log2FC*(-1) }  #***
      if(centerYaxis) {
        maxAbsX=max(abs(range(ANOVAout.all.log2FC,na.rm=TRUE)))
        xRange[[list_element]]=c(-maxAbsX,maxAbsX)
      } else {
        xRange[[list_element]]=range(ANOVAout.all.log2FC,na.rm=TRUE)
      }
    } else {
      xRange[[list_element]]=c(min(as.numeric(df.oneColor[, testIndex + numComp]),na.rm=TRUE), max(as.numeric(df.oneColor[, testIndex + numComp]),na.rm=TRUE))
      yRange[[list_element]]=c(0, max(df.oneColor$negLogP,na.rm=TRUE)*1.01)
      if(centerYaxis) {
        xRange[[list_element]]=c(-abs(max(xRange[[list_element]])),abs(max(xRange[[list_element]])) )
      }
    }

    if(labelTop>0 | "label" %in% colnames(df.oneColor)) {
      options(ggrepel.max.overlaps = Inf)  # Allowing unlimited overlaps improves compatibility of ggrepel labeling in RStudio.
      volcanoweb <- ggplot(data = df.oneColor, aes(x = xdata, y = negLogP, color = color1, text = Symbol, label = label))
    } else {
      volcanoweb <- ggplot(data = df.oneColor, aes(x = xdata, y = negLogP, color = color1, text = Symbol))
    }
    volcanoweb <- volcanoweb +
      scale_colour_manual(values = unique(data.frame(col1 = df.oneColor$color1, col2 = df.oneColor$color2))[order(unique(data.frame(col1 = df.oneColor$color1, col2 = df.oneColor$color2))[, 1]), 2]) + # THIS COLOR(S) IS LOOKED UP ACTIVELY BY PLOTLY IN THE VARIABLE, SO WE'VE USED A LIST ELEMENT THAT IS NEVER CHANGED
      # scale_y_continuous(breaks = seq(0, 8, by = 1))+
      # geom_point(aes(fill=pointColorsVectorListForPlots[[list_element]][,"color2"]), alpha=0.66, size=pointColorsVectorListForPlots[[list_element]]$size, pch=16, color=pointColorsVectorListForPlots[[list_element]][,"color3"]) + #pch=pointColorsVectorListForPlots[[list_element]][,"pch"] not using variable symbol types (21 for outlined circle only)
      geom_point(alpha = 0.66, size = pointColorsVectorListForPlots[[list_element]]$size, pch = 16) + # pch=pointColorsVectorListForPlots[[list_element]][,"pch"] just uses the higher pch code in the web render.
      theme(legend.position = "none") +
      xlim(xRange[[list_element]]) + ylim(yRange[[list_element]]) +
      xlab(xlab.text) +
      ylab(paste0("-log10 p value")) +
      theme(axis.title.x = element_text(size = rel(1.8), angle = 00)) +
      theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) +
      
      geom_hline(yintercept = -log10(sigCutoff), linetype = "dashed", color = "black", size = 1.2) +
       # annotate("text", label=paste0("y=",signif(-log10(sigCutoff),3),"; p=",sigCutoff), x=xRange[[list_element]][1], y=-log10(sigCutoff), vjust = -1, hjust=0) +
      geom_vline(xintercept = cutoff, linetype = "dashed", color = "black", size = 1.2) +
      geom_vline(xintercept = -cutoff, linetype = "dashed", color = "black", size = 1.2) +
      annotate("text", x = min(xRange[[list_element]]) / 2, y = max(yRange[[list_element]]) * .95, size = 5, label = paste0(downregulated.text,": ", bquote(.(length(which(as.numeric(as.character(df.oneColor$threshold1)) == 2)))))) +
      annotate("text", x = max(xRange[[list_element]]) / 2, y = max(yRange[[list_element]]) * .95, size = 5, label = paste0(upregulated.text,": ", bquote(.(length(which(as.numeric(as.character(df.oneColor$threshold1)) == 1)))))) +
      
      theme(
        # axis.text = element_text(size = 14),
        # legend.key = element_rect(fill = "navy"),
        # legend.background = element_rect(fill = "white"),
        # legend.position = c(0.14, 0.80),
        panel.grid.major = element_line(color = "darkgrey", linetype = "dashed"),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white")
      )
    
    if(labelTop>0 | "label" %in% colnames(df.oneColor)) {
      if(!exists("labelSize")) labelSize=4.5
      require(ggrepel,quietly=TRUE)
      volcanoweb <- volcanoweb + geom_label_repel(box.padding = 0.5,label.size=NA, fill=NA, size=labelSize)  #label.size is for border width; size is for font size; fill is for background (behind text) color.
    }

    volcListModColors[[list_element]] <- volcano1
    volcListModColorsWeb[[list_element]] <- volcanoweb
    
    print(volcano1) # prints to active output (separate page)
    rm(volcano1)
    rm(volcanoweb)
    dfListModColors[[list_element]] <- df.oneColor
  } # closes for(eachColor...
} # closes for(testIndex...


## Write Files
if(splitColors) { cat("\n- Creating SplitVolcano folder for all output files...\n"); dir.create(file.path(outputfigs, "/SplitVolcano/")); }

# Print to PDFs, one per color (per comparison, if multiple)
iter <- length(testIndexMasterList) + 1
for (testIndex in testIndexMasterList) {
  iter <- iter - 1
  cat(paste0("\rGenerating PDF volcano for ANOVAout column ", testIndex, " (", comparisonIDs$Comparison[iter], ") ...   "))
  for (eachColor in eachColorSplit) {
    list_element <- paste(comparisonIDs$dfVariable[iter], eachColor, sep = ".")
    df.oneColor <- dfListModColors[[list_element]]
    if(splitColors) {
      file <- paste0(outputfigs,"/SplitVolcano/",outFilePrefix,"Volcano_", eachColor,"-",gsub(" ", "_", comparisonIDs$Comparison[iter]),outFileSuffix,".pdf")
    } else {
      file <- paste0(outputfigs,"/",outFilePrefix,"Volcano_", eachColor,"-",gsub(" ", "_", comparisonIDs$Comparison[iter]),"-",outFileSuffix,".pdf")
    }
    pdf(file = file, colormodel="srgb", height = 8, width = 8)
    par(mfrow = c(1, 1))
    par(mar = c(6, 8.5, 3, 3))
    print(volcListModColorsWeb[[list_element]])
    #     volcListModColors[[list_element]])    #bigspots colors not handled correctly for PDF output of modColored spots:
    dev.off()
  }
}

## Print to html volcano plots (one per module colors and per comparison, if applicable)
if (HTMLout) {
suppressPackageStartupMessages(require(plotly, quietly=TRUE))

iter <- length(testIndexMasterList) + 1
for (testIndex in testIndexMasterList) {
  iter <- iter - 1
  cat(paste0("\rGenerating HTML volcano for ANOVAout column ", testIndex, " (", comparisonIDs$Comparison[iter], ") ...   "))
  for (eachColor in eachColorSplit) {
    list_element <- paste(comparisonIDs$dfVariable[iter], eachColor, sep = ".")
    df.oneColor <- dfListModColors[[list_element]]
    webPlot <- ggplotly(volcListModColorsWeb[[list_element]])
    if(splitColors) {
      tempfilename=paste0(outputfigs,"/SplitVolcano/",outFilePrefix,"HTMLvolcano_",eachColor,"-",gsub(" ", "_", comparisonIDs$Comparison[iter]),"-",outFileSuffix,".html")
    } else {
      tempfilename=paste0(outputfigs,"/",outFilePrefix,"HTMLvolcano_",eachColor,"-",gsub(" ", "_", comparisonIDs$Comparison[iter]),"-",outFileSuffix,".html")
    }
    htmlwidgets::saveWidget(webPlot, tempfilename, selfcontained = TRUE, libdir = "delete.me")
  }
}

} #end if (HTMLout)
cat("\n")

## These variables may be referenced by dependent, downstream pipeline functions like GOparallel, or stackedDEXbar, following creation of ANOVAout data frame and volcanoes, propagating selections.
assign("numComp",numComp, envir=.GlobalEnv)
assign("comparisonIDs",comparisonIDs, envir=.GlobalEnv)
assign("testIndexMasterList",testIndexMasterList, envir=.GlobalEnv)
assign("flip",flip, envir=.GlobalEnv)
assign("FCmin",FCmin, envir=.GlobalEnv)
assign("sigVolcCutoff",sigCutoff, envir=.GlobalEnv)
}



#######################################
## Stacked Barplot Function          ##
#######################################
## - for percents of module members  ##
##   reaching differential abundance ##
#######################################
## 
## requires ggplot2, Cairo, WGCNA, scales packages
##
## ...requires ANOVAout prior creation, and net$MEs and net$colors (latter 2 are blockwiseModules WGCNA function output to slots in the list 'net')
##

DEXpercentStacked <- function(dummyVar="",
## All parameters are set as variables of the global environment, or if not set, fallback to defaults in the function.
#                    selectComps=testIndexMasterList,     # can be set to "ALL" for all pairwise comparisons in ANOVAout, or column indices of ANOVAout where p values of pairwise comparisons are stored.
#                                                 # Defaults to the columns selected for volcanoes.   
#                    flip=c(),                    # p value column numbers in which to swap denominator of pair for x axis range (gene products high in denominator, will be on left)
#                                                 # Already was set for volcanoes and exported from plotVolc() function to global environment.
#                    signifP=0.05,                # p value threshold for counting Differential Abundance points
#                    FCmin=0,                     # 0.25 for 25%, 0 for no threshold (vertical minimum FC threshold dashed lines)
#                    outFilePrefix="4",           # typically "4", or step # in pipeline being run
#                    outFileSuffix=FileBaseName,  # A description of the project, used as a filename suffix
#                    outputfigs=getwd(),          # "drive:/folder/to/location" to save output PDFs and possible HTML files
                    env=.GlobalEnv) {

############################################
## PARAMETER VARIABLES CHECKED / DEFAULTS ##
############################################

suppressPackageStartupMessages(require(ggplot2,quietly=TRUE))
#suppressPackageStartupMessages(require(Cairo,quietly=TRUE))
suppressPackageStartupMessages(require(WGCNA,quietly=TRUE))

if(!exists("corVolc")) { corVolc=FALSE }
if(corVolc & exists("CORout")) { cat("- corVolc=TRUE so plotting volcanoes using correlation statistics in the table stored in variable CORout.\n"); ANOVAout<-CORout; }
if(corVolc & !exists("CORout")) if (!length(dummyVar)==1) { cat("- corVolc=TRUE, but CORout not in memory, using input provided to this function.\n"); ANOVAout=as.data.frame(dummyVar); } else { stop("Variable CORout not found or no input was provided.\nPlease run trait.corStat() function first, and save output to CORout variable or pass its output to this function.\n\n") }

if (!exists("ANOVAout")) if (!length(dummyVar)==1) { cat("- ANOVAout not in memory, using input provided to this function.\n"); ANOVAout=as.data.frame(dummyVar); } else { stop("Variable ANOVAout not found or no input was provided.\nPlease run parANOVA.dex() function first, and save output to ANOVAout variable or pass its output to this function.\n\n") }
if (!ncol(ANOVAout)>3) stop("\nInput or ANOVAout variable contents are not a data (frame) with at least 4 columns. It is not valid output from the parANOVA.dex() function.\n  Please run parANOVA.dex() first.\n\n")

numberOfNonComparisonColumns=length(colnames(ANOVAout)) - if(!corVolc) { length(which(grepl("^diff ",colnames(ANOVAout))))*2 } else { length(which(grepl("^p ",colnames(ANOVAout))))*2 }
numComp <- (length(colnames(ANOVAout)) - numberOfNonComparisonColumns) / 2 # of columns separating comparisons from matched column of log2(diffs), i.e. # of comparisons

if (!exists("NETcolors")) {
  if ("NETcolors" %in% colnames(ANOVAout)) {
    cat("- NETcolors variable not set. But, NETcolors was found in the input ANOVA table. Using this.\n")
    NETcolors=ANOVAout$NETcolors
  } else {
    if ("colors" %in% names(net)) {
      cat("- NETcolors found in net$colors...\n")
      NETcolors=net$colors
    } else {
      stop("\nI could not find a 'NETcolors' column in ANOVAout nor did I find a 'colors' slot of the WGCNA::blockwiseModules() output 'net' list.\nThis function only works if gene products in cleanDat have been clustered into modules.\n\n")
    }
  }
}
if (!length(NETcolors)==nrow(ANOVAout)) { stop("\nNETcolors vector length does not match number of rows in the input ANOVA table.\n"); }
if(!"NETcolors" %in% colnames(ANOVAout)) ANOVAout$NETcolors <- NETcolors
if(exists("testIndexMasterList")) { cat("- Found plotVolc function output variable testIndexMasterList to recall comparisons selected for selectComps. Using the following comparisons:\n"); print(data.frame('Comparisons'=colnames(ANOVAout)[testIndexMasterList])); selectComps=testIndexMasterList; }

if (!exists("selectComps")) { cat("- No comparison p value columns selected in selectComps. Using ALL comparisons.\n"); selectComps="ALL"; }
if (max(selectComps)>numComp+2 | min(selectComps)<3) {
  cat(" - selectComps may not reference valid integer p value column indexes of ANOVAout (or CORout).\n   Output will be for all comparisons or correlation(s).\n")
  selectComps="ALL"
}	
if (selectComps[1]=="ALL" | selectComps[1]=="all" | selectComps[1]=="All") selectComps=c(3:(numComp+2))
selectComps=as.integer(selectComps)

if(corVolc & exists("flip")) { cat("- Plotting using correlation volcano(es) stats. Variable flip will be ignored so positive correlations remain positive.\n"); flip=c(); }
if(!exists("flip")) { cat("- No comparisons selected for flipping numerator and denominator. Variable flip=c().\n"); flip=c(); }
# What columns (column numbers) of ANOVAout (T-test or Tukey p values) do we plot volcanoes for? 
testIndexMasterList <- selectComps

if(corVolc & exists("FCmin")) { cat("- Plotting significant correlation hit counts. Variable FCmin will be ignored.\n"); FCmin=0; }
if (!exists("FCmin")) { cat("- No minimum fold change set in FCmin. Using 0 minimum change.\n"); FCmin=0; }
# log2(1) means NO change minimum to be counted in the volcano bookends; log2(1.25) for 25% FC min.
cutoff <- log2(1+FCmin)

# p value cutoff for Volcano significant hit counting; dashed line at -log10(sigCutoff)
if (exists("sigVolcCutoff")) signifP=sigVolcCutoff
if (!exists("signifP")) { cat("- No minimum significant p value threshold specified. Using signifP < 0.05.\n"); signifP=0.05; }
sigCutoff <- signifP


if (!exists("outFilePrefix")) { outFilePrefix="" } else { if (nchar(outFilePrefix)>0) outFilePrefix=paste0(outFilePrefix,".") }
if (!exists("outFileSuffix")) { if (exists("FileBaseName")) { outFileSuffix=paste0("-",FileBaseName) } else { outFileSuffix="-unspecified_study" }} else { if (nchar(outFileSuffix)>0) outFileSuffix=paste0("-",outFileSuffix) }
if (!exists("outputfigs")) { cat(paste0("- Variable outputfigs not specified. Saving volcano plots to ",getwd()," .\n")); outputfigs=getwd(); }
if (!dir.exists(outputfigs)) { cat(paste0("- Directory ",outputfigs," not found. Attempting to create it.\n")); dir.create(outputfigs); }


#Human DEX Tables Built from ANOVAout dataframe
ANOVAin<-ANOVAout
masterDEXlookup<-selectComps

dexComps<-list()
iter=length(masterDEXlookup)+1
comparisonIDs <- data.frame(dfVariable=rep(NA,length(masterDEXlookup)),Comparison=rep(NA,length(masterDEXlookup)))
for (i in masterDEXlookup) {
 iter=iter-1;
 dexRows<-which(ANOVAin[,i]<0.05) #choose rows where the DEX p<0.05
 comparisonIDs[iter,] <- as.vector( c(paste0("dexTargets.",gsub("-",".",colnames(ANOVAin)[i])),paste0(as.character(gsub("-"," vs ",colnames(ANOVAin)[i])))))
 dexComps[[comparisonIDs[iter,1]]] <- ANOVAin[dexRows,]
 if(!is.na(match(i,flip))) { dexComps[[comparisonIDs[iter,1]]][,i+numComp] <- -1*as.numeric(dexComps[[comparisonIDs[iter,1]]][,i+numComp])
  comparisonIDs[iter,2]<-gsub("(*.*) vs (*.*)","\\2 vs \\1",comparisonIDs[iter,2]) #flip label "vs" in comParisonIDs$Comparison[iter]
 }
}
comparisonIDs #list element names and Logical comparisons for those retrievable Dex measurements in the list elements
ls(dexComps) #list elements are dataframes with the DEX entries for that comparison

#Add the column indexes to retrieve pVals and log2Diff for each comparison in the respective dexComps$ list element
comparisonIDs$pValColIndex<-as.integer(c(rev(masterDEXlookup) ))
comparisonIDs$log2DiffIndex<-as.integer(c((rev(masterDEXlookup)+numComp) ))

#comparisonIDs  # Complete lookup table for Dex Target Comparisons


#precalculate maxUP and maxDN to normalize scale of the plots
maxDN<-0
maxUP<-0
yscaleMax<-0

if("MEs" %in% names(net)) { MEs=net$MEs } else {
  if(!exists("MEs")) {
    cat("- MEs data frame for module eigengenes not found. Attempting to recalculate from cleanDat and net$colors...")
    if(!exists("NETcolors")) if("colors" %in% names(net)) NETcolors=net$colors
    if(!exists("cleanDat") | !exists("NETcolors")) stop("Cannot find NETcolors and/or cleanDat for ME calculation.")
    MEs<-tmpMEs<-data.frame()
    MEList = moduleEigengenes(t(cleanDat), colors = NETcolors)
    MEs = orderMEs(MEList$eigengenes)
    colnames(MEs)<-gsub("ME","",colnames(MEs))  # let's be consistent in case prefix was added, remove it.
    if("grey" %in% colnames(MEs)) MEs[,"grey"] <- NULL
  }
}
colnames(MEs)=gsub("ME","",colnames(MEs))
if("grey" %in% colnames(MEs)) MEs[,"grey"] <- NULL

orderedModules=colnames(MEs)
nModules=length(orderedModules)


cat("\nPerforming calculations for mean differential abundance and fractions of differential gene products per module:\n")
dexCompsStacks<-list()
for (redo in 1:2) { #repeat is necessary to equalize min and max blue and red color scheme across all plots
for (z in 1:length(comparisonIDs$Comparison)) {

dataframeName<-as.character(comparisonIDs$dfVariable[z])
dexTargets<-dexComps[[comparisonIDs$dfVariable[z]]]  #retrieve z'th data frame of dexTargets<-eval(parse(text=paste0("dexComps$",comparisonIDs$dfVariable[z])))

orderedLabels<- cbind(paste("M",seq(1:(nModules+20)),sep=""),labels2colors(c(1:(nModules+20))))
# if you want the modules in order of relatedness from the module relatedness dendrogram:
netcolSizeTable<-table(NETcolors)[which(!names(table(NETcolors))=="grey")]
orderedModules2<-cbind(orderedModules,Size=netcolSizeTable[match(orderedModules,names(netcolSizeTable))])
orderedLabelsByRelatedness<- cbind(Mnum= orderedLabels[ match(orderedModules,orderedLabels[,2]) ,1] ,Color=orderedModules )
orderedLabelsByRelatedness<- cbind( orderedLabelsByRelatedness,Size=orderedModules2[,"Size"] )
## Get fraction occupancy and average log2 diff for each group!
#test of function: length(which(dexTargets$NETcolors=="turquoise" & dexTargets[,comparisonIDs$log2DiffIndex[z]]<0))
downTargets<-sapply(1:nrow(orderedLabelsByRelatedness),function(x) length(which(dexTargets$NETcolors==orderedLabelsByRelatedness[x,"Color"] & dexTargets[,comparisonIDs$log2DiffIndex[z]]<0)))
upTargets<-sapply(1:nrow(orderedLabelsByRelatedness),function(x) length(which(dexTargets$NETcolors==orderedLabelsByRelatedness[x,"Color"] & dexTargets[,comparisonIDs$log2DiffIndex[z]]>0)))
otherTargets<-sapply(1:nrow(orderedLabelsByRelatedness),function(x) length(which(dexTargets$NETcolors==orderedLabelsByRelatedness[x,"Color"] & dexTargets[,comparisonIDs$log2DiffIndex[z]]==0)))
orderedLabelsByRelatedness<-cbind(orderedLabelsByRelatedness,downTargets,upTargets,otherTargets)
fractDown<-as.numeric(orderedLabelsByRelatedness[,"downTargets"])/as.numeric(orderedLabelsByRelatedness[,"Size"])
fractUp<-as.numeric(orderedLabelsByRelatedness[,"upTargets"])/as.numeric(orderedLabelsByRelatedness[,"Size"])
fractOther<-as.numeric(orderedLabelsByRelatedness[,"otherTargets"])/as.numeric(orderedLabelsByRelatedness[,"Size"])
yscaleMax<-max( c(yscaleMax, max(eval(fractDown+fractUp+fractOther))) )
downTargetAvg<-sapply(1:nrow(orderedLabelsByRelatedness),function(x) mean(as.numeric(dexTargets[,comparisonIDs$log2DiffIndex[z]][which(dexTargets$NETcolors==orderedLabelsByRelatedness[x,"Color"] & dexTargets[,comparisonIDs$log2DiffIndex[z]]<0)])))
upTargetAvg<-sapply(1:nrow(orderedLabelsByRelatedness),function(x) mean(as.numeric(dexTargets[,comparisonIDs$log2DiffIndex[z]][which(dexTargets$NETcolors==orderedLabelsByRelatedness[x,"Color"] & dexTargets[,comparisonIDs$log2DiffIndex[z]]>0)])))
otherTargetAvg<-sapply(1:nrow(orderedLabelsByRelatedness),function(x) mean(as.numeric(dexTargets[,comparisonIDs$log2DiffIndex[z]][which(dexTargets$NETcolors==orderedLabelsByRelatedness[x,"Color"] & dexTargets[,comparisonIDs$log2DiffIndex[z]]==0)])))
orderedLabelsByRelatedness<-cbind(orderedLabelsByRelatedness,fractDown,fractUp,fractOther,downTargetAvg,upTargetAvg,otherTargetAvg)

#Colorscale 
#bw<-colorRampPalette(c("#0058CC", "white"))
wb<-colorRampPalette(c("white","#0058CC")) # #0058CC"))
wr<-colorRampPalette(c("white", "#CC3300")) # #CC3300"))
colvecwb<-wb(100)
colvecwr<-wr(100)
maxDN<-max( c(maxDN, max(abs(downTargetAvg),na.rm=T)) )
maxUP<-max( c(maxUP, max(upTargetAvg,na.rm=T)) )
cat(paste0("\rAvg. log2(FC) DEx range: -",round(maxDN,2)," to +",round(maxUP,2)," ..."))

vecDN<- -sapply(1:nrow(orderedLabelsByRelatedness),function(x) round( as.numeric(orderedLabelsByRelatedness[x,"downTargetAvg"])/maxDN*100, 0 ))
vecUP<-sapply(1:nrow(orderedLabelsByRelatedness),function(x) round( as.numeric(orderedLabelsByRelatedness[x,"upTargetAvg"])/maxUP*100, 0 ))
vecOTHER<-sapply(1:nrow(orderedLabelsByRelatedness),function(x) round( as.numeric(orderedLabelsByRelatedness[x,"otherTargetAvg"])*0, 0 )) #*0 insures these will be white.
colvecDN<-sapply(1:nrow(orderedLabelsByRelatedness),function(x) colvecwb[vecDN[x]])
colvecUP<-sapply(1:nrow(orderedLabelsByRelatedness),function(x) colvecwr[vecUP[x]])
colvecOTHER<-sapply(1:nrow(orderedLabelsByRelatedness),function(x) colvecwr[vecOTHER[x]])
colvecDN<-lapply(colvecDN, function(x) if (identical(x,character(0))) {"#FFFFFF" } else { x }) #HANDLE rounded values that were 0 (target averages less than 0.005)
colvecUP<-lapply(colvecUP, function(x) if (identical(x,character(0))) {"#FFFFFF" } else { x })
colvecOTHER<-lapply(colvecOTHER, function(x) if (identical(x,character(0))) {"#FFFFFF" } else { x })

colvec<-c(unlist(colvecDN),unlist(colvecUP),unlist(colvecOTHER))
meltPlotMANUAL<-data.frame(Mnum=rep(as.vector(orderedLabelsByRelatedness[,"Mnum"]),3),variable=c(rep(1,nrow(orderedLabelsByRelatedness)),rep(2,nrow(orderedLabelsByRelatedness)),rep(3,nrow(orderedLabelsByRelatedness))),value=c(as.vector(orderedLabelsByRelatedness[,"fractDown"]),as.vector(orderedLabelsByRelatedness[,"fractUp"]),as.vector(orderedLabelsByRelatedness[,"fractOther"])),sort2=rep(c(1:nrow(orderedLabelsByRelatedness)),3),colvec=as.vector(colvec), meanLog2FC=c(as.vector(orderedLabelsByRelatedness[,"downTargetAvg"]),as.vector(orderedLabelsByRelatedness[,"upTargetAvg"]),as.vector(orderedLabelsByRelatedness[,"otherTargetAvg"])))
rownames(meltPlotMANUAL)<-NULL
meltPlotMANUAL2<-meltPlotMANUAL[order(as.numeric(meltPlotMANUAL$sort2),meltPlotMANUAL[,"variable"], decreasing=FALSE),]
rownames(meltPlotMANUAL2)<-NULL
meltPlotMANUAL2$Mnum<-factor(meltPlotMANUAL2$Mnum,levels=unique(meltPlotMANUAL2$Mnum))
meltPlotMANUAL2$meanLog2FC[meltPlotMANUAL2$meanLog2FC=="NaN"]<- 0
dexCompsStacks[[comparisonIDs$dfVariable[z]]] <- meltPlotMANUAL2
}
} #2 iterations to get min and max colors right having considered all plots min and max values from first iteration.
# note: yscaleMax never needs to be >1

## Plot stacked bars and use common color scale for all plots of pairwise comparison DEP fractions' average log2FC
cat("\n\n- Generating PDF file(s)...\n")
suppressPackageStartupMessages(require(scales,quietly=TRUE))

colorData=data.frame(Mnum=as.character(unique(dexCompsStacks[[comparisonIDs$dfVariable[z]]]$Mnum)),yBlank=c(0), fill=gplots::col2hex(unique(WGCNA::labels2colors(as.numeric(gsub("M","",as.character(unique(dexCompsStacks[[comparisonIDs$dfVariable[z]]]$Mnum))))))), fillName=unique(WGCNA::labels2colors(as.numeric(gsub("M","",as.character(unique(dexCompsStacks[[comparisonIDs$dfVariable[z]]]$Mnum)))))))
plot.inset <- ggplot(colorData, aes(x=Mnum, y=rep(1,nrow(colorData)))) + geom_bar(stat="identity", fill=colorData$fill, color="#000000", aes(fill=fillName)) +
  scale_x_discrete( limits=colorData$Mnum)   + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.3), axis.text.y=element_text(color="#FFFFFF")) + labs(x="", y="") + theme(axis.title = element_text(color="#FFFFFF", face="bold", size=22), axis.text.x = element_text(face="bold", color="#000000", size=14, angle=90), axis.ticks.y=element_blank(), panel.background = element_rect(fill = "transparent",colour = NA), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), plot.background = element_rect(fill = "transparent",colour = NA)) + scale_y_continuous(expand = expansion(mult = c(0, 0)), limits=c(0,1), breaks=c(0), label=c("")) + theme(plot.title=element_text(color="black", face="bold", size=32))
grob.inset <- ggplotGrob(plot.inset)

scalebar.text = if(!corVolc) { "mean log2FC\n" } else { "mean R\n" }
yaxis.text = if(!corVolc) { "DEx, Fraction of Module Members" } else { "Fraction of Module Members with Sig. Correl." }
filename.text = if(!corVolc) { "ANOVA.dex-FractionBar" } else { "corStat.sig-FractionBar" }
for (z in 1:length(comparisonIDs$Comparison)) {
 pdf(file=paste0(outputfigs,"/",outFilePrefix,filename.text,"-",gsub("\\s","_",comparisonIDs$Comparison[z]),outFileSuffix,".pdf"),width=15,height=11.25, onefile=FALSE)  # onefile=FALSE for forcing no blank page 1.
 par(mfrow=c(2,1))
 par(mar=c(5,6,4,2))

 plotTitle.text = if(!corVolc) { as.character(comparisonIDs$Comparison[z]) } else { paste0(as.character(comparisonIDs$Comparison[z]), " sample correlation < ",sigCutoff) }
 
 print(
 ggplot(dexCompsStacks[[comparisonIDs$dfVariable[z]]], aes(x=Mnum, y=as.numeric(value))) + geom_bar(stat="identity", color="#000000", aes(fill=as.numeric(dexCompsStacks[[comparisonIDs$dfVariable[z]]]$meanLog2FC))) + labs(x=" \n ", y=yaxis.text) + theme(axis.title = element_text(color="#000000", face="bold", size=22), axis.ticks.x=element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(face="bold", color="#000000", size=14, angle=0), panel.background = element_rect(fill = "transparent",colour = NA), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), plot.background = element_rect(fill = "transparent",colour = NA)) + scale_y_continuous(limits=c(0,ceiling(yscaleMax*100)/100)) +
               scale_fill_gradient2(name=scalebar.text, limits=c(-maxDN,maxUP), oob = scales::squish, breaks=c(-maxDN*0.99999,-maxDN/2,0,maxUP/2,maxUP*0.99999), labels=c(round(-maxDN,2), round(-maxDN/2,2), 0, round(maxUP/2,2), round(maxUP,2)), minor_breaks=NULL, low="#0058CC", mid="white", high="#CC3300", midpoint=0, space="Lab", na.value="grey50", guide=guide_colorbar(label=TRUE, draw.ulim=TRUE, draw.llim = TRUE, frame.colour = "black", ticks = TRUE, barwidth=30, barheight=1.2, ticks.colour='black', direction='horizontal'), aesthetics="fill") + theme(legend.position="top", legend.text=element_text(size=13), legend.title.align=0) +
               ggtitle(plotTitle.text) + theme(plot.title=element_text(color="black", face="bold", size=32)) +
   annotation_custom(grob = grob.inset, xmin = 0.05-(-3.5/nrow(colorData)+1.8/70*nrow(colorData)), xmax = nrow(colorData)+(0.62+0.3/70*nrow(colorData)), ymin = -(0.11+0.46/nrow(colorData))*yscaleMax, ymax = 0 )
 )

# Prior way of filling bars used precalculated fill; no way to put these on a legend scale.
# ggplot(dexCompsStacks[[comparisonIDs$dfVariable[z]]], aes(x=Mnum, y=as.numeric(value))) + geom_bar(stat="identity", fill=dexCompsStacks[[comparisonIDs$dfVariable[z]]]$colvec, color="#000000") + labs(x=" \n ", y="DEPs, Fraction of Module Members") + theme(axis.title = element_text(color="#000000", face="bold", size=22), axis.ticks.x=element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(face="bold", color="#000000", size=14, angle=0), panel.background = element_rect(fill = "transparent",colour = NA), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), plot.background = element_rect(fill = "transparent",colour = NA)) + scale_y_continuous(limits=c(0,ceiling(yscaleMax*100)/100)) +
#               ggtitle(as.character(comparisonIDs$Comparison[z])) + theme(plot.title=element_text(color="black", face="bold", size=32 rotate=90)) +
#   annotation_custom(grob = grob.inset, xmin = 0.05-(-3.5/nrow(colorData)+1.8/70*nrow(colorData)), xmax = nrow(colorData)+(0.62+0.3/70*nrow(colorData)), ymin = -(0.11+0.46/nrow(colorData))*yscaleMax, ymax = 0 ) )
dev.off()
}

} # close DEXpercentStacked function


## adjustments to swatch positions as inset with annotation_custom() of grob.inset above
#-.55, +.75   @ 23 modules
#-2.35 +1.05  @ 93 modules
#0.3 / 70
#-1.8 / 70
#ymin scaling factor:
#-0.11 @ 93 modules
#-0.2  @ 23 modules
#0.09 / 70

#0.05 to the right (tested @ 23 modules)
