###################################################################################################################
# Cross-species Gene List Fisher's Exact Test - optionally adjusts for symbol lookup inefficiency/loss
# by Eric Dammer & Divya Nandakumar
#-------------------------------------------------------------------------#
# +2/10/19 improved duplicate removal within and across reference lists
# +2/10/19 added toggles and speciescode for biomaRt lookup as parameters
# +8/11/20 fixed calculations to match divya's (swapped moduleList and categories vars, and removed unique() for totProteomeLength
# +5/08/21 added barOption - Divya style barplots
# +1/26/23 Converted to geneListFET function for Levi Wood ALS/FTD network collaboration
# +4/17/23 Added RColorBrewer palette specification to parameter/variable paletteColors (vector of length=# of marker list file inputs),
#          vector character strings must be one of the sequential (first group) or qualitative (third group) palettes shown by:
#          RColorBrewer::display.brewer.all()
# +9/29/25 - Changed error handler to avoid precheck for smaller background (total gene list of inputs or modules) compared to the largest reference list (explicit 2x2 table fisher.text errors shown)
#          - Added labeledHeatmap automatic white text for significant heatmap cells, so text on top of dark colors is readable; helper functions of WGCNA labeledHeatmap now included here
#          - Relative file paths for arguments now ok; illegal characters do not propagate to output file names
# +12/2/25 - Added arguments strictSymmetry, asterisksOnly, maxPcolor, bkgrFileForCategories; deprecated verticalCompression argument
#          - Text and margins (inches) auto-size to fit longest string labels on x and y using helper function at end
#          - PDF closes even on error to avoid multiple lingering open files unable to be viewed after failed function call
#          - maxPcolor pads heatmap scale with white for all values above that p (or FDR) value
#          - bkgrFileForCategories allows a single-column, no-header file input for completion of the gene "universe" or background in case a categoriesFile was specified that does not consider background genes/proteins
#          - strictSymmetry=TRUE dereplicates the background AND all duplicates across different category (column) gene lists. If all reference (row) list genes are within the column lists+background, the statistics
#            of Fisher's exact test should be identical/symmetrical even if row (reference) and column (category) lists are swapped.
#          - Note that only column (module or category) lists plus bkgrFileForCategories genes ever count towards background -- reference list genes not hitting any category/module (column) list do not.
#          - asterisksOnly=TRUE prints *, **, and *** for p (or FDR) <0.05, <0.01, <0.005, but without the explicit printing of p or FDR values in any heatmap cell.
#          - added support for viridisLite palettes: "plasma", "inferno", "magma", "rocket", "cividis", "mako", "turbo" (yellow to cool colors only), and "viridis"; must have viridisLite package installed.
# +12/3/25 - uncoupled heatmapScale from legend scale.  New argument: legendScale="minusLog" (default) or "unlog"
#          - unbranched (no conditional logic and replication of) single call to draw heatmap depending on heatmapScale setting
#          - failsafe for heatmap cells below maxPcolor guarantees first non-white color from paletteColors fills cells with values that otherwise might not get color near the threshold
#-------------------------------------------------------------------------#
# revisited to define fly cell types in Seurat 87 lists 2/10/2019
# Analysis for Laura Volpicelli, mouse a-Syn Bilaterally Injected Brain Regions 2/15/2019
# LFQ-MEGA Cell Type analysis performed with this code, with grey proteins added back in to totProteome, allGenes  4/5/2019 #***## (2 lines)
#=========================================================================#
geneListFET <- function(modulesInMemory=TRUE,categoriesFile=NA,categorySpeciesCode=NA,resortListsDecreasingSize=FALSE,barOption=FALSE,adjustFETforLookupEfficiency=FALSE,allowDuplicates=TRUE,strictSymmetry=FALSE,
                        refDataFiles=NA,speciesCode=NA,refDataDescription="Ref_list(s)-rows-not_described",FileBaseName="FET_to_RefList(s)",paletteColors="YlGnBu", asterisksOnly=TRUE,  #colDataDescription="CategOrModules-cols-not_described",
                        heatmapTitle="Enrichment Heatmap (title not specified)", heatmapScale="minusLogFDR", legendScale="minusLog", maxPcolor=0.25, verticalCompression=1, bkgrFileForCategories=NA, rootdir="./", reproduceHistoricCalc=FALSE, env=.GlobalEnv) {

require(WGCNA,quietly=TRUE)
require(RColorBrewer,quietly=TRUE)
require(biomaRt,quietly=TRUE)

refDataDir<-outputfigs<-outputtabs<-rootdir

# standardize heatmapScale options
if (toupper(heatmapScale)=="FDR") heatmapScale <- "minusLogFDR"
if (toupper(heatmapScale)=="P") heatmapScale <- "p.unadj"
heatmapScale <- match.arg(heatmapScale, c("minusLogFDR","p.unadj"))

# normalize legend scale choice (now uncoupled from heatmapScale)
legendScale <- match.arg(legendScale, c("minusLog","unlog"))

if(!modulesInMemory) {  # Read in Categories as list 
  # old format: 2 column .csv with Symbol and "ClusterID" columns
  #  enumeratedListsDF<-read.csv(file=paste0(refDataDir,"/",categoriesFile),header=TRUE)
  #  enumeratedLists<-list()
  #  for(eachList in unique(enumeratedListsDF[,"ClusterID"])) { enumeratedLists[[as.character(eachList)]] <- enumeratedListsDF[which(enumeratedListsDF$ClusterID==eachList),"GeneSymbol"] }

  # new format: multicolumn .csv with Symbols or UniqueIDs and each cluster's symbols in a separate column with clusterID as column name (in row 1)
	enumeratedLists <- as.list(read.csv(paste0(refDataDir,categoriesFile), stringsAsFactors=FALSE,header=T,check.names=FALSE)) 
	names(enumeratedLists)
		
	#number of entries with no blanks
	length(unlist(lapply(enumeratedLists,function(x) x[!x==''] )))
	#take out blanks from list
	enumeratedLists <- lapply(enumeratedLists,function(x) x[!x==''] )
	#are there symbols duplicated? (yes, if below result is less than above)
	length(unique(unlist(enumeratedLists)))
	
	
	enumeratedLists<-lapply(enumeratedLists,function(x) as.data.frame(do.call("rbind",strsplit(x,"[|]")))[,1] )
	# leave duplicated symbols within each module/category list -
	#enumeratedLists<-lapply(enumeratedLists,unique)
	
	# leave duplicates in the clusters or modules being checked - they should contribute to overlap/enrichment with reference lists more than once if duplicated.
	#if(!allowDuplicates) {
	#  while( length(unique(unlist(enumeratedLists))) < length(unlist(enumeratedLists)) ) {
	#    duplicatedvec<-unique(unlist(enumeratedLists)[which(duplicated(unlist(enumeratedLists)))])
	#    #remove duplicates from any marker list
	#    enumeratedLists<-lapply(enumeratedLists,function(x) { remIndices=as.vector(na.omit(match(duplicatedvec,x))); if (length(remIndices)>0) { x[-remIndices] } else { x }; } )
	#  }
	#}

	if(!is.na(bkgrFileForCategories)) {
	  if(file.exists(bkgrFileForCategories)) {
	    bkgr.1list<-as.list(read.csv(paste0(refDataDir,bkgrFileForCategories), stringsAsFactors=FALSE,header=F,check.names=FALSE))
	    bkgr.vec<-lapply(bkgr.1list,function(x) x[!x==''] )
	    greyToAddToTotProteome<- as.vector( data.frame(do.call("rbind", strsplit(bkgr.vec[[1]], "[|]")))[,1] )  # [[1]] only operate on first column
	  } else { stop("Specified file ",bkgrFileForCategories," specified in optional argument bkgrFileForCategories not found.\nThis should be a single column file with no header, containing all UniqueIDs or symbols expected in any category (heatmap column).") }
	}

} else {  # use WGCNA modules in memory

  # Module lookup table
  nModules<-length(table(net$colors))-1
  modules<-cbind(colnames(as.matrix(table(net$colors))),table(net$colors))
  orderedModules<-cbind(Mnum=paste("M",seq(1:nModules),sep=""),Color=WGCNA::labels2colors(c(1:nModules)))
  modules<-modules[match(as.character(orderedModules[,2]),rownames(modules)),]
  #as.data.frame(cbind(orderedModules,Size=modules))

  # Recalculate Consensus Cohort Eigengenes, i.e. eigenproteins and their relatedness order
  MEs<-data.frame()
  MEList = WGCNA::moduleEigengenes(t(cleanDat), colors = net$colors)
  MEs = orderMEs(MEList$eigengenes)
  net$MEs <- MEs
  colnames(MEs)<-gsub("ME","",colnames(MEs)) #let's be consistent in case prefix was added, remove it.
  rownames(MEs)<-rownames(numericMeta)
  if("grey" %in% colnames(MEs)) MEs[,"grey"] <- NULL

  # Make list of module member gene product official symbols
  enumeratedLists<-sapply( colnames(MEs),function(x) as.vector(data.frame(do.call("rbind",strsplit(rownames(cleanDat),"[|]")))[,1])[which(net$colors==x)] )
  greyToAddToTotProteome<- as.vector( data.frame(do.call("rbind",strsplit(rownames(cleanDat),"[|]")))[,1])[which(net$colors=="grey")] #***##
}
moduleList=enumeratedLists


# ---- build x-axis labels exactly as in your heatmap
if (modulesInMemory) {
  categoryColorSymbols <- paste0("ME", names(moduleList))
  xSymbolsText <- paste0(
    names(moduleList), " ",
    orderedModules[match(names(moduleList), orderedModules[,2]), 1]
  )
} else {
  categoryColorSymbols <- names(moduleList)
  xSymbolsText <- names(moduleList)
}

# choose cex dynamically (reuse your earlier scale_cex or set fixed if you prefer)
scale_cex <- function(chars, n, base = 1.6, per_char = 0.015, per_n = 0.003, mn = 0.6, mx = 1.8) {
  out <- base - per_char * chars - per_n * n
  max(mn, min(mx, out))
}
cexX_global <- scale_cex(max(nchar(xSymbolsText)), length(xSymbolsText))
# for Y we will compute per page, but we also need a worst-case device size:
get_header_names <- function(fp) {
  df <- try(read.csv(fp, stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, nrows = 1), silent = TRUE)
  if (inherits(df, "try-error") || is.null(df)) character(0) else colnames(df)
}
all_y_labels <- unique(unlist(lapply(refDataFiles, function(rdf) get_header_names(paste0(refDataDir, rdf)))))

# worst-case counts across pages
max_rows_overall <- max(1L, length(all_y_labels))
nCols_overall    <- length(xSymbolsText)

# pick a conservative cexY for worst-case sizing
cexY_global <- scale_cex(max(nchar(all_y_labels)), max_rows_overall, base = 1.4)

# compute once: device size that fits the worst page, including margins
sz_overall <- compute_pdf_and_margins(
  xLabs = xSymbolsText,
  yLabs = if (length(all_y_labels)) all_y_labels else "X",
  nCols = nCols_overall,
  nRows = max_rows_overall,
  cexX  = cexX_global,
  cexY  = cexY_global,
  xAngle = 90,
  asterisksOnly = asterisksOnly,
  barOption = barOption
)


# open the PDF (sized for dimensions big enough for the largest page)
refDataDescription.forFilename <- gsub('\\/','.', gsub('\\\\', ".", gsub('\\.\\.','', refDataDescription)))
pdf_path <- paste0(outputfigs,"/",FileBaseName,".Overlap.in.",refDataDescription.forFilename,".pdf")

grDevices::pdf(file = pdf_path, width = sz_overall$width, height = sz_overall$height, onefile = TRUE)
pdf_dev <- grDevices::dev.cur()

# Ensure the PDF device is closed even if an error happens later
on.exit({
  dl <- grDevices::dev.list()
  if (!is.null(dl) && pdf_dev %in% dl) {
    try(grDevices::dev.off(pdf_dev), silent = TRUE)
  }
}, add = TRUE)
############

#***iterating through multiple files (each one a page of output PDF):  
iter=0
for (refDataFile in refDataFiles) {
	iter=iter+1
	this.heatmapScale<-heatmapScale
	
	
	refData <- as.list(read.csv(paste0(refDataDir,refDataFile), stringsAsFactors = FALSE,header=T,check.names=FALSE)) 
	names(refData)
	
	
	#number of entries with no blanks
	length(unlist(lapply(refData,function(x) x[!x==''] )))
	#take out blanks from list
	refData <- lapply(refData,function(x) x[!x==''] )
	#are there symbols duplicated? (yes, if below result is less than above)
	length(unique(unlist(refData)))
	
	##Remove duplicates from all lists if allowDuplicates=FALSE
	#remove duplicated exact symbols within each list regardless:
	if (reproduceHistoricCalc) refData<-lapply(refData,unique)
	refData<-lapply(refData,function(x) as.data.frame(do.call("rbind",strsplit(x,"[|]")))[,1] )
	if (!reproduceHistoricCalc) refData<-lapply(refData,unique)
	
	if(!allowDuplicates) {
	  while( length(unique(unlist(refData))) < length(unlist(refData)) ) {
	    duplicatedvec<-unique(unlist(refData)[which(duplicated(unlist(refData)))])
	    #remove duplicates from any marker list
	    refData<-lapply(refData,function(x) { remIndices=as.vector(na.omit(match(duplicatedvec,x))); if (length(remIndices)>0) { x[-remIndices] } else { x }; } )
	  }
	  duplicateHandling="DupsInRow.REMOVED"
	} else {
	  duplicateHandling="DupsInRow.ALLOWED"
	}
	length(unlist(refData))
	unlist(refData)[which(duplicated(unlist(refData)))]
	refDataMouse<-refData
	
	groupvec<-placeholders<-vector()
	for (i in 1:length(names(refData))) {
		placeholders=c(placeholders,rep(i,length(refData[[i]])))
		groupvec=c(groupvec,rep(names(refData)[i],length(refData[[i]])))
	}
	categoriesData<-data.frame(UniqueID=unlist(refData),Color=labels2colors(placeholders), Annot=groupvec,Mnum=paste0("M",placeholders)) #,row.names=unlist(refData)) #will not work if allowDuplicates==TRUE
	#^holding refData all lists' items, not modules gene symbols
	
	categoriesNameMatcher<-unique(categoriesData[,2:4])
	rownames(categoriesNameMatcher)<-NULL
	
	
	if(!categorySpeciesCode==speciesCode[iter]) {
	  cat(paste0("Converting ",speciesCode[iter]," to ",categorySpeciesCode," for lists in ",refDataFile," ... "))
	
	  this.heatmapTitle=paste0(heatmapTitle," in ",categorySpeciesCode," homologs")
	  #library(biomaRt)
	
	  #human = useEnsembl("genes", dataset = "hsapiens_gene_ensembl", host="https://dec2021.archive.ensembl.org")  #ver=105 equivalent to dec2021
	  #mouse = useEnsembl("genes", dataset = "mmusculus_gene_ensembl", host="https://dec2021.archive.ensembl.org") 
	
	  category.species = useEnsembl("genes", dataset=paste0(categorySpeciesCode,"_gene_ensembl"), host="https://dec2021.archive.ensembl.org")  #ver=105 equivalent to dec2021  ; old code: #useMart("ensembl",dataset=paste0(categorySpeciesCode,"_gene_ensembl"))
	  other = useEnsembl("genes", dataset=paste0(speciesCode[iter],"_gene_ensembl"), host="https://dec2021.archive.ensembl.org")   # old code: #useMart("ensembl",dataset=paste0(speciesCode[iter],"_gene_ensembl"))
	
	  #category species to other species conversion (first column is other species)
	  if(speciesCode[iter]=="hsapiens") {
	    genelist.mouseConv<-getLDS(attributes=c("hgnc_symbol"), filters="hgnc_symbol", values=categoriesData$UniqueID, mart=other, attributesL="external_gene_name",martL = category.species)
	    categoriesData$BiomartFlySymbol <- genelist.mouseConv[match(categoriesData$UniqueID,genelist.mouseConv$HGNC.symbol),"Gene.name"]
	  } else {
	    if(categorySpeciesCode=="hsapiens") {
	      genelist.mouseConv<-getLDS(attributes=c("external_gene_name"), filters="external_gene_name", values=categoriesData$UniqueID, mart=other, attributesL="hgnc_symbol",martL = category.species)
	      categoriesData$BiomartFlySymbol <- genelist.mouseConv[match(categoriesData$UniqueID,genelist.mouseConv$Gene.name),"HGNC.symbol"]
	    } else {
	      genelist.mouseConv<-getLDS(attributes=c("external_gene_name"), filters="external_gene_name", values=categoriesData$UniqueID, mart=other, attributesL="external_gene_name",martL = category.species)
	      categoriesData$BiomartFlySymbol <- genelist.mouseConv[match(categoriesData$UniqueID,genelist.mouseConv$Gene.name),"Gene.name.1"]
	    }
	  }
	} else { this.heatmapTitle=heatmapTitle }
	
	categoriesData.original<-categoriesData
	
	categoriesData.reducedFly<-na.omit(categoriesData)
	categoriesData.reducedFly$MouseID.original<-categoriesData.reducedFly$UniqueID
	if(!categorySpeciesCode==speciesCode[iter]) categoriesData.reducedFly$UniqueID<-categoriesData.reducedFly$BiomartFlySymbol
	categoriesData.reducedFly<-categoriesData.reducedFly[,-5]
	
	refDataFly<-list()
	for (i in unique(categoriesData.reducedFly$Annot)) {
		refDataFly[[i]]<-unique(categoriesData.reducedFly$UniqueID[which(categoriesData.reducedFly$Annot==i)])
	}
	refDataFly.original<-refDataFly
	
	#remove within-list duplicates from any marker list (some homologs map to multiple reference species unique genes)
	refDataFly<-lapply(refDataFly,unique)
	#final check
	length(unlist(refDataFly))
	unlist(refDataFly)[which(duplicated(unlist(refDataFly)))] #any duplicates across lists allowed
	length(unlist(refDataFly))==length(unique(unlist(refDataFly))) #false if duplicates across lists.
	
	refData<-refDataFly
	
	#make data frame of all markers
	mouseSymbolVec<-groupvec<-placeholders<-vector()
	for (i in 1:length(names(refData))) {
		placeholders=c(placeholders,rep(i,length(refData[[i]])))
		groupvec=c(groupvec,rep(names(refData)[i],length(refData[[i]])))
		mouseSymbolVec=c(groupvec,rep(names(refData)[i],length(refData[[i]])))
	}
	categoriesData<-data.frame(UniqueID=unlist(refData),Color=labels2colors(placeholders), Annot=groupvec,Mnum=paste0("M",placeholders))  #,row.names=unlist(refData))
	categoriesData$MouseSymbol=categoriesData.reducedFly$MouseID.original[match(categoriesData$UniqueID,categoriesData.reducedFly$UniqueID)]
	
	categoriesNameMatcher<-unique(categoriesData[,2:4])
	rownames(categoriesNameMatcher)<-NULL
	
	
	
	if(modulesInMemory) {
	  allGenes<- c(unlist(moduleList),greyToAddToTotProteome)  # unique() here decreases significance, totProteomeLength; Not in original code.
	  if (strictSymmetry) allGenes <- unique(allGenes)         # strictSymmetry enforced:  duplicates across all modules (columns) are removed.
	} else {
	  if(!is.na(bkgrFileForCategories)) {
	    allGenes <- unlist(enumeratedLists)
	    if (strictSymmetry) allGenes <- unique(allGenes)       # strictSymmetry enforced: duplicates across all categories (columns) removed.
	    greyToAddToTotProteome <- greyToAddToTotProteome[which(!greyToAddToTotProteome %in% allGenes)]  # only add background if not already in categories. 
	    if (strictSymmetry) greyToAddToTotProteome <- unique(greyToAddToTotProteome)
	    allGenes <- c(allGenes, greyToAddToTotProteome)
	  } else {  # no additional background (grey-like) gene list provided
	    allGenes<- unlist(enumeratedLists)  #                  # here the background is all measured (category/column list) proteins
	    if (strictSymmetry) allGenes <- unique(allGenes)       # strictSymmetry enforced: duplicates across all categories (columns) removed.
	  }
	}
	allGenesNetwork <- as.matrix(allGenes,stringsAsFactors = FALSE) 

	categories <- list()
	categoryNames=names(refData) #reference list names
	for (i in 1:length(categoryNames)) {
		element<-categoryNames[i]
		categories[[element]] <- categoriesData$UniqueID[which(categoriesData$Annot==categoryNames[i])]  #categoriesData$BiomartMouseSymbol[which(categoriesData$colors==modcolors[i])]
	}
	
	##+#+#+#+#+#+#+#+#+#+#+#+#+
	# Final Data Cleaning
	
	nModules <- length(names(moduleList))
	nCategories <- length(names(categories))
	
	for (a in 1:nCategories) {
		categories[[a]] <- unique(categories[[a]][categories[[a]] != ""])
		categories[[a]] <- categories[[a]][!is.na(categories[[a]])]
	}
	for (b in 1:nModules) {
		moduleList[[b]] <- unique(moduleList[[b]][moduleList[[b]] != ""])
	}
	
	if(resortListsDecreasingSize) {
	  categories <- categories[order(sapply(categories,length),decreasing=T)]
	  #only sort lists if 'moduleList' is not a list of WGCNA modules (keep them in relatedness order if they are)
	  if (!modulesInMemory) moduleList <- moduleList[order(sapply(moduleList,length),decreasing=T)]
	} #if FALSE, do not resort lists -- we have them in a precise order already
	
	
	allGenes_cleaned <- na.omit(allGenesNetwork)
	totProteomeLength <- length(allGenes_cleaned)
	## fisher.test below now reports specific contingency table causing error, rather than preventing running here.
#	if(max(sapply(refData,length))>totProteomeLength) {
#	  cat (paste0("One of your reference data lists is larger than the background from the categories (WGCNA or specified categories file--all symbols)!\nUsing the bigger number for Fisher Exact would change all stats. Skipping ",refDataFile,".\n\n"))
#	  next
#	}
	
	### Fisher's Exact Test
	
	cat(paste0("Performing FET for lists [",iter,"] now.\n"))
	
	
	FTpVal <- matrix(,nrow = nModules, ncol = nCategories)
	categoryOverlap <- matrix(,nrow = nModules, ncol = nCategories) 
	numCategoryHitsInDataset <- numCategoryHitsInDataset.UNADJ <- matrix(,nrow = nModules, ncol = nCategories) 
	CategoryHitsInDataset <- list()
	hitLists<-matrix(NA,nrow=nModules,ncol=nCategories) #use a matrix of collapsed (";") gene list strings
	ADJRedundancyAfterLookup=1 #length(unlist(refDataFly.original))/length(unlist(refDataFly)) #bigger than 1
	ADJforCrossSpeciesLookupFailure=nrow(categoriesData)/nrow(categoriesData.original) #less than 1
	totProteomeLength.ADJ <- as.integer(totProteomeLength*ADJforCrossSpeciesLookupFailure*ADJRedundancyAfterLookup)
	RefDataElements<-Categories1<-vector()
	
	for (i in 1:nModules){
		sampleSize <- length(moduleList[[i]])
		RefDataElements=c(RefDataElements,names(moduleList)[i])
		for (j in 1:nCategories){
			if(i==1) { Categories1=c(Categories1,names(categories)[j]) }
			#CategoryHitsInProteome <- categories[[j]] ## If using all of the markers and not just markers in proteome
			CategoryHitsInProteome <- intersect(categories[[j]],allGenesNetwork[,1])
			if (!adjustFETforLookupEfficiency) {
			  ##Unadjusted calculations:
			  numCategoryHitsInProteome <- length(CategoryHitsInProteome) 
			  numNonCategoryHitsInProteome <- totProteomeLength - numCategoryHitsInProteome
			  overlapGenes <- intersect(moduleList[[i]],CategoryHitsInProteome)
			  numOverlap <- length(overlapGenes)
			  otherCategories <- sampleSize - numOverlap
			  notInModule <- numCategoryHitsInProteome - numOverlap
			  notInMod_otherCategories <- totProteomeLength - numCategoryHitsInProteome - otherCategories
			} else {  ##allGenesNetwork has different species Symbols, and categories are also from that full Symbol List, so adjust for comparison to interconverted list overlap
			  #&& adjustments noted
			  numCategoryHitsInProteome <- as.integer(length(CategoryHitsInProteome)*(ADJforCrossSpeciesLookupFailure*ADJRedundancyAfterLookup)) #&& adjusted down for lookup inefficiency
			  numCategoryHitsInProteome.UNADJ <- length(CategoryHitsInProteome)
			  numNonCategoryHitsInProteome <- totProteomeLength.ADJ - numCategoryHitsInProteome #first term is adjusted
			  overlapGenes <- intersect(moduleList[[i]],CategoryHitsInProteome) #does not need adjustment, both subject to lower lookup efficiency
	
			  numOverlap <- length(overlapGenes)
			  otherCategories <- sampleSize - numOverlap
			  notInModule <- numCategoryHitsInProteome - numOverlap #&&using down-adjusted number numCategoryHitsInProteome for lookup efficiency
			  notInMod_otherCategories <- totProteomeLength - numCategoryHitsInProteome - otherCategories #&&first term not adjusted down because this is non-overlap so lookup inefficiency does not apply
									#&& but second term is adjusted because it it the hits subject to lookup efficiency
			}
			hitLists[i,j]<-paste(overlapGenes,collapse=";")
			contingency <- matrix(c(numOverlap,otherCategories,notInModule,notInMod_otherCategories),nrow=2,ncol=2,dimnames=list(c("GenesHit","GenesNotHit"),c("withinCategory","inProteome")))
	#debugging:     if(i==6 & j==3) cat(contingency)
			# variable with presumed explanatory effect should be the row definitions, if known. (can transpose, but no effect on outcome p values)
			FT <- tryCatch(  # now reports specific contingency table causing error, rather than trying to prevent running with precheck above.
		                       fisher.test(contingency, alternative = "greater"),
		                       error = function(e) {
		                         msg <- paste0(
		                           "Fisher test failed at module ", i, ", category ", j, "\n",
		                           "Cause: ", conditionMessage(e), "\n",
		                           "Contingency table:\n",
		                           capture.output(print(contingency))
		                         )
		                         stop(msg, call. = FALSE)
		                       }
		                      )
			FTpVal[i,j] <- FT$p.value
			categoryOverlap[i,j] <- numOverlap
			numCategoryHitsInDataset[i,j] <- numCategoryHitsInProteome
			numCategoryHitsInDataset.UNADJ[i,j] <- if(adjustFETforLookupEfficiency) { numCategoryHitsInProteome.UNADJ } else { numCategoryHitsInProteome }
			if (i==1){
				CategoryHitsInDataset[[j]] <- array(CategoryHitsInProteome)
			}		
		}
	}
	
	

	rownames(FTpVal) <- RefDataElements
	colnames(FTpVal) <- Categories1
	rownames(categoryOverlap) <- RefDataElements
	colnames(categoryOverlap) <- Categories1
	colnames(numCategoryHitsInDataset) <- Categories1
	rownames(numCategoryHitsInDataset) <- RefDataElements
	names(CategoryHitsInDataset) <- Categories1
	rownames(hitLists) <- RefDataElements
	colnames(hitLists) <- Categories1
	
	
	#### Format Data for Plotting ########
	
	NegLogUncorr <- -log10(FTpVal)
	rownames(NegLogUncorr) <- rownames(FTpVal)
	colnames(NegLogUncorr) <- colnames(FTpVal)
	NegLogUncorr <- as.matrix(NegLogUncorr)
	
	nCategories = ncol(FTpVal)
	nModules = nrow(FTpVal)
	
	FisherspVal <- unlist(FTpVal)
	adjustedPVal <- p.adjust(FisherspVal, method = "fdr", n=length(FisherspVal))
	adjustedPval <- matrix(adjustedPVal,nrow=nModules,ncol=nCategories)
	rownames(adjustedPval) <- rownames(FTpVal)
	colnames(adjustedPval) <- colnames(FTpVal)
	NegLogCorr <- -log10(adjustedPval)
	
	## Transpose above stats and hits matrices
	categoryOverlap<-t(categoryOverlap)
	numCategoryHitsInDataset<-t(numCategoryHitsInDataset)
	numCategoryHitsInDataset.UNADJ<-t(numCategoryHitsInDataset.UNADJ)
	CategoryHitsInDataset<-t(CategoryHitsInDataset)
	hitLists<-t(hitLists)
	NegLogUncorr<-t(NegLogUncorr)
	NegLogCorr<-t(NegLogCorr)
	adjustedPval<-t(adjustedPval)
	FTpVal<-t(FTpVal)
	
	
	##Make sure colors are in correct (WGCNA) order before changing to numbered modules!
	if(modulesInMemory) {
	  orderedLabels<-cbind(paste("M",seq(1:nCategories),sep=""),labels2colors(c(1:nCategories)))
	} else {
	  orderedLabels<- cbind(paste("M",seq(1:nModules),sep=""),labels2colors(c(1:nModules))) #these go from M1 to M(# of reference lists)
	}
	
	#if you want the modules in order of relatedness from the module relatedness dendrogram:
	if(!modulesInMemory) {
	  orderedLabelsByRelatedness<- orderedLabels #(this is chronol. order)
	  if (!length(na.omit(match(orderedLabelsByRelatedness[,2],RefDataElements)))==nrow(orderedLabelsByRelatedness)) orderedLabelsByRelatedness[,2]<- RefDataElements; dummyColors=orderedLabelsByRelatedness[,2]; # our category/cluster names on categoriesFile row 1 are not WGCNA colors.
	  NegLogUncorr<-NegLogUncorr[,match(orderedLabelsByRelatedness[,2],colnames(NegLogUncorr))]
	  NegLogCorr<-NegLogCorr[,match(orderedLabelsByRelatedness[,2],colnames(NegLogCorr))]
	  adjustedPval<-adjustedPval[,match(orderedLabelsByRelatedness[,2],colnames(adjustedPval))]
	  FTpVal<-FTpVal[,match(orderedLabelsByRelatedness[,2],colnames(FTpVal))]
	} else {
	  orderedLabelsByRelatedness<- cbind( orderedLabels[ match(gsub("ME","",colnames(MEs)),orderedLabels[,2]) ,1] ,gsub("ME","",colnames(MEs)) )
	
	  NegLogUncorr<-NegLogUncorr[,match(orderedLabelsByRelatedness[,2],colnames(NegLogUncorr))]
	  NegLogCorr<-NegLogCorr[,match(orderedLabelsByRelatedness[,2],colnames(NegLogCorr))]
	  adjustedPval<-adjustedPval[,match(orderedLabelsByRelatedness[,2],colnames(adjustedPval))]
	  FTpVal<-FTpVal[,match(orderedLabelsByRelatedness[,2],colnames(FTpVal))]
	}
	xlabels <- orderedLabelsByRelatedness[,1]
	

	### Write p Values to a table/file
	#rownames(hitLists)<-categoriesNameMatcher$Annot[match(rownames(categoryOverlap),categoriesNameMatcher$Annot)]
	#rownames(FTpVal)<-categoriesNameMatcher$Annot[match(rownames(FTpVal),categoriesNameMatcher$Annot)]
	#rownames(adjustedPval)<-categoriesNameMatcher$Annot[match(rownames(adjustedPval),categoriesNameMatcher$Annot)]
	#rownames(categoryOverlap)<-categoriesNameMatcher$Annot[match(rownames(categoryOverlap),categoriesNameMatcher$Annot)]
	
	outputData <- rbind("FET pValue", FTpVal,"FDR corrected",adjustedPval,"Overlap",categoryOverlap,"CategoryHitsInDataSet(ADJ)",numCategoryHitsInDataset,"CategoryHitsInDataSet(Unadj)",numCategoryHitsInDataset.UNADJ,"OverlappedGeneLists",hitLists)
        refDataFile.forFilename<-gsub('\\/','.',gsub('\\\\', ".", gsub('\\.\\.','',refDataFile)))  # avoid relative path invalid characters inclusion in csv filename
	#write.csv(outputData,file = paste0(outputtabs,"/",FileBaseName,".Overlap.in.",refDataFile.forFilename,"-",duplicateHandling,"-hitListStats.csv"))
	# --- write hitListStats with fallback to short unique name if open/write fails ---
	primary_csv <- file.path(
	  outputtabs,
	  paste0(FileBaseName, ".Overlap.in.", refDataFile.forFilename, "-", duplicateHandling, "-hitListStats.csv")
	)
	
	.try_write_csv <- function(dat, path) {
	  ok <- tryCatch({
	    utils::write.csv(dat, file = path)
	    TRUE
	  }, error = function(e) {
	    FALSE
	  })
	  ok
	}
	
	if (suppressWarnings(.try_write_csv(outputData, primary_csv))) {
	  cat(sprintf("- Wrote hitListStats to: %s\n", primary_csv))
	} else {
	  cat("x Primary CSV write failed (likely path too long). Trying short fallback name...\n")
	  i <- 1L
	  wrote <- FALSE
	  repeat {
	    fallback_csv <- file.path(outputtabs, paste0(i, ".hitListStats.csv"))
	    if (!file.exists(fallback_csv)) {
	      if (suppressWarnings(.try_write_csv(outputData, fallback_csv))) {
	        cat(sprintf("+ Fallback write succeeded: %s\n", fallback_csv))
	        wrote <- TRUE
	        break
	      }
	    }
	    i <- i + 1L
	    if (i > 10000L) {
	      stop("x Failed to write hitListStats.csv after 10000 attempts (both primary and fallback).")
	    }
	  }
	  if (!wrote) stop("x Could not write hitListStats.csv (short filename fallback failed).")
	}

	#auto-check if all FET (BH) calculations = 1, then switch to p value visualization
		if(length(dim(adjustedPval))<2) {  #*** handle single input list
	  adjustedPval<-matrix(adjustedPval,byrow=TRUE,nrow=1)
	  FTpVal<-matrix(FTpVal,byrow=TRUE,nrow=1)
	  NegLogCorr<-matrix(NegLogCorr,byrow=TRUE,nrow=1)
	  NegLogUncorr<-matrix(NegLogUncorr,byrow=TRUE,nrow=1)
	}
	if(mean(rowMeans(adjustedPval,na.rm=T),na.rm=T)==1) { this.heatmapScale<-"p.unadj"; addText="-fallback b/c no FDR values lower than 100%"; } else { this.heatmapScale<-heatmapScale; addText=""; }
	
	## Use the text function with the FDR filter in labeledHeatmap to add asterisks, e.g. * 
	 txtMat <- adjustedPval
	 txtMat[adjustedPval>=0.05] <- ""
	  txtMat[adjustedPval <0.05&adjustedPval >0.01] <- "*"
	  txtMat[adjustedPval <0.01&adjustedPval >0.005] <- "**"
	  txtMat[adjustedPval <0.005] <- "***"
	
	  txtMat1 <- signif(adjustedPval,2)
	  txtMat1[adjustedPval>maxPcolor] <- ""
	
	  
	  if (!asterisksOnly) { textMatrix1 = paste( txtMat1, txtMat , sep = ' ') } else { textMatrix1 = txtMat }
	  textMatrix1= matrix(textMatrix1,ncol=ncol(adjustedPval),nrow=nrow(adjustedPval))
	
	  #for textMatrix of p.unadj
	 txtMat <- FTpVal
	 txtMat[FTpVal>=0.05] <- ""
	  txtMat[FTpVal <0.05&FTpVal >0.01] <- "*"
	  txtMat[FTpVal <0.01&FTpVal >0.005] <- "**"
	  txtMat[FTpVal <0.005] <- "***"
	
	  txtMat.p.unadj <- signif(FTpVal,2)
	  txtMat.p.unadj[FTpVal>maxPcolor] <- ""
	
	  if (!asterisksOnly) { textMatrix.p.unadj = paste( txtMat.p.unadj, txtMat , sep = ' ') } else { textMatrix.p.unadj = txtMat }
	  textMatrix.p.unadj= matrix(textMatrix.p.unadj,ncol=ncol(FTpVal),nrow=nrow(FTpVal))
	
	
	## Plotting
	if(!barOption) {
		if(exists("colvec")) suppressWarnings(rm(colvec))
		
		#RColorBrewer::display.brewer.all()
		if(iter>length(paletteColors)) {
		   cat(paste0("  - paletteColors specified being recycled for additional heatmaps for inputs after #",length(paletteColors),"\n"))
		   paletteColors<-c(paletteColors,rep(paletteColors,ceiling(length(refDataFiles)/length(paletteColors))))
		}
		if(!paletteColors[iter] %in% c("YlOrRd","YlOrBr","YlGnBu","YlGn","Reds","RdPu","Purples","PuRd","PuBuGn","PuBu","OrRd","Oranges","Greys","Greens","GnBu","BuPu","BnGn","Blues",
		                           "Spectral","RdYlGn","RdYlBu","RdGy","RdBu","PuOr","PRGn","PiYG","BrBG")) {
		   if(!paletteColors[iter] %in% c("plasma","inferno","magma","rocket","cividis","mako","turbo","viridis")) {
		     cat(paste0("  - paletteColors specified as '",paletteColors[iter],"' is not in RColorBrewer::display.brewer.all() groups 1 or 3; nor is it provided by viridisLite package.\n    Using palette 'YlGnBu' (yellow, green, blue)...\n"))
		     paletteColors[iter]="YlGnBu"
		   } else {
		     require(viridisLite,quietly=TRUE)
		   }
		}
		if(paletteColors[iter] %in% c("YlOrRd","YlOrBr","YlGnBu","YlGn","Reds","RdPu","Purples","PuRd","PuBuGn","PuBu","OrRd","Oranges","Greys","Greens","GnBu","BuPu","BnGn","Blues")) {
		   paletteLength=9
		   outOfParkColor=brewer.pal(paletteLength,paletteColors[iter])[paletteLength]
		   colvec<- brewer.pal(paletteLength,paletteColors[iter])[1:6]
		} else {
		   if(paletteColors[iter] %in% c("plasma","inferno","magma","rocket","cividis","mako","turbo","viridis")) {  # viridisLite palettes
		     if(paletteColors[iter]=="turbo") {
		       paletteLength=50
		       fullPalette=rev(do.call(paletteColors[iter], as.list(paletteLength)))
		       outOfParkColor=fullPalette[paletteLength]
		       colvec<- c("#FFFF37FF", fullPalette[c(22,24,26,28,30,32,34,36,38,40,42,44)]) #45-49 skipped for outOfPark effect; starts with yellow
		     } else {  # not turbo (2-sided intense ends of palette)
		       paletteLength=15
		       fullPalette=rev(do.call(paletteColors[iter], as.list(paletteLength)))
		       outOfParkColor=fullPalette[paletteLength]
		       colvec<- fullPalette[1:(paletteLength-4)]
		     }
		   } else {  # display.bewer.all() group 3; take half of the shifting gradient
		     paletteLength=11
		     outOfParkColor=brewer.pal(paletteLength,paletteColors[iter])[paletteLength]
		     colvec<- rev(brewer.pal(paletteLength,paletteColors[iter])[1:6])  #rev so we take the left side of the palette color swatches
		   }
		}
		   
		colvecBase<- vector()
		for (k in 1:(length(colvec)-1)) {
		   gradations <- if (k<4) { 6 } else { 25 }
		   temp<-colorRampPalette(c(colvec[k],colvec[k+1]))
		   colvecBase<-c(colvecBase, temp(gradations))
		}
		
		temp2<-colorRampPalette(c(colvecBase[length(colvecBase)], outOfParkColor)) ## grade to outOfParkColor at top of scale
		colvecBase<-c(colvecBase, temp2(gradations))
		
                ## helpers to generate data-aware palettes for the two legend scales (robust to edge cases)
                make_cols_neglog <- function(zmax, thrNeg, eps = 1e-9) {
                  if(thrNeg==0) return(c("#FFFFFF",colvecBase))
                  cols <- colvecBase
                  nC   <- length(cols)
                  if (!is.finite(zmax) || zmax <= 0) return(rep("white", nC))
                  thrNeg <- max(0, thrNeg)
                  ## Case A: nothing exceeds threshold -> all white (old logic equivalent to zmax <= thrNeg)
                  if (zmax <= thrNeg + eps) return(rep("white", nC))
                  ## Case B: threshold ~0 -> no whites
                  if (thrNeg <= eps) return(cols)
                  ## General case: left-pad with white up to thrNeg
                  frac <- thrNeg / zmax
                  ## Keep strictly below 1 to avoid division by zero
                  frac <- max(0, min(1 - eps, frac))
                  nW <- ceiling((frac / (1 - frac)) * nC)
                  nW <- max(0L, nW)  #min(nW, 5L * nC))  # safety clamp
                  c(rep("white", nW), cols)
                }
                make_cols_unlog <- function(pmin, pmax, maxPcolor, eps = 1e-9) {
                  cols <- rev(colvecBase)
                  nC   <- length(cols)
                  if (!is.finite(pmin)) pmin <- 0
                  if (!is.finite(pmax) || pmax <= 0) pmax <- 1
                  rng <- max(pmax - pmin, .Machine$double.eps)
                  ## Case A: everything above threshold -> all white (min > maxPcolor)
                  if (pmin >= maxPcolor - eps) return(rep("white", nC))
                  ## Case B: nothing above threshold -> no whites (max <= maxPcolor)
                  if (pmax <= maxPcolor + eps) return(cols)
                  ## General case: append whites at the TOP end for p >= maxPcolor
                  frac <- (pmax - maxPcolor) / rng
                  frac <- max(0, min(1 - eps, frac))
                  nW <- ceiling((frac / (1 - frac)) * nC)
                  nW <- max(0L, nW)  #min(nW, 5L * nC))  # safety clamp
                  c(cols, rep("white", nW))
                }

		## Build palette for minusLog (left padding = whites for small −logP)
		cols_neglog <- colvecBase
		allWhite_log <- FALSE
		if (maxPcolor == 1) {
		  cols_neglog <- c("#FFFFFF", cols_neglog)
		} else {
		  maxNeg <- max(NegLogCorr[is.finite(NegLogCorr)], na.rm = TRUE)
		  thrNeg <- -log10(maxPcolor)
		  if (!is.finite(maxNeg) || maxNeg <= 0) { maxNeg <- thrNeg }
		  if (maxNeg < thrNeg) {
		    cols_neglog <- rep("white", length(cols_neglog))
		    allWhite_log <- TRUE
		  } else {
		    whiteProp_left <- max(0, min(1, thrNeg / maxNeg))
		    nWhite <- if (whiteProp_left >= 1 - 1e-12) length(cols_neglog)
		              else ceiling((whiteProp_left / (1 - whiteProp_left)) * length(cols_neglog))
		    cols_neglog <- c(rep("white", nWhite), cols_neglog)
		  }
		}
		## Build palette for unlogged p (top padding = whites for p >= maxPcolor)
		cols_unlog <- rev(colvecBase)
		allWhite_p  <- FALSE
		whiteProp_top <- max(0, min(1, 1 - maxPcolor))
		if (whiteProp_top >= 0.999) {
		  cols_unlog <- rep("white", length(cols_unlog))
		  allWhite_p <- TRUE
		} else if (whiteProp_top > 0) {
		  nWhite <- ceiling((whiteProp_top / (1 - whiteProp_top)) * length(cols_unlog))
		  cols_unlog <- c(cols_unlog, rep("white", nWhite))   ## append to TOP end
		}

		if (modulesInMemory) { categoryColorSymbols=paste0("ME",names(moduleList)) } else { if(!length(na.omit(match(orderedLabelsByRelatedness[,2],rownames(NegLogUncorr))))==nrow(orderedLabelsByRelatedness)) { categoryColorSymbols=dummyColors } else { categoryColorSymbols=names(moduleList) } }
		xSymbolsText= ifelse ( rep(modulesInMemory,length(names(moduleList))), paste0(names(moduleList)," ",orderedModules[match(names(moduleList),orderedModules[,2]),1]), names(moduleList) )


		# --- per-page margins in inches (x uses same labels; y = names(categories) this page)
		cexY_page <- scale_cex(max(nchar(names(categories))), length(names(categories)), base = 1.4)

		sz_page <- compute_pdf_and_margins(
		  xLabs = xSymbolsText,
		  yLabs = names(categories),
		  nCols = length(xSymbolsText),
		  nRows = length(names(categories)),
		  cexX  = cexX_global,
		  cexY  = cexY_page,
		  xAngle = 90,
		  asterisksOnly = asterisksOnly,
		  barOption = barOption
		)

		# Set margins in inches to prevent 'figure margins too large'
		par(mfrow = c(verticalCompression, 1))
		par(mai = sz_page$mai)    # <-- replaces previous par(mar = ...)


		# --- *X* dynamic cex for axis labels ---
		nX <- length(xSymbolsText)
		nY <- length(names(categories))
		maxY <- if (nY) max(nchar(names(categories)), na.rm = TRUE) else 0
		maxX <- if (nX) max(nchar(xSymbolsText),      na.rm = TRUE) else 0

		scale_cex <- function(chars, n, base = 1.6, per_char = 0.015, per_n = 0.003, mn = 0.6, mx = 1.8) {
		  out <- base - per_char * chars - per_n * n
		  max(mn, min(mx, out))
		}
		cexX <- scale_cex(maxX, nX, base = 1.5)     # for x labels (bottom)
		cexY <- scale_cex(maxY, nY, base = 1.4)  # for y labels (left)

          ## -------- Unified plotting: choose base P-matrix, transform by legendScale, and use the same
          ##          palette & limits for BOTH the heatmap and the legend.
          baseP <- if (identical(this.heatmapScale, "minusLogFDR")) adjustedPval else FTpVal
          overlayText <- if (identical(this.heatmapScale, "minusLogFDR")) textMatrix1 else textMatrix.p.unadj
  
          if (identical(legendScale, "minusLog")) {
            ## plot in −log10 space
            Z <- -log10(pmax(baseP, 1e-323))
            zlim_use <- c(0, max(Z[is.finite(Z)], na.rm = TRUE))
            cols_use <- make_cols_neglog(zmax = zlim_use[2], thrNeg = -log10(maxPcolor))
            legend_label <- if (identical(this.heatmapScale,"minusLogFDR"))
                               as.expression(bquote(-log[10]~FDR))
                             else as.expression(bquote(-log[10]~P))
          } else {
            ## plot in unlogged p space
            Z <- baseP
            zlim_use <- c(suppressWarnings(min(Z, na.rm = TRUE)), min(1,max(Z[is.finite(Z)], na.rm=TRUE)))
            if (!is.finite(zlim_use[1])) zlim_use[1] <- 0
            cols_use <- make_cols_unlog(pmin = zlim_use[1], pmax = zlim_use[2], maxPcolor = maxPcolor)
            legend_label <- if (identical(this.heatmapScale,"minusLogFDR")) "FDR" else "P Value"
          }
  
          ## Build color matrix exactly as labeledHeatmap would, then repair any
          ## white cells that are actually below the unlogged threshold.
          cm_raw <- WGCNA::numbers2colors(Z, colors = cols_use, signed = FALSE, lim = zlim_use)
          dim(cm_raw) <- dim(Z)

          ## Force any p <= maxPcolor to be colored,
          ## even if binning nudged them into the white band at the top.
          eps <- max(.Machine$double.eps * 32, 1e-12,
                     (abs(zlim_use[2] - zlim_use[1]) / max(1L, length(cols_use))) / 2048)
          ## First non-white body color from the legend palette actually used:
          body_cols <- cols_use[!(toupper(cols_use) %in% c("WHITE", "#FFFFFF"))]
          if (!length(body_cols)) body_cols <- rev(colvecBase)  # fallback
          first_body <- body_cols[1L]

          ## Base p-values (works for either heatmapScale)
          baseP_here <- if (identical(this.heatmapScale, "minusLogFDR")) adjustedPval else FTpVal

          white_mask <- toupper(cm_raw) %in% c("WHITE", "#FFFFFF")
          fix_mask   <- white_mask & is.finite(baseP_here) & (baseP_here <= (maxPcolor + eps))
          if (any(fix_mask, na.rm = TRUE)) cm_raw[fix_mask] <- first_body

          colorMatrix <- cm_raw

          ## Brightness-aware text color on the *final* matrix:
          hex2brightness <- function(hex) { rgb <- col2rgb(hex); (0.299*rgb[1,]+0.587*rgb[2,]+0.114*rgb[3,])/255 }
          brightness <- matrix(hex2brightness(colorMatrix), nrow = nrow(Z))
          textColors <- ifelse(brightness < 0.5, "white", "black"); dim(textColors) <- dim(Z)
  
          labeledHeatmap.txtMatCols(
            Matrix = Z,
            yLabels = names(categories),
            xLabels = categoryColorSymbols,
            xLabelsAngle = 90,
            xSymbols = xSymbolsText,
            xColorLabels = FALSE,
            colors = cols_use,
            colorMatrix = colorMatrix,  # precise colors avoiding any white if p<maxPcolor
            allWhite = if(length(which(cols_use=="white"))==length(cols_use)) { TRUE } else { FALSE },
            colvecBase = colvecBase,
            textMatrix = overlayText,
            txtMatCols = textColors,
            setStdMargins = FALSE,
            cex.text = if (asterisksOnly) 0.6 else 0.43,
            cex.lab.x = cexX,
            cex.lab.y = cexY,
            verticalSeparator.x = c(rep(c(1:length(names(moduleList))), nrow(orderedLabelsByRelatedness))),
            verticalSeparator.col = 1,
            verticalSeparator.lty = 1,
            verticalSeparator.lwd = 1,
            verticalSeparator.ext = 0,
            horizontalSeparator.y = c(rep(c(1:length(names(categories))), nrow(orderedLabelsByRelatedness))),
            horizontalSeparator.col = 1,
            horizontalSeparator.lty = 1,
            horizontalSeparator.lwd = 1,
            horizontalSeparator.ext = 0,
            zlim = zlim_use,
            main = paste0("Enrichment of ", this.heatmapTitle, "\n",
                          "Row marker list(s) file: ", refDataFile, " - Gene Symbol Enr. (", duplicateHandling, ")\n",
                          if(!is.na(categoriesFile)) { paste0("Column (categories) file: ",categoriesFile) }, 
                          if (!is.na(categoriesFile) & !is.na(bkgrFileForCategories)) { " | " }, 
                          if(!is.na(bkgrFileForCategories)) { paste0("Added bkgr: ", bkgrFileForCategories) }, 
                          if(strictSymmetry) { " (Redundancy ACROSS column lists REMOVED)\n" } else { "\n" },
                          if (identical(this.heatmapScale,"minusLogFDR")) {
                            "Heatmap: FDR-corrected p values"
                          } else {
                            "Heatmap: Fisher Exact p value, Uncorrected"
                          },
                          if (!asterisksOnly) " (values shown)" else "",
                          "\n(values colored when below ", maxPcolor, ") ",if(this.heatmapScale=="minusLogFDR") { "FDR" } else { "p" },"<0.05*, <0.01**, <0.005***"),
            cex.main = 1.3,
            this.heatmapScale = this.heatmapScale,
            legendScale = legendScale,
            ## no legend overrides; legend uses the SAME zlim and SAME palette as the heatmap
            legendLabelOverride = legend_label
          )

	} else {  #if barOption==TRUE:  PLOT BAR PLOTS FOR EACH REFERENCE LIST

		# --- per-page margins in inches (x uses same labels; y = names(categories) this page)
		cexY_page <- scale_cex(max(nchar(names(categories))), length(names(categories)), base = 1.4)

		sz_page <- compute_pdf_and_margins(
		  xLabs = xSymbolsText,
		  yLabs = names(categories),
		  nCols = length(xSymbolsText),
		  nRows = 20,  # fixed height;  length(names(categories)),
		  cexX  = cexX_global,
		  cexY  = cexY_page,
		  xAngle = 90,
		  asterisksOnly = asterisksOnly,
		  barOption = barOption
		)

		# Set margins in inches to prevent 'figure margins too large'
		par(mfrow = c(verticalCompression, 1))
		par(mai = sz_page$mai)    # <-- replaces previous par(mar = ...)
		#par(mar=c(15,7,4.5,1))

		# --- *X* dynamic cex for axis labels ---
		nX <- length(xSymbolsText)
#		nY <- length(names(categories))
#		maxY <- if (nY) max(nchar(names(categories)), na.rm = TRUE) else 0
		maxX <- if (nX) max(nchar(xSymbolsText),      na.rm = TRUE) else 0

		scale_cex <- function(chars, n, base = 1.6, per_char = 0.015, per_n = 0.003, mn = 0.6, mx = 1.8) {
		  out <- base - per_char * chars - per_n * n
		  max(mn, min(mx, out))
		}
		cexX <- scale_cex(maxX, nX, base = 1.25)     # for x labels (bottom)
#		cexY <- scale_cex(maxY, nY, base = 1.4)  # for y labels (left)

	
		moduleColors= if (modulesInMemory) { names(moduleList) } else { "bisque4" }  # if (modulesInMemory), expect colors for names(moduleList)
		xSymbolsText= ifelse ( rep(modulesInMemory,length(names(moduleList))), paste0(names(moduleList)," ",orderedModules[match(names(moduleList),orderedModules[,2]),1]), names(moduleList) )
		if (this.heatmapScale=="p.unadj") {
		
			for( i in 1:nrow(NegLogUncorr)) {
				plotting <- NegLogUncorr[i,]
				#*** added to handle inf values
				plotting[which(!is.na(plotting) & !is.finite(plotting))]<- -log10(1e-323)   #minimum double-precision number // vs. <- max(plotting[is.finite(plotting)], na.rm=TRUE)
				cellType <- rownames(NegLogUncorr)[i]
				barplot(plotting,main = cellType, ylab="",cex.names=cexX, las=2, cex.main=1.7, cex.axis=1.4, legend.text=F,col=moduleColors,names.arg=xSymbolsText)
				mtext(side=2, line=2.6, "-log(p)\n(Uncorrected)", col="black", font=2, cex=1.5)
				abline(h=1.3,col="red")
			}
		}
		
		if (this.heatmapScale=="minusLogFDR") {
			for( i in 1:nrow(NegLogCorr)) {
				plotting <- NegLogCorr[i,]
				#*** added to handle inf values
				plotting[which(!is.na(plotting) & !is.finite(plotting))]<- -log10(1e-323)   #minimum double-precision number // vs. <- max(plotting[is.finite(plotting)], na.rm=TRUE)
				cellType <- rownames(NegLogCorr)[i]
				barplot(plotting,main = cellType, ylab="",cex.names=cexX, las=2, cex.main=1.7, cex.axis=1.4, legend.text=F,col=moduleColors,names.arg=xSymbolsText)
				mtext(side=2, line=2.6, "-log(FDR)\n(Benjamini-Hochberg Correction)", col="black", font=2, cex=1.5)
				abline(h=1.3,col="red")
			}
		}
	
	} #ends if(!barOption)
	
	
	#+#+#+#+#+#+#+#+#+#+#+#+#+
} #end for(refDataFile ...
# Optional explicit PDF close (safe if on.exit already fired)
dl <- grDevices::dev.list()
if (!is.null(dl) && pdf_dev %in% dl) grDevices::dev.off(pdf_dev)
}



labeledHeatmap.txtMatCols <- function (Matrix, xLabels, yLabels = NULL, xSymbols = NULL, ySymbols = NULL, 
    colorLabels = NULL, xColorLabels = FALSE, yColorLabels = FALSE, 
    checkColorsValid = TRUE, invertColors = FALSE, setStdMargins = TRUE, 
    xLabelsPosition = "bottom", xLabelsAngle = 45, xLabelsAdj = 1, 
    yLabelsPosition = "left", xColorWidth = 2 * strheight("M"), 
    yColorWidth = 2 * strwidth("M"), xColorOffset = strheight("M")/3, 
    yColorOffset = strwidth("M")/3, colorMatrix = NULL, 
    colors = NULL, allWhite=FALSE, colvecBase= NULL, naColor = "grey", textMatrix = NULL, 
    cex.text = NULL, textAdj = c(0.5, 0.5), cex.lab = NULL, cex.lab.x = cex.lab, 
    cex.lab.y = cex.lab, colors.lab.x = 1, colors.lab.y = 1, 
    font.lab.x = 1, font.lab.y = 1, bg.lab.x = NULL, bg.lab.y = NULL, 
    x.adj.lab.y = 1, plotLegend = TRUE, keepLegendSpace = plotLegend, 
    verticalSeparator.x = NULL, verticalSeparator.col = 1, verticalSeparator.lty = 1, 
    verticalSeparator.lwd = 1, verticalSeparator.ext = 0, verticalSeparator.interval = 0, 
    horizontalSeparator.y = NULL, horizontalSeparator.col = 1, 
    horizontalSeparator.lty = 1, horizontalSeparator.lwd = 1, 
    horizontalSeparator.ext = 0, horizontalSeparator.interval = 0, 
    showRows = NULL, showCols = NULL, txtMatCols = "black", cex.main=1.5,
    this.heatmapScale = NA, legendScale = "minusLog",
    legendColorsOverride = NULL, legendLimOverride = NULL, legendLabelOverride = NULL, legendTickTransform = NULL, ...)
{
    textFnc = match.fun("text")
    if (!is.null(colorLabels)) {
        xColorLabels = colorLabels
        yColorLabels = colorLabels
    }
    if (is.null(yLabels) & (!is.null(xLabels)) & (dim(Matrix)[1] == 
        dim(Matrix)[2])) 
        yLabels = xLabels
    nCols = ncol(Matrix)
    nRows = nrow(Matrix)
    if (length(xLabels) != nCols) 
        stop("Length of 'xLabels' must equal the number of columns in 'Matrix.'")
    if (length(yLabels) != nRows) 
        stop("Length of 'yLabels' must equal the number of rows in 'Matrix.'")
    if (is.null(showRows)) 
        showRows = c(1:nRows)
    if (is.null(showCols)) 
        showCols = c(1:nCols)
    nShowCols = length(showCols)
    nShowRows = length(showRows)
    if (nShowCols == 0) 
        stop("'showCols' is empty.")
    if (nShowRows == 0) 
        stop("'showRows' is empty.")
    if (checkColorsValid) {
        xValidColors = !is.na(match(substring(xLabels, 3), colors()))
        yValidColors = !is.na(match(substring(yLabels, 3), colors()))
    }
    else {
        xValidColors = rep(TRUE, length(xLabels))
        yValidColors = rep(TRUE, length(yLabels))
    }
    if (sum(xValidColors) > 0) 
        xColorLabInd = xValidColors[showCols]
    if (sum(!xValidColors) > 0) 
        xTextLabInd = !xValidColors[showCols]
    if (sum(yValidColors) > 0) 
        yColorLabInd = yValidColors[showRows]
    if (sum(!yValidColors) > 0) 
        yTextLabInd = !yValidColors[showRows]
    if (setStdMargins) {
        if (xColorLabels & yColorLabels) {
            par(mar = c(2, 2, 3, 5) + 0.2)
        }
        else {
            par(mar = c(7, 7, 3, 5) + 0.2)
        }
    }
    xLabels.show = xLabels[showCols]
    yLabels.show = yLabels[showRows]
    if (!is.null(xSymbols)) {
        if (length(xSymbols) != nCols) 
            stop("When 'xSymbols' are given, their length must equal the number of columns in 'Matrix.'")
        xSymbols.show = xSymbols[showCols]
    }
    else xSymbols.show = NULL
    if (!is.null(ySymbols)) {
        if (length(ySymbols) != nRows) 
            stop("When 'ySymbols' are given, their length must equal the number of rows in 'Matrix.'")
        ySymbols.show = ySymbols[showRows]
    }
    else ySymbols.show = NULL
    xLabPos = charmatch(xLabelsPosition, c("bottom", "top"))
    if (is.na(xLabPos)) 
        stop("Argument 'xLabelsPosition' must be (a unique abbreviation of) 'bottom', 'top'")
    yLabPos = charmatch(yLabelsPosition, c("left", "right"))
    if (is.na(yLabPos)) 
        stop("Argument 'yLabelsPosition' must be (a unique abbreviation of) 'left', 'right'")
    if (is.null(colors)) 
        colors = heat.colors(30)
    if (invertColors) 
        colors = rev(colors)
    
    # pick legend label: override > legendScale > legacy behavior
    legend_label_local <- if (!is.null(legendLabelOverride)) {
      legendLabelOverride
    } else if (legendScale == "minusLog") {
      if (this.heatmapScale == "minusLogFDR") as.expression(bquote(-log[10]~FDR)) else as.expression(bquote(-log[10]~p))
    } else {
      if (this.heatmapScale == "minusLogFDR") "FDR P Value" else "P Value"
    }
    labPos = .heatmapWithLegend(Matrix[showRows, showCols, drop = FALSE], 
        signed = FALSE, colorMatrix = colorMatrix, colors = colors, 
        naColor = naColor, cex.legendLabel = cex.main, plotLegend = plotLegend, 
        keepLegendSpace = keepLegendSpace,
        allWhite=allWhite, colvecBase=colvecBase, legendScale=legendScale,  # in case we want to add color to the legend beyond an all-white range plotted
        #legendLabel = if (this.heatmapScale=="minusLogFDR") { as.expression(bquote(-log[10]~FDR)) } else { "P Value" },
        legendLabel = legend_label_local,
        legendSpan = if(dim(Matrix)[1]<4) { "plot+bottom" } else { "plot" }, 
        legendAlign = if(dim(Matrix)[1]<4) { "span" } else { "heatmapCenter" }, 
        legendLabelOverride = legend_label_local, ...)
#        legendColorsOverride = legendColorsOverride,
#        legendLimOverride = legendLimOverride,
#        legendLabelOverride = legendLabelOverride, 
#        tickLabTransform = legendTickTransform, ...) 

    plotbox = labPos$box
    xmin = plotbox[1]
    xmax = plotbox[2]
    ymin = plotbox[3]
    yrange = plotbox[4] - ymin
    ymax = plotbox[4]
    xrange = xmax - xmin
    xLeft = labPos$xLeft
    xRight = labPos$xRight
    yTop = labPos$yTop
    yBot = labPos$yBot
    xspacing = labPos$xMid[2] - labPos$xMid[1]
    yspacing = abs(labPos$yMid[2] - labPos$yMid[1])
    offsetx = .extend(xColorOffset, nCols)[showCols]
    offsety = .extend(yColorOffset, nRows)[showRows]
    xColW = xColorWidth
    yColW = yColorWidth
    textOffsetY = strheight("M") * cos(xLabelsAngle/180 * 
        pi)
    if (any(xValidColors)) 
        offsetx = offsetx + xColW
    if (any(yValidColors)) 
        offsety = offsety + yColW
    extension.left = par("mai")[2] * par("cxy")[1]/par("cin")[1]
    extension.right = par("mai")[4] * par("cxy")[1]/par("cin")[1]
    extension.bottom = par("mai")[1] * par("cxy")[2]/par("cin")[2] - 
        offsetx
    extension.top = par("mai")[3] * par("cxy")[2]/par("cin")[2] - 
        offsetx
    figureBox = par("usr")
    figXrange = figureBox[2] - figureBox[1]
    figYrange = figureBox[4] - figureBox[3]
    if (!is.null(bg.lab.x)) {
        bg.lab.x = .extend(bg.lab.x, nCols)[showCols]
        if (xLabPos == 1) {
            y0 = ymin
            ext = extension.bottom
            sign = 1
        }
        else {
            y0 = ymax
            ext = extension.top
            sign = -1
        }
        figureDims = par("pin")
        angle = xLabelsAngle/180 * pi
        ratio = figureDims[1]/figureDims[2] * figYrange/figXrange
        ext.x = -sign * ext * 1/tan(angle)/ratio
        ext.y = sign * ext * sign(sin(angle))
        offset = offsetx + textOffsetY
        for (cc in 1:nShowCols) polygon(x = c(xLeft[cc], xLeft[cc], 
            xLeft[cc] + ext.x, xRight[cc] + ext.x, xRight[cc], 
            xRight[cc]), y = c(y0, y0 - sign * offset[cc], y0 - 
            sign * offset[cc] - ext.y, y0 - sign * offset[cc] - 
            ext.y, y0 - sign * offset[cc], y0), border = bg.lab.x[cc], 
            col = bg.lab.x[cc], xpd = TRUE)
    }
    if (!is.null(bg.lab.y)) {
        bg.lab.y = .extend(bg.lab.y, nRows)
        reverseRows = TRUE
        if (reverseRows) 
            bg.lab.y = rev(bg.lab.y)
        bg.lab.y = bg.lab.y[showRows]
        if (yLabPos == 1) {
            xl = xmin - extension.left
            xr = xmin
        }
        else {
            xl = xmax
            xr = xmax + extension.right
        }
        for (r in 1:nShowRows) rect(xl, yBot[r], xr, yTop[r], 
            col = bg.lab.y[r], border = bg.lab.y[r], xpd = TRUE)
    }
    colors.lab.x = .extend(colors.lab.x, nCols)[showCols]
    font.lab.x = .extend(font.lab.x, nCols)[showCols]
    if (sum(!xValidColors) > 0) {
        xLabYPos = if (xLabPos == 1) 
            ymin - offsetx - textOffsetY
        else ymax + offsetx + textOffsetY
        if (is.null(cex.lab)) 
            cex.lab = 1
        mapply(textFnc, x = labPos$xMid[xTextLabInd], y = xLabYPos, 
            labels = xLabels.show[xTextLabInd], col = colors.lab.x[xTextLabInd], 
            font = font.lab.x[xTextLabInd], MoreArgs = list(srt = xLabelsAngle, 
                adj = xLabelsAdj, xpd = TRUE, cex = cex.lab.x))
    }
    if (sum(xValidColors) > 0) {
        baseY = if (xLabPos == 1) 
            ymin - offsetx
        else ymax + offsetx
        deltaY = if (xLabPos == 1) 
            xColW
        else -xColW
        rect(xleft = labPos$xMid[xColorLabInd] - xspacing/2, 
            ybottom = baseY[xColorLabInd], xright = labPos$xMid[xColorLabInd] + 
                xspacing/2, ytop = baseY[xColorLabInd] + deltaY, 
            density = -1, col = substring(xLabels.show[xColorLabInd], 
                3), border = substring(xLabels.show[xColorLabInd], 
                3), xpd = TRUE)
        if (!is.null(xSymbols)) 
            mapply(textFnc, x = labPos$xMid[xColorLabInd], y = baseY[xColorLabInd] - 
                textOffsetY - sign(deltaY) * strwidth("M")/3, 
                labels = xSymbols.show[xColorLabInd], col = colors.lab.x[xColorLabInd], 
                font = font.lab.x[xColorLabInd], MoreArgs = list(adj = xLabelsAdj, 
                  xpd = TRUE, srt = xLabelsAngle, cex = cex.lab.x))
    }
    x.adj.lab.y = .extend(x.adj.lab.y, nRows)[showRows]
    if (yLabPos == 1) {
        marginWidth = par("mai")[2]/par("pin")[1] * 
            xrange
    }
    else {
        marginWidth = par("mai")[4]/par("pin")[1] * 
            xrange
    }
    xSpaceForYLabels = marginWidth - 2 * strwidth("M")/3 - 
        ifelse(yValidColors[showRows], yColW, 0)
    xPosOfYLabels.relative = xSpaceForYLabels * (1 - x.adj.lab.y) + 
        offsety
    colors.lab.y = .extend(colors.lab.y, nRows)[showRows]
    font.lab.y = .extend(font.lab.y, nRows)[showRows]
    if (sum(!yValidColors) > 0) {
        if (is.null(cex.lab)) 
            cex.lab = 1
        if (yLabPos == 1) {
            x = xmin - strwidth("M")/3 - xPosOfYLabels.relative[yTextLabInd]
            adj = x.adj.lab.y[yTextLabInd]
        }
        else {
            x = xmax + strwidth("M")/3 + xPosOfYLabels.relative[yTextLabInd]
            adj = 1 - x.adj.lab.y[yTextLabInd]
        }
        mapply(textFnc, y = labPos$yMid[yTextLabInd], labels = yLabels.show[yTextLabInd], 
            adj = lapply(adj, c, 0.5), x = x, col = colors.lab.y[yTextLabInd], 
            font = font.lab.y[yTextLabInd], MoreArgs = list(srt = 0, 
                xpd = TRUE, cex = cex.lab.y))
    }
    if (sum(yValidColors) > 0) {
        if (yLabPos == 1) {
            xl = xmin - offsety
            xr = xmin - offsety + yColW
            xtext = xmin - strwidth("M")/3 - xPosOfYLabels.relative[yColorLabInd]
            adj = x.adj.lab.y[yColorLabInd]
        }
        else {
            xl = xmax + offsety - yColW
            xr = xmax + offsety
            xtext = xmin + strwidth("M")/3 + xPosOfYLabels.relative[yColorLabInd]
            adj = 1 - x.adj.lab.y[yColorLabInd]
        }
        rect(xleft = xl[yColorLabInd], ybottom = rev(labPos$yMid[yColorLabInd]) - 
            yspacing/2, xright = xr[yColorLabInd], ytop = rev(labPos$yMid[yColorLabInd]) + 
            yspacing/2, density = -1, col = substring(rev(yLabels.show[yColorLabInd]), 
            3), border = substring(rev(yLabels.show[yColorLabInd]), 
            3), xpd = TRUE)
        if (!is.null(ySymbols)) 
            mapply(textFnc, y = labPos$yMid[yColorLabInd], labels = ySymbols.show[yColorLabInd], 
                adj = lapply(adj, c, 0.5), x = xtext, col = colors.lab.y[yColorLabInd], 
                font = font.lab.y[yColorLabInd], MoreArgs = list(srt = 0, 
                  xpd = TRUE, cex = cex.lab.y))
    }
    showCols.ext = c(if (1 %in% showCols) 0 else NULL, showCols)
    showCols.shift = if (0 %in% showCols.ext) 
        1
    else 0
    if (length(verticalSeparator.x) > 0) {
        if (any(verticalSeparator.x < 0 | verticalSeparator.x > 
            nCols)) 
            stop("If given. 'verticalSeparator.x' must all be between 0 and the number of columns.")
        colSepShowIndex = which(verticalSeparator.x %in% showCols.ext)
        verticalSeparator.x.show = .restrictIndex(verticalSeparator.x, 
            showCols.ext) - showCols.shift
    }
    else if (verticalSeparator.interval > 0) {
        verticalSeparator.x.show = verticalSeparator.x = seq(from = verticalSeparator.interval, 
            by = verticalSeparator.interval, length.out = floor(length(showCols)/verticalSeparator.interval))
        colSepShowIndex = 1:length(verticalSeparator.x)
    }
    else verticalSeparator.x.show = NULL
    if (length(verticalSeparator.x.show) > 0) {
        nLines = length(verticalSeparator.x)
        vs.col = .extend(verticalSeparator.col, nLines)[colSepShowIndex]
        vs.lty = .extend(verticalSeparator.lty, nLines)[colSepShowIndex]
        vs.lwd = .extend(verticalSeparator.lwd, nLines)[colSepShowIndex]
        vs.ext = .extend(verticalSeparator.ext, nLines)[colSepShowIndex]
        x.lines = ifelse(verticalSeparator.x.show > 0, labPos$xRight[verticalSeparator.x.show], 
            labPos$xLeft[1])
        nLines.show = length(verticalSeparator.x.show)
        for (l in 1:nLines.show) lines(rep(x.lines[l], 2), c(ymin, 
            ymax), col = vs.col[l], lty = vs.lty[l], lwd = vs.lwd[l])
        angle = xLabelsAngle/180 * pi
        if (angle == 0) 
            angle = pi/2
        if (xLabelsPosition == "bottom") {
            sign = 1
            y0 = ymin
            ext = extension.bottom
        }
        else {
            sign = -1
            y0 = ymax
            ext = extension.top
        }
        figureDims = par("pin")
        ratio = figureDims[1]/figureDims[2] * figYrange/figXrange
        ext.x = -sign * ext * 1/tan(angle)/ratio
        ext.y = sign * ext * sign(sin(angle))
        offset = offsetx + textOffsetY
        for (l in 1:nLines.show) lines(c(x.lines[l], x.lines[l], 
            x.lines[l] + vs.ext[l] * ext.x[l]), c(y0, y0 - sign * 
            offset[l], y0 - sign * offset[l] - vs.ext[l] * ext.y[l]), 
            col = vs.col[l], lty = vs.lty[l], lwd = vs.lwd[l], 
            xpd = TRUE)
    }
    showRows.ext = c(if (1 %in% showRows) 0 else NULL, showRows)
    showRows.shift = if (0 %in% showRows.ext) 
        1
    else 0
    if (length(horizontalSeparator.y) > 0) {
        if (any(horizontalSeparator.y < 0 | horizontalSeparator.y > 
            nRows)) 
            stop("If given. 'horizontalSeparator.y' must all be between 0 and the number of rows.")
        rowSepShowIndex = which(horizontalSeparator.y %in% showRows.ext)
        horizontalSeparator.y.show = .restrictIndex(horizontalSeparator.y, 
            showRows.ext) - showRows.shift
    }
    else if (horizontalSeparator.interval > 0) {
        horizontalSeparator.y.show = horizontalSeparator.y = seq(from = horizontalSeparator.interval, 
            by = horizontalSeparator.interval, length.out = floor(length(showRows)/horizontalSeparator.interval))
        rowSepShowIndex = 1:length(horizontalSeparator.y)
    }
    else horizontalSeparator.y.show = NULL
    if (length(horizontalSeparator.y.show) > 0) {
        reverseRows = TRUE
        if (reverseRows) {
            horizontalSeparator.y.show = nShowRows - horizontalSeparator.y.show + 
                1
            y.lines = ifelse(horizontalSeparator.y.show <= nShowRows, 
                labPos$yBot[horizontalSeparator.y.show], labPos$yTop[nShowRows])
        }
        else {
            y.lines = ifelse(horizontalSeparator.y.show > 0, 
                labPos$yBot[horizontalSeparator.y.show], labPos$yTop[1])
        }
        nLines = length(horizontalSeparator.y)
        vs.col = .extend(horizontalSeparator.col, nLines)[rowSepShowIndex]
        vs.lty = .extend(horizontalSeparator.lty, nLines)[rowSepShowIndex]
        vs.lwd = .extend(horizontalSeparator.lwd, nLines)[rowSepShowIndex]
        vs.ext = .extend(horizontalSeparator.ext, nLines)[rowSepShowIndex]
        nLines.show = length(horizontalSeparator.y.show)
        for (l in 1:nLines.show) {
            if (yLabPos == 1) {
                xl = xmin - vs.ext[l] * extension.left
                xr = xmax
            }
            else {
                xl = xmin
                xr = xmax + vs.ext[l] * extension.right
            }
            lines(c(xl, xr), rep(y.lines[l], 2), col = vs.col[l], 
                lty = vs.lty[l], lwd = vs.lwd[l], xpd = TRUE)
        }
    }
    if (!is.null(textMatrix)) {
        if (is.null(cex.text)) 
            cex.text = par("cex")
        if (is.null(dim(textMatrix))) 
            if (length(textMatrix) == prod(dim(Matrix))) 
                dim(textMatrix) = dim(Matrix)
        if (is.null(dim(txtMatCols))) {
            txtMatCols=rep(txtMatCols[1],dim(textMatrix)[1]*dim(textMatrix)[2])
            dim(txtMatCols) = dim(txtMatrix)
        }
        if (!isTRUE(all.equal(dim(textMatrix), dim(Matrix)))) 
            stop("labeledHeatmap: textMatrix was given, but has dimensions incompatible with Matrix.")
        for (rw in 1:nShowRows) for (cl in 1:nShowCols) {
            text(labPos$xMid[cl], labPos$yMid[rw], as.character(textMatrix[showRows[rw], 
                showCols[cl] ]), col=txtMatCols[rw,cl], xpd = TRUE, cex = cex.text, adj = textAdj)
        }
    }
    axis(1, labels = FALSE, tick = FALSE)
    axis(2, labels = FALSE, tick = FALSE)
    axis(3, labels = FALSE, tick = FALSE)
    axis(4, labels = FALSE, tick = FALSE)
    invisible(labPos)
}


.heatmapWithLegend = function(data, signed, 
                     colorMatrix = NULL,
                     colors, naColor = "grey", zlim = NULL, 
                     reverseRows = TRUE,
                     plotLegend = TRUE,
                     keepLegendSpace = plotLegend,
                     cex.legendAxis = 1, 
                     legendShrink = 0.94,
                     legendPosition = 0.5, ## center; 1 means at the top, 0 means at the bottom
                     legendLabel = "",
                     cex.legendLabel = 1,
                     ## The following arguments are now in inches
                     legendSpace = 0.5 + (as.character(legendLabel)!="") * 1.5*
                            strheight("M",units = "inch", cex = cex.legendLabel),   
                     legendWidth = 0.13,
                     legendGap = 0.09,
                     maxLegendSize = 4,
                     legendLengthGap = 0.15,
                     frame = TRUE,
                     frameTicks = FALSE, tickLen = 0.09,
                     tickLabelAngle = -90,  #previously oriented vertically: 0
                     legendSpan = c("plot","figure","plot+bottom","plot+top"),
                     legendAlign = c("heatmapCenter","span"), # can now center legend on heatmap rows (used with >4 rows in heatmap); spans to bottom with 4 or fewer rows
                     allWhite = FALSE,
                     colvecBase = NULL,
                     legendScale = NULL,
                     legendColorsOverride = NULL,
                     legendLimOverride    = NULL,
                     legendLabelOverride  = NULL,
                     tickLabTransform = NULL,
                     ...)
{
 
  if (length(naColor)==0) naColor = 0;  ### Means transparent (as opposed to white) color.
  # apply legend overrides *without* altering the heatmap color mapping
  legend_colors_final <- if (is.null(legendColorsOverride)) colors else legendColorsOverride
  legend_lim_final    <- if (is.null(legendLimOverride))    zlim   else legendLimOverride
  legend_label_final  <- if (is.null(legendLabelOverride))  ""     else legendLabelOverride

  data = as.matrix(data); nCols = ncol(data); nRows = nrow(data);
  if (is.null(zlim)) 
  {
    zlim = range(data, na.rm = TRUE);
    if (signed) zlim = c(-max(abs(zlim)), max(abs(zlim)));
  }

  barplot(1, col = "white", border = "white", axisnames = FALSE,
                  axes = FALSE, ...);

  pin = par("pin");
  box = par("usr");
  xminAll = box[1]; 
  xmaxAll = box[2]; 
  yminAll = box[3]; 
  ymaxAll = box[4]; 

  legendSpace.usr = legendSpace/pin[1] * (xmaxAll-xminAll);
  legendWidth.usr = legendWidth/pin[1] * (xmaxAll-xminAll);
  legendGap.usr = legendGap/pin[1] * (xmaxAll-xminAll);
  tickLen.usr = tickLen/pin[1] * (xmaxAll-xminAll);
  maxLegendSize.usr = maxLegendSize/pin[2] * (ymaxAll-yminAll);
  legendLengthGap.usr = legendLengthGap/pin[2] * (ymaxAll-yminAll)

  if (!keepLegendSpace && !plotLegend)
  {
     legendSpace.usr = 0;
     legendWidth.usr = 0;
     legendGap.usr = 0;
  }

  ymin = yminAll; 
  ymax = ymaxAll; 
  xmin = xminAll; 
  xmax = xmaxAll - legendSpace.usr;
  if (xmax < xmin) stop("'legendSpace is too large, not enough space for the heatmap."); 

  xStep = (xmax - xmin)/nCols; 
  xLeft = xmin + c(0:(nCols-1)) * xStep;
  xRight = xLeft + xStep; 
  xMid = (xLeft + xRight)/2;

  yStep = (ymax - ymin)/nRows; yBot  = ymin + c(0:(nRows-1)) * yStep;
  yTop  = yBot + yStep; yMid = c(yTop+ yBot)/2;

  
  if (isTRUE(allWhite)) {
    colorMatrix <- matrix("white", nrow = nRows, ncol = nCols)
  } else {
    if (is.null(colorMatrix)) colorMatrix = numbers2colors(data, signed, colors = colors, lim = zlim, naColor = naColor)
    dim(colorMatrix) = dim(data);
    if (reverseRows) colorMatrix = .reverseRows(colorMatrix);
  }

  for (c in 1:nCols)  # draw heatmap cells
  {
    rect(xleft = rep(xLeft[c], nRows), xright = rep(xRight[c], nRows),
         ybottom = yBot, ytop = yTop, col = ifelse(colorMatrix[, c]==0, 0, colorMatrix[, c]), 
                border = ifelse(colorMatrix[, c]==0, 0, colorMatrix[, c]));
    ## Note: the ifelse seems superfluous here but it essentially converts a potentially character "0" to the number 0
    ## which the plotting system should understand as transparent color.
  }

  if (frame) lines( c(xmin, xmax, xmax, xmin, xmin), c(ymin, ymin, ymax, ymax, ymin) );

  if (plotLegend)
  {
      # Now plot the legend.
#      legendSize.usr = legendShrink * (ymaxAll - yminAll);
#      if (legendSize.usr > maxLegendSize.usr) legendSize.usr = maxLegendSize.usr
#      if (legendLengthGap.usr > 0.5*(ymaxAll - yminAll)*(1-legendShrink)) 
#          legendLengthGap.usr = 0.5*(ymaxAll - yminAll)*(1-legendShrink);
#      y0 = yminAll + legendLengthGap.usr;
#      y1 = ymaxAll - legendLengthGap.usr;
#      movementRange = (y1-y0 - legendSize.usr);
#      if (movementRange < -1e-10) {browser(".heatmapWithLegend: movementRange is negative."); movementRange = 0;}
#      ymin.leg = y0 + legendPosition * movementRange;
#      ymax.leg = y0 + legendPosition * movementRange + legendSize.usr
#      legendPosition = .plotColorLegend(xmin = xmaxAll - (legendSpace.usr - legendGap.usr),
#                       xmax = xmaxAll - (legendSpace.usr - legendGap.usr - legendWidth.usr),
#                       ymin = ymin.leg,
#                       ymax =  ymax.leg,
#                       lim = zlim,
#                       colors = colors,
#                       tickLen.usr = tickLen.usr,
#                       cex.axis = cex.legendAxis,
#                       lab = legendLabel,
#                       cex.lab = cex.legendLabel,
#                       tickLabelAngle = tickLabelAngle
#                       );

    legendSpan <- match.arg(legendSpan)
    legendAlign <- match.arg(legendAlign)

    pin <- par("pin")                    # inner plot size [inches]
    box <- par("usr")                    # plot region coordinates
    xminAll <- box[1]; xmaxAll <- box[2]
    yminAll <- box[3]; ymaxAll <- box[4]
    yRange  <- ymaxAll - yminAll

    mai <- par("mai")                    # inches
    bottom_usr <- mai[1] / pin[2] * yRange
    top_usr    <- mai[3] / pin[2] * yRange

    # space (usr) on the right for legend band
    legendSpace.usr <- legendSpace / pin[1] * (xmaxAll - xminAll)
    legendWidth.usr <- legendWidth / pin[1] * (xmaxAll - xminAll)
    legendGap.usr   <- legendGap   / pin[1] * (xmaxAll - xminAll)
    tickLen.usr     <- tickLen     / pin[1] * (xmaxAll - xminAll)
    maxLegendSize.usr   <- maxLegendSize    / pin[2] * yRange
    legendLengthGap.usr <- legendLengthGap  / pin[2] * yRange

    if (!keepLegendSpace && !plotLegend) {
      legendSpace.usr <- legendWidth.usr <- legendGap.usr <- 0
    }

    # Heatmap plot box (no margins)
    ymin <- yminAll
    ymax <- ymaxAll

    # Choose legend vertical span in usr
    span_min <- switch(legendSpan,
      "plot"        = ymin,
      "figure"      = yminAll - bottom_usr,
      "plot+bottom" = yminAll - bottom_usr,
      "plot+top"    = yminAll
    )
    span_max <- switch(legendSpan,
      "plot"        = ymax,
      "figure"      = ymaxAll + top_usr,
      "plot+bottom" = ymaxAll,
      "plot+top"    = ymaxAll + top_usr
    )
    spanRange <- span_max - span_min

    # Legend height
    legendSize.usr <- legendShrink * spanRange
    if (legendSize.usr > maxLegendSize.usr) legendSize.usr <- maxLegendSize.usr

    # --- Vertical positioning
    if (legendAlign == "heatmapCenter") {
      # center on the heatmap rows (plot region), then clamp into span
      y_center <- (ymin + ymax) / 2
      ymin.leg <- y_center - legendSize.usr / 2
      ymax.leg <- y_center + legendSize.usr / 2

      lowBound  <- span_min + legendLengthGap.usr
      highBound <- span_max - legendLengthGap.usr
      if (ymin.leg < lowBound) {
        shift <- lowBound - ymin.leg
        ymin.leg <- ymin.leg + shift;  ymax.leg <- ymax.leg + shift
      }
      if (ymax.leg > highBound) {
        shift <- ymax.leg - highBound
        ymin.leg <- ymin.leg - shift;  ymax.leg <- highBound
      }
    } else {
      # prior behavior: position within the span by legendPosition
      y0 <- span_min + legendLengthGap.usr
      y1 <- span_max - legendLengthGap.usr
      movementRange <- max(0, (y1 - y0 - legendSize.usr))
      ymin.leg <- y0 + legendPosition * movementRange
      ymax.leg <- ymin.leg + legendSize.usr
    }

    # Draw legend (ticks horizontal because tickLabelAngle = -90 => angle 0 for vertical legend)
    # Use overrides if provided; otherwise use heatmap palette and limits
    legend_colors_use <- if (is.null(legendColorsOverride)) colors else legendColorsOverride
    legend_lim_use    <- if (is.null(legendLimOverride))    zlim   else legendLimOverride
    if (isTRUE(allWhite)) {
      ## Build a legend that shows the (all-white) heatmap range first,
      ## then appends the body palette for values ABOVE that range.
      ## Use the current legend palette length as the size of the white band.
      base_cols <- tryCatch(get("colvecBase", inherits = TRUE), error = function(e) NULL)
      if (is.null(base_cols) || !length(base_cols)) {
        base_cols <- if (is.matrix(colors)) as.vector(colors[, 1]) else colors
      }
      nWhiteBand <- if (is.matrix(legend_colors_use)) nrow(legend_colors_use) else length(legend_colors_use)
      ## If legend is in −log10 space, extend the upper limit so non-white colors are shown.
      is_minuslog <- legendScale=="minusLog"  #!is.null(tickLabTransform) ||
                     #(is.language(legend_label_final) &&
                     # grepl("-log", paste(deparse(legend_label_final), collapse = "")))
      if (is.null(legend_lim_use) || length(legend_lim_use) != 2L || any(!is.finite(legend_lim_use))) {
        legend_lim_use <- if (!is.null(zlim)) zlim else c(0, 1)
      }
      low  <- ifelse(is.finite(legend_lim_use[1]), legend_lim_use[1], 0)  # should always be finite
      high <- ifelse(is.finite(legend_lim_use[2]), legend_lim_use[2], 1)
      if (is_minuslog) {
        legend_colors_use <- c(rep("white", nWhiteBand), base_cols)
        legend_lim_use <- c(low, max(high * 2, low + (high - low) * 2))
      } else {  # unlogged legend scale
        base_cols100<-grDevices::colorRampPalette(base_cols)(100)
        nWhiteBand.percent <- 100 - floor(range(data, na.rm=T)[1] * 100) 
        legend_colors_use <- c(rev(base_cols100), rep("white", ceiling(nWhiteBand.percent/(1-(nWhiteBand.percent/100)))) )
        legend_lim_use <- c(0, 1)
      }
    }
    legendPosition <- .plotColorLegend(
      xmin = xmaxAll - (legendSpace.usr - legendGap.usr),
      xmax = xmaxAll - (legendSpace.usr - legendGap.usr - legendWidth.usr),
      ymin = ymin.leg, ymax = ymax.leg,
      #lim = zlim, colors = colors, tickLen.usr = tickLen.usr,
      lim = legend_lim_use,
      colors = legend_colors_use, tickLen.usr = tickLen.usr,
      cex.axis = cex.legendAxis, lab = legend_label_final, # legendLabel,
      cex.lab = cex.legendLabel, tickLabelAngle = tickLabelAngle
    )
    
  } else legendPosition = NULL

  invisible(list(xMid = xMid, yMid = if (reverseRows) rev(yMid) else yMid, 
       box = c(xmin, xmax, ymin, ymax), xLeft = xLeft, xRight = xRight,
       yTop = yTop, yBot = yBot,
       legendPosition = legendPosition));
  
}


.reverseRows = function(Matrix)
{
  ind = seq(from=dim(Matrix)[1], to=1, by=-1);
  Matrix[ind,, drop = FALSE];
  #Matrix
}

.extend = function(x, n)
{
  nRep = ceiling(n/length(x));
  rep(x, nRep)[1:n];
}

# Adapt a numeric index to a subset
# Aim: if 'index' is a numeric index of special entries of a vector,
#    create a new index that references 'subset' elements of the vector  
.restrictIndex = function(index, subset)
{
  out = match(index, subset);
  out[!is.na(out)];
}


.autoTicks = function(min, max, maxTicks = 6, tickPos = c(1,2,5))
{
  if (max < min) { x = max; max = min; min = x }
  range = max - min;
  if (range==0) return(max);
  tick0 = range/(maxTicks+1-1e-6)
  maxTick = max(tickPos);
  # Ticks can only be multiples of tickPos
  mult = 1;
  if (tick0 < maxTick/10)
  {
     while (tick0 < maxTick/10) {tick0 = 10*tick0; mult = mult*10; }
  } else
     while (tick0 >=maxTick ) {tick0 = tick0/10; mult = mult/10; }

  ind = sum(tick0 > tickPos) + 1;
  tickStep = tickPos[ind] / mult;

  lowTick = min/tickStep;
  if (floor(lowTick)!=lowTick) lowTick = lowTick + 1;
  lowTick = floor(lowTick);

  ticks = tickStep * (lowTick:(lowTick + maxTicks+1));
  ticks = ticks[ticks <= max];
  ticks;
}

.plotStandaloneLegend = function(
                            colors,
                            lim,
                            ## These dimensions are in inches
                            tickLen = 0.09,
                            tickGap = 0.04,
                            minBarWidth = 0.09,
                            maxBarWidth = Inf,
                            mar = c(0.5, 0.2, 0.5, 0.1),
                            lab = "",
                            horizontal = FALSE,
                            ...)
{
  par(mar = mar);
  plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "");
  box = par("usr");
  if (horizontal) box.eff = box[c(3,4,1,2)] else box.eff = box;
  tickVal = .autoTicks(lim[1], lim[2]);
  pin = par("pin");
  pin.eff = if (horizontal) pin[c(2,1)] else pin;
  wrange = box.eff[2] - box.eff[1];
  tickLen.usr = tickLen/pin.eff[1] * wrange
  tickGap.usr = tickGap/pin.eff[1] * wrange
  minBarWidth.usr = minBarWidth/pin.eff[1] * wrange
  maxBarWidth.usr = maxBarWidth/pin.eff[1] * wrange
  sizeFnc = if (horizontal) strheight else strwidth;
  maxTickWidth = max(sizeFnc(tickVal));
  if (maxTickWidth + tickLen.usr + tickGap.usr > box.eff[2]-box.eff[1]-minBarWidth.usr) 
     warning("Some tick labels will be truncated.");
  haveLab = length(lab) > 0
  if (haveLab && is.character(lab)) haveLab = lab!="";
  width = max(box.eff[2]-box.eff[1]-maxTickWidth - tickLen.usr - tickGap.usr- haveLab * 3*sizeFnc("M"), minBarWidth.usr);
  if (width > maxBarWidth.usr) width = maxBarWidth.usr;
  .plotColorLegend(box[1], if (horizontal) box[2] else box[1] + width,
                   if (horizontal) box[4]-width else box[3], box[4], 
                   colors = colors,
                   lim = lim,
                   tickLen.usr = tickLen.usr, horizontal = horizontal,
                   tickGap.usr = tickGap.usr, lab = lab, ...);
}

if (FALSE)
{
   source("~/Work/RLibs/WGCNA/R/heatmapWithLegend.R")
   .plotStandaloneLegend(colors = blueWhiteRed(10), lim = c(-25, 25))
   d = matrix(rnorm(100), 10, 10);
   par(mar = c(2,2,2,0));
   
   .heatmapWithLegend(d,
                     signed = TRUE,
                     colors = blueWhiteRed(20), 
                     plotLegend = TRUE,
                     cex.legendAxis = 1,
                     legendShrink = 0.94,
                     legendLabel = "",
                     cex.legendLabel = 1)
                     ## The following arguments are now in inches
                     #legendSpace = 0.5 + (legendLabel!="") * 1.5*strheight("M",units = "inch", cex = cex.legendLabel),
                     #legendWidth = 0.13,
                     #legendGap = 0.09,
                     #frame = TRUE,
                     #frameTicks = FALSE, tickLen = 0.09);

}

.plotColorLegend = function(xmin, xmax, ymin, ymax,
                            # colors can be a vector or a matrix (in which case a matrix of colors will be plotted)
                            colors,
                            horizontal = FALSE,
### FIXME: it would be good if these could respect settings in par("mgp")
                            tickLen.usr = 0.5* (if (horizontal) strheight("M") else strwidth("M")),
                            tickGap.usr = 0.5 * (if (horizontal) strheight("M") else strwidth("M")),
                            lim, cex.axis = 1, tickLabelAngle = if (horizontal) 0 else -90,
                            lab = "", cex.lab = 1, labAngle = 0, 
                            labGap = 0.6 * (if (horizontal) strheight("M") else strwidth("M")),
                            tickLabTransform = NULL
                            )
{
  tickVal = .autoTicks(lim[1], lim[2]);
  tickLab = if (is.null(tickLabTransform)) tickVal else tickLabTransform(tickVal)
  nTicks = length(tickVal);

  if (horizontal) {
    lmin = xmin; lmax = xmax; 
    tmin = ymin; tmax = ymax;
  } else {
    tmin = xmin; tmax = xmax; 
    lmin = ymin; lmax = ymax;
  }
  tickPos = (tickVal - lim[1]) / (lim[2] - lim[1]) * (lmax - lmin) + lmin;
  pin = par("pin");
  box = par("usr");
  asp = pin[2]/pin[1] * ( box[2]-box[1])/(box[4] - box[3]);
  # Ticks:
  
  if (horizontal) {
    angle0 = 0;
    angle = angle0 + tickLabelAngle;
    if (angle==0) adj = c(0.5, 1) else adj = c(1, 0.5);
    for (t in 1:nTicks) 
      lines(c(tickPos[t], tickPos[t]), c(ymin, ymin - tickLen.usr), xpd = TRUE);
    text(tickPos, rep(ymin - tickLen.usr - tickGap.usr), tickVal, adj = adj, cex = cex.axis,
           xpd = TRUE, srt = angle);
    tickLabelWidth = if (angle==0) max(strheight(tickVal)) else max(strwidth(tickVal))/asp;
  } else {
    angle0 = 90;
    angle = angle0 + tickLabelAngle;
    if (angle==0) adj = c(0, 0.5) else adj = c(0.5, 1);
    for (t in 1:nTicks) 
      lines(c(xmax, xmax + tickLen.usr), c(tickPos[t], tickPos[t]), xpd = TRUE);
    text(rep(xmax + tickLen.usr + tickGap.usr), tickPos, tickVal, adj = adj, cex = cex.axis,
         xpd = TRUE, srt = angle);
    tickLabelWidth = if (angle==0) max(strwidth(tickVal)) else max(strheight(tickVal)) * asp;
  }
  # Fill with color:
  colors = as.matrix(colors);
  nColumns = ncol(colors);
  nColors = nrow(colors);
  bl = (lmax-lmin)/nColors * (0:(nColors-1)) + lmin;
  tl = (lmax-lmin)/nColors * (1:nColors) + lmin;
  wi.all = tmax - tmin;
  wi1 = wi.all/nColumns
  if (horizontal) {
    for (col in 1:nColumns)
      rect(xleft = bl, xright = tl,
         ybottom = rep(tmin + (col-1) * wi1, nColors), ytop = rep(tmin + wi1*col, nColors), 
            col = colors[, col], border = colors[, col], xpd = TRUE);
  } else {
    for (col in 1:nColumns)
       rect(xleft = rep(tmin + (col-1) * wi1, nColors), xright = rep(tmin + wi1*col, nColors),
          ybottom = bl, ytop = tl, col = colors[, col], border = colors[, col], xpd = TRUE);
  }
  # frame for the legend
  lines(c(xmin, xmax, xmax, xmin, xmin), c(ymin, ymin, ymax, ymax, ymin), xpd = TRUE );

  if (nColumns > 1) for (col in 2:nColumns) 
    if (horizontal) lines(c(xmin, xmax), c(tmin + (col-1) * wi1, tmin + (col-1) * wi1)) else 
                    lines(c(tmin + (col-1) * wi1, tmin + (col-1) * wi1), c(ymin, ymax));
  # Axis label
  if (length(lab)>0 && as.character(lab) != "")
  {
    if (horizontal)
    {
      y = ymin - tickLen.usr - tickGap.usr - tickLabelWidth - labGap;
      x = (xmin + xmax)/2;
      adj = if (labAngle==0) c(0.5, 1) else c(1, 0.5)
      angle = labAngle;
      text(x, y, lab, cex = cex.lab, srt = labAngle, xpd = TRUE, adj = adj);
    } else {
      y = (ymin + ymax)/2;
      x = xmax + tickLen.usr + tickGap.usr + tickLabelWidth + labGap;
      adj = if (labAngle==0) c(0.5, 1) else c(0, 0.5);
      angle = labAngle+90;
      text(x, y, lab, cex = cex.lab, srt = labAngle+90, xpd = TRUE, adj = adj);
    }
    height = strheight(lab);
    if (!horizontal) height = height * asp;
    labelInfo = list(x = x, y = y, angle = angle, adj = adj,
                     space.usr = height, gap.usr = labGap);
  } else labelInfo = list(space.usr = 0, gap.usr = 0);
  #### FIXME: also include a component named box that gives the outer coordinates of the area used by the legend, to the
  ###best approximation. Maybe include the padding around the color bar.
  invisible(list(bar = list(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, 
                            space.usr = tmax - tmin),
       ticks = list(length.usr = tickLen.usr, gap.usr = tickGap.usr, labelSpace.usr = tickLabelWidth),
       label = labelInfo));
}



.boxDimensionsForHeatmapWithLegend = function(
                     data,
                     plotLegend = TRUE,
                     keepLegendSpace = plotLegend,
                     cex.legend = 1,
                     legendShrink = 0.94,
                     ## The following arguments are now in inches
                     legendSpace = 0.5,
                     legendWidth = 0.13,
                     legendGap = 0.09, 
                     startTempPlot = TRUE,
                     plotDevice = "pdf",
                     plotDeviceOptions = list(),
                     width = 7, height = 7,...)
{
  data = as.matrix(data); nCols = ncol(data); nRows = nrow(data);

  if (startTempPlot)
  {
    if (!is.null(plotDevice))
    {
      if (plotDevice == "x11") 
      {
        do.call(match.fun(plotDevice), c(list(width = width, height = height), plotDeviceOptions));
        on.exit(dev.off());
      } else {
        file = tempfile();
        do.call(match.fun(plotDevice), c(list(file = file, width = width, height = height), plotDeviceOptions))
        on.exit({ dev.off(); unlink(file)});
      }
      par(mar = par("mar"));
    }
    barplot(1, col = "white", border = "white", axisnames = FALSE,
                  axes = FALSE, ...);
  }
  pin = par("pin");
  box = par("usr");
  xminAll = box[1];
  xmaxAll = box[2];
  yminAll = box[3];
  ymaxAll = box[4];

  legendSpace.usr = legendSpace/pin[1] * (xmaxAll-xminAll);
  legendWidth.usr = legendWidth/pin[1] * (xmaxAll-xminAll);
  legendGap.usr = legendGap/pin[1] * (xmaxAll-xminAll);

  if (!keepLegendSpace && !plotLegend)
  {
     legendSpace.usr = 0;
     legendWidth.usr = 0;
     legendGap.usr = 0;
  }

  ymin = yminAll;
  ymax = ymaxAll;
  xmin = xminAll;
  xmax = xmaxAll - legendSpace.usr;
  if (xmax < xmin) stop("'legendSpace is too large, not enough space for the heatmap.");
  xStep = (xmax - xmin)/nCols;
  xLeft = xmin + c(0:(nCols-1)) * xStep;
  xRight = xLeft + xStep;
  xMid = (xLeft + xRight)/2;

  yStep = (ymax - ymin)/nRows; yBot  = ymin + c(0:(nRows-1)) * yStep;
  yTop  = yBot + yStep; yMid = c(yTop+ yBot)/2;

  list(xMin = xmin, xMax = xmax, yMin = ymin, yMax = ymax,
       xLeft = xLeft, xRight = xMid, xMid = xMid,
       yTop = yTop, yMid = yMid, yBottom = yBot);
}


# --- measure a vector of labels in inches on a temp device
.measure_text_inches <- function(labels, angle = 0, cex = 1, fam = "", which = c("width","height")) {
  which <- match.arg(which)
  if (length(labels) == 0L) return(0)

  tf <- tempfile(fileext = ".pdf")

  # Remember the previously active device (if any)
  prev_dev <- try(grDevices::dev.cur(), silent = TRUE)

  # Open a throwaway PDF device just to get font metrics
  grDevices::pdf(tf, width = 7, height = 7)
  tmp_dev <- grDevices::dev.cur()

  op <- graphics::par(family = fam)

  # Always restore/cleanup the specific temp device we opened
  on.exit({
    try(graphics::par(op), silent = TRUE)
    try(grDevices::dev.off(tmp_dev), silent = TRUE)
    try(unlink(tf), silent = TRUE)
    # Restore previous device if it still exists
    dl <- grDevices::dev.list()
    if (!inherits(prev_dev, "try-error") && !is.null(dl) && prev_dev %in% dl) {
      try(grDevices::dev.set(prev_dev), silent = TRUE)
    }
  }, add = FALSE)

  if (which == "width") {
    max(strwidth(labels, units = "inches", cex = cex, srt = angle))
  } else {
    max(strheight(labels, units = "inches", cex = cex, srt = angle))
  }
}

# --- compute required PDF size and par(mai) given labels + grid size
compute_pdf_and_margins <- function(
  xLabs, yLabs,
  nCols, nRows,
  cexX = 1, cexY = 1,
  xAngle = 90,
  # legend / padding
  legend_space_in = 0.60,
  pad_in         = 0.33,
  title_in       = 1.40,
  right_pad_in   = 0.25,
  # cell-size constraints
  min_cell_w_in  = 0.20,
  min_cell_h_in  = 0.16,
  max_cell_w_in  = 0.40,    # NEW: cap width per column (optional)
  max_cell_h_in  = 0.24,    # NEW: cap height per row   (requested)
  # plot-area floors (heatmap area only; margins are added on top)
  min_total_plot_w_in = 6.0,
  min_total_plot_h_in = 3.0,
  # absolute hard caps on *total* device (optional)
  max_total_width_in  = 60,
  max_total_height_in = 60,
  asterisksOnly=FALSE,
  barOption=FALSE
) {
  if(!asterisksOnly) min_cell_w_in = min_cell_w_in * 3.3
  if(barOption) {
    min_total_plot_h_in=min_total_plot_h_in * 2.2
    nRows = 20
  }

  # label extents in inches; for x at 90°, vertical need is the *width* of text
  x_label_vert_in <- .measure_text_inches(xLabs, angle = xAngle, cex = cexX, which = "width")
  y_label_horiz_in <- .measure_text_inches(yLabs, angle = 0,      cex = cexY, which = "width")

  # margins in inches
  bottom_in <- pad_in + x_label_vert_in + 0.25
  left_in   <- pad_in + y_label_horiz_in + 0.25
  top_in    <- title_in
  right_in  <- legend_space_in + right_pad_in

  # ---- plot area height/width (heatmap grid only) with min & max per-cell clamps
  plot_w_min <- max(min_total_plot_w_in, nCols * min_cell_w_in)
  plot_w_max <- if (is.finite(max_cell_w_in)) nCols * max_cell_w_in else Inf
  plot_w_in  <- if (is.finite(plot_w_max)) min(plot_w_min, plot_w_max) else plot_w_min

  plot_h_min <- max(min_total_plot_h_in, nRows * min_cell_h_in)
  plot_h_max <- if (is.finite(max_cell_h_in)) nRows * max_cell_h_in else Inf
  plot_h_in  <- if (is.finite(plot_h_max)) min(plot_h_min, plot_h_max) else plot_h_min

  # total device size
  width_in  <- left_in + plot_w_in + right_in
  height_in <- bottom_in + plot_h_in + top_in

  # clamp to sane bounds
  width_in  <- max(8,  min(max_total_width_in,  width_in))
  height_in <- max(2.2,  min(max_total_height_in, height_in))

  list(width = width_in, height = height_in, mai = c(bottom_in, left_in, top_in, right_in))
}
