#!/usr/bin/env $CODEBASE/Rscript
#
# Define constants and methods used in the analysis of the ATAC-seq data.
#
# run("atac", "init")


###################################
### ATAC-SEQ-SPECIFIC LIBARIES ###
###################################

if(lib$projectName!="lch") run("common")
setCurrentAnalysis("")
loadLibraries(c("fastcluster", "cluster","data.table", "ggplot2", "LOLA", "simpleCache", "GenomicRanges")) #"ChIPpeakAnno", "org.Hs.eg.db","DiffBind", "limma", 




############################
### FUNCTION DEFINITIONS ###
############################

summarizeReads <- function(regs, sample.names, file.paths, useDuplicates=F) {
	loadLibraries(c("Rsamtools", "GenomicAlignments"))
	
	rblapply(1:length(sample.names), function(i) {
		msg(sample.names[i])
		so <- as.data.table(assay(summarizeOverlaps(features = regs, reads = BamFileList(file.paths[i], yieldSize = 5000000), ignore.strand = TRUE, singleEnd = TRUE, fragments = FALSE, param = ScanBamParam(flag = scanBamFlag(isDuplicate = useDuplicates)))))
		setnames(so, c("count"))
		so[,regNum:=1:nrow(so)]
		so[,sampleName:=sample.names[i]]
		so
	}, "sampleNum")
}
loadUCSCRepeats <- function() { 
	loadLibrary("LOLA")
	stopifnot(requireNamespace("LOLA"))
	repeats <- LOLA::getRegionSet(lib$resourceDir("regions/LOLACore/", genomeVer), collections="ucsc_features", "rmsk.bed")[[1]]
	return(repeats)
}
# define a color as a weighted mixture of pre-defined colors (used to get colors for regulatory modules by cell types in which they are accessible):
makeColor <- function(x, annot=dA, grpBy="sample_group", cols=colorPalettes[[grpBy]]) {
	grps <- annot[,get(grpBy)]
	x <- sapply(split(x, f=grps), mean)
	x <- x-min(x)
	x <- (x/max(x)*10)^2	
	grps <- unique(grps)
	alpha(lib$blendColors(unlist(sapply(1:length(grps), function(i) rep(cols[grps[i]], round(x[i]))))), 1)
}
writeFASTA <- function(sequences, sequenceFile) {
	msg("\t* create input file: ", sequenceFile)
	tmp <- paste0(rep(">seq", length(sequences)*2), floor(1:(length(sequences)*2) / 2 +1))
	tmp[seq(2,length(tmp),by=2)] <- as.character(sequences)
	write.table(tmp, file=sequenceFile, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
	sequenceFile
}
readFIMOResult <- function(fimoResultFile, motifFile=paste0(Sys.getenv("RESOURCES"),"/motifs/motif_databases_feb2017/JASPAR/JASPAR_CORE_2016_vertebrates.meme"), nSeqs=NULL, motifPThresh=0.05) {
	cmd <- paste("grep 'MOTIF '",motifFile)
	msg(cmd)
	motifNames <- tryCatch(lib$dt2namedVec(fread(cmd, select=c(2,3), header=FALSE), "V2"), error=function(e) {
		tmp <- fread(cmd, select=c(2), header=FALSE)
		structure(tmp[,V2], names=tmp[,V2])
	})
	
	dt <- fread(fimoResultFile, select=c(1,3,8))
		
	if(is.null(nSeqs)) nSeqs <- dt[,length(unique(sequence_name))]
	dt[, motifName:=motifNames[motif_id]]

	dt[`p-value`<=motifPThresh, .(seq=sequence_name, motifName, motifId=motif_id)]
}
runFIMO <- function(sequences, outFile="fimo_out.txt", motifP=1e-4, motifFile=paste0(Sys.getenv("RESOURCES"),"/motifs/motif_databases_feb2017/JASPAR/JASPAR_CORE_2016_vertebrates.meme"), params="--no-qvalue --text --bgfile motif-file", fimoExec= toolsDir("meme/bin/fimo"), inputIsFile=length(sequences)==1) {

	dir.create(dirname(outFile), recursive=TRUE, showWarnings=FALSE)

	if(!inputIsFile) {
		sequenceFile <- writeFASTA(sequences, paste0(outFile, "_input.fasta"))
	}
	else {
		msg("\t* use existing input file: ", sequences)
		sequenceFile <- sequences
	}

	cmd <- paste(fimoExec,paste("--thresh",motifP), params,motifFile,sequenceFile,">",outFile)
	msg("\t* run external tool: $ ", cmd)
	system(cmd)

	if(!inputIsFile) {
		msg("\t* remove temporary input file: ", sequenceFile)
		file.remove(sequenceFile)
	}
	
	outFile
}
augmentLolaRes <- function(lolaResAll, qThresh=0.05, qvalCol="qvalMean", pvalCol="pvalMean", orCol="oddsRatioMean") {
	lolaResAll[,log2odds:=log2(get(orCol))]
	
	# harmonize synonyms:
	lolaResAll[,term:=antibody]
	lolaResAll[is.na(term),term:=gsub("_\\(.+\\)$","",gsub("GSM\\d+_","",gsub("Human_","",gsub("wgEncode.wg","",gsub(".(bed|narrowPeak)","",filename)))))]
	lolaResAll[collection=="sheffield_dnase", term:=paste0("DNase #",gsub(".bed","",filename), " (", sapply(strsplit(description,";"),function(x) paste(substr(x,1,3),collapse=";")), ")")]
	lolaResAll[,term:=gsub("^(EGFP|C)\\W","",gsub("_\\(.+\\)$","",toupper(term)))]
	lolaResAll[,term:=gsub("EP300", "P300", term)]
	lolaResAll[,term:=gsub("PU\\.?1", "SPI1", term)]
	lolaResAll[,term:=gsub("[−-]", "", term)]
	lolaResAll[,term:=gsub("POL(II|LL)", "POL2", term)]
	lolaResAll[,term:=gsub("POLIII", "POL3", term)]
	lolaResAll[antibody%in%c("GR","NR3C1"), antibody:="NR3C1"]
	lolaResAll[grepl("^E?P300",antibody,ignore.case=T), antibody:="EP300"]
	lolaResAll[antibody%in%c("GABP"), antibody:="GABPA"]
	lolaResAll[term=="P300",term:="EP300"]
	lolaResAll[term=="ERALPHA_A",term:="ESR1"]
	lolaResAll[term=="TCF7L2_C9B9",term:="TCF7L2"]

	# kick out Pol2:
	lolaResAll <- lolaResAll[!grepl("^POL",antibody,ignore.case=T),]

	lolaResAll
}
curateRegionDb <- function(regionDB) {
	capFirst <- function(str) paste0(toupper(substr(str, 1, 1)), substr(str, 2, nchar(str)))
	regionDB$regionAnno[,cellType:=capFirst(gsub("(cell|progenitor|precursor|phage|cyte|blast)s","\\1", tolower(cellType), perl=TRUE))]
	regionDB
}
matchLenient <- function(x, y) {
	z <- x[sapply(x, function(xx) {
		sum(sapply(y, function(yy) {
					grepl(paste0("^",gsub("\\d+[ABCDEL]?$","",xx)), yy) || grepl(paste0("^",gsub("\\d+[ABCDEL]?$","",yy)), xx)
		})>0)>0
	})]
	setdiff(z, paste0("SP",1:10))
}



#########################################
### SCRNA-SEQ-SPECIFIC CONFIG OPTIONS ###
#########################################

dA <- loadAnnot("atac")[use==T,]
setkey(dA,sample_name)
sampleToPatient <- lib$dt2namedVec(dA, "sample_name", "pat_id")

lfcThreshDEG <- log2(2)
qThreshDEG <- 1
clustK <- 6
useDuplicates <- FALSE

dA[,peaks_file:=sprintf(dataDir("%s/peaks_redo_q0.1/%s_peaks.narrowPeak"), sample_name, sample_name)]
dA[,reads_file:=sprintf(dataDir("%s/mapped/%s.trimmed.bowtie2.filtered.bam"), sample_name, sample_name)]

setCurrentAnalysis(paste0(cfg$ver,"/atac"))

# keep lists of genes that are expressed in LCH cell subsets:
nCCS <- 7
exprDir <- paste0(cfg$ver, "/10x/traj_cca_",nCCS)
clusterExpr <- list(all=lib$dtToDf(fread(baseResultsDir(exprDir,"/cluster_prototypes_scaled.csv"))))
for(patId in gsub("cluster_prototypes_scaled_(.+).csv","\\1",list.files(baseResultsDir(exprDir), pattern="cluster_prototypes_scaled_.*.csv"))) {
	clusterExpr[[patId]] <- lib$dtToDf(fread(baseResultsDir(exprDir,"/cluster_prototypes_scaled_",patId,".csv")))
}
minExprThresh <- as.numeric(quantile(as.matrix(clusterExpr$all),0.9))
clusterIsExpr <- sapply(apply(clusterExpr$all>=minExprThresh,2,which), function(i) rownames(clusterExpr$all)[i])
allExprGenes <- sort(unique(unlist(clusterIsExpr)))
clusterIsExprE <- sapply(apply(clusterExpr$E>=minExprThresh,2,which), function(i) rownames(clusterExpr$E)[i])
allExprGenesE <- sort(unique(unlist(clusterIsExprE)))

# define order in which to list LCH subsets and regulatory modules:
# (N.B. this uses in parts an outdated nomenclature that we tidied up for the manuscript;
# e.g. a5 = accessible in ("old") 5 = accessible in LCH-C10 = accessible in LCH-S1;
# apologies for the mess!)
moduleOrder <- c("a5","a2","a8","i5","i2","i8")
cellOrder <- c("C5", "C2", "C8") #"C1", 
atacOrder <- paste_("LCH", c("C5_3", "C5_4", "C2_3", "C2_4", "C8_3", "C8_4")) 
cOrder <- fread(baseResultsDir(exprDir,"/clusters.csv"))$cluster
	
cOrderSel <- paste0("LCH-C",c(10,7,12,5,0,13))
modDict <- list(
	C5=c("LCH-C10"),
	C2=c("LCH-C12"),
	C8=c("LCH-C5")
)

cOrderSel <- cOrder[cOrder%in%unlist(modDict)]

subsetMarkers <- fread(baseResultsDir(exprDir,"/cluster_markers_traj_cluster_id_sel_strat_w7.csv"))
degs <- list(
	C5=subsetMarkers[cluster%in%modDict$C5,sort(unique(gene))],
	C2=subsetMarkers[cluster%in%modDict$C2,sort(unique(gene))],
	C8=subsetMarkers[cluster%in%modDict$C8,sort(unique(gene))]
)
allDegs <- unique(as.character(unlist(degs)))

moduleClusters <- list(
	a2 = c("C2"), a5 = c("C5"), a8 = c("C8"),
	i2 = c("C5","C8"), i5 = c("C2","C8"), i8 = c("C2","C5")
)

moduleNames <- c("a5"="a1", "i5"="i1", "a2"="a12", "i2"="i12", "a8"="a11", "i8"="i11")