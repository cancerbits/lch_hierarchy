# load libraries:
setwd(paste0(Sys.getenv("CODEBASE"),"/lch"))
lib$projectName <- "lch"
loadLibraries(c("data.table", "ggplot2", "simpleCache", "digest", "tidyr", "ggrepel","plyr", "dplyr", "pryr", "pheatmap", "BiocManager", "devtools", "viridis"))
tryCatch({library("Cairo")}, error=function(e) print(e))


############################
### FUNCTION DEFINITIONS ###
############################

loadAnnot <- function(type=gsub("^(.*/)(.+)$","\\2",lib$getCurrentAnalysis())) {
	f <- list.files(metaDir(), pattern=paste_("samples",type), full.names=T)
	f <- f[sapply(f, file.size)>0] # drop empty files
	dA <- rblapply(f, fread, "annot_file")
	dA <- merge(dA, fread(metaDir("pat_info.csv")), by="pat_id", all.x=T, all.y=F, cartesian=T)
	dA[,disease_extent_detail:=disease_extent]
	dA[,disease_extent:=gsub("^.*(SS|MS).*$", "\\1", disease_extent)]
	dA[!sex%in%c("M","F"),sex:=NA]
	for(curCol in c("treated","genetics_v600e")) dA[get(curCol)=="" | get(curCol)=="N.d." | get(curCol)=="n/a", paste(curCol):=NA]
	setkey(dA, "sample_name")
	dA
}
alpha <- function(colour, alpha) {
	x <- scales::alpha(colour, alpha)
	names(x) <- names(colour)
	x
}
percent <- function(dt, byTerm, byGroup) {
	y <- as.data.table(dt)[,.N,by=c(byGroup)]
	y <- lib$dt2namedVec(y, byGroup)
	x <- as.data.table(dt)[,.(.N, perc=.N/y[get(byGroup)]),by=c(byTerm, byGroup)]
	x
}
myHclust <- function(x, dist.meth="euclidean", clust.meth="complete") {
	dists <- dist(x, method=dist.meth)
	dists[is.na(dists)] <- max(dists,na.rm=T)
	hclust(dists, method=clust.meth)
}
fwrite <- function(x, ...) data.table::fwrite(tryCatch(as.data.table(x, keep.rownames=T), error=function(e) as.data.table(as.matrix(x), keep.rownames=T)), ...)
resDir <- lib$resourceDir
toolsDir <- function(...) {
	paste0(Sys.getenv("TOOLS"), "/",...)
}
writeCombos <- function(combos, f=resultsDir("combos.csv")) {
	dt <- data.table(sample_name=dA$sample_name, sapply(combos, function(x) dA$sample_name%in%x))
	fwrite(dt, file=f)
}
readCombos <- function(f=resultsDir("combos.csv")) {
	dt <- fread(f)	
	apply(dt[,-1], 2, function(i) dt[i,sample_name])
}
sanitizeEnrichrLabels <- function(txt) {
	sapply(txt, function(x) {
		x <-  gsub(" (Mus musculus|Mouse)", " (mouse)", gsub(" (Homo sapiens|Human)", " (human)", gsub("_", " ", x), ignore.case=TRUE), ignore.case=TRUE)
		x <- gsub(" \\(NAFLD\\)", "", x)
		x <- gsub("^MP\\d+\\s+", "", x, perl=TRUE)
		x <- gsub("\\s+(hsa|WP)\\d+$", "", x, perl=TRUE)
		x <- gsub("\\s+\\(GO:.+\\)$", "", x, perl=TRUE)
		x <- gsub(" \\d+ ChIP-.+ ", " ", x, perl=TRUE)
		x <- gsub("positive", "pos.", x, perl=TRUE)
		x <- gsub("negative", "neg.", x, perl=TRUE)
		x <- gsub("regulation of", "reg. of", x, perl=TRUE)
		x <- gsub("involved in ", "in ", x, perl=TRUE)
		x <- gsub("ligase activity", "ligase", x, perl=TRUE)		
		x <- gsub("(GSE\\d+) sample \\d+", "\\1", x, perl=TRUE)
		x <- gsub("UMLS ", "", x)
		x <- gsub("\\s+", " ", x, perl=TRUE)
		x <- gsub("in DMSO-Rat-Primary rat ", "DMSO-Rat ", x)
		x <- gsub("\\_(\\w|\\d){8}-(\\w|\\d){4}-(\\w|\\d){4}-(\\w|\\d){4}-(\\w|\\d){12}", "",x)
		x <- gsub(" (\\w+)*-\\w+-\\w+$","",x)
		x <- gsub(" .human.$","",x)
		x <- gsub(" GSE.+$","",x)
		x <- gsub(" (activity|kinase ARCHS4 coexpression)", "", gsub("(modification) process", "\\1", gsub("(tissue)? development", "dev.", gsub("(catabo|metabo)lic process","\\1ism",gsub("biosynthetic process","biosynthesis", gsub("([Ss]ignaling|regulatory) [Pp]athway", "pathway", gsub("\\s+$","",gsub(".human. (\\w+-?)*$","",x))))))))
		x <- gsub("GTEX-[^\\s]+ (.+ (fe)?male).+$","GTEX \\1",x)
		lib$capFirst(substr(x, 1, 64))
	})
}
runEnrichr <- function(geneLists, dbs=c("GO_Biological_Process_2018")) {
	enrichrRes <- rblapply(names(geneLists), function(curGrp) { #
		msg(curGrp)
		curGenes <- geneLists[[curGrp]]			

		res <- enrichr(curGenes, databases=dbs)
		rbindlist(res, "database", use.names=T, fill=T)
	}, "grp")

	enrichrRes[, n.hits:=as.numeric(gsub("^(\\d+)/(\\d+)$","\\1",Overlap))]
	enrichrRes[, n.total:=as.numeric(gsub("^(\\d+)/(\\d+)$","\\2",Overlap))]
	
	enrichrRes[, original.term:=Term]
	enrichrRes[, Term:=NULL]
	enrichrRes[, term:=sanitizeEnrichrLabels(original.term)]
	enrichrRes[, uniq.term:=paste_(database,original.term)]

	enrichrRes
}



#################
### CONSTANTS ###
#################

options(bitmapType='cairo') # to enable PNG plots without X11

ver <- ""

colorPalettes <- sapply(split(fread(metaDir("colors.csv")),by="category"), lib$dt2namedVec, nameCol="term", valCol="value", simplify=F)
#colorPalettes$tsne <- rev(brewer.pal(8,"Spectral"))
colorPalettes$tsne <- c("#DDDDDD",brewer.pal(5,"YlGnBu")[2:5])	

colorPalettes$heatmap_expression <- rev(brewer.pal(9,"PiYG"))
colorPalettes$heatmap_chromatin <- rev(brewer.pal(9,"PiYG"))
colorPalettes$heatmap_enrichment <- inferno(9)

colorPalettes$cell_type$DC <- "#ff7f00"
colorPalettes$cell_type$MacMono <- "#984ea3"
colorPalettes$cell_type$TNKcell <- "#377eb8"
colorPalettes$cell_type$Bcell <- "#4daf4a"
colorPalettes$cell_type <- unlist(colorPalettes$cell_type)

# add colors for all annotation categories which have not been explicitly defined yet:
dAtmp <- loadAnnot()
for(annotType in setdiff(colnames(dAtmp),names(colorPalettes))) {
	lvls <- dAtmp[,sort(unique(get(annotType)))]
	lvls[is.na(lvls) | lvls==""] <- "NA"
	colorPalettes[[annotType]] <- lib$getCategoryColors(lvls)
}
rm(dAtmp,annotType)

markerGenes <- fread(resDir("marker_genes.csv"))

enrichrDBs <- as.character(read.csv(metaDir("enrichr_dbs.csv"))[,1]) 
enrichrDBs <- enrichrDBs[!grepl("#",enrichrDBs)]

cfg <- list(
	ver="tsne5",
	#sel_clustering=list(name="res.0.8", res=0.8, pruning=1/15, reduction="tsne", perplexity=30, tsne_n_pcs=5, nPcMax=6),
	sel_clustering=list(name="res.0.2", res=0.2, pruning=1/15, reduction="tsne", perplexity=30, tsne_n_pcs=5, nPcMax=6),
	thresh_max_mito = 0.1,
	thresh_min_cells = 3,
	thresh_min_genes = 1000,
	thresh_max_genes = Inf,
	thresh_min_ncells_total = 1000,
	thresh_min_ngenes_total = 10000,
	regress_vars=c("nUMI", "perc_mito", "G1", "G2M", "S", "pat_id", "biopsy", "sex", "age"),
	deg_min_pct = 0.25, 
	deg_thresh_use = 0.25,
	deg_thresh_lfc = 0,
	deg_thresh_fdr = 0.005,
	lch_marker_thresh = 1,
	lch_clust_thresh = 1/4,
	lch_nonclust_thresh = 0.6,
	lch_neg_list = c("CD19","CD3D","CD27","IL32","CD7","NKG7","CD163"), # neg for CD163: https://www.sciencedirect.com/science/article/pii/B9780323497145000107
	lch_pos_list = c("CD1A","CD207")
)

options("RCACHE.DIR"=baseResultsDir(cfg$ver,"/RCache/"))


minPlotCatVals <- 2
maxPlotCatVals <- 32

patientGroupings <- c("pat_id", "biopsy", "age", "disease_extent")

genomeVer <- "hg38" #"GRCh38"

# a lookup table to translate between names in an outdated version of the LCH cell
# subset nomenclature (still used in parts of the code) and the one used in the paper:
finalDict <- c(
	C5="LCH-S1",
	C2="LCH-S12",
	C8="LCH-S11",
	"LCH-C0"="LCH-S13",
	"LCH-C1"="LCH-S7",
	"LCH-C2"="LCH-S4",
	"LCH-C3"="LCH-S3",
	"LCH-C4"="LCH-S9",
	"LCH-C5"="LCH-S11",
	"LCH-C6"="LCH-S10",
	"LCH-C7"="LCH-S2",
	"LCH-C8"="LCH-S5",	
	"LCH-C9"="LCH-S8",
	"LCH-C10"="LCH-S1",
	"LCH-C11"="LCH-S6",
	"LCH-C12"="LCH-S12",	
	"LCH-C13"="LCH-S14",
	"a2"="a12",
	"a5"="a1",
	"a8"="a11",
	"i2"="i12",
	"i5"="i1",
	"i8"="i11"
)
