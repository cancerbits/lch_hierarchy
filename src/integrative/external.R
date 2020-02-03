#!/usr/bin/env $CODEBASE/Rscript
#
# Analyze our gene signatures in external datasets.
#
# run("integrative", "external")

setCurrentAnalysis(paste0(cfg$ver,"/integrative"))
loadLibraries(c("Biobase", "GEOquery", "affy", "data.table","limma"))


# download and cache external data:
externalDataList <- fread(metaDir("external_data.csv"))
simpleCache("geo_data", {
	# get raw microarray data from GEO and quantile-normalize them:

	geoData <- sapply(externalDataList[,unique(series_accession)], function(acc) { 
		# doesn't work anymore (GEO changed its URLs):
		# getGEO(acc, GSEMatrix =TRUE, AnnotGPL=TRUE)
		
		# do it by hand instead:
		AnnotGPL <- TRUE
		GSEMatrix <- TRUE
		getGPL <- TRUE
		stub <- gsub("\\d{1,3}$", "nnn", acc, perl=TRUE)
		f <- resultsDir(acc,"_geo.matrix")
		gz <- paste0(f,".gz")
		if(!file.exists(f)) {
			download.file(sprintf("ftp://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/%s_series_matrix.txt.gz", stub, acc, acc), destfile = gz, mode = "wb", method = getOption("download.file.method.GEOquery"))
			gunzip(gz, destname=f, overwrite=T)
		}
		gse <- GEOquery:::parseGSEMatrix(f, destdir=resultsDir(), AnnotGPL=AnnotGPL, getGPL = getGPL)$eset
		if(nrow(exprs(gse))==0) gse <- NULL
		return(gse)
	},  simplify=FALSE)

	missing <- sapply(geoData, is.null)
	msg("missing: ", paste(names(which(missing)), collapse=", "))
	geoData <- geoData[!missing]

	geoDataNorm <- rblapply(geoData, function(x) { 
		qnorm <- normalizeBetweenArrays(exprs(x), "quantile") 
		tmp <- as.data.table(fData(x))[,.(ID, `Gene symbol`)]
		setnames(tmp, c("ID", "geneSymbol"))
		setkey(tmp,ID)
		tmp <- melt(cbind(tmp[rownames(qnorm),], qnorm), id.vars=c("ID","geneSymbol"))
		tmp$sampleName <- pData(x)[as.character(tmp$variable),"title"]

		if(tmp[,max(value)] >= 32) tmp[,value:=log2(value)]

		tmp
	}, "gse")
	geoDataNorm[,dataType:="Microarray log2(qnorm)"]


	# get processed data for missing experiments (that is, RNA-seq data): 

	geoDataNorm <- rbind(geoDataNorm, rblapply(names(which(missing)), function(gse) {
		# the below isn't really generic, it assumes there's one column with a gene name and the rest
		# are sample-wise gene expression values! ... might need adjusting if more samples are added

		f <- externalDataList[series_accession==gse,unique(file_link)]
		download.file(f, lib$plotDir("tmpfile"))
		if(grepl("gz$", f)) {
			f <- paste("zcat",lib$plotDir("tmpfile"))
		} else {
			f <- lib$plotDir("tmpfile")
		}
		
		dt <- melt(fread(f))
		dt <- dt[,.(ID=Gene, geneSymbol=Gene, variable, value, sampleName=variable, dataType=externalDataList[series_accession==gse,unique(data_type)])]
		dt
	}, "gse"))


	# write outputs:

	fwrite(geoDataNorm, file=resultsDir("geo_data.csv"))
	
	geoDataNorm
}, assignToVar="geoData", recreate=F)



# load selected gene signatures:
pT <- cfg$deg_thresh_fdr
fcT <- cfg$deg_thresh_lfc
allGenesets <-  getMarkerLists(pT, fcT)[c("traj_cluster_id_sel_strat_w7","LCH_vs_cTypes_strat_overlap")]
for(curClust in paste0("LCH-C",c(10,0,5,12,13))) allGenesets[[paste_( "traj_cluster_id_sel_strat_w7", curClust)]] <- list(genesets=allGenesets[[ "traj_cluster_id_sel_strat_w7"]]$genesets[sel==T & p_val_adj<=cfg$deg_thresh_fdr & (avg_logFC)>=0.25 & cluster==curClust,])


# generate plots:
n <- "traj_cluster_id_sel_strat_w7"
for(n in names(allGenesets)) {
	clusterMarkers <- allGenesets[[n]]$genesets[p_val_adj<=pT & abs(avg_logFC)>=fcT,]
	markerFor <- clusterMarkers[, .(cluster=paste(cluster,collapse="/"),p_val_adj=min(p_val_adj), avg_logFC=mean(avg_logFC)), by=gene][order(cluster),.(gene,cluster,p_val_adj,avg_logFC)]    	
	setkey(markerFor,gene)
	selDegMarkers <- markerFor[order(cluster,p_val_adj),gene]
	
	curSeries <- "GSE74442"
	for(curSeries in geoData[,unique(gse)]) {

		plotData <- geoData[gse==curSeries & geneSymbol%in%selDegMarkers, ]
		msgF("%s: %s (%d data points)", n, curSeries, nrow(plotData))
			
		if(curSeries=="GSE16395") plotData[, sampleName:=gsub("TonsilPool.? ?CD", "Tonsil CD", (gsub("[-−]", " ", sampleName)))] # do some cosmetic renaming
		if(curSeries=="GSE35457") plotData[, sampleName:=gsub("dc1c","cd1c",gsub("bd-","blood-",sampleName))] # do some cosmetic renaming

		hmData <- (lib$dtToDf(dcast(plotData[,.(m=mean(value)),by=.(geneSymbol,sampleName)], geneSymbol~sampleName, value.var="m")))
		hmData <- hmData[sort(rownames(hmData)),]
		
		hcCol <- myHclust(t(hmData))
		hcRow <- myHclust(hmData)

		doScale <- TRUE

		hmDataX <- hmData[selDegMarkers,]
		rownames(hmDataX) <- selDegMarkers
		
		if(doScale) {
			hmDataX <- t(apply(hmDataX, 1, scale))
			colnames(hmDataX) <- colnames(hmData)
		}
		
		annot <- lib$dtToDf(markerFor[rownames(hmDataX),.(gene,cluster)])
		
		lib$pdfPlot(paste_("hm_degs", n, curSeries, doScale), 16, 50)
		pheatmap(hmDataX, gaps_row=cumsum(as.data.table(annot[rownames(hmDataX),,drop=F])[,.N,by=cluster][,N]), cellwidth=9, cellheight=9, annotation_row=annot[rownames(hmDataX),,drop=F], col=colorPalettes$heatmap_expression, treeheight_row=0, cluster_cols=hcCol, cluster_rows=F, scale="none", border_color="white", main=curSeries)
		dev.off()
		
		lib$pdfPlot(paste_("hm_degs", n, curSeries, doScale, "small"), 16, 6)
		pheatmap(hmDataX, show_rownames=F, cellwidth=9, gaps_row=cumsum(as.data.table(annot[rownames(hmDataX),,drop=F])[,.N,by=cluster][,N]), annotation_row=annot[rownames(hmDataX),,drop=F], col=colorPalettes$heatmap_expression, treeheight_row=0, cluster_cols=hcCol, cluster_rows=F, scale="none", border_color="white", main=curSeries)
		dev.off()

		lib$pdfPlot(paste_("hm_degs", n, curSeries, doScale, "sort"), 16, 50) 
		hmDataX <- hmDataX[sort(rownames(hmDataX)),]
		pheatmap(hmDataX, gaps_row=cumsum(as.data.table(annot[rownames(hmDataX),,drop=F])[,.N,by=cluster][,N]), cellwidth=9, cellheight=9, annotation_row=annot[rownames(hmDataX),,drop=F], col=colorPalettes$heatmap_expression, treeheight_row=0, cluster_cols=hcCol, cluster_rows=F, scale="none", border_color="white", main=curSeries)
		dev.off()
		
		pheatmap(hmDataX, file=resultsDir(paste_("hm_degs", n, curSeries, doScale, "clust.pdf")), cluster_rows=hcRow, cellwidth=9, cellheight=9, annotation_row=annot[rownames(hmDataX),,drop=F], col=colorPalettes$heatmap_expression, treeheight_row=0, cluster_cols=hcCol,scale="none", border_color="white", main=curSeries)
	}
	
}
