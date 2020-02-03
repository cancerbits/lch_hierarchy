#!/usr/bin/env $CODEBASE/Rscript
#
# Define constants and methods used in the analysis of the 10x scRNA-seq data.
#
# run("single_cell", "init")

###################################
### SCRNA-SEQ-SPECIFIC LIBARIES ###
###################################

if(lib$projectName!="lch") run("common")
setCurrentAnalysis("")
loadLibraries(c("fpc","doMC", "Rtsne", "tsne"))
loadLibraries(c("monocle", "DelayedArray", "DelayedMatrixStats", "org.Hs.eg.db", "org.Mm.eg.db")) 
if(!"Seurat"%in%installed.packages()) devtools::install_version(package = 'Seurat', version = package_version('2.3.4'))
loadLibrary("Seurat")




############################
### FUNCTION DEFINITIONS ###
############################

scaleMinMax <- function(x) {
	x <- x-min(x)
	x/max(x)
}
getMarkerLists <- function(pThresh=cfg$deg_thresh_fdr, fcThresh=0) {
	allGenesets <- list()

	simpleCache(paste_("clusterMarkers", "w", "traj_cluster_id_sel_strat_w7", allComboStr), assignToVar="clusterMarkers")
	clusterMarkers[, sel:=p_val_adj<=pThresh & avg_logFC>fcThresh]
	allGenesets[["traj_cluster_id_sel_strat_w7"]] <- list(data=NULL, genesets=clusterMarkers[sel==T,])

	simpleCache(paste_("scData", "cca", nCCS, allComboStr,"ext"), assignToVar="scData")

	simpleCache(paste_("clusterMarkers", "w", "LCH_vs_cTypes_strat", allComboStr), assignToVar="clusterMarkers")
	clusterMarkers[, sel:=p_val_adj<=pThresh & avg_logFC>fcThresh]
	allGenesets[["LCH_vs_cTypes_strat"]] <- list(data=scData, groupBy="cell_type", bg="LCH", genesets=clusterMarkers[p_val_adj<=pThresh,])
	allGenesets[["LCH_vs_cTypes_strat_overlap"]] <- list(data=scData, groupBy="cell_type", genesets=clusterMarkers[p_val_adj<=pThresh & avg_logFC>fcThresh, .(.N, avg_logFC=mean(avg_logFC), p_val_adj=max(p_val_adj)), by=gene][N>=3,.(gene,avg_logFC,p_val_adj,cluster="LCH")])
	allGenesets[["LCH_vs_cTypes_strat_rev"]] <- list(data=NULL, groupBy=NULL, genesets=clusterMarkers[p_val_adj<=pT & avg_logFC<fcThresh,])

	allGenesets
}
doStratifiedDifferentialAnalysis <- function(scDataMerged, groupBy="cluster_id", stratBy="sex", sampleSize=NULL, nSamples=25, sampleDeThresh=0.9) {
	scDataMerged <- SetAllIdent(scDataMerged, groupBy)
		
	stratBy <- "sex"
	strata <- unique(scDataMerged@meta.data[,stratBy])
	
	clusterMarkersTmp <- rblapply(strata, function(stratum) {
		msgF("\t* %s", stratum)
		scDataCur <- SubsetData(scDataMerged, cells.use=scDataMerged@meta.data[,stratBy]==stratum)
		if(!is.null(sampleSize)) {
			clusterMarkersSamp <- rblapply(1:nSamples, function(n) {
				msg("\t\t",n)
				clusterMarkers <- FindAllMarkers(object = scDataCur, return.thresh=1, random.seed=n, max.cells.per.ident=sampleSize, only.pos = T)	
				clusterMarkers <- as.data.table(clusterMarkers)
			}, "samp")
			clusterMarkers <- clusterMarkersSamp[,.(.N, p_val=max(p_val), p_val_adj=max(p_val_adj), avg_logFC=mean(avg_logFC)),by=.(gene,cluster)]
			clusterMarkers <- clusterMarkers[N>=nSamples*sampleDeThresh,]		
		clusterMarkers			
		} else {
			clusterMarkers <- FindAllMarkers(object = scDataCur, return.thresh=1, only.pos = T)				
		}
		clusterMarkers <- as.data.table(clusterMarkers)
		clusterMarkers
	},"stratum")
	
	nPerClust <- lib$dt2namedVec(as.data.table(scDataMerged@meta.data)[,melt(table(get(stratBy))),by=groupBy][value>=20, length(unique(Var1)), by=.(ID=get(groupBy))])
	
	clusterMarkers <- clusterMarkersTmp[,.(.N, p_val=max(p_val), p_val_adj=max(p_val_adj), avg_logFC=mean(avg_logFC)),by=.(gene,cluster)]
	clusterMarkers <- clusterMarkers[N>=nPerClust[cluster],]
		
	clusterMarkers[, rnkFC:=rank(-avg_logFC)]
	clusterMarkers[, rnkP:=rank(p_val)]
	clusterMarkers[, rnkCombined:=rank(pmax(rnkFC,rnkP),ties.method="random")]	
	clusterMarkers[, sel:=avg_logFC>=log2(1.5) & p_val_adj<=cfg$deg_thresh_fdr
	
	clusterMarkers
}
plotTsneMetadata <- function(scData,  sub=NULL, suffix="", cols=colnames(scData@meta.data), forceColors=NULL) {
	msg("t-SNE plots for metadata:")
	#p <- TSNEPlot(object = scData)
	#gg(p, "tsne", 8, 6, type="pdf", sub=sub)
	plotData <- data.table(scData@dr$tsne@cell.embeddings, scData@meta.data)
	for(n in intersect(colnames(scData@meta.data),cols)) {
		nVals <- length(unique(scData@meta.data[,n]))
		isNum <- is.numeric(scData@meta.data[,n]) 
		msgF("\t* %s (%s)", n, ifelse(isNum, "numeric", sprintf("categorical, %d values", nVals)))
		if(isNum | (nVals<=maxPlotCatVals & nVals>=minPlotCatVals)) {
			p <- ggplot(plotData[sample(nrow(plotData)),], aes(x=tSNE_1, y=tSNE_2, color=get(n))) + geom_point(alpha=0.7, size=1) + defTheme(topLegend=T, noLegendTitle=T) + xlab("t-SNE 1") + ylab("t-SNE 2") + ggtitle(n) #TSNEPlot(object = scData, group.by=n) 
			if(!isNum) {
				p <- p + geom_label(aes(label=lbl), data=as.data.table(plotData)[,.(lbl=unique(get(n)),tSNE_1=median(tSNE_1), tSNE_2=median(tSNE_2)), by=get(n)], color="black")
			}
			if(!is.null(forceColors)) {
				if(isNum) {
					p <- p + scale_color_gradientn(colours=alpha(forceColors, 0.6)) 
				} else {
					p <- p + scale_color_manual(values=alpha(forceColors, 0.6)) 
				}
			} else {
				if(n%in%names(colorPalettes)) p <- p + scale_color_manual(values=alpha(colorPalettes[[n]], 0.6))
				if(isNum) p <- p + scale_color_gradientn(colours=colorPalettes$tsne)
			}
			gg(p, paste0(paste_("tsne",lib$escapeStringInFileName(n)), suffix), 8, 6, type="pdf", sub=sub) 
		}
	}
}
makeTsneMetadataPanel <- function(scData, n, cols=c("pat_id"), forceColors=NULL, ncols=3) {	
	cols <- intersect(colnames(scData@meta.data),cols)
	plotData <- melt(data.table(scData@dr$tsne@cell.embeddings, scData@meta.data), id.vars=colnames(scData@dr$tsne@cell.embeddings), measure.vars=cols)
	plotData[, variable:=factor(variable,levels=cols)]
	isNum <- is.numeric(plotData$value) 
	p <- ggplot(plotData[sample(nrow(plotData)),], aes(x=tSNE_1, y=tSNE_2, color=value)) + geom_point(alpha=0.7, size=1) + defTheme(topLegend=T, noLegendTitle=T) + xlab("t-SNE 1") + ylab("t-SNE 2") + ggtitle(n) + facet_wrap(~variable, ncol=ncols)
	if(!is.null(forceColors)) {
			if(isNum) {
				p <- p + scale_color_gradientn(colours=alpha(forceColors, 0.6)) 
			} else {
				p <- p + scale_color_manual(values=alpha(forceColors, 0.6)) 
			}
	} else {
			if(n%in%names(colorPalettes)) p <- p + scale_color_manual(values=alpha(colorPalettes[[n]], 0.6))
			if(isNum) p <- p + scale_color_gradientn(colours=colorPalettes$tsne)
	}
	
	p
}
plotTsneMetadataPanel <- function(scData, n, suffix="", sub=NULL, cols=c("pat_id"), forceColors=NULL, ncols=3) {	
	p <- makeTsneMetadataPanel(scData, n, cols=cols, forceColors=forceColors, ncols=ncols)
	nrows <- ceiling(length(cols) / ncols)		
	gg(p, paste0(paste_("tsne",lib$escapeStringInFileName(n)),suffix), min(ncols,length(cols))*3.5, nrows*3, type="pdf", addBase=T, sub=sub) 
}
plotTsneGenePanel <- function(scData, n, genes, plotsPerPage=9) {
	markerGeneSymbols <- unique(intersect(genes,rownames(scData@data)))
	suppressWarnings(paginatedMarkers <- split(markerGeneSymbols, ceiling(1:length(markerGeneSymbols) / plotsPerPage)))
	msgF("\t* %s (n = %d, %d pages)", n, length(markerGeneSymbols), length(paginatedMarkers))
	forEach(paginatedMarkers, function(markerGenePage) {
		if(length(markerGenePage)>0) FeaturePlot(object = scData, no.legend=F, pt.size=4.5, vector.friendly=T, features.plot = markerGenePage, cols.use = colorPalettes$tsne, reduction.use = "tsne")
	})
}
plotTsneGenePanels <- function(scData, sub=NULL, suffix="", w=12, h=9, selMarkersX=NULL, cols=colnames(scData@meta.data), plotsPerPage=9, steps=c("cutom","generic")) {
	msg("t-SNE gene panels:")
	
	allMarkerLists <- list()
	if("generic" %in% steps) allMarkerLists[["generic"]] = markerGenes[order(cell_type), symbol]
	if(!is.null(selMarkersX) & "custom"%in%steps) {
		allMarkerLists <- c(list(custom=selMarkersX), allMarkerLists)
	}
	
	lib$pdfPlot(paste0("tsne_panel_genes",suffix), w, h, sub=sub)
	forEach(names(allMarkerLists), function(n) {
		plotTsneGenePanel(scData, n, allMarkerLists[[n]], plotsPerPage=plotsPerPage)
	})
	dev.off()
}
callLCH <- function(scData) {
	msg("Call LCH cells...")
	d <- as.matrix(scData@data)[,rownames(scData@meta.data)	]
	
	thresh <- cfg$lch_marker_thresh
	n <- c(cfg$lch_neg_list,cfg$lch_pos_list)
	if(grepl("q\\d+",thresh)) {
		thresh <- gsub("q","", thresh)
		thresh <- (as.numeric(thresh)/100)
		thresh <- sapply(n, function(x) as.numeric(quantile(d[x,], thresh)))
	} else {
		thresh <- structure(rep(thresh, length(n)), names=n)
	}
	
	scData@meta.data$neg_lch <- rowSums(sapply(cfg$lch_neg_list, function(x) d[x,]>thresh[x]))>0
	scData@meta.data$pos_lch <- rowSums(sapply(cfg$lch_pos_list, function(x) d[x,]>thresh[x]))>=length(cfg$lch_pos_list)
		
	scData@meta.data$lch_pos_perc <- lib$dt2namedVec(as.data.table(scData@meta.data)[,.(
		perc=round(sum(pos_lch)/length(pos_lch)*100,1)
	),by=cluster_id],"cluster_id")[scData@meta.data$cluster_id]
	scData@meta.data$lch_neg_perc <- lib$dt2namedVec(as.data.table(scData@meta.data)[,.(
		perc=round(sum(neg_lch)/length(neg_lch)*100,1)
	),by=cluster_id],"cluster_id")[scData@meta.data$cluster_id]
	
	scData@meta.data$is_lch_cell <- scData@meta.data$lch_pos_perc/100>=cfg$lch_clust_thresh & scData@meta.data$lch_neg_perc/100<cfg$lch_nonclust_thresh
	msgF("\t* found %d LCH cells (%.1f%%)", sum(scData@meta.data$is_lch_cell), sum(scData@meta.data$is_lch_cell)/length(scData@meta.data$is_lch_cell)*100)
	
	scData
}
callCellType <- function(scData, lbl, posList, negList) {
	msgF("Call %s cells...", lbl)
	d <- as.matrix(scData@data)[,rownames(scData@meta.data)	]
	
	thresh <- cfg$lch_marker_thresh #as.numeric(quantile(d, cfg$lch_marker_quant))
	n <- c(posList, negList)
	if(grepl("q\\d+",thresh)) {
		thresh <- gsub("q","", thresh)
		thresh <- (as.numeric(thresh)/100)
		thresh <- sapply(n, function(x) as.numeric(quantile(d[x,], thresh)))
	} else {
		thresh <- structure(rep(thresh, length(n)), names=n)
	}
	
	scData@meta.data[,paste_("neg",lbl)] <- rowSums(sapply(negList, function(x) d[x,]>thresh[x]))>0
	scData@meta.data[,paste_("pos",lbl)] <- rowSums(sapply(posList, function(x) d[x,]>thresh[x]))>=length(posList)
		
	scData@meta.data[,paste_("pos","perc", lbl)] <- lib$dt2namedVec(as.data.table(scData@meta.data)[,.(
		perc=round(sum(get(paste_("pos",lbl)))/length(get(paste_("pos",lbl)))*100,1)
	),by=cluster_id],"cluster_id")[scData@meta.data$cluster_id]
	scData@meta.data[,paste_("neg","perc", lbl)] <- lib$dt2namedVec(as.data.table(scData@meta.data)[,.(
		perc=round(sum(get(paste_("neg",lbl)))/length(get(paste_("neg",lbl)))*100,1)
	),by=cluster_id],"cluster_id")[scData@meta.data$cluster_id]
	
	scData@meta.data[,paste_("is",lbl)] <- scData@meta.data[,paste_("pos","perc", lbl)]/100>=cfg$lch_clust_thresh & scData@meta.data[,paste_("neg","perc", lbl)]/100<cfg$lch_nonclust_thresh
	msgF("\t* found %d %s cells (%.1f%%)", sum(scData@meta.data[,paste_("is",lbl)]), lbl, sum(scData@meta.data[,paste_("is",lbl)])/length(scData@meta.data[,paste_("is",lbl)])*100)
	
	scData
}
process10xData <- function(scData, s, scoreCellCyle=FALSE, doClustering=TRUE, doDimRed=TRUE, doScaling=TRUE, callLCH=TRUE, doVarGenes=doScaling, selClustering=cfg$sel_clustering, regressVars=cfg$regress_vars, aovVars=intersect(colnames(dA),regressVars)) {
	if(doScaling) {
		msg("Normalize data...")	
		scData <- NormalizeData(object = scData, normalization.method = "LogNormalize", scale.factor = 1000000)
	}
	
	# score cell cycle stage:
	if(scoreCellCyle) { # if FALSE, assume this has been done previously
		loadLibrary("scran") 
		
		msg("Calculate cell cycle scores...")	
		ccMarkers <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
		anno <- AnnotationDbi::select(org.Hs.eg.db, keys=rownames(scData@data), keytype="SYMBOL", column="ENSEMBL")
		ensembl <- anno$ENSEMBL[match(rownames(scData@data), anno$SYMBOL)]
		assignments <- cyclone(as.matrix(scData@data), ccMarkers, gene.names=ensembl, BPPARAM=MulticoreParam(workers=8), verbose=TRUE) #
		for(n in colnames(assignments$score)) {
			x <- as.numeric(assignments$score[,n])
			names(x) <- colnames(scData@data)
			scData <- AddMetaData(scData, x, n)
		}
	}
	
	# find variable genes, clusters, and perform dimensionality reduction:
	if(doScaling) {
		msg("Scale data and regress confounders...")		
		rvars <- intersect(regressVars, names(which(apply(scData@meta.data,2,function(x) length(unique(x)))>1)))
		#if(selClustering$tsne_n_pcs=="aov") rvars <- setdiff(rvars, aovVars)
		scData <- ScaleData(object=scData, vars.to.regress=rvars)
	}
	if(doVarGenes) {
		msg("Find variable genes...")	
		scData <- FindVariableGenes(object = scData, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
	}
	
	whichPcs <- 1:5
	
	if(doDimRed) {
		msg("Perform PCA...")	
		lib$pdfPlot("snn_graph", 5, 5, sub=s)
		scData <- RunPCA(object = scData, plot.SNN=T, save.SNN=T, pc.genes = scData@var.genes,do.print = F, pcs.print = 1:5, genes.print = 5)
		dev.off()
		scData <- ProjectDim(object = scData)
		
		nPcs <- selClustering$tsne_n_pcs
		nPcMax <- min(selClustering$nPcMax,ncol(scData@dr$pca@gene.loadings))
		
		if(nPcs=="aov") {
			padj <- linkPCAtoAnnotation(scData@dr$pca, scData@meta.data[rownames(scData@dr$pca@cell.embeddings), aovVars], do.plot=F)
			whichPcs <- as.numeric(which(rowSums(padj <= cfg$deg_thresh_fdr) == 0))
			whichPcs <- whichPcs[whichPcs<=nPcMax]
		} else if(nPcs=="jackstraw") {
			msg("Run JackStraw...")	
			scData <- JackStraw(object = scData, num.pc=nPcMax, num.replicate = 100)
			lib$pdfPlot("jackstraw", 9, 9, sub=s)
			scData <- JackStrawPlot(scData, PCs=1:nPcMax)
			dev.off()
			whichPcs <- 1:max(2,sum(scData@dr$pca@jackstraw@overall.p.values[,2]<=cfg$deg_thresh_fdr))
			print(whichPcs)
		} else {
			whichPcs <- 1:nPcs
		}
	}
	
	if(length(whichPcs)<2) {
		msg("<2 PCs selected, add")
		whichPcs <- c(whichPcs, setdiff(1:3, whichPcs)[1:(2-length(whichPcs))])
		print(length(whichPcs))
	}
	
	if(doClustering) {
		msg("Run t-SNE...")	
		#nc <- max(2,min(nc,nPcs))
		scData <- RunTSNE(object = scData, dims.use = whichPcs, do.fast = TRUE, perplexity=selClustering$perplexity)
		msg("Find clusters...")	
		if(cfg$sel_clustering$reduction=="tsne") whichPcs <- 1:2
		scData <- FindClusters(object = scData, force.recalc=T, reduction.type=selClustering$reduction, dims.use = whichPcs, prune.SNN=selClustering$pruning, resolution=selClustering$res, print.output = 0, save.SNN = TRUE)	
		scData@meta.data$cluster_id <- paste0("c",scData@meta.data[, selClustering$name])	
	}

	# figure out which cells are LCH cells (by looking for clusters with a certain percentage of CD1A/CD207-double-positive cells):
	if(callLCH) {
		scData <- callLCH(scData)
	}
	
	scData
}

#### GSEA plots
# https://github.com/PeeperLab/Rtoolbox/blob/master/R/ReplotGSEA.R
## Script by Thomas Kuilman
## path argument: path to output folder of analysis (e.g. PATH/my_analysis.GseaPreranked.1470948568349)
## gene.set argument: name of the gene set (e.g. V$AP1_Q2).
## It is used in a grep command, so multiple matching is possible.
## Also, R regular expressions can be handled, e.g. "IL2[0-9]$"
## Leading "V$" from gene set names are stripped to allow using the grep command.
## In case of multiple grep matches a warning is given and the first option is plotted.
## class.name: the name of the class / variable to which genes have been correlated (e.g. drug-treatment)
##
## * adapted to allow for non [0;1] ranking metric, custom colors, custom axis title
replotGSEA <- function(path, gene.set, class.name, negCol="blue", posCol="red", rnkMetricTitle="Ranking metric") {
  
  if(missing(path)) {
    stop("Path argument is required")
  }
  if (!file.exists(path)) {
    stop("The path folder could not be found. Please change the path")
  }
  if(missing(gene.set)) {
    stop("Gene set argument is required")
  }

  ## Load .rnk data
  path.rnk <- list.files(path = file.path(path, "edb"),
                         pattern = ".rnk$", full.names = TRUE)
  gsea.rnk <- read.delim(file = path.rnk, header = FALSE)
  colnames(gsea.rnk) <- c("hgnc.symbol", "metric")
  
  ## Load .edb data
  path.edb <- list.files(path = file.path(path, "edb"),
                         pattern = ".edb$", full.names = TRUE)
  gsea.edb <- read.delim(file = path.edb,
                         header = FALSE, stringsAsFactors = FALSE)
  gsea.edb <- unlist(gsea.edb)
  gsea.metric <- gsea.edb[grep("METRIC=", gsea.edb)]
  gsea.metric <- unlist(strsplit(gsea.metric, " "))
  gsea.metric <- gsea.metric[grep("METRIC=", gsea.metric)]
  gsea.metric <- gsub("METRIC=", "", gsea.metric)
  gsea.edb <- gsea.edb[grep("<DTG", gsea.edb)]
  
  # Select the right gene set
  if (length(gsea.edb) == 0) {
    stop(paste("The gene set name was not found, please provide",
               "a correct name"))
  }
  if (length(grep(paste0(gsub(".\\$(.*$)", "\\1", gene.set), " "), gsea.edb)) > 1) {
    warning(paste("More than 1 gene set matched the gene.set",
                  "argument; the first match is plotted"))
  }
  gsea.edb <- gsea.edb[grep(paste0(gsub(".\\$(.*$)", "\\1", gene.set), " "), gsea.edb)[1]]
  
  # Get template name
  gsea.edb <- gsub(".*TEMPLATE=(.*)", "\\1", gsea.edb)
  gsea.edb <- unlist(strsplit(gsea.edb, " "))
  gsea.template <- gsea.edb[1]
  
  # Get gene set name
  gsea.gene.set <- gsea.edb[2]
  gsea.gene.set <- gsub("GENESET=gene_sets.gmt#", "", gsea.gene.set)
  
  # Get enrichment score
  gsea.enrichment.score <- gsea.edb[3]
  gsea.enrichment.score <- gsub("ES=", "", gsea.enrichment.score)
  
  # Get gene set name
  gsea.normalized.enrichment.score <- gsea.edb[4]
  gsea.normalized.enrichment.score <- gsub("NES=", "",
                                           gsea.normalized.enrichment.score)
  
  # Get nominal p-value
  gsea.p.value <- gsea.edb[5]
  gsea.p.value <- gsub("NP=", "", gsea.p.value)
  gsea.p.value <- as.numeric(gsea.p.value)
  
  # Get FDR
  gsea.fdr <- gsea.edb[6]
  gsea.fdr <- gsub("FDR=", "", gsea.fdr)
  gsea.fdr <- as.numeric(gsea.fdr)
  
  # Get hit indices
  gsea.edb <- gsea.edb[grep("HIT_INDICES=", gsea.edb):length(gsea.edb)]
  gsea.hit.indices <- gsea.edb[seq_len(grep("ES_PROFILE=", gsea.edb) - 1)]
  gsea.hit.indices <- gsub("HIT_INDICES=", "", gsea.hit.indices)
  gsea.hit.indices <- as.integer(gsea.hit.indices)
  
  # Get ES profile
  gsea.edb <- gsea.edb[grep("ES_PROFILE=", gsea.edb):length(gsea.edb)]
  gsea.es.profile <- gsea.edb[seq_len(grep("RANK_AT_ES=", gsea.edb) - 1)]
  gsea.es.profile <- gsub("ES_PROFILE=", "", gsea.es.profile)
  gsea.es.profile <- as.numeric(gsea.es.profile)
  
  
  ## Create GSEA plot
  # Save default for resetting
  def.par <- par(no.readonly = TRUE)
  
  # Create a new device of appropriate size
  #dev.new(width = 3, height = 3)
  
  # Create a division of the device
  gsea.layout <- layout(matrix(c(1, 2, 3, 4)), heights = c(1.7, 0.5, 0.2, 2))
  layout.show(gsea.layout)
  
  # Create plots
  par(mar = c(0, 5, 2, 2))
  plot(c(1, gsea.hit.indices, length(gsea.rnk$metric)),
       c(0, gsea.es.profile, 0), type = "l", col = "red", lwd = 1.5, xaxt = "n",
       xaxs = "i", xlab = "", ylab = "Enrichment score (ES)",
       main = list(gsea.gene.set, font = 1, cex = 1),
       panel.first = {
          abline(h = seq(round(min(gsea.es.profile), digits = 1),
                 max(gsea.es.profile), 0.1),
                 col = "gray95", lty = 2)
         abline(h = 0, col = "gray50", lty = 2)
       })
  plot.coordinates <- par("usr")
  if(gsea.enrichment.score < 0) {
    text(length(gsea.rnk$metric) * 0.01, plot.coordinates[3] * 0.98,
         paste("Nominal p-value:", gsea.p.value, "\nFDR:", gsea.fdr, "\nES:",
               gsea.enrichment.score, "\nNormalized ES:",
               gsea.normalized.enrichment.score), adj = c(0, 0))
  } else {
    text(length(gsea.rnk$metric) * 0.99, plot.coordinates[4] - ((plot.coordinates[4] - plot.coordinates[3]) * 0.03),
         paste("Nominal p-value:", gsea.p.value, "\nFDR:", gsea.fdr, "\nES:",
               gsea.enrichment.score, "\nNormalized ES:",
               gsea.normalized.enrichment.score, "\n"), adj = c(1, 1))
  }
  
  par(mar = c(0, 5, 0, 2))
  plot(0, type = "n", xaxt = "n", xaxs = "i", xlab = "", yaxt = "n",
       ylab = "", xlim = c(1, length(gsea.rnk$metric)))
  abline(v = gsea.hit.indices, lwd = 0.75)
  
  par(mar = c(0, 5, 0, 2))

	
	abscap <- function(x, cap=0.99) {
		thresh <- quantile(abs(x), cap)
		i <- abs(x)>=thresh
		x[i] <- sign(x[i]) * thresh
		x
	}

  rank.colors <- abscap(gsea.rnk$metric)
  rank.colors <-  rank.colors + max(abs( rank.colors)) # min(gsea.rnk$metric)
  rank.colors <- rank.colors / max(rank.colors)
  rank.colors <- ceiling(rank.colors * 255 + 1)
  rank.colors <- colorRampPalette(c(negCol, "white", posCol))(256)[rank.colors]
  # Use rle to prevent too many objects
  rank.colors <- rle(rank.colors)
  barplot(matrix(rank.colors$lengths), col = rank.colors$values, border = NA, horiz = TRUE, xaxt = "n", xlim = c(1, length(gsea.rnk$metric)))
  box()
  text(length(gsea.rnk$metric) / 2, 0.7,
       labels = ifelse(!missing(class.name), class.name, gsea.template))
  text(length(gsea.rnk$metric) * 0.01, 0.7, "Positive", adj = c(0, NA))
  text(length(gsea.rnk$metric) * 0.99, 0.7, "Negative", adj = c(1, NA))
  
  par(mar = c(5, 5, 0, 2))
  rank.metric <- rle(round(gsea.rnk$metric, digits = 2))
	
  yb <- range(rank.metric$values)

  plot(gsea.rnk$metric, type = "n", xaxs = "i",
	     xlab = "Rank in ordered gene list", xlim = c(0, length(gsea.rnk$metric)),
	     ylim = yb, yaxs = "i",
	     ylab = if(gsea.metric == "None") {rnkMetricTitle} else {gsea.metric},
	     panel.first = abline(h = seq(-0.5, 0.5, 0.5), col = "gray95", lty = 2))

  barplot(rank.metric$values, col = "lightgrey", lwd = 0.1, xaxs = "i",
       xlab = "Rank in ordered gene list", xlim = c(0, length(gsea.rnk$metric)),
       ylim = yb, yaxs = "i", width = rank.metric$lengths, border = NA,
       ylab = ifelse(gsea.metric == "None", rnkMetricTitle, gsea.metric), space = 0, add = TRUE)
  box()
  
  # Reset to default
  par(def.par)

}
getEntropyData <- function(selCombos, recluster=F) {
	dx <- sapply(selCombos, function(combo) {
		comboStr <- paste_(paste_(combo,collapse="."),"lch_only")
		cacheName <- paste_("scData", "SCENTSLICE3", "cca", nCCS, comboStr)
		if(recluster) cacheName <- paste_(cacheName,"recluster")
		
		sc <- NULL
		if(paste0(cacheName,".RData")%in%listCaches()) {
			simpleCache(cacheName, assignToVar="sc", nofail=T)  			
			
			sc <- list(combo=combo, n=comboStr, data=sc)			
		} else {
			msgF("no trajectories for: %s", comboStr)
		}
		
		return(sc)	
	}, simplify=F)
	dx[!sapply(dx,is.null)]
}
makeEntropyPlotPanel <- function(x, grpBy="cluster_id", renameClusts=F, forceMetaColors=NULL, metaCols=c("pat_id"), selMarkers=c("CD1A","CD207"), thm=NULL) {
	pTmp <- makeEntropyPlot(x$data, grpBy=grpBy, renameClusts=renameClusts)
	
	pEnt <- pTmp$p	
	if(!is.null(thm) & thm%in%colnames(dA)) {
		curAnnot <- unique(dA[x$combo,get(thm)])
		pEnt <- pEnt + ggtitle(sprintf("%s = %s\n(%s))", thm, paste0(curAnnot,collapse=", "), paste0(x$combo,collapse=".")))
	} else {
		pEnt <- pEnt + ggtitle(x$n)
	}
	
	o <- pTmp$pdata[,.(e=median(entropy)),by=.(clustId=get(grpBy))][order(-e), clustId]
	pBoxEnt <- ggplot(pTmp$pdata, aes(x=factor(get(grpBy), levels=o), y=entropy)) + geom_boxplot(aes(fill=get(grpBy))) + defTheme(flipX=T) + xlab(NULL) + scale_fill_manual(values=pTmp$cols, guide=F)
			
	# metadata:
	pMeta <- makeTsneMetadataPanel(x$data, "void", cols=metaCols, forceColors=forceMetaColors, ncols=ceiling(sqrt(length(metaCols))))
	
	# genes:
	colorPalettes$tsne <- c("#EEEEEE",brewer.pal(5,"YlGnBu")[2:5])
	plTmp <- FeaturePlot(object = x$data, pt.size=4.5, vector.friendly=T, features.plot = intersect(rownames(x$data@data), selMarkers), no.legend=T, cols.use = colorPalettes$tsne, reduction.use = "tsne", do.return=T)
	
	plTmp <- lapply(plTmp[selMarkers], function(p) {
		if(is.null(p)) p <- ggplot()
		return(p)
	})
	#print(length(plTmp))
	
	plist <- c(list(pEnt=pEnt, pBoxEnt=pBoxEnt, meta=pMeta), plTmp)
	
	res <- list(plotlist=plist, pdata=pTmp$pdata)
	
	return(res)
}
makeEntropyPlot <- function(sc, showLegend=F, mar=0.1, ps=2.25, cols=NULL, grpBy="cluster_id", renameClusts=F) {
	loadLibraries(c("ggrepel","plyr","viridis"))
	
	pData <- cbind(as.data.table(sc@meta.data[,c("entropy",grpBy)],keep.rownames="cell"), sc@dr$tsne@cell.embeddings)
	pData[,c("x","y"):=.(tSNE_1, tSNE_2)]
	
	hullData <- rblapply(pData[,unique(get(grpBy))], function(s) {
		d <- pData[get(grpBy)==s,]
		r1 <- d[, quantile(x, c(mar,1-mar))]
		r2 <- d[, quantile(y, c(mar,1-mar))]
		d <- d[x>=r1[1] & x<=r1[2] & y>=r2[1] & y<=r2[2],]
		if(nrow(d) < 5) d <- NULL
		d
	}, grpBy)		
	hullData <- as.data.table(ddply(hullData, grpBy, lib$findHull, dimName1="x", dimName2="y"))

	centerData <- pData[, .(x=median(x), y=median(y)), by=c(grpBy)]
	
	curColors <- colorPalettes[[grpBy]]
	if(is.null(curColors)) curColors <- lib$getCategoryColors(pData[,unique(get(grpBy))])
	
	if(renameClusts) {
		n <- pData[,.(e=median(entropy)),by=c(grpBy)][order(-e), .(newID=sprintf("c%02d", 1:.N),ID=get(grpBy))]	
		pData[, (grpBy):=as.character(lib$dt2namedVec(n)[get(grpBy)])]
		hullData[, (grpBy):=as.character(lib$dt2namedVec(n)[get(grpBy)])]
		curColors <- structure(rev(viridis(nrow(n))), names=n$newID)
	}
	if(!is.null(cols)) {
		curColors <- cols
	}
	p <- ggplot(data=centerData, aes(x=x, y=y, color=get(grpBy))) + labs(x="t-SNE1", y="t-SNE2") + defTheme(noLegendTitle=T, topLegend=T)
	p <- p + geom_point(size=ps, shape=16, data=pData)
	p <- p + geom_polygon(aes(fill=get(grpBy)), data = hullData, linetype=1, size=0.25, alpha=0.5)
	p <- p + geom_point(size=ps, shape=1, color="black") + geom_text_repel(aes_string(label=grpBy), color="black")
	p <- p + scale_color_manual(values=curColors, guide=showLegend) + scale_fill_manual(values=curColors, guide=showLegend)
	
	return(list(p=p,pdata=pData,cols=curColors))
}




#########################################
### SCRNA-SEQ-SPECIFIC CONFIG OPTIONS ###
#########################################

dA <- loadAnnot("10x")[disease=="LCH" & use==1,]
setkey(dA,sample_name)
sampleToPatient <- lib$dt2namedVec(dA, "sample_name", "pat_id")

selMarkersX <- unique(c(cfg$lch_pos_list, cfg$lch_neg_list, "ZBTB45", "CCND1", "MKI67", "CD300A", "CLEC9A", "ID2", "BATF3", "IRF8", "RUNX1", "FLT3", "HMMR", "MYBL2", "LAMP3", "JDP2", "CCL2", "IL7R", "MMP12", "GPR183", "GPR84", "HLA-DQB2", "CD300A", "ATP1B3"))

setCurrentAnalysis(paste0(cfg$ver,"/10x"))