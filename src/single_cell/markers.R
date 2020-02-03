#!/usr/bin/env $CODEBASE/Rscript
#
# Perform differential analysis between cell types (LCH vs other immune cells) and LCH cell subsets.
#
# run("single_cell", "markers")

# LCH cell subsets vs. each other:
simpleCache(paste_("clusterMarkers", "w", "traj_cluster_id_sel_strat_w7", allComboStr), {
	msg("Finding LCH subset markers...")
	simpleCache(paste_("scrnaseq_pat", "SCENT3", "cca", "stage1", nCCS, allComboStr), assignToVar="scDataMerged", reload=T)
	scDataMerged@meta.data$cluster_id <- paste0("LCH-C",scDataMerged@meta.data$cluster_id)
	msg("cells total = ", ncol(scDataMerged@scale.data))
	
	selClusts <- paste0("LCH-C",c(10,7,12,5,0,13))
	scDataMerged <- SubsetData(scDataMerged, cells.use=scDataMerged@meta.data$cluster_id %in% selClusts)
	msg("cells selected = ", ncol(scDataMerged@scale.data))
	
	# stratify differential comparisons by sample sex to avoid biases results due to differing composition 
	# of female and male cells in compared groups:
	clusterMarkers <- doStratifiedDifferentialAnalysis(scDataMerged, groupBy="cluster_id", stratBy="sex")
		
	clusterMarkers[sel==T,.(cluster,gene)]
	
		
	#clusterMarkers[,.(cluster,gene)]
	clusterMarkers[sel==T & cluster=="LCH-C12",.(cluster,gene)]
	
	clusterMarkers
}, recreate=F, noload=T)

# LCH cells vs. non-LCH immune cells:
simpleCache(paste_("clusterMarkers", "w", "LCH_vs_cTypes_strat", allComboStr), {
	msg("Finding LCH markers...")
	simpleCache(paste_("scData", "cca", nCCS, allComboStr, "ext"), assignToVar="scDataMerged", reload=T)

	cTypes <- c("Bcell"="CD19","TNKcell"="CD3D",MacMono="CD14","DC"="IL3RA")

	scDataMerged <- SubsetData(scDataMerged, cells.use=!is.na(scDataMerged@meta.data$cell_type))
	
	scDataMerged <- SetAllIdent(scDataMerged, "cell_type")
	
	
	# compare cell types but stratify differential comparisons by sample sex to avoid biases results due 
	# to differing composition of female and male cells in compared groups:
	
	stratBy <- "sex"
	groupBy <- "cell_type"
	strata <- unique(scDataMerged@meta.data[,stratBy])
		
	nPerClust <- lib$dt2namedVec(as.data.table(scDataMerged@meta.data)[,melt(table(get(stratBy))),by=groupBy][value>=20, length(unique(Var1)), by=.(ID=get(groupBy))])
	
	clusterMarkersTmp <- rblapply(strata, function(stratum) {
		scDataCur <- SubsetData(scDataMerged, cells.use=scDataMerged@meta.data[,stratBy]==stratum)
					
		clusterMarkers <- rblapply(names(cTypes), function(ctBg) {	
			msgF("\t* %s %s", stratum, ctBg)
			clusterMarkers <- FindMarkers(object = scDataCur, ident.1="LCH", ident.2=ctBg, only.pos = F)
			clusterMarkers <- as.data.table(clusterMarkers, keep.rownames="gene")
			clusterMarkers
		}, "cluster")
		
		clusterMarkers
	},"stratum")
	
	clusterMarkers <- clusterMarkersTmp[,.(.N, p_val=max(p_val), p_val_adj=max(p_val_adj), avg_logFC=mean(avg_logFC)),by=.(gene,cluster)]
	clusterMarkers <- clusterMarkers[N>=nPerClust[cluster],]
			

	clusterMarkers[, sel:=p_val_adj<=cfg$deg_thresh_fdr & abs(avg_logFC)>=log2(1.5)]
	
	clusterMarkers
}, recreate=F, noload=T)




# generate feature plots, heatmaps, etc. with selected marker genes:

if(doAllPlots) {
	msg("Generating plots for LCH markers...")

	subName <- paste_("traj","cca",nCCS)

	simpleCache(paste_("scrnaseq_pat", "SCENT3", "cca", "stage1", nCCS, allComboStr), assignToVar="scDataMerged")
	scDataMerged@meta.data$cluster_id <- paste0("LCH-C",scDataMerged@meta.data$cluster_id )
	
	markerVer <- "traj_cluster_id_sel_strat_w7"
	simpleCache(paste_("clusterMarkers", "w", markerVer, allComboStr), assignToVar="clusterMarkers")

	o <- fread(resultsDir(subName,"/clusters.csv"))$cluster
	
	clusterMarkers[, sel:=p_val_adj<=cfg$deg_thresh_fdr & avg_logFC>=cfg$deg_thresh_lfc]
	clusterMarkers[, cluster:=factor(gsub("_","-",cluster), levels=o)]
	
	clusterMarkers[,.N,by=cluster]
	
	
	if(!file.exists(resultsDir(subName,"/cluster_prototypes.csv"))) {
		clusterPrototypesScaled <- sapply(o, function(curClust) {
			rowMeans(as.matrix(scDataMerged@scale.data[,scDataMerged@meta.data$cluster_id==curClust]))
		})
		fwrite(clusterPrototypesScaled, file=resultsDir(subName,"/cluster_prototypes_scaled.csv"))
		
		for(patId in unique(scDataMerged@meta.data[,"pat_id"])) {
			clusterPrototypesScaled <- sapply(o, function(curClust) {
				rowMeans(as.matrix(scDataMerged@scale.data[,scDataMerged@meta.data$pat_id==patId & scDataMerged@meta.data$cluster_id==curClust]))
			})
			fwrite(clusterPrototypesScaled, file=resultsDir(subName,"/cluster_prototypes_scaled_",patId,".csv"))
		}
	}
	
	fwrite(clusterMarkers[sel==TRUE,], file=resultsDir(subName,"/cluster_markers_",markerVer,".csv"))
		
	markerFor <- clusterMarkers[sel==T, .(cluster=paste(cluster,collapse="/")), by=gene][order(cluster),.(gene,cluster)]    	
		
	m <- as.data.table(scDataMerged@meta.data, keep.rownames=T)	
	setkey(m, rn)
	
	m[, cluster_id:=factor(cluster_id, levels=o)]
	
	m[,.(id=finalDict[as.character(cluster_id)],.N),by=.(cluster_id, pat_id)][order(id,pat_id),.(id,pat_id,N)]#[,.N,by=id]
	
	hmData <- abscap(t(apply(as.matrix(scDataMerged@scale.data[markerFor$gene,]),1,scale)),0.95)
	colnames(hmData) <- colnames(scDataMerged@scale.data)
	
	highlightGenes <- c("IFI6","IFIT3","IFIT1","CD83","HLA-DQB2","LAMP3","CCR7","HLA-DQA2","CXCR4","HMMR","AURKB","MKI67","AURKA","IL22RA2","CLEC9A","BATF3","IRF8","S100A3","MMP9","NRP2","CD68","MYBL2","HSF2","JDP2","BATF3","FLT3","SPI1","EZH2","HDAC2","RUNX1","IRF8","BCL6","FOS","STAT3","JAK2","MAP2K1","BRAF","CBFB","CEBPB","CDK1","FANCG","CHEK1","MMP12","CCL2","CHD2","MAP2K6","MAP4K3","MAPK9","MAP3K13","MAP4K1","PLEKHM1",selMarkersX)
	
	lib$pdfPlot(paste_("hm_markers_all",markerVer), 15, 12, sub=subName)
	i <- m[order(cluster_id), rn]
	hmD <- hmData[,i]				
	lbls <- rownames(hmD)
	lbls[!lbls%in%highlightGenes] <- ""
	pheatmap(hmD, 
		main="all",
		col=colorPalettes$heatmap_expression, 
		annotation_row=lib$dtToDf(markerFor[,.(gene,geneset=cluster)]), 
		cluster_rows=F, cluster_cols=F, 
		gaps_row=cumsum(table(markerFor$cluster)), 
		gaps_col=cumsum(table(m[i,"cluster_id"])), 
		annotation_col=lib$dtToDf(m[i,.(rn,cluster=cluster_id,sample_name=group)]),
		treeheight_row=10, treeheight_col=10, 
		show_rownames=T, 
		fontsize_row=10,
		show_colnames=F,
		labels_row=lbls,
		annotation_color=colorPalettes, 
		scale="none"
	)
	dev.off()
	
	selClusts <- paste0("LCH-C",c(10,7,12,5,0,13))
	i <- m[cluster_id%in%selClusts,][order(factor(cluster_id, levels=selClusts)), rn]
	markerForSel <- clusterMarkers[sel==T & cluster%in%selClusts,][order(cluster),.(k=make.unique(gene),gene,cluster)] #[, .(cluster=paste(cluster,collapse="/")), by=gene]    	
	markerForSel[, cluster:=factor(cluster, levels=selClusts)]
	markerForSel <- markerForSel[order(cluster),]
	hmD <- hmData[markerForSel$gene,i]
	rownames(hmD) <- markerForSel$k
	lbls <- markerForSel$gene
	lbls[!lbls%in%highlightGenes] <- ""
	
	# plot the heatmap 3 times: once with just a few rows (PDF), once with just a few columns (PDF), and once completely (PNG)
	# --> this is a hack to enable high-res import into Inkscape
	
	pheatmap(hmD, 
		main="all",
		col=colorPalettes$heatmap_expression, 
		annotation_row=lib$dtToDf(markerForSel[,.(k, cluster)]), 
		cluster_rows=F, cluster_cols=F, 
		gaps_row=markerForSel[,.N,by=cluster][,cumsum(N)], 
		gaps_col= m[i,.N,by=cluster_id][,cumsum(N)], 
		annotation_col=lib$dtToDf(m[i,.(rn,cluster=cluster_id,pat_id=dA[group,pat_id])]),
		show_rownames=T, 
		fontsize_row=10,
		show_colnames=F,
		labels_row=lbls,
		annotation_color=colorPalettes, 
		scale="none",
		file=lib$resultsDir(subName,"/hm_markers_sel_",markerVer,".png")
	)
	pheatmap(hmD[,1:3], 
		main="col",
		col=colorPalettes$heatmap_expression, 
		annotation_row=lib$dtToDf(markerForSel[,.(k, cluster)]), 
		cluster_rows=F, cluster_cols=F, 
		gaps_row=markerForSel[,.N,by=cluster][,cumsum(N)], 
		show_rownames=T, 
		fontsize_row=10,
		show_colnames=F,
		labels_row=lbls,
		annotation_color=colorPalettes, 
		scale="none",
		file=lib$resultsDir(subName,"/hm_markers_sel_",markerVer,"_c.pdf")
	)
	pheatmap(hmD[1:3,], 
		main="row",
		col=colorPalettes$heatmap_expression, 
		gaps_col= m[i,.N,by=cluster_id][,cumsum(N)], 
		annotation_col=lib$dtToDf(m[i,.(rn,cluster=cluster_id,pat_id=dA[group,pat_id])]),
		cluster_rows=F, cluster_cols=F, 
		show_rownames=T, 
		fontsize_row=10,
		show_colnames=F,
		labels_row=lbls,
		annotation_color=colorPalettes, 
		scale="none",
		file=lib$resultsDir(subName,"/hm_markers_sel_",markerVer,"_r.pdf")
	)
				
	subName <- paste_("traj","cca",nCCS)

	scData <- scDataMerged
	
	plotTsneMetadata(scData, sub=subName, suffix="_traj5")
	
	fwrite(as.matrix(scDataMerged@data), file=resultsDir("traj_cca_7/d_expr_lch.csv"))
	fwrite(as.matrix(scDataMerged@meta.data), file=resultsDir("traj_cca_7/d_meta_lch.csv"))
		
	colorPalettes$tsne <- c("#DDDDDD",brewer.pal(5,"YlGnBu")[2:5])
	
	plotTsneGenePanels(scData, suffix="_x", w=16, h=12, plotsPerPage=16, selMarkersX=c("JDP2", "FOS", "CBFB", "CEBPB", "STAT3", "HSF2", "MYBL2", "BATF3", "SPI1", "EZH2", "RUNX1", "BCL6", "IRF8"), sub=subName, steps="custom")
	plotTsneGenePanels(scData, suffix="_facs", w=16, h=12, plotsPerPage=16, selMarkersX=c("CD1A", "CD207", "PTPRC", "CLEC9A", "HMMR", "LAMP3", "CD300A", "IL2RG", "ANPEP"), sub=subName, steps="custom")
	plotTsneGenePanels(scData, suffix="_trans", w=16, h=12, plotsPerPage=16, selMarkersX=c("CD1A","CD207",sort(c("HMMR","MYBL2","AURKB","CD300A","CLEC9A","LAMP3","JDP2","ANPEP","CXCR4","MMP12","CD68","CD83","IRF8","EZH2"))), sub=subName, steps="custom")
	
	selColors <- (unlist(lapply(c("pat_id","biopsy","sex","disease_extent"), function(curCol) {
		colorPalettes[[curCol]][unique(scData@meta.data[,curCol])]
	})))
	
	p <- plotTsneMetadataPanel(scData, n="sel", suffix="_trans_meta", forceColors=selColors, cols=c("pat_id","biopsy","sex","disease_extent"), sub=subName, ncols=2)
	
		
	exprDataCap <- scData@scale.data
	
	markerProfiles <- rblapply(split(clusterMarkers,by="cluster"), function(curMarkers) {
		gs <- curMarkers$gene
		data.table(barcode=colnames(exprDataCap), means=apply(exprDataCap[gs,,drop=F], 2, mean))				
	},"cluster")
	pData <- cbind(markerProfiles, scData@dr$tsne@cell.embeddings[markerProfiles$barcode,], patient=scData@meta.data[markerProfiles$barcode,c("pat_id")])
	pData[, cluster:=factor(cluster,levels=paste0("LCH-C",c(10,7,12,5,0,13)))]
	
	cols <- c(colorPalettes$tsne[1],brewer.pal(4,"YlOrBr"),"black")
	
	
	p <- ggplot(pData[means>0,], aes(x=tSNE_1, y=tSNE_2, color=means))  + scale_color_gradientn(colors=cols, guide=guide_colorbar(title="Mean scaled UMI count")) + defTheme(noLegendTitle=TRUE) + xlab("t-SNE 1") + ylab("t-SNE 2") + facet_wrap(~cluster,nrow=1) + coord_cartesian(xlim=range(pData$tSNE_1), ylim=range(pData$tSNE_2))
	gg(p + geom_point(size=0.5), "all_marker_profiles", 13, 2.5, type="pdf", addBase=T, sub=subName)
	gg(p, "all_marker_profiles_nodata", 13, 2.5, type="pdf", addBase=T, sub=subName)

	p <- ggplot(pData, aes(x=tSNE_1, y=tSNE_2, color=means)) + geom_point(size=0.25) + scale_color_gradientn(colors=cols, guide=guide_colorbar(title="Mean scaled UMI count")) + defTheme(noLegendTitle=TRUE) + xlab("t-SNE 1") + ylab("t-SNE 2") + facet_grid(cluster~patient) + coord_cartesian(xlim=range(pData$tSNE_1), ylim=range(pData$tSNE_2))
	gg(p, "all_marker_profiles_split", 12, 7, type="pdf", addBase=T, sub=subName)#



	exprDataCap <- scData@data
	
	indivGenes <- c("CD1A","CD207", "MKI67","HMMR","CLEC9A")
	markerProfiles <- rblapply(split(data.table(gene=indivGenes,cluster=indivGenes),by="cluster"), function(curMarkers) {
		gs <- curMarkers$gene
		data.table(barcode=colnames(exprDataCap), means=apply(exprDataCap[gs,,drop=F], 2, mean))				
	},"cluster")
	pData <- cbind(markerProfiles, scData@dr$tsne@cell.embeddings[markerProfiles$barcode,], patient=scData@meta.data[markerProfiles$barcode,c("pat_id")])
	pData[, cluster:=factor(cluster,levels=indivGenes)]
		
	p <- ggplot(pData, aes(x=tSNE_1, y=tSNE_2, color=means))  + scale_color_gradientn(colors=colorPalettes$tsne, guide=guide_colorbar(title="Scaled UMI count")) + defTheme(noLegendTitle=TRUE) + xlab("t-SNE 1") + ylab("t-SNE 2") + facet_wrap(~cluster,nrow=1) + coord_cartesian(xlim=range(pData$tSNE_1), ylim=range(pData$tSNE_2))
	gg(p + geom_point(size=0.5), "all_marker_profiles_indivGenes", 13, 2.5, type="pdf", addBase=T, sub=subName)#c("white","#F0F0F0",rev(heat.colors(3)))
	gg(p, "all_marker_profiles_indivGenes_nodata", 13, 2.5, type="pdf", addBase=T, sub=subName)#c("white","#F0F0F0",rev(heat.colors(3)))
	
	p <- ggplot(pData, aes(x=tSNE_1, y=tSNE_2, color=means)) + geom_point(size=0.25) + scale_color_gradientn(colors=colorPalettes$tsne, guide=guide_colorbar(title="Scaled UMI count")) + defTheme(noLegendTitle=TRUE) + xlab("t-SNE 1") + ylab("t-SNE 2") + facet_grid(cluster~patient) + coord_cartesian(xlim=range(pData$tSNE_1), ylim=range(pData$tSNE_2))
	gg(p, "all_marker_profiles_indivGenes_split", 12, 7, type="pdf", addBase=T, sub=subName)#



	exprDataCap <- scData@scale.data
		
	cols <- c(colorPalettes$tsne[1],brewer.pal(4,"YlOrBr"),"black") # "darkred",
		
	
	facsMarkers <- data.table(gene=c("CD1A","CD300A","HMMR", "CD1A","CLEC9A", "CD1A","LAMP3"), subset=c("C5","C5","C5", "C2","C2", "C8","C8"))
	exprDataCap <- scData@scale.data

	markerProfiles <- rblapply(split(facsMarkers,by="subset"), function(curMarkers) {
		gs <- curMarkers$gene
		data.table(barcode=colnames(exprDataCap), means=apply(exprDataCap[gs,,drop=F], 2, mean))				
	},"cluster")
	pData <- cbind(markerProfiles, scData@dr$tsne@cell.embeddings[markerProfiles$barcode,], patient=scData@meta.data[markerProfiles$barcode,c("pat_id")])
	
	cols <- c(colorPalettes$tsne[1],brewer.pal(4,"YlOrBr"),"black") # "darkred",
	
	p <- ggplot(pData[means>0,], aes(x=tSNE_1, y=tSNE_2, color=means))  + scale_color_gradientn(colors=cols, guide=guide_colorbar(title="Mean scaled UMI count")) + defTheme(noLegendTitle=TRUE) + xlab("t-SNE 1") + ylab("t-SNE 2") + facet_wrap(~cluster,ncol=1) + coord_cartesian(xlim=range(pData$tSNE_1), ylim=range(pData$tSNE_2))
	gg(p + geom_point(size=0.5), "facs_marker_profiles", 3, 6, type="pdf", addBase=T, sub=subName)#c("white","#F0F0F0",rev(heat.colors(3)))
	gg(p, "facs_marker_profiles_nodata", 3, 6, type="pdf", addBase=T, sub=subName)#c("white","#F0F0F0",rev(heat.colors(3)))

	p <- ggplot(pData, aes(x=tSNE_1, y=tSNE_2, color=means)) + geom_point(size=0.25) + scale_color_gradientn(colors=cols, guide=guide_colorbar(title="Mean scaled UMI count")) + defTheme(noLegendTitle=TRUE) + xlab("t-SNE 1") + ylab("t-SNE 2") + facet_grid(cluster~patient) + coord_cartesian(xlim=range(pData$tSNE_1), ylim=range(pData$tSNE_2))
	gg(p, "facs_marker_profiles_split", 12, 7, type="pdf", addBase=T, sub=subName)#

}

# ROC curves for discrimination of LCH vs. non-LCH cells based on LCH gene signature vs. all genes vs. random:
if(doAllPlots) {
	run("single_cell",c("init","load"))

	loadLibraries(c("glmnet","ROCR","caret"))

	pT <- cfg$deg_thresh_fdr
	fcT <- cfg$deg_thresh_lfc
	signatureGenes <-  getMarkerLists(pT, fcT)[["LCH_vs_cTypes_strat_overlap"]]

	X <- as.matrix(signatureGenes$data@data)
	y <- (as.factor(!is.na(signatureGenes$data@meta.data$cell_type) & signatureGenes$data@meta.data$cell_type=="LCH"))

	geneSels <- list(all=1:nrow(X), signature=signatureGenes$geneset$gene)

	for(rnd in 1:10) {
		rndGenes <- sample(rownames(signatureGenes$data@data), length(signatureGenes$geneset$gene))
		geneSels[[paste_("rnd", rnd)]] <- rndGenes
	}

	n <- 10

	allRes <- list()

	for(i in 1:8) {

		j <- 1:ncol(X) %% n == i
		
		for(g in names(geneSels)) {
			msgF("cv %d, sel %s", i, g)
			simpleCache(paste_("glmnet",i,g), {
				lassoMod <- cv.glmnet(x = t(X[geneSels[[g]],!j]), y = y[!j], family = 'binomial', type.measure = 'auc')
				res <- predict(lassoMod, type="response", newx = t(X[geneSels[[g]],j]), s = 'lambda.min')
				pred <- prediction(res, y[j])

				perf <- performance(pred,"tpr","fpr")
				auc <- performance(pred,"auc")
				
				list(i=i, g=g, perf=perf, auc=auc, mod=lassoMod)
			}, assignToVar="tmp")
			
			allRes[[paste_(g, i)]] <- tmp
		}
	}

	pCols <- c("all"="black", "signature"="red", "rnd"="grey")

	lib$pdfPlot("roc", 6, 6)
	plot(c(0,1),c(0,1),col = "gray", lty = 2, type="l", xlab="False positive rate", ylab="True positive rate")
	pData <- rbindlist(lapply(allRes, function(res) {
		#msg(gsub("_\\d+","",res$g))
		plot(res$perf, col=pCols[gsub("_\\d+","",res$g)], add=T)
		data.table(auc=as.numeric(res$auc@y.values), i=res$i, g=res$g)
	}))
	legend("bottomright", c(pData[g=="all",sprintf("All genes (AUC = %.2f)", mean(auc))],pData[g=="signature",sprintf("LCH gene signature (AUC = %.2f)", mean(auc))], pData[grepl("rnd",g),sprintf("Random gene signature (AUC = %.2f)", mean(auc))]), col=c("black","red"), lty=1)
	dev.off()

	pData <- rbindlist(lapply(allRes, function(res) {
		data.table(x=res$perf@x.values[[1]], y=res$perf@y.values[[1]], auc=as.numeric(res$auc@y.values), i=res$i, g=res$g)
	}))

	pCols2 <- pCols
	names(pCols2) <- c(pData[g=="all",sprintf("All genes (AUC = %.2f)", mean(auc))],pData[g=="signature",sprintf("LCH gene signature (AUC = %.2f)", mean(auc))], pData[grepl("rnd",g),sprintf("Random gene signature (AUC = %.2f)", mean(auc))])

	pData[, type:=factor(gsub("_\\d+","",g), levels=names(pCols), labels=names(pCols2))]
	p <- ggplot(pData, aes(x=round(x*100)/100, y=y, group=g, color=type, fill=type)) + stat_summary(geom="line") + defTheme(topLegend=T, noLegendTitle=T) + scale_color_manual(values=pCols2) + scale_fill_manual(values=alpha(pCols2,0.25)) + xlab("False positive rate") + ylab("True positive rate")
	gg(p, "roc2", 4, 4, type="pdf")+ stat_summary(geom="errorbar") 

}