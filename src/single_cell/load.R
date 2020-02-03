#!/usr/bin/env $CODEBASE/Rscript
#
# Load in and process 10x scRNA-seq data. Store results in caches.
#
# run("single_cell", "load")

# load pre-processed (CellRanger) scRNA-seq data using Seurat -- patient by patient:
simpleCache(paste_("scrnaseq", paste0(sort(dA$sample_name),collapse="."), "complete"), {
	allData <- sapply(dA$sample_name, function(s) {
		msgF("==== %s ====", s)
		
		dir.create(resultsDir(s), showWarnings=F)	
			
		simpleCache(paste_("scrnaseq_pat", s, "filt"), {

			scData <- Read10X(dataDir(s, "/filtered_gene_bc_matrices/GRCh38/")) 
			scData <- CreateSeuratObject(raw.data = scData, min.cells = cfg$thresh_min_cells, min.genes = cfg$thresh_min_genes, project = s)
			
			mitoGenes <- grep(pattern = "^MT-", x = rownames(x = scData@data), value = TRUE)
			percMito <- Matrix::colSums(scData@raw.data[mitoGenes, ])/Matrix::colSums(scData@raw.data)
			scData <- AddMetaData(scData, metadata = percMito, col.name = "perc_mito")
			gg(VlnPlot(scData, features.plot = c("nGene", "nUMI", "perc_mito"), nCol = 3), "ngene_umi_mito", 9, 3, type="pdf", sub=s)		
			
			if(nrow(scData@raw.data)>=cfg$thresh_min_cells) {
				scData <- FilterCells(object = scData, subset.names = c("nGene", "perc_mito"), low.thresholds = c(cfg$thresh_min_genes, -Inf), high.thresholds = c(cfg$thresh_max_genes, cfg$thresh_max_mito))
			}		
			scData <- RenameCells(scData, add.cell.id=paste0(s,"_"))
						
			scData
		}, assignToVar="scData", recreate=F, noload=F)
		
		scData
	})
	
	qc <- data.table(
		sample_name=names(allData),
		ngenes=sapply(allData, function(d) nrow(d@raw.data)),
		ncells=sapply(allData, function(d) ncol(d@raw.data))
	)
	fwrite(qc, file=plotDir("qc_info.csv"))
	
	allData
}, assignToVar="allData", recreate=F, noload=T)

qc <- fread(plotDir("qc_info.csv"))
qcOk <- qc[ngenes>=cfg$thresh_min_ngenes_total & ncells>=cfg$thresh_min_ncells_total, sample_name]
msgF("%d samples loaded. %d samples passing quality control.", nrow(qc), length(qcOk))
dA <- dA[qcOk,]


# consider merging the datasets at different levels:
# per patient = no merging,
# per annotation, e.g. same biopsy tissue,
# or all samples together
combosToAnalyze <- c(
	as.list(dA[,unique(sample_name)]),
	# # all leave-one-out combinations:
	# combn(dA$sample_name, nrow(dA)-1, FUN = function(i) sort(dA$sample_name[i]), simplify = F),	
	# all datasets:
	list(sort(dA$sample_name))
)
for(patGrp in patientGroupings) {
	lvls <- setdiff(dA[,sort(unique(get(patGrp)))],c("","N/A","-"))
	
	combosToAnalyze <- c(combosToAnalyze, lapply(lvls, function(b) dA[get(patGrp)==b, sort(sample_name)]))
}
combosToAnalyze <- unique(combosToAnalyze)
combosToAnalyze <- combosToAnalyze[sapply(combosToAnalyze,length)>0]
writeCombos(combosToAnalyze)

allCombo <- dA[,sort(sample_name)]
allComboStr <- paste_("merged", paste_(allCombo,collapse="."))
allComboStrF <- paste_(allComboStr,"lch_only")	
dir.create(resultsDir(allComboStr), showWarnings=F)
dir.create(resultsDir(allComboStr,"_prelim"), showWarnings=F)
dir.create(resultsDir(allComboStrF), showWarnings=F)	
msgF("==== %s ====", allComboStr)	

nCCS <- length(allCombo)

# combine datasets using CCA + run t-SNE projection and clustering:
simpleCache(paste_("scData", "cca", nCCS, allComboStr), {
	msgF("Load data for %d samples...", length(allCombo))
	
	simpleCache(paste_("scrnaseq", paste0(sort(dA$sample_name),collapse="."), "complete"), assignToVar="allData")

	scData <- allData
	while(length(scData)>1) {
		scData <- c(MergeSeurat(scData[[1]],scData[[2]], scale.factor = 1000000),scData[-(1:2)])
	}		
	scData <- scData[[1]]
	scData <- process10xData(scData, scoreCellCyle=T, doClustering=F, doDimRed=F, doScaling=T, callLCH=F, doVarGenes=T)

	scDataSplit <- lapply(allCombo, function(s) {
		scDataX <- SubsetData(object = scData, cells.use=scData@meta.data$orig.ident==s)
		scDataX@meta.data[,"group"] <- s
		scDataX
	})
	
	varG <- scData@var.genes

	msgF("CCA analysis with %d components", nCCS)
	
	selClustering <- cfg$sel_clustering
	
	scDataCur <- RunMultiCCA(scDataSplit, genes.use=varG, num.ccs=nCCS)
	
	lib$pdfPlot(paste_("tsne_cca_dims", nCCS), 12, 12, sub=allComboStr)
	DimHeatmap(object = scDataCur, reduction.type = "cca", cells.use = 500, dim.use = 1:nCCS, do.balanced = TRUE)
	dev.off()
	
	scDataCur <- AlignSubspace(scDataCur, reduction.type="cca", grouping.var="group", dims.align=1:nCCS)
	
	scDataCur <- RunTSNE(scDataCur, reduction.use="cca.aligned", dims.use = 1:min(10,nCCS), do.fast=T, check_duplicates=F)	
	scDataCur <- FindClusters(scDataCur, reduction.type="tsne", dims.use=1:2, save.SNN=T, resolution=selClustering$res, force.recalc=T)
	scDataCur@meta.data[,"cluster_id"] <- scDataCur@meta.data[,selClustering$name]
	scDataCur <- callLCH(scDataCur)
	
	p1 <- TSNEPlot(scDataCur, do.return=T, pt.size=0.5, group.by="group")
	p2 <- TSNEPlot(scDataCur, do.return=T, pt.size=0.5, group.by="cluster_id")
	p <- plot_grid(p1, p2)
	gg(p, paste_("tsne_cca", nCCS), 12, 4, type="pdf", sub=allComboStr)
	  
	plotTsneGenePanels(scDataCur, selMarkersX=unique(c("CD1A","CD207","THBD","CD19","CD14","CD3D","IL3RA","CLEC9A","HMMR","LAMP3","MMP12",cfg$lch_neg_list)), sub=allComboStr, suffix=paste_("_cca", nCCS),steps="custom")

	cols <- setdiff(setdiff(colnames(dA), colnames(scDataCur@meta.data)), c("annot_file","sample_name","use"))
	scDataCur@meta.data[,cols] <- as.matrix(dA[as.character(scDataCur@meta.data$group), cols, with=F])
	p <- makeTsneMetadataPanel(scDataCur, "CCA", cols=c("cluster_id","pat_id","is_lch_cell","sex","biopsy","disease_extent"), ncols=4)
	gg(p, paste_("tsne_meta", nCCS), 12, 8, type="pdf", sub=allComboStr)

	scDataCur
}, assignToVar="scDataCur", recreate=F, noload=T)


# subset on LCH cells, re-run integration, t-SNE projection, and clustering:
simpleCache(paste_("scrnaseq_pat", "SCENT3", "cca", "stage1", nCCS, allComboStr), {

	subName <- paste_("traj","cca",nCCS)

	simpleCache(paste_("scData", "cca", nCCS, allComboStr), assignToVar="scDataMerged")

	varGOld <- scDataMerged@var.genes
	scDataFilt <- SubsetData(object = scDataMerged, cells.use=scDataMerged@meta.data$is_lch_cell)
	scDataFilt <- process10xData(scDataFilt, scoreCellCyle=T, doClustering=F, doDimRed=F, doScaling=T, callLCH=F, doVarGenes=T)
	varG <- scDataFilt@var.genes

	scDataSplit <- lapply(allCombo, function(s) {
		scDataX <- SubsetData(object = scDataFilt, cells.use=scDataFilt@meta.data$is_lch_cell & scDataFilt@meta.data$group==s )
		scDataX
	})
		
	
	selClustering <- cfg$sel_clustering

	scDataCur <- RunMultiCCA(scDataSplit, genes.use=varG, num.ccs=nCCS)

	dir.create(resultsDir(subName), showWarnings=F)

	lib$pdfPlot(paste_("tsne_cca_dims", nCCS), 12, 12, sub=subName)
	DimHeatmap(object = scDataCur, reduction.type = "cca", cells.use = 500, dim.use = 1:nCCS, do.balanced = TRUE)
	dev.off()

	scDataCur <- AlignSubspace(scDataCur, reduction.type="cca", grouping.var="group", dims.align=1:nCCS)

	scDataCur <- RunTSNE(scDataCur, reduction.use="cca.aligned", dims.use = 1:min(10,nCCS), do.fast=T, check_duplicates=F)

	scDataCur <- FindClusters(scDataCur, reduction.type="tsne", dims.use=1:2, save.SNN=T, resolution=selClustering$res, force.recalc=T)
	scDataCur@meta.data[,"cluster_id"] <- scDataCur@meta.data[,selClustering$name]

	p1 <- TSNEPlot(scDataCur, do.return=T, pt.size=0.5, group.by="group")
	p2 <- TSNEPlot(scDataCur, do.return=T, pt.size=0.5, group.by="cluster_id")
	p <- plot_grid(p1, p2)
	gg(p, paste_("tsne_cca", nCCS), 12, 4, type="pdf", sub=subName)
	  
	plotTsneGenePanels(scDataCur, selMarkersX=unique(c("CD1A","CD207","CLEC9A","HMMR","LAMP3","MMP12",cfg$lch_neg_list)), sub=subName, suffix=paste_("_cca", nCCS),steps="custom")

	cols <- setdiff(setdiff(colnames(dA), colnames(scDataCur@meta.data)), c("annot_file","sample_name","use"))
	scDataCur@meta.data[,cols] <- as.matrix(dA[as.character(scDataCur@meta.data$group), cols, with=F])
	p <- makeTsneMetadataPanel(scDataCur, "CCA", cols=c("cluster_id","pat_id","is_lch_cell","sex","biopsy","disease_extent"), ncols=4)
	gg(p, paste_("tsne_meta", nCCS), 12, 8, type="pdf", sub=subName)
	
	scDataCur
}, assignToVar="scData", noload=T)

# call cell types by marker gene expression and generate extended cache:
simpleCache(paste_("scData", "cca", nCCS, allComboStr,"ext"), {
	simpleCache(paste_("scData", "cca", nCCS, allComboStr), assignToVar="scDataMerged")

	cTypes <- c("Bcell"="CD19","TNKcell"="CD3D",MacMono="CD14","DC"="IL3RA")

	scDataMerged@meta.data$cell_type <- NA

	for(ct in names(cTypes)) {
		scDataMerged <- callCellType(scDataMerged, ct, posList=cTypes[[ct]], negList=c("CD1A","CD207"))
		scDataMerged@meta.data[scDataMerged@meta.data[,paste_("is",ct)],"cell_type"] <- ct
	}	
	scDataMerged@meta.data[scDataMerged@meta.data$is_lch_cell,"cell_type"] <- "LCH"
	
	plotTsneMetadata(scDataMerged, sub=allComboStr, cols=c("cell_type",paste_("is",names(cTypes))))	
	
	as.data.table(scDataMerged@meta.data)[, .N, by=.(cell_type,pat_id)][order(cell_type,pat_id),]
	
	scDataMerged
}, noload=T)



# plot features, clusters, etc.
if(doAllPlots) {
	msg("Generating plots for complete dataset...")
	
	simpleCache(paste_("scData", "cca", nCCS, allComboStr), assignToVar="scDataMerged")	
	
	msg("Number of cells per dataset:")
	print(as.data.table(scDataMerged@meta.data)[,.(total=.N, lch=sum(is_lch_cell)),by=pat_id][order(pat_id),])
	
	tmp <- colorPalettes$tsne
	colorPalettes$tsne <- rev(brewer.pal(8,"Spectral"))	
	
	scDataMerged@meta.data$nUMIminMax <- scaleMinMax(scDataMerged@meta.data$nUMI)
	scDataMerged@meta.data$nGeneMinMax <- scaleMinMax(scDataMerged@meta.data$nGene)
	scDataMerged@meta.data$perc_mito_minMax <- scaleMinMax(scDataMerged@meta.data$perc_mito)
	
	p <- plotTsneMetadataPanel(scDataMerged, n="sel", suffix="_selMeta2", cols=c("nUMIminMax","nGeneMinMax","perc_mito_minMax","G1","G2M","S"), sub=allComboStr, ncols=3)
	
	plotTsneMetadata(scDataMerged, sub=allComboStr)	
	
	colorPalettes$tsne <- tmp
		
	plotTsneGenePanels(scDataMerged, suffix="_fig1e", selMarkersX=c("CD1A","CD207","MAP2K1","CD33","BRAF","CD3D","IL32","FOXP3","CD19","CD79A","IL3RA","CLEC4A","TCF4","CD14","CD163","CD1C"), sub=allComboStr, steps="custom")	
	
	cnames <- intersect(colnames(scDataMerged@meta.data), colnames(dA))
	scDataMerged@meta.data[,cnames] <- dA[as.character(scDataMerged@meta.data$group), cnames, with=F]
	p <- ggplot(as.data.table(scDataMerged@meta.data)[,.N,by=.(cell_type=ifelse(is_lch_cell,"LCH","non_LCH"),pat_id)], aes(x=pat_id, y=N, fill=cell_type)) + geom_bar(stat="identity", width=0.5) + defTheme() + scale_fill_manual(values=colorPalettes$cell_type) + xlab(NULL) + ylab("Cell number")
	gg(p, "bars_cell_num", 6, 1.5, type="pdf", sub=allComboStr, expand0=T)
	
}	