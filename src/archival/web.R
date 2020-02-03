#!/usr/bin/env $CODEBASE/Rscript
#
# Generate additional files for the supplementary website.
#
# run("archival", "web")


loadLibrary("R.utils")


### scRNA-seq data ###

run("single_cell", c("init","load"))
setCurrentAnalysis("archival")

simpleCache(paste_("scData", "cca", nCCS, allComboStr, "ext"), assignToVar="scDataMerged")

dMeta <- as.data.table(scDataMerged@meta.data, keep.rownames="cell_id")[, .(cell_id, barcode="placeholder", pat_id, perc_mito, G1, S, G2M, biopsy, sex, disease_extent, cluster_id, cell_type)]
dMeta[cell_type=="TNKcell", cell_type:="Tcell"]
dMeta[, barcode:=gsub("^.+__(.+)$","\\1",cell_id)]
simpleCache(paste_("scrnaseq_pat", "SCENT3", "cca", "stage1", nCCS, allComboStr), assignToVar="scDataSub")
dMeta[, lch_subset_id:=structure(finalDict[paste0("LCH-C",scDataSub@meta.data[,"cluster_id"])], names=rownames(scDataSub@meta.data))[as.character(cell_id)]]
dMeta[, cell_id:=paste_(pat_id, barcode)]
dMeta <- cbind(dMeta, scDataMerged@dr$tsne@cell.embeddings)

dRaw <- as.matrix(scDataMerged@raw.data)
dNorm <- round(as.matrix(scDataMerged@data),3)
colnames(dRaw) <- colnames(dNorm) <- dMeta$cell_id
dRaw <- as.data.table(dRaw, keep.rownames="gene_symbol")
dNorm <- as.data.table(dNorm, keep.rownames="gene_symbol")


fwrite(dMeta, file=resultsDir("lch_10x_meta.csv"))
fwrite(dRaw, file=resultsDir("lch_10x_raw.csv"))
fwrite(dNorm, file=resultsDir("lch_10x_norm.csv"))

gzip(resultsDir("lch_10x_meta.csv"),destname=resultsDir("lch_10x_meta.csv.gz"), overwrite=T)
gzip(resultsDir("lch_10x_raw.csv"),destname=resultsDir("lch_10x_raw.csv.gz"), overwrite=T)
gzip(resultsDir("lch_10x_norm.csv"),destname=resultsDir("lch_10x_norm.csv.gz"), overwrite=T)


pT <- cfg$deg_thresh_fdr
fcT <- cfg$deg_thresh_lfc
allGenesets <-  getMarkerLists(pT, fcT)

fwrite(allGenesets[["traj_cluster_id_sel_strat_w7"]]$genesets[p_val_adj<=pT & avg_logFC>=fcT, .(subset=finalDict[as.character(cluster)], gene, lfc=avg_logFC, pval=p_val, padj=p_val_adj)][order(subset,-lfc),], file=resultsDir("lch_deg_subsets.csv"))
fwrite(allGenesets[["LCH_vs_cTypes_strat"]]$genesets[p_val_adj<=pT,.(comparison=paste("LCH vs.", gsub("TNK","T",cluster)), gene, lfc=avg_logFC, pval=p_val, padj=p_val_adj)][order(comparison,-lfc),], file=resultsDir("lch_deg_vsnonlch.csv"))
fwrite(allGenesets[["LCH_vs_cTypes_strat_overlap"]]$genesets[p_val_adj<=pT & avg_logFC>=fcT,.(gene, lfc=avg_logFC, padj=p_val_adj)][order(-lfc),], file=resultsDir("lch_deg_signature.csv"))

selGenes <- c("CD1A","CD207","MAP2K1","ARAF","BRAF","FOXP3","MMP9","CD1E","MYC","CD33","CD3D","IL32","CD19","CD79A","IL3RA","CLEC4C","TCF4","CD14","CD163","HLA-DQB2","TACSTD2","CCL4","CCL3","LYZ","S100A8","CST3","S100A4","NFKBIA","IGHG1","IGKC","IRF7","IRF8")
forEach(selGenes, function(g) {
	lib$pngPlot(paste_("tsne","all", g), 6, 4, res=300)
	FeaturePlot(scDataMerged, g, no.legend = FALSE, vector.friendly = FALSE,  cols.use = colorPalettes$tsne)
	dev.off()
})

selGenes <- c("MKI67","AURKA","AURKB","CD83","LAMP3","CXCR4","CCR7","CLEC9A","BATF3","IRF8","MMP9","CD68","CD63","ANPEP","TUBB","TUBA1B","CD1A","CD207","BRAF","CDK1","IL22RA2","CD300A","IL7R","CXCR4","HLA-DQA2","IFI6","IFIT3","MMP12","S100A3","JDP2","MMP9","CD68","DNMT1","CENPA","SMYD3","KDM2B","PHF19","MXD3","TCF7L2","FOSL1","STAT3","BCL11A","NFATC1","MTA3","IRF4","GTF2B","CDK9","IRF1","KDM5C","RBBP5","KDM5B","HMGN3","MYB","MXI1","CDK7","CBX3","STAT2","BCL3","STAT5A","SOX4","PRDM1","SPIB","TFDP1","MITF","IRF8","REL","RELB","JUNB","ETV3","FOSL2","SPI1","KDM1A","MYBL2","BATF3","JDP2","KDM4A","BCL6","STAT1","STAT3","EP300")
forEach(selGenes, function(g) {
	lib$pngPlot(paste_("tsne","lch", g), 6, 4, res=300)
	FeaturePlot(scDataSub, g, no.legend = FALSE, vector.friendly = FALSE,  cols.use = colorPalettes$tsne)
	dev.off()
})




### ATAC-seq data ###

run("atac", c("init","load"))
setCurrentAnalysis("archival")

dA[, old_name:=sample_name]
for(curClust in c("C2","C5","C8")) {
	dA[, sample_name:=gsub(paste_("LCH",curClust),finalDict[[curClust]],sample_name)]	
	dA[, sample_group:=gsub(paste_("LCH",curClust),finalDict[[curClust]],sample_group)]	
}
dA[, sample_name:=gsub("_3$","_1",sample_name)]
dA[, sample_name:=gsub("_4$","_2",sample_name)]
setkey(dA, sample_name)

fwrite(dA[,.(sample_name, sample_group, flowcell, lane, organism, `library`)], file=resultsDir("lch_atac_meta.csv"))
gzip(resultsDir("lch_atac_meta.csv"),destname=resultsDir("lch_atac_meta.csv.gz"), overwrite=T)

forEach(dA$sample_name, function(s) {
	fwrite(dA[s, fread(peaks_file, sel=1:3, col.names=c("chrom","start","end"))], file=resultsDir("lch_atac_",s,"_peaks.bed"))
	gzip(resultsDir("lch_atac_",s,"_peaks.bed"),destname=resultsDir("lch_atac_",s,"_peaks.bed.gz"), overwrite=T)
})

simpleCache("atac_peaks_annot", {
	peaksDt
}, assignToVar="peaksDt", reload=T, recreate=F)

fwrite(peaksDt[, .(chrom=unique(chrom), start=unique(start), end=unique(end), module_names=paste(sort(finalDict[as.character(clusterId)]),collapse=",")), by=.(region_id=rid)], file=resultsDir("lch_atac_merged_peaks.csv"))
gzip(resultsDir("lch_atac_merged_peaks.csv"),destname=resultsDir("lch_atac_merged_peaks.csv.gz"), overwrite=T)

setkey(dA, old_name)

simpleCache("atac_counts", assignToVar="counts")	
dCounts <- lib$dtToDf(dcast(counts, paste0("r",regNum)~sampleName, value.var="count"))
#dCounts[counts[order(regNum), paste0("r",unique(regNum))],]
colnames(dCounts) <- dA[colnames(dCounts), sample_name]
dCounts <- as.data.table(dCounts[counts[order(regNum), paste0("r",unique(regNum))],], keep.rownames="region_id")

fwrite(dCounts, file=resultsDir("lch_atac_raw.csv"))
gzip(resultsDir("lch_atac_raw.csv"),destname=resultsDir("lch_atac_raw.csv.gz"), overwrite=T)



simpleCache("atac_deseq2", assignToVar="tmp", recreate=F)
dNorm <- round(tmp$dNorm,3)
colnames(dNorm) <- dA[colnames(dNorm), sample_name]
dNorm <- as.data.table(dNorm[counts[order(regNum), paste0("r",unique(regNum))],], keep.rownames="region_id")
fwrite(dNorm, file=resultsDir("lch_atac_norm.csv"))
gzip(resultsDir("lch_atac_norm.csv"),destname=resultsDir("lch_atac_norm.csv.gz"), overwrite=T)
