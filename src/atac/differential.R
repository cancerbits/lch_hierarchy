#!/usr/bin/env $CODEBASE/Rscript
#
# Perform differential analysis of ATAC-seq data to define regulatory modules.
#
# run("atac", "differential")

# use DESeq2 to perform differential accessibility analysis:
simpleCache("atac_deseq2", {
	loadLibrary("DESeq2")
	
	simpleCache("atac_counts", assignToVar="counts")	
	dCounts <- lib$dtToDf(dcast(counts, paste0("r",regNum)~sampleName, value.var="count"))
	
	grpVar <- "sample_group"
	dds <- DESeqDataSetFromMatrix(countData = dCounts, colData = lib$dtToDf(dA[colnames(dCounts),c("sample_name",grpVar),with=F]), design= as.formula(sprintf("~ %s", grpVar)))
	dds <- DESeq(dds, test="Wald", fitType="local", betaPrior=T)

	grps <- dA[,unique(get(grpVar))]
	allRes <- rblapply(grps, function(A) { 
			msg(A)
			cont <- list(paste0(grpVar,A), paste0(grpVar,setdiff(grps,A)))
			return(cbind(as.data.table(results(dds, contrast=cont, pAdjustMethod="fdr", format="DataFrame", tidy=TRUE)),B="Rest"))
	}, "A")
	allRes[, rid:=row]
	allRes[, row:=NULL]
 	
	allRes[, cmp:=paste_(A,"v",gsub("LCH_","",B))]
	allRes[, dir:=sign(log2FoldChange)]

	d <- log2(DESeq2::counts(dds, normalized=TRUE, replaced=TRUE)+1)

	write.table(d, file=resultsDir("atac_peak_intensities.csv"), sep=",", col.names=NA, row.names=T)
		
	lib$pdfPlot("boxplots_counts_norm", 6, 15)
	par(mfrow=c(2,1))
	boxplot(log2(dCounts+1), las=3, main="counts")
	boxplot(d, las=3, main="deseq2")
	dev.off()

	p <- ggplot(allRes, aes(log2FoldChange, color=A)) + geom_freqpoly() + defTheme(noLegendTitle=T) + scale_color_manual(values=colorPalettes$sample_group)
	gg(p, "hist_lfc", 5, 4, type="pdf")

	return(list(dds=dds, allRes=allRes, dNorm=d))
}, assignToVar="tmp", recreate=F)

allRes <- tmp$allRes
d <- tmp$dNorm
rm(tmp)

# define significant enrichments:
selI <- allRes[,which(abs(log2FoldChange)>=lfcThreshDEG & padj<=qThreshDEG)]
selIds <- allRes[selI,sort(unique(rid))]
print(allRes[selI,.N,by=cmp])

corTable <- cor(d)
colnames(corTable) <- rownames(corTable) <- gsub("atac_", "", colnames(corTable))
lib$pdfPlot("cortable_deqseq2", 8.3, 8)
pheatmap(corTable, border_color="white", col=colorRampPalette(rev(brewer.pal(7, "Spectral")))(19), breaks=seq(0.1,1,length.out=20), display_numbers=T, number_color="black")
dev.off()

corTable <- cor(d[selIds,])
colnames(corTable) <- rownames(corTable) <- gsub("atac_", "", colnames(corTable))
lib$pdfPlot("cortable_deqseq2_degsonly", 8.3, 8)
pheatmap(corTable, border_color="white", col=colorRampPalette(rev(brewer.pal(7, "Spectral")))(19), breaks=seq(0.1,1,length.out=20), display_numbers=T, number_color="black")
dev.off()


simpleCache("atac_peaks", assignToVar="tmp")

peaks <- tmp$peaks
peaksDt <- tmp$peaksDt
rm(tmp)

# get chromatin accessibility scores for all differentially accessible regions:
hmD <- d[selIds, ]
hmAnnot <- lib$dtToDf(dcast(allRes[selI,], rid~cmp, value.var="dir"))
hmAnnot[is.na(hmAnnot)] <- 0
hmAnnotCol <- sapply(allRes[,unique(cmp)], function(x) {
	a <- gsub("(LCH_C.)_._(.+)","\\1",x)
	c("-1"="black","0"="white","1"=colorPalettes$sample_group[a])
}, simplify=F) 
hmDscaled <- t(apply(hmD, 1, scale))
dimnames(hmDscaled) <- dimnames(hmD)

# define region modules as those differentially more/less accessible in one LCH subset than in the others:
allRes[,clustId:=factor(paste0(ifelse(dir>0,"a","i"), gsub("LCH_C","",A)), levels=moduleOrder)]
allRes[,A:=factor(A, levels=paste_("LCH",cellOrder))]
allRes[,unique(clustId)]
clusterIds <- lib$dt2namedVec(allRes[selI,][order(A, -dir, -log2FoldChange),.(ID=rid, clustId)])

hmD <- hmD[names(clusterIds),]
hmDscaled <- hmDscaled[names(clusterIds),]
hmAnnot <- hmAnnot[names(clusterIds),]

# N.B. a few peaks occur multiple times since regulatory modules are NOT define mutually exclusively
# this has a minor impact on the clustering (which is only used for visualization). Make the row IDs unique:
uniqRid <- make.unique(rownames(hmD))
rownames(hmD) <- rownames(hmDscaled) <- rownames(hmAnnot) <- uniqRid
hmAnnot$clusterId <- clusterIds
setkey(peaksDt, rid)

# add cluster annotations to peak list
noClustRid <- setdiff(peaksDt$rid, names(clusterIds))
peaksDt <- rbind(
	cbind(k=uniqRid, peaksDt[names(clusterIds),], clusterId=clusterIds),
	cbind(k=noClustRid, peaksDt[noClustRid,], clusterId=NA)
)
setkey(peaksDt, k )

# save annotated peaks:
simpleCache("atac_peaks_annot", {
	peaksDt
}, assignToVar="peaksDt", reload=T, recreate=F)

# assign colors to each module:
regionModuleColors <- apply(t(sapply(split(as.data.table(hmD), f=hmAnnot[rownames(hmD),"clusterId"]), colMeans)), 1, makeColor)
hmAnnotCol$clusterId <- regionModuleColors
hmCols <- colorRampPalette(c("black", "steelblue",  "#ffffbf", "goldenrod","firebrick2"))(32)


lib$pdfPlot("phm_degs_clusters", 4, 4)
hmMain <- pheatmap(hmD, gaps_row=cumsum(as.data.table(hmAnnot)[,.N,by=clusterId]$N), col=hmCols, cluster_rows=F, cutree_row=clustK, border_color="white", treeheight_row=10, treeheight_col=10, show_rownames=FALSE, annotation_row=hmAnnot, annotation_colors=hmAnnotCol, scale="none") #cluster_rows=hm$tree_row, cluster_cols=hm$tree_col, col=c("white","#555555"), col=colorRampPalette(rev(brewer.pal(7, "RdBu")))(19), breaks=seq(-3.5, 3.5, length.out=20),[,c("clusterId"),drop=F]
dev.off()

lib$pdfPlot("phm_degs_clusters_scaled", 4, 4)
hmMain <- pheatmap(hmDscaled, gaps_row=cumsum(as.data.table(hmAnnot)[,.N,by=clusterId]$N), col=hmCols, cluster_rows=F, cutree_row=clustK, border_color="white", treeheight_row=10, treeheight_col=10, show_rownames=FALSE, annotation_row=hmAnnot, annotation_colors=hmAnnotCol, scale="none") 
dev.off()



fwrite(peaksDt[!is.na(clusterId),.("Region ID"=rid, Chromosome=chrom, Start=start, End=end, "Region module"=moduleNames[clusterId], "LCH-S1 R1"=LCH_C5_3, "LCH-S1 R2"=LCH_C5_4, "LCH-S12 R1"=LCH_C2_3, "LCH-S12 R2"=LCH_C2_4, "LCH-S11 R1"=LCH_C8_3, "LCH-S11 R2"=LCH_C8_4)], file=resultsDir("differential_peaks.csv"))
fwrite(cbind(
	peaksDt[,.(region_id=rid, chrom, start, end)], 
	as.data.table(d[peaksDt$rid,])[, .(LCH_S1_3=LCH_C5_3, LCH_S1_4=LCH_C5_4, LCH_S12_3=LCH_C2_3, LCH_S12_4=LCH_C2_4, LCH_S11_3=LCH_C8_3, LCH_S11_4=LCH_C8_4)]
), file=resultsDir("atac_peaks_geo.csv"))
