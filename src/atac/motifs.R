#!/usr/bin/env $CODEBASE/Rscript
#
# Scan peaks for matches to known motifs of DNA-binding transcription factors.
#
# run("atac", "motifs")

simpleCache("atac_peaks", assignToVar="peaks")
peaks <- peaks$peaksDt

loadLibrary("BSgenome.Hsapiens.UCSC.hg38")
gen <- BSgenome.Hsapiens.UCSC.hg38

motifDbs <- list(
	"hocomoco_full_v11"= list(meme=lib$resourceDir("motifs/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme"), anno=lib$resourceDir("motifs/HOCOMOCOv11_full_annotation_HUMAN_mono.tsv"))
)


#### RUN MOTIF SEARCH WITH FIMO  ####

dir.create(resultsDir("fimo_out/"), showWarnings=F)

msg("Export sequences...")
seqFile <- resultsDir("fimo_out/peaks.fa")
if(!file.exists(seqFile)) {
	seqs <- getSeq(gen, lib$dtToGr(peaks))
	seqFile <- writeFASTA(seqs, seqFile)
} else {								
	msg("\t* sequence file exists -- skipped: ", seqFile)
}
for(motifDb in names(motifDbs)) {
	motifFile <- motifDbs[[motifDb]]$meme

	if(file.exists(motifFile)) {
		msg("Scan for motif occurrences in ", motifDb, " (", motifFile, ")...")
		
		f <- resultsDir("fimo_out/fimo_", motifDb, "_t0.001.txt")
		if(!file.exists(f)) {			
			f <- runFIMO(seqFile, outFile=f, motifFile=motifFile, inputIsFile=TRUE, fimoExec=)
		} else {								
			msg("\t* output file exists -- skipped: ",f)
		}
		f
	}
	else {
		msg("ERROR: Unknown motif database: ", motifDb, " (", motifFile, ")...")			
	}
}



#### ANALYZE RESULTS OF MOTIF SEARCH ####

motifDbSel <- "hocomoco_full_v11"
annoMotif <- fread(motifDbs[[motifDbSel]]$anno, key="Model")
lchClusts <- paste0("C", c(2,5,8)) # N.B. these IDs refer to an outdated nomenclature for LCH subsets (see common.R)
selMotifDbs <- c(motifDbSel)
motifPThresh <- 1e-4  # default on FIMO website


simpleCache("atac_peaks_annot", assignToVar="peaksDt", reload=T)

# read in FIMO results and merge with peak annotations:
fimo <- rbindlist(sapply(selMotifDbs, function(motifDb) {
	motifFile <- motifDbs[[motifDb]]$meme
	resultsFile <- resultsDir("fimo_out/fimo_", motifDb, "_t0.001.txt")
	data.table(readFIMOResult(resultsFile, motifFile=motifFile, motifPThresh=motifPThresh), file=resultsFile)	
}, simplify=FALSE), idcol="motifDb")
fimo[,file:=NULL]
fimo[,motifName:=NULL]
fimo <- cbind(fimo, rid=peaks[as.numeric(gsub("seq","",fimo$seq)), rid])
fimo <- merge(fimo, peaksDt, by="rid", allow.cartesian=T)


# calculate motif enrichments using Fisher's exact test:

nMod <- lib$dt2namedVec(peaksDt[!is.na(clusterId),.N,by=.(ID=clusterId)])
n <- peaksDt[, length(unique(rid))]

fimoWide <- lib$dtToDf(dcast(fimo, rid~motifId))
motifEnrich <- fimo[!is.na(clusterId), .(hitsMod=length(unique(rid))), by=.(clusterId, motifId)]
motifEnrich[, motifName:=annoMotif[motifId, `Transcription factor`]]
motifEnrich[, hitsBg:=as.numeric(colSums(fimoWide>0)[motifId])]
motifEnrich[, totalMod:=as.numeric(nMod[as.character(clusterId)])]
motifEnrich[, percMod:=hitsMod/totalMod]
motifEnrich[, percBg:=hitsBg/n]
motifEnrich[, percOnlyBg:=(hitsBg-hitsMod)/(n-totalMod)]
motifEnrich[, logOdds:=log2( percMod / percOnlyBg )]

motifEnrich[, pval:=apply(motifEnrich,1,function(x) {
	hitsMod <- as.numeric(x["hitsMod"])
	hitsBg <- as.numeric(x["hitsBg"])
	totalMod <- as.numeric(x["totalMod"])
	m <- matrix(c(hitsMod,totalMod-hitsMod,hitsBg-hitsMod,n-totalMod-hitsBg+hitsMod),nrow=2,byrow=T)
	fisher.test(m)$p.value
})]
motifEnrich[, padj:=p.adjust(pval,method="fdr")]

# select motifs by p-value, odds ratio, and only those that belong to TFs which are expressed:
filtGenes <- subsetMarkers[,unique(gene)]
motifsSel <- motifEnrich[padj<=0.05 & logOdds>=log2(1.5) & annoMotif[motifId, `Transcription factor`]%in%filtGenes,.(clusterId, motifId, motifName, logOdds)][order(clusterId),]
motifsUnsel <- motifEnrich[padj<=0.05 & logOdds>=log2(1.5) & !(annoMotif[motifId, `Transcription factor`]%in%filtGenes),.(clusterId, motifId, motifName, logOdds)][order(clusterId),]
tmp <- motifEnrich[percMod>=0.05 & percMod>percBg & motifId%in%c(motifsSel$motifId,motifsUnsel$motifId),.("Region module"=moduleNames[as.character(clusterId)], "Motif ID"=motifId, "Transcription factor"=motifName, "% hits in region module"=percMod, "% hits in background"=percBg, "P-value"=pval, "Adjusted p-value"=padj)]
fwrite(tmp, file=resultsDir("motif_enrich.csv"))
motifsSelIds <- motifsSel[, unique(motifId)]

lib$pdfPlot("phm_motif_enrichment", 12, 5)
hmD <- data.matrix(lib$dtToDf(dcast(motifEnrich[motifId%in%motifsSelIds, ], clusterId~motifName, fun.aggregate=max, value.var="logOdds")))
hmD[!is.finite(hmD)] <- NA
pheatmap(hmD, main="enrichment of motif occurrence relative to background", cluster_rows=F, cluster_cols=myHclust(t(hmD)), border_color="white", cellheight=8, cellwidth=10, col=colorRampPalette(c(brewer.pal(5,"PuRd")))(7), treeheight_column=10, treeheight_row=10,  number_format="%s", number_color="black", display_numbers=F) 
hmD <- data.matrix(lib$dtToDf(dcast(motifEnrich[motifId%in%motifsSelIds, ], clusterId~motifName, fun.aggregate=max, value.var="percMod")))
hmD[!is.finite(hmD)] <- NA
pheatmap(hmD*100, main="percent of peaks in module", cluster_rows=F, cluster_cols=myHclust(t(hmD)), border_color="white", cellheight=18, cellwidth=18, col=colorRampPalette(c(brewer.pal(5,"PuRd")))(7), treeheight_column=10, treeheight_row=10,  number_format="%.1f", number_color="black", display_numbers=T) 
dev.off()
	
	


# display p-values, odds ratios, and accessibility scores for selected motifs in heatmaps:

simpleCache("atac_deseq2", assignToVar="tmp")
d <- tmp$dNorm 
hmD <- lib$dtToDf(dcast(motifEnrich[motifId%in%motifsSelIds, ], clusterId~motifName, fun.aggregate=max, value.var="logOdds"))
hmDpval <- lib$dtToDf(dcast(motifEnrich[motifId%in%motifsSelIds, ], clusterId~motifName, fun.aggregate=function(x) min(x[is.finite(x)])[1], value.var="padj"))
hmDp <- lib$pToSig(hmDpval)
dimnames(hmDp) <- dimnames(hmDpval)

ub <- 4
m <- melt(as.data.table(fimoWide[,motifsSelIds],keep.rownames="rid"), id.vars="rid", variable.name="motifId")[value>0,]
m[, motifTf:=annoMotif[as.character(motifId), `Transcription factor`]]
hmAtac <- merge(m[,.(rid,motifTf)], as.data.table(d,keep.rownames="rid"), by="rid")
hmAtac <- melt(hmAtac, measure.vars=atacOrder)
hmAtac <- lib$dtToDf(dcast(hmAtac, motifTf~variable, value.var="value", fun.aggregate=mean))#

lib$pdfPlot("phm_net_motifs", 12, 5)
hm <- pheatmap(hmD[moduleOrder,], cluster_rows=F, cluster_cols=T, border_color="white", cellheight=8, cellwidth=10,  breaks=seq(0,ub,length.out=ub*3), col=colorRampPalette(c(brewer.pal(5,"PuRd")))(ub*3-1), treeheight_column=10, treeheight_row=10,  number_format="%s", number_color="black", display_numbers=(hmDp[moduleOrder,])) 
pheatmap(t(hmAtac)[atacOrder,colnames(hmD)], cluster_rows=F,  cluster_cols=hm$tree_col, scale="none", border_color="white", cellheight=8, cellwidth=10,  breaks=seq(2.6,4.3,length.out=21), col=colorRampPalette(c("white",brewer.pal(9,"YlOrRd"),"black"))(20), treeheight_column=0, treeheight_row=10,  number_color="black")	
for(patId in names(clusterExpr)) {	
	hmExpr <- clusterExpr[[patId]][colnames(hmD),]
	rownames(hmExpr) <- colnames(hmD)
	hmExpr[hmExpr>2] <- 2
	pheatmap(t(hmExpr)[cOrderSel,], main=sprintf("scRNA-seq: %s", patId), cluster_rows=F, cluster_cols=hm$tree_col, border_color="white", cellheight=8, cellwidth=10, col=colorRampPalette(c("white",brewer.pal(5,"YlGnBu")[2:5],"#000033"))(32), treeheight_column=10, treeheight_row=10) 
}
dev.off()

selTfs <- colnames(hmD)

tmp <- lib$dtToDf(dcast(m[motifTf%in%selTfs,], rid~motifTf))>0
write.table(tmp, file=lib$plotDir("motif_peak_overlaps.csv"), sep=",", col.names=NA, row.names=T, quote=F) 















