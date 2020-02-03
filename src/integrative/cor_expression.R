#!/usr/bin/env $CODEBASE/Rscript
#
# Correlate gene expression differences between LCH subsets 
# with differences in chromatin accessibility.
#
# run("integrative", "cor_expression")


loadLibraries(c("Seurat","ggrepel"))

# ATAC-seq data:
simpleCache("atac_peaks_annot", assignToVar="peaksDt", reload=T)
simpleCache("atac_deseq2", assignToVar="tmp")
dAtac <- tmp$dNorm
dAtac <- cbind(peaksDt[,.(k, rid, geneSymbols, distToNextGene, clusterId)], dAtac[peaksDt$rid,])
dAtac <- unnest(dAtac, geneSymbols)

# expression data:
allDegs <- as.character(unlist(degs))
dRna <- fread(baseResultsDir(exprDir,"/d_expr_lch.csv"))
meta <- lib$dt2namedVec(fread(baseResultsDir(exprDir,"/d_meta_lch.csv"))[!is.na(cluster_id),.(ID=rn, cluster_id)]) #fread(baseResultsDir(exprDir,"/d_meta_lch.csv"))


# multiple peaks per gene, collapse per gene:
dRnaAtacX <- merge(dRna, dAtac, by.x="rn", by.y="geneSymbols", allow.cartesian=TRUE)
dRnaAtacX[,g:=rn]
cols1 <- c("k","g","rid","rn","clusterId","geneSymbols","distToNextGene")
cols2 <- atacOrder
cols3 <- setdiff(colnames(dRnaAtacX), c(cols1,cols2))
# there's only one unique expression value per gene, actually
simpleCache("rna_atac_epragg", {
	dExprAgg <- sapply(split(dRnaAtacX, by="g"), function(x) {
		x[1,cols3,with=F]
	})
	dExprAgg
}, assignToVar="dExprAgg")

# collapse ATAC-seq per gene by taking the mean across all associated peaks:
# (cache the result)
simpleCache("rna_atac_atacagg", {
	dAtacAgg <- list(
		expMean=sapply(split(dRnaAtacX, by="g"), function(x) {
			apply(x[,cols2,with=F],2,function (x) return(log2(mean(2^as.numeric(x) - 1) + 1)))
		})
	)
	dAtacAgg
}, assignToVar="dAtacAgg")


cors <- data.table()
aggMeth <-"expMean"

showGenes <- fread("http://humantfs.ccbr.utoronto.ca/download/v_1.01/TF_names_v_1.01.txt", header=F)$V1

aggMeth <- names(dAtacAgg)[1]
for(aggMeth in names(dAtacAgg)) {
	msg(aggMeth)
	tmp <- t(as.data.frame(dExprAgg))

	dRnaAtac <- data.table(g=rownames(tmp), tmp, t(dAtacAgg[[aggMeth]]))
	
	allClusts <- paste0("C", c(2, 5, 8))

	bg <- "C5"
	clusts <- setdiff(allClusts,bg)

	exprMeanFun <- Seurat::ExpMean

	iBg <- names(meta)[which(meta%in%modDict[[bg]])]
	rnaBg <- apply(data.matrix(dRnaAtac[,iBg,with=F]),1,exprMeanFun)
	atacBg <- rowMeans(data.matrix(as.data.frame(dRnaAtac[, grep(bg,atacOrder,value=T), with=F])))

	r <- cor(rnaBg, atacBg, method="pearson")
	s <- cor(rnaBg, atacBg, method="spearman")
	cors <- rbind(data.table(aggMeth=aggMeth, dataset=bg, pearson=r, spearman=s), cors)

	allPlotData <- data.table()

	for(fg in clusts) {
		fgX <- fg
		msg("\t", fgX, " vs. ", bg)
	
		fgXX <- modDict[[fgX]]
		iFg <- names(meta)[which(meta%in%fgXX)]
		rnaFg <- apply(data.matrix(dRnaAtac[,iFg,with=F]),1,exprMeanFun)
		atacFg <- rowMeans(data.matrix(as.data.frame(dRnaAtac[, grep(fg,atacOrder,value=T), with=F])))
		
		r <- cor(rnaFg, atacFg, method="pearson")
		s <- cor(rnaFg, atacFg, method="spearman")
		cors <- rbind(data.table(aggMeth=aggMeth, dataset=fgX, pearson=r, spearman=s), cors)

		plotData <- data.table(rnaDiff=rnaFg-rnaBg, atacDiff=atacFg-atacBg, gene=dRnaAtac$g)
		
		plotData[,rnaRnk:=rank(-abs(rnaDiff))]
		plotData[,atacRnk:=rank(-abs(atacDiff))]
		plotData[,rnk:=rank((rnaRnk+atacRnk)/2)]

		plotData[,sel:=rnk<=20]
		plotData[,show:=sel | gene%in%allDegs]

		plotData[,foreground:=fgX]

		allPlotData <- rbind(plotData, allPlotData)
	}

	allPlotData[, deg:="-"]
	forEach(c(bg, clusts), function(x) {
		allPlotData[gene%in%degs[[x]], deg:=x]
	})
	allPlotData[gene%in%degs[[fg]], deg:=fg] # just in case a gene is in multiple lists, give preference to the current foreground
	
	atacT <- 0.5
	rnaT <- 0.25
	
	allPlotData[,alwaysShow:=(tolower(gene)%in%showGenes & (abs(atacDiff)>=atacT | abs(rnaDiff)>=rnaT))]
	
	allPlotData[,isTf:=ifelse(gene%in%showGenes,"TF","NoTF")]
	p <- ggplot(allPlotData, aes(x=rnaDiff, y=atacDiff)) + defTheme(topLegend=T) + 
			geom_hline(yintercept=0, linetype=2) + geom_vline(xintercept=0, linetype=2) +
			geom_point(shape=16, size=0.25, alpha=0.5, data=allPlotData[(show==TRUE | alwaysShow==T),]) + # & isTf!="TF"
			geom_point(aes(color=deg, shape=isTf), data=allPlotData[(show==TRUE | alwaysShow==T) & isTf=="TF",]) +
			geom_text_repel(aes(label=gene), data=allPlotData[gene%in%allDegs & gene%in%showGenes & (abs(rnaDiff)>=rnaT | abs(atacDiff)>=atacT) & show==TRUE , ], min.segment.length=unit(0.000001, "lines")) +
			xlab(sprintf("scRNA-seq: X - %s (Scaled UMI, log2)", bg)) + ylab(sprintf("ATAC-seq: X - %s (Norm. RPM, log2)",  bg))  +
			coord_cartesian(xlim=c(-1.5,1.8), ylim=c(-1.8,1.8)) + scale_shape_manual(values=c("TF"=1, "NoTF"=16)) + #scale_color_manual(values=c("TF"="red", "NoTF"=alpha("black",0.5))) + 
			facet_wrap(~foreground, nrow=1, scales="fixed")  + scale_color_manual(values=c(colorPalettes$transition_state, "-"="black"))
	gg(p, paste_("scatter","all",bg, "agg", aggMeth), 5, 2.8, type="pdf", addBase=T)

}

