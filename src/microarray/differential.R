#!/usr/bin/env $CODEBASE/Rscript
#
# Perform differential analysis of microarray data to define cell type signatures.
#
# run("microarray", "differential")


# perform exploratory analysis of the data:

mdsData <- as.data.table(cmdscale(dist(t(dWide)), 2))
mdsData[,sampleName:=colnames(dWide)]
mdsData[,sampleGroup:=gsub("_?(\\d+|[A-Z])$","",sampleName)]
p <- ggplot(mdsData, aes(x=V1, y=V2, color=sampleGroup)) + geom_point(size=3) + defTheme() + defTheme() + xlab("MDS1") + ylab("MDS2") + scale_color_manual(values=colorPalettes$cell_type)
gg(p, "mds", 3.9, 8/3, type="pdf")

mod <- model.matrix(~0+sampleGroup, data=mdsData[sampleGroup%in%c("LC","CD14","CD1c","pDC"),.(sampleGroup)])
colnames(mod) <- gsub("sampleGroup","",colnames(mod))

degData <- dWide[,mdsData[sampleGroup%in%c("LC","CD14","CD1c","pDC"),sampleName]]

grpMeans <- sapply(split(as.data.frame(t(degData)), f=as.factor(mdsData[sampleGroup%in%c("LC","CD14","CD1c","pDC"),sampleGroup])), colMeans, na.rm=T)

lib$pdfPlot("deg_data", 6, 4)
boxplot(degData)
dev.off()


# use limma for differential expression analysis:
# (to define signatures of epidermal Langerhans cells [LC], CD14+ monocytes, CD1c dendritic cells, and pDCs compared to each other)

simpleCache("limma_fit", {
	fit <- lmFit(degData,mod)
	contDif <- makeContrasts(contrasts=list(LC="LC-(CD14+CD1c+pDC)/3", CD14="CD14-(LC+CD1c+pDC)/3", CD1c="CD1c-(LC+CD14+pDC)/3", pDC="pDC-(LC+CD14+CD1c)/3"),levels=mod)
	colnames(contDif) <- gsub("^(.+)-.+$","\\1",colnames(contDif))
	fit <- contrasts.fit(fit, contDif)
	fit <- eBayes(fit, robust=T)
	fit
}, recreate=F, assignToVar="fit")

deRes <- decideTests(fit, method="global", adjust.method="fdr", lfc=0, p.value=0.05)
print(summary(deRes))

deResX <- data.table(gene=rownames(dWide), grpMeans, dWide)

for(curCol in colnames(mod)) {
	tmp <- topTable(fit, coef=curCol, lfc=0, adjust.method="none", n=Inf)

	deResX[,paste_(curCol,"lfc"):=tmp[deResX$gene, "logFC"]]
	deResX[,paste_(curCol,"pval"):=tmp[deResX$gene, "P.Value"]]
}

tmp <- lib$padjustMatrix(as.matrix(deResX[,paste_(colnames(mod),"pval"), with=FALSE]))
colnames(tmp) <-  paste_(colnames(mod),"padj")
deResX <- cbind(deResX, tmp)
for(curCol in colnames(mod)) {
	deResX[, paste_(curCol,"deg"):=0]
	deResX[get(paste_(curCol,"padj"))<=0.05 & abs(get(paste_(curCol,"lfc")))>=log2(2), paste_(curCol,"deg"):=sign(get(paste_(curCol,"lfc")))]

	deResX[get(paste_(curCol,"deg"))!=0, paste_(curCol,"_down_rnk"):=rank((rank(get(paste_(curCol,"pval")))+rank(get(paste_(curCol,"lfc"))))/2, ties.method="random")]
	deResX[get(paste_(curCol,"deg"))!=0, paste_(curCol,"_up_rnk"):=rank((rank(get(paste_(curCol,"pval")))+rank(-get(paste_(curCol,"lfc"))))/2, ties.method="random")]
	deResX[get(paste_(curCol,"deg"))!=0, paste_(curCol,"_abs_rnk"):=rank((rank(get(paste_(curCol,"pval")))+rank(-abs(get(paste_(curCol,"lfc")))))/2, ties.method="random")]
}

# save lists of marker genes for each cell type of interest:
markerLists <- rblapply(colnames(mod), function(cType) {
	rblapply(c("-1","1"), function(dir) {
		deResX[get(paste_(cType,"deg"))==dir, .(gene)]
	}, "dir") 
}, "cellType")
markerLists[, geneSet:=paste_(cellType, dir)]

markerLists <- rbindlist(list(allDEGs=markerLists[,.(geneSet,gene)]), idcol="listType")

fwrite(markerLists, file=resultsDir("marker_lists.csv"))
