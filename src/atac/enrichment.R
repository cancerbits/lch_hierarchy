#!/usr/bin/env $CODEBASE/Rscript
#
# Functional enrichment analysis using Enrichr and LOLA.
#
# run("atac", "enrichment")



##############################################################################
####             Enrichr: genes linked to regulatory modules             #####
##############################################################################

# install_github("wjawaid/enrichR")
loadLibrary("enrichR")

# get genes linked to regulatory modules:
moduleGenes <- sapply(as.character(unique(clusterIds)), function(cid) {
	genes <- sort(unique(unlist(peaksDt[!is.na(clusterId) & clusterId==cid & !is.na(geneSymbols),geneSymbols]))) 
	write.table(genes, file=lib$plotDir("cluster-",cid,"-genes.tsv"), sep="\t", col.names=F, row.names=F, quote=F)
	genes
}, simplify=FALSE)
allGenes <- sort(unique(c(unlist(peaksDt$geneSymbol)))) #,dM$targetName
write.table(allGenes, file=resultsDir("allexpr-genes.csv"), sep=",", col.names=F, row.names=F, quote=F)

enrichrDBs <- as.character(read.csv("metadata/enrichr_dbs.csv")[,1]) 
enrichrDBs <- enrichrDBs[!grepl("#",enrichrDBs)]
enrichrLists <- moduleGenes

# use Enrichr API to perform enrichment analysis:
simpleCache(paste0(paste_("enrichr_res_x",digest(enrichrDBs)),""), { 
	runEnrichr(enrichrLists, dbs=enrichrDBs)
}, recreate=F, reload=T, assignToVar="e")
e[, grp:=factor(grp,levels=moduleOrder)]
e[, combscore:= Combined.Score]
e[, pval:= P.value]
e[, padj:= Adjusted.P.value]

dbSels <- list(
	"sel"=c("ARCHS4_Tissues","Jensen_DISEASES","BioCarta_2016","GO_Biological_Process_2018","KEGG_2016","NCI-Nature_2016","MSigDB_Oncogenic_Signatures")
)
for(dbSel in names(dbSels)) {
	selEnrichrRes <- e[database%in%dbSels[[dbSel]],]
	selEnrichrRes[,database:=as.factor(database)]

	selEnrichrRes[,rnk:=NULL]
	forEach(selEnrichrRes[,unique(grp)], function(curCmp) {
		forEach(selEnrichrRes[,unique(database)], function(db) {
			selEnrichrRes[grp==curCmp & database==db & combscore>0, rnk:=rank(pval,ties.method="random")]
		})
	})
	
	fwrite(selEnrichrRes[rnk<=50,][order(grp,database,rnk),.("Group"=grp, "Database"=database, "Term"=term, "Hits"=n.hits, "Total list size"=n.total, "P-value"=pval, "Corrected p-value"=padj, "Combined score"=combscore, "Full term"=original.term, "Genes in overlap"=Genes)], file=resultsDir("enrichr_atac_",dbSel,".csv"))
	
	plotData <- selEnrichrRes[ is.finite(rnk) & rnk<=3,]
	mx <- plotData[,ceiling(max(combscore))*1.25]
	p <- ggplot(plotData, aes(x=reorder(substr(paste0(as.numeric(database), "::",grp,":: ",term),0,48),combscore), y=combscore)) + geom_bar(stat="identity", position="dodge", width=0.7) + defTheme(topLegend=T) + coord_flip() + scale_x_discrete(labels=function(x) { gsub("^\\d+::.+::(.+)$","\\1",x) }) + ylab("Enrichr combined score") + xlab(NULL) + facet_wrap(database~grp, scales="free", drop=F, dir="v", ncol=plotData[,length(unique(database))]) + geom_text(aes(label=lib$pToSig(padj),y=combscore*1.05, hjust=0, vjust=0.75)) + scale_y_continuous(expand=c(0.1,0)) #limits=c(0,mx), 
	gg(p, paste_("enrichr_atac_bars", dbSel), length(plotData[,unique(database)])*4, length(plotData[,unique(grp)])*1.1, addBase=T, type="pdf")
	p <- ggplot(plotData, aes(x=reorder(substr(paste0(as.numeric(database), "::",grp,":: ",term),0,48),combscore), y=combscore)) + geom_bar(stat="identity", position="dodge", width=0.7) + defTheme(topLegend=T) + coord_flip() + scale_x_discrete(labels=function(x) { gsub("^\\d+::.+::(.+)$","\\1",x) }) + ylab("Enrichr combined score") + xlab(NULL) + facet_wrap(database~grp, scales="free_y", drop=F, dir="v", ncol=plotData[,length(unique(database))]) + geom_text(aes(label=lib$pToSig(padj),y=combscore*1.05, hjust=0, vjust=0.75)) + scale_y_continuous(expand=c(0.1,0)) #limits=c(0,mx), 
	gg(p, paste_("enrichr_atac_bars_fixed", dbSel), length(plotData[,unique(database)])*4, length(plotData[,unique(grp)])*1.1, addBase=T, type="pdf")
}







##############################################################################
####                  LOLA: overlaps with ChIP-seq data                  #####
##############################################################################

loadLibrary("LOLA")

selCollections <- c("codex","encode_tfbs","cistrome_cistrome") #,"sheffield_dnase"[,unique(collection)] #lolaResAll[,setdiff(unique(collection), c("encode_tfbs", "codex", "sheffield_dnase"))], "cistrome_epigenome", "ucsc_features"
qThresh <- cfg$deg_thresh_fdr
percThresh <- 1/3*100
effCol <- "log2odds"
effThresh <- log2(2)
tfCollections <- c("codex","encode_tfbs","cistrome_cistrome")
selCollections <- tfCollections
filtBy <- "padj"

oddsRatioLimit <- 8

lolaUserSets <- sapply(split(peaksDt, peaksDt[,clusterId]), lib$dtToGr) 
lolaUserSets[["rest"]] <- lib$dtToGr(peaksDt[is.na(clusterId),])

regionDB <- loadRegionDB(lib$resourceDir("regions/LOLACore/", genomeVer))
regionDB <- curateRegionDb(regionDB)

simpleCache("lola", {
	selUserSets <- GRangesList(lolaUserSets) 
	univ <- buildRestrictedUniverse(selUserSets)
	runLOLA(selUserSets, univ, regionDB)
}, assignToVar="lolaResAll", recreate=F)


lolaResAll <- lolaResAll[collection%in%tfCollections,]
lolaResAll[, pval:=10^(-pValueLog)] # older versions of LOLA used log not log10 -- need to make sure this fits the version used (i.e. use exp(-p) for older versions instead)
lolaResAll[, padj:=p.adjust(pval, method="fdr")]
lolaResAll <- augmentLolaRes(lolaResAll, qThresh=qThresh, orCol="oddsRatio", qvalCol="padj", pvalCol="pval")

# rename some TFs to be compatible with the naming used for the scRNA-seq data:
# GR NR3C1
# UBF UBTF
# COREST RCOR1
# TBLR1 TBL1X
# NRSF REST
# INI1 SMARCB1
# PAX5C20 PAX5
# SIN3AK20 SIN3A
# PLU1 KDM5B 
# HAE2F1 E2F1

lolaResAll[term=="GR", term:="NR3C1"]
lolaResAll[term=="UBF", term:="UBTF"]
lolaResAll[term=="COREST", term:="RCOR1"]
lolaResAll[term=="TBLR1", term:="TBL1X"]
lolaResAll[term=="NRSF", term:="REST"]
lolaResAll[term=="INI1", term:="SMARCB1"]
lolaResAll[term=="PAX5C20", term:="PAX5"]
lolaResAll[term=="SIN3AK20", term:="SIN3A"]
lolaResAll[term=="PLU1", term:="KDM5B"]
lolaResAll[term=="HAE2F1", term:="E2F1"]

# filter to only genes that are expressed:
lolaResAll <- lolaResAll[term%in%as.character(unlist(clusterIsExprE)),]

# define significant enrichments:
lolaResAll[,sig:=!grepl("rest",userSet) & get(filtBy)<=qThresh & get(effCol)>=effThresh & collection%in%selCollections] #
lolaResAll[sig==T,.N,by=userSet]
sigI <-  lolaResAll[,which(sig)]
selFiles <- lolaResAll[sigI,unique(filename)]
selTerms <- lolaResAll[filename%in%selFiles,][,.(s=absmax(10*log10(get(filtBy)))),by=term][order(s),term] 

filtGenes <- subsetMarkers[,unique(gene)]
selTermsFocus <- matchLenient(selTerms,filtGenes) #intersect(selTerms, allDegs) #matchLenient(selTerms,allDegs) #intersect(matchLenient(selTerms,allDegs), lolaResAll[sig & grepl("a",userSet) , unique(term)]) # grepl("a",userSet)   #lolaResAll[sig & term%in%allDegs, unique(term)] # grepl("a",userSet) 
selTermsUnfocus <- setdiff(selTerms, selTermsFocus) #lolaResAll[sig & !(term%in%allDegs), unique(term)]  #grepl("i",userSet)

fwrite(lolaResAll[term%in%selTerms,][order(moduleNames[userSet],pval), .(`Regulatory region module`=moduleNames[userSet], collection, term, term%in%selTermsFocus, description, cellType, tissue, antibody, filename, a=support, b, c, d, log2odds, pval, padj)], file=resultsDir("lola_res_selTerms.csv"))

lolaResAll[,colValBin:="black"]
lolaResAll[filename%in%selFiles,colValBin:="red"]



# plot heatmaps of all significant LOLA results:

hmData <- as.data.frame(dcast(lolaResAll[term%in%selTerms,], userSet~term, value.var=effCol, fun.aggregate=function(x) max(x[is.finite(x)])[1]))
rownames(hmData) <- hmData[,1]
hmData <- as.matrix(hmData[,-1])
hmData[!is.finite(hmData)] <- min(hmData[is.finite(hmData)],na.rm=T)
hmDataP <- as.data.frame(dcast(lolaResAll[term%in%selTerms,], userSet~term, value.var=filtBy, fun.aggregate=function(x) min(x[is.finite(x)])[1]))
rownames(hmDataP) <- hmDataP[,1]
hmDataP <- as.matrix(hmDataP[,-1])
hmDataP <- matrix(lib$pToSig(hmDataP), ncol=ncol(hmData), dimnames=dimnames(hmData))
rownames(hmData) <- rownames(hmDataP) <- gsub("cluster_(.+)_n_.+$","\\1",rownames(hmData))
ub <- 4
hmData[hmData > ub] <- ub
hmData[hmData < -ub] <- -ub

lib$pdfPlot(paste_("lola","PHM","tfs"), 25, 5)
selR <- moduleOrder
colSels <- list(all=selTerms, unfocus=selTermsUnfocus, focus=selTermsFocus)
for(colSel in names(colSels)) {
	selC <- colSels[[colSel]]
	pheatmap((hmData[selR,selC]), labels_row=selR, border_color="white", cellheight=10, cellwidth=10, cluster_rows=F, cluster_cols=T, breaks=seq(0,ub,length.out=ub*3), col=colorRampPalette(c(brewer.pal(5,"PuRd")))(ub*3-1), treeheight_column=10, treeheight_row=10,  number_format="%s", number_color="black", display_numbers=(hmDataP[selR,selC]), main=colSel)
}
dev.off()


# get overlaps of enriched LOLA peak lists with ATAC-seq peaks:

selFilesTmp <- lolaResAll[gsub("cluster_(.+)_n_.+$","\\1",userSet)%in%selR & sig & term%in%selTermsFocus,unique(filename)]

if(!any(class(peaks)=="GRanges")) peaks <- lib$dtToGr(peaks)
xUserSets <- GRangesList(all=peaks)
simpleCache(paste_("lola_overlaps",digest(selFilesTmp)), {
	enrichOverlaps <- rblapply(selFilesTmp, function(f) {
		msg(f)
		x <- unique(lolaResAll[filename==f, .(collection,filename)])
		tmp <- rblapply(names(xUserSets), function(uSet) {
			msg("\t",uSet)
			lib$grToDt(extractEnrichmentOverlaps(data.table(x,userSet=uSet), xUserSets, regionDB))
		}, "uSet")
		cbind(tmp, collection=x$collection)
	}, "f")
	enrichOverlaps
}, assignToVar="enrichOverlaps")
enrichOverlapsM <- merge(unique(lolaResAll[,.(collection, filename, term)]),unique(merge(enrichOverlaps, peaksDt, by=c("chrom","start","end"))[,.(rid,uSet,f,collection,clusterId)]), by.x=c("collection","filename"), by.y=c("collection","f"))
enrichOverlapsM <- merge(enrichOverlapsM, as.data.table(melt(d)), by.x="rid", by.y="Var1",allow.cartesian=T, all.x=T, all.y=F)
enrichOverlapsM[,sampleGroup:=dA[as.character(Var2),sample_group]]
enrichOverlapsW1<- lib$dtToDf(dcast(enrichOverlapsM, term~Var2, value.var="value", fun.aggregate=mean))#[!is.na(clusterId),]
enrichOverlapsW2 <- lib$dtToDf(dcast(enrichOverlapsM[!is.na(clusterId),], term~Var2, value.var="value", fun.aggregate=mean))#

tmp <- (lib$dtToDf(dcast(enrichOverlapsM, rid~term, value.var="value"))>0)
write.table(tmp, file=lib$plotDir("chipseq_peak_overlaps.csv"), sep=",", col.names=NA, row.names=T, quote=F) 


# summarize ATAC-seq signals in enriched peak sets:

cDirSub <- resultsDir("rcache_atac_sigs")
dir.create(cDirSub, showWarnings=F)

meanAcc <- sapply(selFilesTmp, function(f) {	
	simpleCache(paste_("counts",gsub(".bed","",f)), {
		r <- getRegionSet(regionDB, f, collections=selCollections)
		counts <- summarizeReads(r, dA$sample_name, dA$reads_file, useDuplicates=useDuplicates)
		counts
	}, assignToVar="counts", reload=T)
	
	counts
})
selC <- colSels$focus

lib$pdfPlot(paste_("lola","PHM","tfs", "acc"), 12, 5)
hm <- pheatmap(t(enrichOverlapsW2)[atacOrder,selC], cluster_rows=F,  scale="none",border_color="white", cellheight=6, cellwidth=10,  cluster_cols=T, col=colorRampPalette(c("white",brewer.pal(9,"YlOrRd"),"black"))(20), treeheight_column=0, treeheight_row=10,  number_color="black")	
pheatmap((hmData[moduleOrder,selC]), cluster_rows=F, cluster_cols=hm$tree_col, border_color="white", cellheight=8, cellwidth=10,  breaks=seq(0,ub,length.out=ub*3), col=colorRampPalette(c(brewer.pal(5,"PuRd")))(ub*3-1), treeheight_column=10, treeheight_row=10,  number_format="%s", number_color="black", display_numbers=(hmDataP[moduleOrder,selC]))
for(patId in names(clusterExpr)) {	
	hmExpr <- clusterExpr[[patId]][selC,]
	hmExpr[hmExpr>2] <- 2
	pheatmap(t(hmExpr)[cOrderSel,selC], main=sprintf("scRNA-seq: %s", patId), cluster_rows=F, cluster_cols=hm$tree_col, border_color="white", cellheight=8, cellwidth=10, col=colorRampPalette(c("white",brewer.pal(5,"YlGnBu")[2:5],"#000033"))(32), treeheight_column=10, treeheight_row=10) 
}
dev.off()
