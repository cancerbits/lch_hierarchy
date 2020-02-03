#!/usr/bin/env $CODEBASE/Rscript
#
# Functional enrichment analysis using Enrichr and GSEA.
#
# run("single_cell", "enrichment")

# install_github("wjawaid/enrichR")
loadLibraries(c("ggrepel"))
subName <- "enrichment"
dir.create(resultsDir(subName), showWarnings=F)




###############################################################################
####             Enrichr: (LCH vs. NON-LCH) and (LCH subsets)             #####
###############################################################################

pT <- cfg$deg_thresh_fdr
fcT <- cfg$deg_thresh_lfc
allGenesets <-  getMarkerLists(pT, fcT)

cTypes <- c("Bcell"="CD19","TNKcell"="CD3D",MacMono="CD14","DC"="IL3RA")

n <- "LCH_vs_cTypes_strat_overlap"
n <- "traj_cluster_id_sel_strat_w7"

for(n in names(allGenesets)) {
		
	msg(n)
	
	fwrite(allGenesets[[n]]$genesets, file=resultsDir(subName, "/genes_", n, ".csv"))
	
	clusterMarkers <- allGenesets[[n]]$genesets
	
	# scatterplots between groups:
	if(!is.null(allGenesets[[n]]$data) & !is.null(allGenesets[[n]]$groupBy)) {
		scData <- allGenesets[[n]]$data
		grpVals <- scData@meta.data[,allGenesets[[n]]$groupBy]
		grps <- as.character(na.omit(setdiff(unique(grpVals), "NA")))
	
		N <- 12
		bg <- allGenesets[[n]]$bg
	
		for(fg in setdiff(grps,bg)) {
			msgF("\t\t%s vs %s", fg, paste(bg,collapse="+"))
		
			if(is.null(bg)) bg <- setdiff(grps,fg)
			
			plotData <- data.table(gene=rownames(scData@data), x=apply(as.matrix(scData@data[,!is.na(grpVals) & grpVals==fg]), 1, ExpMean), y=apply(as.matrix(scData@data[,!is.na(grpVals) & grpVals%in%bg]), 1, ExpMean))
			plotData <- merge(plotData,clusterMarkers[cluster==fg,],by="gene",all=T)
			plotData[is.na(cluster), cluster:="-"]
			plotData[cluster!="-" & avg_logFC>0, cluster:=paste(bg,collapse="+")]
			
			topN <- clusterMarkers[cluster==fg & avg_logFC>0,] %>% group_by(cluster) %>% top_n(N, avg_logFC)
			bottomN <- clusterMarkers[cluster==fg & avg_logFC<0,] %>% group_by(cluster) %>% top_n(N, -avg_logFC)
      
			p <- ggplot(plotData, aes(x=x, y=y)) + stat_density2d(aes(fill=(..density..^0.1)), color=NA, geom = "tile", contour = FALSE) + scale_fill_continuous(low = "white", high = "darkgrey") + geom_point(aes(color=as.factor(cluster)), alpha=0.5, data=plotData[cluster!="-",]) + defTheme(topLegend=TRUE, noLegendTitle=TRUE) + geom_text_repel(aes(label=gene), size=2, color="black", data=plotData[gene%in%c(bottomN$gene,topN$gene),], min.segment.length=0) + xlab(paste("Average UMI count (log):",fg))  + ylab(paste("Average UMI count (log):",paste(bg,collapse="+")))
			p <- p + scale_color_manual(values=colorPalettes$cell_type)
			
			gg(p, paste_("scatter", n, fg), 3, 3.3, type="pdf", sub=subName) #expand0=T, 
		}
	}
	
	# use Enrichr API:
	simpleCache(paste_("clusterMarkers", n, "enrichr", digest(enrichrDBs), allComboStr, fcT, pT), { 
		loadLibrary("enrichR")
		
		enrichrLists <- sapply(split(clusterMarkers, by="cluster"), function(x) x$gene, simplify=F)
		#[avg_logFC>fcT,]  # this was used in the code, but I removed it because filtering is already done in getMarkerLists() and >0 doesn't hold for the reverse marker selection
		
		runEnrichr(enrichrLists, dbs=enrichrDBs)
	}, recreate=F, reload=T, assignToVar="e")
	
	e[, combscore:= Combined.Score]
	e[, pval:= P.value]
	e[, padj:= Adjusted.P.value]

	# plot selected results:
	
	dbSels <- list(
		"few2"=c("ARCHS4_Tissues","Jensen_DISEASES","BioCarta_2016","GO_Biological_Process_2018","KEGG_2016","NCI-Nature_2016","MSigDB_Oncogenic_Signatures"),
		"archs4"=c("ARCHS4_Tissues"),
		"gobp"=c("GO_Biological_Process_2018"),
		"kegg"=c("KEGG_2016"),
		"nci"=c("NCI-Nature_2016")
	)

	for(dbSel in names(dbSels)) {
			
		selEnrichrRes <- e[database%in%dbSels[[dbSel]],]
		selEnrichrRes[,database:=as.factor(database)]

		suppressWarnings(selEnrichrRes[,rnk:=NULL])
		forEach(selEnrichrRes[,unique(grp)], function(curCmp) {
			forEach(selEnrichrRes[,unique(database)], function(db) {
				selEnrichrRes[grp==curCmp & database==db, rnk:=rank(-combscore,ties.method="random")]
			})
		})

		fwrite(selEnrichrRes[rnk<=50,][order(grp,database,rnk),.("Group"=grp, "Database"=database, "Term"=term, "Hits"=n.hits, "Total list size"=n.total, "P-value"=pval, "Corrected p-value"=padj, "Combined score"=combscore, "Full term"=original.term, "Genes in overlap"=Genes)], file=resultsDir(subName,"/enrichr_",n,"_",dbSel,".csv"))
			
		plotData <- selEnrichrRes[ is.finite(rnk) & rnk<=3,]
		mx <- plotData[,ceiling(max(combscore))*1.25]
		p <- ggplot(plotData, aes(x=reorder(substr(paste0(as.numeric(database), "::",grp,":: ",term),0,48),combscore), y=combscore)) + geom_bar(stat="identity", position="dodge", width=0.7) + defTheme(topLegend=T) + coord_flip() + scale_x_discrete(labels=function(x) { gsub("^\\d+::.+::(.+)$","\\1",x) }) + ylab("Enrichr combined score") + xlab(NULL) + facet_wrap(database~grp, scales="free", drop=F, dir="v", ncol=plotData[,length(unique(database))]) + geom_text(aes(label=lib$pToSig(padj),y=combscore*1.05, hjust=0, vjust=0.75)) + scale_y_continuous(expand=c(0.1,0)) #limits=c(0,mx), 
		gg(p, paste_("enrichr_bars", n, dbSel), length(plotData[,unique(database)])*4, length(plotData[,unique(grp)])*1.1, sub=subName, addBase=T, type="pdf")
		p <- ggplot(plotData, aes(x=reorder(substr(paste0(as.numeric(database), "::",grp,":: ",term),0,48),combscore), y=combscore)) + geom_bar(stat="identity", position="dodge", width=0.7) + defTheme(topLegend=T) + coord_flip() + scale_x_discrete(labels=function(x) { gsub("^\\d+::.+::(.+)$","\\1",x) }) + ylab("Enrichr combined score") + xlab(NULL) + facet_wrap(database~grp, scales="free_y", drop=F, dir="v", ncol=plotData[,length(unique(database))]) + geom_text(aes(label=lib$pToSig(padj),y=combscore*1.05, hjust=0, vjust=0.75)) + scale_y_continuous(expand=c(0.1,0)) #limits=c(0,mx), 
		gg(p, paste_("enrichr_bars_fixed", n, dbSel), length(plotData[,unique(database)])*4, length(plotData[,unique(grp)])*1.1, sub=subName, addBase=T, type="pdf")	
	}
}




###############################################################################
####                        GSEA: LCH vs. NON-LCH                         #####
###############################################################################

simpleCache(paste_("scData", "cca", nCCS, allComboStr, "ext"), assignToVar="scData")
m <- scData@meta.data[,"is_lch_cell"]==TRUE
exprDataBase <- as.matrix(scData@data)
	
exprData <- t(apply(exprDataBase,1,scale))
colnames(exprData) <- colnames(scData@data)

exprDataCap <- abscap(exprData, 0.95)

nTerms <- 30

# calculate ranking criteria for GSEA:
simpleCache(paste_("clusterMarkers", "ranks", "gsea", allComboStr), { 
	simpleCache(paste_("scData", "cca", nCCS, allComboStr, "ext"), assignToVar="scData")

	exprDataBase <- as.matrix(scData@data)
		
	dAll <- rblapply(c("LCH",names(cTypes)), function(x) {
		msgF("\t* %s",x)
	
		m <- !is.na(scData@meta.data[,"cell_type"]) & scData@meta.data[,"cell_type"]==x

		# log fold change:
		d <- data.table(gene=rownames(scData@data), x=apply(exprDataBase[,m], 1, ExpMean), y=apply(exprDataBase[,!m], 1, ExpMean))
		d[,fc:=x-y]

		# Wilcoxon test statistic:
		wtmp <- apply(exprDataBase, 1, function(x) {
			wilcox.test(x[m],x[!m])$statistic
		})
		d[,wstat:=wtmp]
		d[,wstatS:=scale(wstat)]
		
		d
	}, "cell_type")
	
	dAll
}, recreate=F, reload=T, assignToVar="dAll")
fwrite(dAll, file=resultsDir(subName,"/all_ranks.csv"))

for(ct in dAll[,unique(cell_type)]) {
	d <- dAll[cell_type==ct,]
	setkey(d,gene)

	for(rnkBy in c("wstatS")) {
		n <- paste_("rnk", ct, rnkBy)
		gseaResDir <- resultsDir(subName,"/gsea_", n)
		
		msgF("==== GSEA %s ====", n)
		
		if(!dir.exists(gseaResDir) | length(list.files(gseaResDir, pattern="GseaPreranked"))==0) {
			rnkFile <- resultsDir(subName,"/", n, ".rnk")
			write.table(d[order(get(rnkBy)),.(gene,get(rnkBy))], file=rnkFile, sep="\t", quote=F, col.names=F, row.names=F)
			cmd <- sprintf("%s -cp %s -Xmx%s xtools.gsea.GseaPreranked -gmx %s -norm %s -nperm 1000 -rnk %s -scoring_scheme %s -rpt_label %s -create_svgs true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report true -out %s -gui false", toolsDir("java8"), toolsDir("gsea-3.0.jar"), "1024m", "gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/h.all.v6.2.symbols.gmt", "meandiv", rnkFile, "weighted", n, gseaResDir)
			system(cmd)
		} else {
			msgF("\tOutput directory exists (%s), skipping GSEA run...", gseaResDir)
		}

		gseaResult <- list.files(gseaResDir, pattern="GseaPreranked")[1]
		
		if(length(gseaResult)!=1) stop(paste("Something unexpected. !=1 directory with GSEA results for", n))
		
		loadLibrary("qusage")
		geneSets <- read.gmt(paste0(gseaResDir,"/",gseaResult,"/edb/gene_sets.gmt"))

		for(enrichDir in c("pos")) {
			selTerms <- fread(list.files(paste0(gseaResDir,"/",gseaResult), full.name=T, pattern=paste0("gsea_report_for_na_",enrichDir,"_.*xls")))	

			msg("\t barplot")
			plotData <- selTerms[1:min(nrow(selTerms),nTerms),.(NAME=sprintf("%44s",lib$capFirst(tolower(gsub("_"," ",gsub("HALLMARK_","",NAME))))),NES,`FDR q-val`)]
			if(nrow(plotData)<nTerms) plotData <- rbind(plotData, data.table(NAME=sapply(1:(nTerms-nrow(plotData)), function(i) paste0(rep(" ",i), collapse="")), NES=0, `FDR q-val`=1))
			plotData[,resultSet:=n]
			plotData[,rnk:=1:nrow(plotData)]
			mx <- plotData[,round(max(NES*1.1),1)]
			p <- ggplot(plotData, aes(x=reorder(NAME,NES), y=NES)) + geom_bar(stat="identity", width=0.7) + defTheme(topLegend=T) + coord_flip() + scale_x_discrete(labels=function(x) { gsub("^\\d+::.+::(.+)$","\\1",x) }) + ylab("NES") + xlab(NULL) + geom_text(aes(label=lib$pToSig(`FDR q-val`),y=NES+0.05), hjust=0, vjust=0.75) + scale_y_continuous(limits=c(0,mx), expand=c(0,0)) + ggtitle(paste(n,enrichDir))
			gg(p, paste_("gsea_bars", enrichDir, n), 7, 5, type="pdf", sub=subName)

			selTerms <- selTerms[order(NES),][`FDR q-val`<=0.05,NAME]

			if(length(selTerms)>0) {					
				lib$pdfPlot(paste_("gsea_gsea", enrichDir, n), 4, 4, sub=subName)
				for(curTerm in selTerms) {
					msg("\t GSEA plot: ", curTerm)
					replotGSEA(paste0(gseaResDir,"/",gseaResult), curTerm, "LCH - non-LCH", rnkMetricTitle=sprintf("Ranking metric (%s)",rnkBy), negCol="#00bebe", posCol="#bc4676")						
				}
				dev.off()	
			}
		}
	}
}
