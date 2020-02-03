#!/usr/bin/env $CODEBASE/Rscript
#
# Constants and functions for integrative analysis.
#
# run("integrative", "init")

##############################################
### INTEGRATIVE-ANALYSIS-SPECIFIC LIBARIES ###
##############################################

if(lib$projectName!="lch") run("common")
setCurrentAnalysis("")


getEnrichrPathwayGenes <- function(sel, dbCol="database", termCol="category") {
	pathwayGenes <- rblapply(sel[,unique(database)], function(curDb) {
		fread(paste0("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=", curDb), sep="~", header=F)	
	}, dbCol)
	pathwayGenes[, paste(termCol):=gsub("^(.+)\t\t.+$","\\1",V1)]
	
	enrichrGenes <- rblapply( intersect(pathwayGenes[,get(termCol)],sel[,get(termCol)]), function(curTerm) {
		data.table(gene=pathwayGenes[get(termCol)==curTerm, unique(unlist(strsplit(gsub(",\\d+\\.\\d+","",gsub("^.+\t\t(.+)$","\\1",V1)), "\t")))])
	}, termCol)

	enrichrGenes
}

plotNetworkWithExpr <- function(n, node.w, edge.w, edge.show=structure(rep(T,length(edge.w)), names=E(n)$id), l=layout_nicely(n), node.col="grey", lab="network", s=12, label.scale=0.6) {
	require("igraph")

	msg(lab)	

	i <- as.character(E(n)$id)
	ew <- as.numeric(edge.w[i])
	es <- as.numeric(edge.show[i])
	i <- as.character(V(n)$name)
	nw <- pmax(0,as.numeric(node.w[i]))+0.005
	nc <- alpha(V(n)$col,pmin(1,nw/max(nw)))
	data.table(i, nw, nc)


	plot.igraph(n, 
		edge.arrow.size=.4, edge.curved=.1, 
		vertex.frame.color="grey", 
		#vertex.frame.color=NA, 

		vertex.size=nw * 40,
		vertex.shape=V(n)$shape,
		vertex.label.cex=nw * label.scale,
		vertex.label.color="black",
		vertex.color=nc,
		vertex.label=V(n)$lbl,

		edge.lty=ifelse(grepl("i",E(n)$type), 2, ifelse(es==0, 0, 1)),
		edge.width=es*ew,	
		edge.color=alpha(eCols[E(n)$type],pmax(0,pmin(es * (ew/max(ew))^2,1))), #es*sqrt(ew)*0.5

		layout=l*s,
		xlim=c(-s,s), ylim=c(-s,s),
		rescale=F,
		main=lab
	)
}
preparePlotData <- function(expr.data, meta.data, edge.weights, node.names, fun.size, sel.cells, cluster.col="cluster_id", fun.agg.atac=sum, node.size.max=6, min.node.size.quant=0.95, node.colors=nCols, col.name.prefix="", min.atac.thresh=quantile(edge.weights,0.25)) { #
	if(is.data.table(meta.data)) meta.data <- lib$dtToDf(meta.data)
	
	sel.cells <- intersect(sel.cells, colnames(expr.data))
	curExpr <- rowMeans(expr.data[,sel.cells])

	edgeW <- rowMeans(edge.weights[,meta.data[sel.cells,cluster.col]])

	toW <- sapply(split(edgeW,gsub(".+_(.+)_.+","\\1",names(edgeW))),fun.agg.atac)
	fromW <- sapply(split(edgeW,gsub("(.+)_.+_.+","\\1",names(edgeW))),fun.agg.atac)

	i <- node.names%in%names(curExpr)
	nSizeTmp <- fromW
	nSize <- nSizeTmp
	ns <- node.names[i]
	nSize[ns] <- fun.size(nSize[ns], curExpr[ns])
	ns <- node.names[!i]
	nSize[!is.finite(nSize)] <- 0
	
	nSizeNoCap <- nSize
	nSize <- pmin(nSize,node.size.max)

	curNodeCol <- rowMeans(col2rgb(node.colors[paste0(col.name.prefix,meta.data[sel.cells,cluster.col])]))
	curNodeCol <- rgb(curNodeCol[1],curNodeCol[2],curNodeCol[3],maxColorValue=255)

	list(nSize=nSize, nSizeNoCap=nSizeNoCap, edgeW=edgeW, curNodeCol=curNodeCol)
}

####################################################
### INTEGRATIVE-ANALYSIS-SPECIFIC CONFIG OPTIONS ###
####################################################

dA <- loadAnnot("10x")

setCurrentAnalysis(paste0(cfg$ver,"/integrative"))