#!/usr/bin/env $CODEBASE/Rscript
#
# Build regulatory networks from combined ATAC-seq and scRNA-seq data
# (specific to LCH-S1, LCH-S11 and LCH-S12.
#
# run("integrative", "networks")

loadLibrary("igraph")
loadLibrary("visNetwork")
loadLibrary("igraph")
loadLibrary("shiny")
loadLibrary("htmlwidgets")
loadLibrary("knitr")


# some config:
exclGrps <- c("-")
nCols <- colorPalettes$cluster
eCols <- c("-"="#EEEEEE", "i1"="#a2514b", "i2"="#578a68", "i5"="#534463", "i8"="#b5a43f", "a1"="#329959", "a2"="#950444", "a5"="#fcfc33", "a8"="#2b338e")
nSizeMax <- 6
sizeFun <- function(n, ex) {
	x <- pmax(0,ex-exprThresh)
	x^2 * n^2
}
withAnim <- F
#lo <- "layout_with_fr"
lo <- "layout_with_mds"
curClust <- "C5"
selPatients <- "E"


# load ATAC-seq data:
chipOverlap <- fread(lib$resultsDir("../atac/chipseq_peak_overlaps.csv")) # check ok
motifOverlap <- fread(lib$resultsDir("../atac/motif_peak_overlaps.csv")) # check ok
simpleCache("atac_peaks_annot", assignToVar="peaksDt", reload=T)
simpleCache("atac_deseq2", assignToVar="tmp")
peakIntens <- tmp$dNorm
maxIntens <- max(peakIntens)
# mean accessibility per LCH subset:
peaksDt <- merge(merge(cbind(peaksDt,sapply(paste0("C",c(2,5,8)), function(curClust) {
	rowMeans(peakIntens[peaksDt$rid,grep(curClust, colnames(peakIntens))])
})), chipOverlap, by.x="rid", by.y="V1", all.x=T), motifOverlap, by.x="rid", by.y="V1", all.x=T)
# indicate which peaks have matches for TF motifs and/or ChIP-seq data:
for(tf in gsub(".x$", "", grep("\\.x$",colnames(peaksDt),value=T))) { # for TFs that occurred in motifs and LOLA
	peaksDt[, (tf):= get(paste0(tf,".x")) | get(paste0(tf,".y"))]
}

# gene expression data per LCH subset:
clusterExprS <- clusterExpr$E[,cOrderSel]
clusterExprS <- clusterExprS / apply(clusterExprS,1,max)
dExpr <- lib$dtToDf(fread(baseResultsDir(exprDir,"/d_expr_lch.csv")))
exprThresh <- as.numeric(quantile(as.matrix(dExpr),0.5))

# expression metadata:
dMeta <- fread(baseResultsDir(exprDir,"/d_meta_lch.csv"))

setCurrentAnalysis(paste0(cfg$ver,"/networks"))









# stretch out into long table of connections from TFs (`from`) to targets (`geneSymbols`)
allTfs <- setdiff(sort(unique(c(colnames(chipOverlap),colnames(motifOverlap)))),"V1")
simpleCache("net_dt", {
	dt <- melt(peaksDt, measure.vars=allTfs, variable.name="from")[value==T,]
	dt[,value:=NULL]
	dt <- unnest(dt, geneSymbols)
	dt
}, assignToVar="dt")

# get LCH-subset-specific network connections by keeping only nodes that correspond to expressed genes, 
# and peaks which were called peaks in the ATAC-seq data of that subset:
dtNet <- rblapply(cellOrder, function(curClust) {
	exprG <- as.character(unlist(clusterIsExprE[modDict[[curClust]]]))
	as.data.table(dt)[from%in%exprG & geneSymbols%in%exprG & (get(sprintf("LCH_%s_3",curClust))==T | get(sprintf("LCH_%s_4",curClust))==T),.(from,to=geneSymbols,rid,grp=clusterId, w=get(curClust))]
}, "net")
dtNet[is.na(grp),grp:="-"] # set the ID of connections that are NOT differentially accessible (i.e., not part of any of the regulatory modules)




# build network graph:
net <- unique(dtNet[!grp%in%exclGrps ,.(from=as.character(from),to=as.character(to),type=grp,w=1)])
allN <- unique(c(net[,from],net[,to]))
net <- rbind(
	net
)
allNx <- unique(c(net[,from],net[,to])) 
net[,id:=paste_(from,to,type)]
vs <- data.table(id=allNx, lbl=allNx)
n <- graph_from_data_frame(d=as.data.frame(net), vertices=as.data.frame(vs), directed=F) #
n <- simplify(n, remove.multiple = F, remove.loops = T) 

V(n)$shape <- "none"
V(n)$shape[V(n)$name%in%allTfs] <- "circle"
V(n)$col <- "#AAAAAA"
V(n)$col[V(n)$name%in%allTfs] <- "#CCCCCC"
i <- V(n)$name%in%allDegs
V(n)$col[i] <- sapply(V(n)$name[i], function(g) {
	lib$blendColors(nCols[subsetMarkers[gene==g, cluster]])
})

# get edge weights from ATAC-seq strengths in measured data:
edgeWeights <- sapply(dtNet[,unique(net)], function(curClust) {	
	lib$dt2namedVec(dtNet[!grp%in%exclGrps & net==curClust,.(from=as.character(from),to=as.character(to),type=grp,w)][,.(ID=paste_(from,to,type),w)])[E(n)$id]
})
edgeWeights[is.na(edgeWeights)] <- 0
rownames(edgeWeights) <- E(n)$id
# infer others by similarity of expression profiles to one of the measured ones:
edgeWeights <- cbind(edgeWeights, sapply(setdiff(cOrder, colnames(edgeWeights)), function(curClust) {
	i <- as.numeric(round(exp((sapply(colnames(edgeWeights), function(ref) {
		clusts <- modDict[[ref]]
		cor(clusterExpr$E[,curClust], apply(clusterExpr$E[,clusts,drop=F],1,Seurat::ExpMean))
	})-min(cor(clusterExpr$E)))*100)))
	rowMeans(t(t(edgeWeights) * i)/sum(i))
}))
edgeWeights <- edgeWeights/5

atacThresh <- as.numeric(quantile(edgeWeights, 0.8))
fAgg <- function(x) sqrt(sum(2^(x-atacThresh)))

# get layout:
l <- get(lo)(n) 
l <- norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)



msg("Making dynamic network visualizations...")

# http://kateto.net/network-visualization


# reverse y-coordinate in layout to make visNetwork display it in the same way as plotted:
l2 <- l
l2[,2] <- -l2[,2]
	
	
selSubsets <- names(modDict)
nr <- 1
nc <- ceiling(length(selSubsets)/nr)

for(curClust in selSubsets) {
	msg(curClust)
	
	curCells <- dMeta[cluster_id%in%modDict[[curClust]] & pat_id%in%selPatients,rn]
	pData <- preparePlotData(dExpr, dMeta, edgeWeights, allNx, sizeFun, curCells, node.size.max=5, fun.agg.atac=fAgg )
	sel <- names(sort(-pData$nSize)[1:10])

	i <- as.character(E(n)$id)
	ew <- as.numeric(pData$edgeW[i])
	es <- as.numeric(T)
	i <- as.character(V(n)$name)
	nw <- pmax(0,as.numeric(pData$nSize[i])/3)+0.005
	nc <- alpha(V(n)$col,pmin(1,nw/max(nw)*8))

	nn <- toVisNetworkData(n, idToLabel = TRUE)
	nn$nodes <- as.data.frame(cbind(nn$nodes,color.background=nc))
	nn$nodes$shape <- "circle"
	nn$nodes$font <- lapply(nw, function(s) list(size=round(s*40)+4))
	nn$edges <- cbind(nn$edges,
		dashes=grepl("i",E(n)$type),
		color=alpha(eCols[E(n)$type],pmax(0.1,pmin(es * (ew/max(ew))^2,1)))
	)

	i <- sapply(alpha(nCols,1), function(x) x%in% unique(alpha(nn$nodes$color.background,1)))
	lnodes <- data.frame(
		label=as.character(finalDict[names(nCols)[i]]),
		shape = c("ellipse"), 
		color = nCols[i],
		title = "Marker genes:"
	)
	rownames(lnodes) <- finalDict[rownames(lnodes)]

	i <- sort(unique(E(n)$type))
	ledges <- data.frame(
		color = eCols[i],
		label = finalDict[i],
		lty = ifelse(grepl("i",i), 2, 1)
	)

	htmlFile <- resultsDir("visnet_",finalDict[curClust],".html")

	nn$edges$type <- as.character(finalDict[nn$edges$type])
	
	# use visNetwork to write HTML pages:
	visNetwork(nodes=nn$nodes, edges=nn$edges, main=finalDict[[curClust]], footer="<hr style=\"margin-top: 1em; margin-bottom: 1em;\"/><p>Instructions: Use the <b>mouse wheel</b> to zoom in and out. <b>Click and drag on individual nodes</b> in the network to move them. <b>Click and drag on the network pane</b> to shift the view and move the entire network. <b>Select a node</b> (by clicking it or using the dropdown menu) to highlight only this node and those immediately connected to it.</p><p>Node sizes in the network are proportional to the expression levels of the respective genes and to their centrality in the network (\"node importance\"). Transcription factors have been linked to their putative targets by genomic proximity in 2D and 3D. Nodes corresponding to genes that are differentially expressed in individual LCH subsets are coloured according to the respective subset (or gray otherwise). Edge colours indicate regulatory modules.</p></div>") %>% 
		visIgraphLayout(randomSeed=0, physics=F, smooth=F, type="full", layout="layout.norm", layoutMatrix=l2) %>%
		visNodes(borderWidth=0.1, color=list(border="black")) %>% 
		visOptions(highlightNearest = list(enabled = T, labelOnly=F, degree = 1, hover = F, hideColor = "rgba(0,0,0,0)", algorithm="all"), nodesIdSelection = list(enabled=TRUE, values=sort(nn$nodes$lbl), main="Highlight gene:")) %>%
		visLegend(addNodes = lnodes, addEdges = ledges, useGroups = FALSE, stepY=60, ncol = 2) %>% 
		visInteraction(navigationButtons = TRUE)%>%
		visSave(file = htmlFile) 
		
		
	# hack the HTML code to customize the visualization:
	htmlCode <- readChar(htmlFile, file.info(htmlFile)$size)
	htmlLinks <- paste(sprintf("<a href=\"visnet_%s.html\">%s</a>", finalDict[selSubsets], finalDict[selSubsets]), collapse=" | ")
	htmlCodeX <- gsub("(<body[^>]*>)",paste0("\\1","<div style=\"text-align:center;\">Choose LCH subset to visualise regulatory network: ", htmlLinks, "</div><div id=\"htmlwidget_container\"><hr style=\"margin-bottom: 2em;\"/><div class=\"container-fluid\">"), htmlCode)
	htmlCodeX <- gsub("<title>.+</title>",sprintf("<title>%s</title>", finalDict[curClust]), htmlCodeX)
	htmlCodeX <- gsub("</body>","<script type=\"text/javascript\" src=\"https://ajax.googleapis.com/ajax/libs/jquery/1.11.1/jquery.min.js\"></script><script type=\"text/javascript\" src=\"https://maxcdn.bootstrapcdn.com/bootstrap/3.3.0/js/bootstrap.min.js\"></script></body>", htmlCodeX) #<script type=\"text/javascript\"> $(function(){ $(\".dropdown\").each(function(i,e) { var options = $(e).find('option'); var arr = options.map(function(_, o) { return { t: $(o).text(), v: o.value }; }).get(); arr.sort(function(o1, o2) { if(o1.t==\"Highlight gene:\") return -1; if(o2.t==\"Highlight gene:\") return 1; return o1.t > o2.t ? 1 : o1.t < o2.t ? -1 : 0; }); options.each(function(i, o) { o.value = arr[i].v; $(o).text(arr[i].t); }); }); }); </script>
	htmlCodeX <- gsub("<head>","<head><meta http-equiv=\"X-UA-Compatible\" content=\"IE=edge\"><meta name=\"viewport\" content=\"width=device-width, initial-scale=1\"><meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\"><link rel=\"stylesheet\" href=\"https://maxcdn.bootstrapcdn.com/bootstrap/3.3.0/css/bootstrap.min.css\"><link rel=\"stylesheet\" href=\"https://maxcdn.bootstrapcdn.com/bootstrap/3.3.0/css/bootstrap-theme.min.css\"><link href=\"inc/style.css\" rel=\"stylesheet\" type=\"text/css\"><!-- HTML5 shim and Respond.js for IE8 support of HTML5 elements and media queries --><!-- WARNING: Respond.js doesn't work if you view the page via file:// --><!--[if lt IE 9]><script src=\"https://oss.maxcdn.com/html5shiv/3.7.2/html5shiv.min.js\"></script><script src=\"https://oss.maxcdn.com/respond/1.4.2/respond.min.js\"></script><![endif]-->", htmlCodeX)
	htmlCodeX <- gsub("\"idselection\":\\{\"enabled\":true\\}","\"idselection\":{\"enabled\":true,\"useLabels\":true, \"main\":\"Highlight gene:\"}", htmlCodeX)
	htmlCodeX <- gsub("class=\"visNetwork html-widget\" style=\"[^\"]+\"","class=\"visNetwork html-widget\" style=\"width:1100px;height:500px;margin-left:auto;margin-right:auto;\"", htmlCodeX)
	htmlCodeX <- gsub("font-family:Georgia, Times New Roman, Times, serif;","", htmlCodeX)
	
	
	for(x in selSubsets) htmlCodeX <- gsub(paste0("LCH-",x), finalDict[curClust], htmlCodeX)
	
	writeChar(htmlCodeX, htmlFile)
}
