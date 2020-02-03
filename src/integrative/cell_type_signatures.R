#!/usr/bin/env $CODEBASE/Rscript
#
# Interrogate microarray-based cell type signatures in LCH single-cell data.
#
# run("integrative", "cell_type_signatures")

setCurrentAnalysis(paste0(cfg$ver,"/integrative"))

cellEntropies <- fread(baseResultsDir(cfg$ver, "/10x/traj_cca_",nCCS,"/cell_entropies.csv"))
cellExpression <- fread(baseResultsDir(cfg$ver, "/10x/traj_cca_",nCCS,"/cell_expr.csv"))
markerLists <- fread(baseResultsDir(cfg$ver, "/microarray/marker_lists.csv"))
o <- fread(baseResultsDir(cfg$ver, "/10x/traj_cca_",nCCS,"/clusters.csv"))$cluster
	
lt <- "allDEGs"     
markers <- markerLists[listType==lt,]

m <- merge(cellExpression[,.(cell=Var2,gene=Var1,e=value)], markers[,.(gene,geneSet)], by="gene", all.x=F, allow.cartesian=T)
m <- merge(m, cellEntropies[,.(rn,entropy,pat_id,group,cluster_id)], by.x="cell", by.y="rn")

pData <- m[,.(e=median(e), entropy=mean(entropy)), by=.(cell,group,cluster_id,pat_id,geneSet)]
pData[, cluster_id:=factor(cluster_id, levels=o)]

p <- ggplot(pData[grepl("up|_1",geneSet),], aes(x=cluster_id, y=e, fill=cluster_id)) + geom_boxplot() + facet_wrap(~geneSet, scales="free") + defTheme(flipX=T) + scale_fill_manual(values=colorPalettes$cluster) + xlab(NULL) + ylab("Mean norm. UMI count")
gg(p, paste_("celltype_box",lt), 10, 6, type="pdf", addBase=T)
