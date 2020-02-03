#!/usr/bin/env $CODEBASE/Rscript
#
# Calculate and display the single-cell transcriptional entropy of cells as a measure for
# relative differentiated-ness (transcriptional promiscuity being a hallmark of stem cells). 
# run("single_cell", "entropy")

#install_github("ChenWeiyan/LandSCENT")
loadLibrary("LandSCENT") 
loadLibrary("corpcor")
loadLibraries(c("AnnotationDbi","org.Hs.eg.db")) # for ID mapping



#############################################################################################
# replace some LandSCENT functions to make it possible to work with our pre-defined clusters:
# (only minor changes; original code by Weiyan Chen & Andrew E. Teschendorff: https://github.com/ChenWeiyan/LandSCENT)
#############################################################################################
CompSRanaPRL <- function (idx, exp.m, adj.m, local = TRUE, maxSR = NULL) {
	if(idx%%100==1) cat("\tsample ", idx)
	exp.v <- exp.m[, idx]
	sumexp.v <- as.vector(adj.m %*% matrix(exp.v, ncol = 1))
	invP.v <- exp.v * sumexp.v
	nf <- sum(invP.v)
	invP.v <- invP.v/nf
	p.m <- t(t(adj.m) * exp.v)/sumexp.v
	S.v <- apply(p.m, 1, LandSCENT:::CompS)
	SR <- sum(invP.v * S.v)
	if (is.null(maxSR) == FALSE) {
		SR <- SR/maxSR
	}
	if (local) {
		NS.v <- apply(p.m, 1, LandSCENT:::CompNS)
	}
	else {
		NS.v <- NULL
	}
	return(list(sr = SR, inv = invP.v, s = S.v, ns = NS.v))
}
CompSRana <- function (Integration.l, local = FALSE, mc.cores = 1) {
    Integration.l <-  LandSCENT:::CompMaxSR(Integration.l)
    maxSR <- Integration.l$maxSR
    idx.l <- as.list(seq_len(ncol(Integration.l$expMC)))
	out.l <- lapply(1:length(idx.l), function(x) NULL)
	j <- 1
	# The SR calculation with mclapply returns NULL's sometimes when one parallel thread fails 
	# (I don't know why, but suspect it's related to the OS killing processes when it's overloaded),
	# so we just call the process repeatedly until all values have been filled:
    while(sum(sapply(out.l, is.null))>0) {
		i <- sapply(out.l, is.null)
		msgF("iteration %d, SRs left to calculate = %d", j, sum(i))
		out.l[i] <- mclapply(idx.l[i], LandSCENT:::CompSRanaPRL, exp.m = Integration.l$expMC, adj.m = Integration.l$adjMC, local = local, maxSR = maxSR, mc.cores = max(1,mc.cores-j*2))
		j <- j+1
	}
	
    SR.v <- sapply(out.l, function(v) return(v[[1]]))
    invP.v <- sapply(out.l, function(v) return(v[[2]]))
    S.v <- sapply(out.l, function(v) return(v[[3]]))
    NS.v <- sapply(out.l, function(v) return(v[[4]]))
    Integration.l$SR <- SR.v
    Integration.l$inv <- invP.v
    Integration.l$s <- S.v
    Integration.l$ns <- NS.v
    if (!is.null(Integration.l$data.sce)) {
        colData(Integration.l$data.sce)$SR <- SR.v
    }
    else if (!is.null(Integration.l$data.cds)) {
        pData(Integration.l$data.cds)$SR <- SR.v
    }
    return(Integration.l)
}
InferLandmark <- function(Integration.l, pheno.v = NULL, pctG = 0.01, reduceMethod = c("PCA", "tSNE","manual"), red = null, clusterMethod = c("PAM", "dbscan","manual"), clust = null, k_pam = 9, eps_dbscan = 10, minPts_dbscan = 5, pctLM = 0.05, pcorTH = 0.1) {
    reduceMethod <- match.arg(reduceMethod)
    clusterMethod <- match.arg(clusterMethod)
    
    if (!is.null(Integration.l$data.sce)) {
        if (is.null(pheno.v)) {
            pheno.v <- 
                SingleCellExperiment::colData(Integration.l$data.sce)$phenoInfo
        }
        if (is.null(pheno.v)) {
            warning("No phenotype information, make sure it was stored as name of phenoInfo!")
        }
    }else if (!is.null(Integration.l$data.cds)) {
        if (is.null(pheno.v)) {
            pheno.v <- 
                Biobase::pData(Integration.l$data.cds)$phenoInfo
        }
        if (is.null(pheno.v)) {
            warning("No phenotype information, make sure it was stored as name of phenoInfo!")
        }
    }
    exp.m <- Integration.l$expMC
    sr.v <- Integration.l$SR
    ordpotS.v <- Integration.l$potencyState
    nPS <- length(levels(as.factor(ordpotS.v)))
    
    ### set an integer for gene selection
    ntop <- floor(pctG*nrow(exp.m))
    
    ### now cluster cells independently of SR
    ### select genes over which to cluster
    print("Using RMT to estimate number of significant components of variation in scRNA-Seq data");
    tmp.m <- exp.m - rowMeans(exp.m);
    rmt.o <- isva::EstDimRMT(tmp.m, plot = FALSE);
    svd.o <- svd(tmp.m);
    tmpG2.v <- vector();
    print(paste("Number of significant components=",rmt.o$dim,sep=""));
    for(cp in seq_len(rmt.o$dim)){
        tmp.s <- sort(abs(svd.o$u[,cp]),decreasing=TRUE,index.return=TRUE);
        tmpG2.v <- union(tmpG2.v,rownames(exp.m)[tmp.s$ix[seq_len(ntop)]]);
    }
    selGcl.v <- tmpG2.v;
    
    ### do dimension reduction
    if (reduceMethod == "tSNE") {
        print("Do dimension reduction via tSNE")
        irlba_res <- irlba::prcomp_irlba(t(exp.m), n = rmt.o$dim
                                         , center = TRUE)
        irlba_pca_res <- irlba_res$x
        topDim_pca <- irlba_pca_res
        tsne_res <- Rtsne::Rtsne(as.matrix(topDim_pca), dims = 2, 
                                 pca = FALSE)
        coordinates <- tsne_res$Y[, 1:2]
        Integration.l$coordinates <- coordinates
    } else if(reduceMethod == "PCA") {
        print("Do dimension reduction via PCA")
        coordinates <- svd.o$v[, 1:2]
        Integration.l$coordinates <- coordinates
    } else {
		coordinates <- red
        Integration.l$coordinates <- coordinates
	}
    
	map.idx <- match(selGcl.v,rownames(exp.m));
	
    ### now perform clustering of all cells over the selected genes
    if (clusterMethod == "dbscan") {
        print("Identifying co-expression clusters via dbscan")
        dbsc.o <- dbscan::dbscan(coordinates, eps = eps_dbscan, minPts = minPts_dbscan)
        clust.idx <- dbsc.o$cluster
        names(clust.idx) <- colnames(Integration.l$expMC)
        k.opt <- length(unique(as.factor(dbsc.o$cluster)))
        print(paste("Inferred ",k.opt," clusters",sep=""))
        psclID.v <- paste("PS",ordpotS.v,"-CL",clust.idx,sep="")
    } else if (clusterMethod == "PAM") {
        print("Identifying co-expression clusters via PAM");
        
        distP.o <- as.dist( 0.5*(1-cor(exp.m[map.idx,])) );
        asw.v <- vector();
        for(k in 2:k_pam){
            pam.o <- pam(distP.o,k,stand=FALSE);
            asw.v[k-1] <- pam.o$silinfo$avg.width
        }
        k.opt <- which.max(asw.v)+1;
        pam.o <- pam(distP.o,k=k.opt,stand=FALSE)
        clust.idx <- pam.o$cluster
        print(paste("Inferred ",k.opt," clusters",sep=""))
        psclID.v <- paste("PS",ordpotS.v,"-CL",clust.idx,sep="")
    } else {
		clust.idx <- as.numeric(as.factor(clust))
		k.opt <- length(unique(clust.idx))
        psclID.v <- paste("PS",ordpotS.v,"-CL",clust.idx,sep="")
	}
    
    ### identify landmark clusters
    print("Now identifying landmarks (potency co-expression clusters)")
    distPSCL.m <- table(paste("CL",clust.idx,sep=""),
                        paste("PS",ordpotS.v,sep=""))
    sizePSCL.v <- as.vector(distPSCL.m)
    namePSCL.v <- vector()
    ci <- 1
    for(ps in seq_len(nPS)){
        for(cl in seq_len(k.opt)){
            namePSCL.v[ci] <- paste("PS",ps,"-CL",cl,sep="")
            ci <- ci+1
        }
    }
    names(sizePSCL.v) <- namePSCL.v
    ldmkCL.idx <- which(sizePSCL.v > pctLM*ncol(exp.m))
    print(paste("Identified ",length(ldmkCL.idx)," Landmarks",sep=""))
    
    ### distribution of phenotypes among LMs
    if(!is.null(pheno.v)){
        tab.m <- table(pheno.v,psclID.v)
        tmp.idx <- match(names(sizePSCL.v)[ldmkCL.idx],colnames(tab.m))
        distPHLM.m <- tab.m[,tmp.idx];
    }
    else {
        distPHLM.m <- NULL;
    }
    
    ### medoids
    print("Constructing expression medoids of landmarks");
    med.m <- matrix(0,nrow=length(selGcl.v),ncol=nPS*k.opt);
    srPSCL.v <- vector();
    ci <- 1;
    for(ps in seq_len(nPS)){
        for(cl in seq_len(k.opt)){
            tmpS.idx <- 
                intersect(which(ordpotS.v==ps),which(clust.idx==cl));
            med.m[,ci] <- apply(matrix(exp.m[map.idx,tmpS.idx],
                                       nrow=length(map.idx)),1,median);
            srPSCL.v[ci] <- mean(sr.v[tmpS.idx]);
            ci <- ci+1;
        }
    }
    names(srPSCL.v) <- namePSCL.v;
    srLM.v <- srPSCL.v[ldmkCL.idx];
    medLM.m <- med.m[,ldmkCL.idx];
    colnames(medLM.m) <- namePSCL.v[ldmkCL.idx];
    rownames(medLM.m) <- selGcl.v;
    
    ### now project each cell onto two nearest landmarks
    print("Inferring dependencies/trajectories/transitions between landmarks");
    cellLM2.v <- vector(); cellLM.v <- vector();
    for(c in seq_len(ncol(exp.m))){
        distCellLM.v <- 0.5*(1-as.vector(cor(exp.m[map.idx,c],medLM.m)));
        tmp.s <- sort(distCellLM.v,decreasing=FALSE,index.return=TRUE);
        cellLM2.v[c] <- paste("LM",tmp.s$ix[1],"-LM",tmp.s$ix[2],sep="");
        cellLM.v[c] <- colnames(medLM.m)[tmp.s$ix[1]];
    }
    adjLM.m <- matrix(0,nrow=ncol(medLM.m),ncol=ncol(medLM.m));
    rownames(adjLM.m) <- colnames(medLM.m);
    colnames(adjLM.m) <- colnames(medLM.m);
    for(lm1 in seq_len(ncol(medLM.m))){
        for(lm2 in seq_len(ncol(medLM.m))){
            adjLM.m[lm1,lm2] <- 
                length(which(cellLM2.v==paste("LM",lm1,"-LM",lm2,sep="")));
        }
    }
    sadjLM.m <- adjLM.m + t(adjLM.m);
    corLM.m <- cor(medLM.m);
    pcorLM.m <- cor2pcor(corLM.m);
    rownames(pcorLM.m) <- rownames(corLM.m);
    colnames(pcorLM.m) <- rownames(corLM.m);  
    netLM.m <- pcorLM.m;
    diag(netLM.m) <- 0;
    netLM.m[pcorLM.m < pcorTH] <- 0;
    netLM.m[pcorLM.m > pcorTH] <- 1;    
    
    Integration.l$InferLandmark.l <- list(cl=clust.idx,pscl=psclID.v,
                                          distPSCL=distPSCL.m,medLM=medLM.m,
                                          srPSCL=srPSCL.v,srLM=srLM.v,
                                          distPHLM=distPHLM.m,cellLM=cellLM.v,
                                          cellLM2=cellLM2.v,adj=sadjLM.m,
                                          pcorLM=pcorLM.m,netLM=netLM.m)
    
    return(Integration.l)
}
# /end functions
################




# calculate the entropy of each cell in the combined dataset:
simpleCache(paste_("scrnaseq_pat", "SCENT3", allComboStr), {		
	simpleCache(paste_("scData", "cca", nCCS, allComboStr), assignToVar="scDataMerged")
	
	# load pre-defined PPI and mapped Entrez IDs to gene symbols:
	netName <- "net13Jun12.m"
	data(list=netName)
	net <- get(netName)
	rm(list=netName)
	rownames(net) <- colnames(net) <- mapIds(org.Hs.eg.db, keys = rownames(net), keytype = "ENTREZID", column = "SYMBOL", multiVals = "first")
	
	# harmonize available genes in network and dataset:
	commonIds <- intersect(rownames(net), rownames(scDataMerged@data))
	m <- scDataMerged@scale.data[commonIds,]
	n <- net[commonIds, commonIds]
	
	# need to make sure minimum value is >0 (~0.1 suggested by LandSCENT)
	m <- m - min(m) + 0.1
	
	# integrate network and dataset:
	mn <- DoIntegPPI(exp.m = m, ppiA.m = n)
	
	# calculate signaling entropy:
	sigEnt <- CompSRana(mn, local=F, mc.cores=max(4, detectCores()))
	
	sigEnt$SR <- sapply(sigEnt$SR, function(x) mean(unlist(x)), simplify=T)
	names(sigEnt$SR) <- colnames(sigEnt$expMC)
		
	data.table(cell_id=names(sigEnt$SR), sr=sigEnt$SR, maxSR=sigEnt$maxSR)	
}, assignToVar="sigEntData", noload=F)	
setkey(sigEntData, cell_id)
	
	

subName <- paste_("traj","cca",nCCS)

# for all sample combinations to analyze (e.g., all merged, by disease extent, etc. 
# -- only "across all samples" and "per patient" used in the paper), ...
combosToAnalyze <- combosToAnalyze[sapply(combosToAnalyze,length)%in%c(1,nrow(dA))]
forEach(combosToAnalyze, function(combo) {
		
	# filter down to only LCH-cells:
	comboStr <- paste_(combo,collapse=".")
	comboStrF <- paste_(comboStr,"lch_only")	
	
	dir.create(resultsDir(comboStrF), showWarnings=F)
	
	# filter to LCH cells and combine with entropy values while maintaining existing projection and clustering:
	cacheName <- paste_("scData", "SCENTSLICE3", "cca", nCCS, comboStrF)
	simpleCache(cacheName, {
		msgF("==== %s ====", comboStrF)	
		simpleCache(paste_("scrnaseq_pat", "SCENT3", "cca", "stage1", nCCS, allComboStr), assignToVar="scData", reload=T)

		selClustering <- cfg$sel_clustering
	
		scDataCurCombo <- SubsetData(object = scData, cells.use=sapply(scData@meta.data$group, function(x) x%in%combo))
		
		sigEnt <- sigEntData[colnames(scDataCurCombo@data),]
		scDataCurCombo@meta.data[,"cluster_id"] <- scDataCurCombo@meta.data[,selClustering$name]
		scDataCurCombo@meta.data[,"entropy"] <- sigEntData[colnames(scDataCurCombo@data),]$sr
	
		return(scDataCurCombo)
	}, recreate=F, noload=T)	

	
	# re-cluster with only filtered cells:		(used in supplementary figure)
	cacheName <- paste_("scData", "SCENTSLICE3", "cca", nCCS, comboStrF, "recluster")
	simpleCache(cacheName, {
		msgF("==== %s ====", comboStrF)	
		simpleCache(paste_("scrnaseq_pat", "SCENT3", "cca", "stage1", nCCS, allComboStr), assignToVar="scData", reload=T)

		selClustering <- cfg$sel_clustering
		selClustering$res <- 0.8 # parameter from standard Seurat workflow here
		selClustering$name <- "res.0.8"
	
		scDataCurCombo <- SubsetData(object = scData, cells.use=sapply(scData@meta.data$group, function(x) x%in%combo))		
		scDataCurCombo <- FindVariableGenes(object = scDataCurCombo, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
		scDataCurCombo <- RunPCA(object = scDataCurCombo, plot.SNN=F, save.SNN=T, pc.genes = scDataCurCombo@var.genes, do.print = F)
		scDataCurCombo <- ProjectDim(object = scDataCurCombo)
		scDataCurCombo <- RunTSNE(scDataCurCombo, do.fast=T, check_duplicates=F)
		scDataCurCombo <- FindClusters(scDataCurCombo, reduction.type="tsne", dims.use=1:2, save.SNN=T, resolution=selClustering$res, force.recalc=T)
		scDataCurCombo@meta.data[,"cluster_id"] <- scDataCurCombo@meta.data[,selClustering$name]
		scDataCurCombo@meta.data[,"entropy"] <- sigEntData[colnames(scDataCurCombo@data),]$sr
		
		return(scDataCurCombo)
	}, assignToVar="sc", recreate=F, noload=T)	
})






subName <- paste_("traj","cca",nCCS)
dir.create(resultsDir(subName), showWarnings=F)
selMarkersY <- unique(c("HMMR", "MKI67", "CD300A", "LAMP3", "MMP12", "JDP2", "CCL2", "CLEC9A", "ATP1B3", "CD1A", "CD207", "LAMP1", "ANPEP", "CD33", "CD1E", "CSF1R", "CD86", "CD9", "CCND1", "MAP2K1", "CST3", "CD44", "CD1C"))



# plot projections, entropies, marker genes for analysis on merged dataset across all patients:
patientGroupings <- c("disease") # "pat_id",
thm <- patientGroupings[1]
lvls <- dA[,sort(unique(get(thm)))]
if(length(lvls)>1) stop("sanity check failed: this will only work for a single combo")
nc <- length(lvls) # 4
metaCols <- c("pat_id", "biopsy", "sex", "disease_extent")
nr <- length(selMarkersY) + 3 # feature plots plus title and entropies and metadata
l <- matrix(1:(nr*nc), ncol = nc, nrow = nr, byrow=F)
combos <- lapply(lvls, function(b) dA[get(thm)==b, sample_name])	
combos <- combos[sapply(combos, length)>0]
names(combos) <- sapply(combos, function(combo) paste_(combo,collapse="."))

# retrieve cached data:
d <- getEntropyData(combos)
	
if(length(d)==1) {
	plotList <- sapply(d, makeEntropyPlotPanel, thm=thm, selMarkers=selMarkersY, metaCols=metaCols, simplify=F)
	
	lib$pdfPlot(paste_("traj5", thm), 4*nc, nr*3.8, sub=subName)
	lib$multiplot(plotlist=unlist(lapply(plotList, function(x) x$plotlist), recursive=F), layout=l)
	dev.off()
	
	x <- d[[1]]
	
	x$data@meta.data$cluster_id <- paste0("LCH-C", x$data@meta.data$cluster_id)

	pData <- cbind(as.data.table(x$data@meta.data[,c("entropy","cluster_id")],keep.rownames=T), x$data@dr$tsne@cell.embeddings)
	pData[,c("x","y"):=.(tSNE_1, tSNE_2)]
	pDataCenters <- pData[, .(x=median(x), y=median(y)), by=cluster_id]
	o <- as.character(pData[,.(e=median(entropy)),by=cluster_id][order(-e), cluster_id])
	curColors <- structure(rev(viridis(length(o))), names=o)
	
	if(all.equal(sort(names(colorPalettes$cluster)), sort(names(curColors)))) {
		curColors <- colorPalettes$cluster
	}
	
	if(!file.exists(resultsDir(subName,"/cell_expr.csv"))) {
		fwrite(melt(x$data@scale.data), file=resultsDir(subName,"/cell_expr.csv"))
		fwrite(data.table(cluster=o, col=curColors[o], lib$dt2namedVec(pData[,.(e=median(entropy)),by=.(ID=cluster_id)])[o]), file=resultsDir(subName,"/clusters.csv"))
		fwrite(cbind(pData,x$data@meta.data[pData$rn,c("group","pat_id")]), file=resultsDir(subName,"/cell_entropies.csv"))
	}
	
	pClust <- ggplot(data=pData, aes(x=x, y=y, color=cluster_id)) + labs(x="t-SNE1", y="t-SNE2") + defTheme(noLegendTitle=T, topLegend=T) + geom_point(size=0.5, shape=16) + geom_text_repel(aes(label=cluster_id), color="black", data=pDataCenters)
	pClust <- pClust + scale_color_manual(values=curColors, guide=F) 
	gg(pClust, paste_("traj5", "clust", thm, x$n), 4, 3.8, type="pdf", sub=subName)
	
	pBox <- ggplot(pData, aes(x=factor(cluster_id, levels=o), y=entropy)) + geom_boxplot(aes(fill=cluster_id)) + defTheme(flipX=T) + xlab(NULL) + scale_fill_manual(values=curColors, guide=F)
	gg(pBox, paste_("traj5", "box", thm, x$n), 4, 3, type="pdf", sub=subName)		

	curData <- x$data@scale.data[,pData$rn]
	cors <- apply(curData, 1, function(x) {
		suppressWarnings(cor(x, pData$entropy, method="pearson"))
	})
	cors <- sort(na.omit(cors))
	
	fwrite(as.data.table(melt(cors), keep.rownames="gene")[order(-value),], file=resultsDir(subName,"/cors.csv"))
	 
	selByCor <- c(names(rev(cors))[1:3],"CD1A")
	
	pData2 <- melt(data.table(e=rank(pData$entropy), t(curData[selByCor,])), id.vars="e")
	pData2[, v:=sprintf("%s (r = %.2f)", variable, cors[as.character(variable)])]
	
	p <- ggplot(pData2, aes(x=-e, y=value, fill=v, color=v)) + geom_smooth() + defTheme(topLegend=T, noLegendTitle=T) + scale_color_brewer(type="qual", palette="Dark2") + scale_fill_brewer(type="qual", palette="Dark2")
	gg(p, "cor_to_ent", 5, 2.8, type="pdf", sub=subName)
	
	m <- as.data.table(x$data@meta.data, keep.rownames=T)	
	setkey(m, rn)
	
	m[, cluster_id:=factor(cluster_id, levels=o)]
		
	pData <- melt(m[,.(cluster_id, pat_id, sex, age, biopsy, disease_extent)], id.vars="cluster_id")
	
	p <- ggplot(pData, aes(x=cluster_id, fill=value)) + geom_bar(position="fill")  + defTheme() + facet_wrap(~variable,ncol=1)
	gg(p, "meta_distr", 5, 5, type="pdf", sub=subName)
	
	for(n in pData[,unique(variable)]) {
		p <- ggplot(pData[variable==n,], aes(x=cluster_id, fill=value)) + geom_bar(position="fill")  + defTheme(topLegend=T, flipX=T) + scale_fill_manual(values=colorPalettes[[n]], guide=guide_legend(nrow=1)) + xlab(NULL) + ylab(NULL)
		gg(p, paste_("meta_distr",n), 2.5, 2, type="pdf", sub=subName)
	}
}

# write entropy values to file:
cell2cluster <- lib$dt2namedVec(fread(resultsDir(subName,"/cell_entropies.csv"))[,.(ID=rn, cluster_id)])




# plot projections, entropies, marker genes for analysis re-clustered per patient (Supplementary Fig. S5):

patientGroupings <- c("pat_id")
for(thm in patientGroupings) {
	lvls <- dA[,sort(unique(get(thm)))]
	nc <- length(lvls) # 4
	metaCols <- c("cluster") # "pat_id", "biopsy"
		
	combos <- lapply(lvls, function(b) dA[get(thm)==b, sample_name])	
	combos <- combos[sapply(combos, length)>0]
	names(combos) <- sapply(combos, function(combo) paste_(combo,collapse="."))
	
	d <- getEntropyData(combos, recluster=T)
	d <- sapply(d, function(x) {
		x$data@meta.data[,"cluster"] <- cell2cluster[rownames(x$data@meta.data)]
		x
	}, simplify=F)
		
	if(length(d)>0) {		
		selMarkersZ <- c("CD1A","CD207","HMMR","MKI67","AURKA","CLEC9A","BATF3","LAMP3","MMP12")
	
		plotList <- sapply(d, makeEntropyPlotPanel, renameClusts=T, thm=thm, selMarkers=selMarkersZ, forceMetaColors=colorPalettes$cluster, metaCols=metaCols, simplify=F)
		tryCatch(dev.off(),error=function(e) {})
		
		nr <- 3
		l <- matrix(1:(nr*nc), ncol = nc, nrow = nr, byrow=F)

		lib$pdfPlot(paste_("traj5",thm,"recluster","nogenes"), 4*nc, nr*3.8, sub=subName)
		lib$multiplot(plotlist=unlist(lapply(plotList, function(x) { 
			x$plotlist[1:3] 
		}), recursive=F), layout=l)
		dev.off()
		
		nr <- 1
		l <- matrix(1:(nr*nc), ncol = nc, nrow = nr, byrow=F)

		lib$pdfPlot(paste_("traj5",thm,"recluster","tsne"), 4*nc, nr*3.8, sub=subName)
		lib$multiplot(plotlist=unlist(lapply(plotList, function(x) { 
			x$plotlist[1] 
		}), recursive=F), layout=l)
		dev.off()
		lib$pdfPlot(paste_("traj5",thm,"recluster","box"), 4*nc, nr*1.8, sub=subName)
		lib$multiplot(plotlist=unlist(lapply(plotList, function(x) { 
			x$plotlist[2] 
		}), recursive=F), layout=l)
		dev.off()
		lib$pdfPlot(paste_("traj5",thm,"recluster","meta"), 4*nc, nr*3.8, sub=subName)
		lib$multiplot(plotlist=unlist(lapply(plotList, function(x) { 
			x$plotlist[3] 
		}), recursive=F), layout=l)
		dev.off()
		
		nr <- length(selMarkersZ)
		l <- matrix(1:(nr*nc), ncol = nc, nrow = nr, byrow=F)

		lib$pdfPlot(paste_("traj5",thm,"recluster","genes"), 4*nc, nr*3.8, sub=subName)
		lib$multiplot(plotlist=unlist(lapply(plotList, function(x) { 
			lapply(x$plotlist[-(1:3)], function(p) { p + xlab(NULL) + ylab(NULL) +
			  theme(axis.title=element_blank(),
					axis.text=element_blank(),
					axis.ticks=element_blank())
			})
		}), recursive=F), layout=l)
		dev.off()
	}
}
