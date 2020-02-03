#!/usr/bin/env $CODEBASE/Rscript
#
# Load ATAC-seq data and auxiliary data required for their analysis into local caches.
#
# run("atac", "load")

# load coordinates of repeats:
simpleCache("repeats", loadUCSCRepeats(), noload=T)

# load pcHi-C data from Javierre et al. (2016):
simpleCache("hic_javierre", {
	loadLibrary("rtracklayer")
	
	msg("Reading in pcHi-C data...")
	f <- "PCHiC_peak_matrix_cutoff5.tsv"
	dBait <- fread(lib$resourceDir("regions/enhancer_promoter_interactions/javierre_2016/hg38/", f, ".cut.bait.bed"))
	setnames(dBait, c("chrom","start","end","info"))
	dOe <- fread(lib$resourceDir("regions/enhancer_promoter_interactions/javierre_2016/hg38/", f, ".cut.oe.bed"))
	setnames(dOe, c("chrom","start","end","info"))
	
	msg("Merging baits and 'other ends'...")
	dM <- merge(dBait, dOe, by="info", all=F, suffixes=c("Bait", "Oe"))
	tmp <- t(dM[,strsplit(`info`,"___")])
	dM <- cbind(dM[,2:7], tmp)
	rm(tmp, dBait, dOe, f)
	gc()
	
	msg("Relabel bait coordinates with annotation used for 10x analysis...")
	fp <- lib$resourceDir("10X_Genomics/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf")
	geneAnnot10x <- import.gff(fp)
	baitGr <- with(dM, GRanges(gsub("^chr","",chromBait), IRanges(startBait, endBait)))
	o <- findOverlaps(baitGr, geneAnnot10x)
	newBaitNames <- unique(data.table(q=queryHits(o), gname=geneAnnot10x$gene_name[subjectHits(o)]))
	newBaitNames <- newBaitNames[,.(gname=paste(gname,collapse=";")),by=q]
	dM[as.numeric(newBaitNames$q),targetName:=newBaitNames$gname]

	fwrite(dM, file=resultsDir("hic_dM.csv"))
	
	dM
}, noload=T)

# define peaks:
simpleCache("atac_peaks", {
	msg("Load auxiliary data...")

	simpleCache("hic_javierre", assignToVar="dM")
	simpleCache("repeats")

	# only use interactions for expressed genes:	
	dM <- dM[targetName%in%allExprGenes,]

	# extend contact regions by 1kb in either direction:
	extension <- 1000
	baitGr <- with(dM, GRanges(chromBait, IRanges(startBait-extension, endBait+extension)))
	oeGr <- with(dM, GRanges(chromOe, IRanges(startOe-extension, endOe+extension)))

	msg("Read and process peaks...")
	
	# read in all peaks called in each sample (and throw out those on the mitochondrial and weird parts of the genome):
	umPeaksDt <- rblapply(dA[,peaks_file], fread, "f")[!grepl("_|M",V1)]	
	umPeaks <- with(umPeaksDt, GRanges(V1, IRanges(V2,V3)))
	
	# write / plot peak summary:
	peakNum <- rblapply(dA[,peaks_file], fread, "f")[,.N,by=f]
	fwrite(peakNum, file=resultsDir("peaks.csv"))
	plotData <- merge(peakNum, dA, by.x="f", by.y="peaks_file")
	p <- ggplot(plotData, aes(x=reorder(sample_name,-N), y=N, fill=gsub("LCH_","",sample_group))) + geom_bar(stat="identity") + defTheme(flipX=T) + scale_fill_manual(values=colorPalettes$transition_state)
	gg(p, "peak_num", 4, 4, type="pdf", expand0=TRUE)

	# merge overlapping peaks:
	peaks <- reduce(umPeaks)

	# remove blacklisted regions (see https://sites.google.com/site/anshulkundaje/projects/blacklists):
	blacklist <- with(fread("wget -O - http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg38-human/hg38.blacklist.bed.gz | zcat"), GRanges(V1,IRanges(V2,V3)))
	overlapBlacklist <- countOverlaps(peaks, blacklist)>0
	overlapRepeats <- countOverlaps(peaks, repeats)>0
	peaks <- peaks[!overlapRepeats & !overlapBlacklist]
	peakLen <- end(peaks)-start(peaks)

	# annotate each merged peak with the matching source peaks:
	tmp <- sapply(dA$sample_name, function(s) {
		(countOverlaps(peaks, umPeaks[umPeaksDt[,f==dA[s,peaks_file]],])>0)*1
	})
	
	# pack all peak info in a big table:
	peaksDt <- data.table(chrom=as.character(seqnames(peaks)), start=start(peaks), end=end(peaks), tmp)
	peaksDt[,rid:=paste0("r", 1:nrow(peaksDt))]

	msg("Annotate peaks...")
	
	# annotate with genes by checking overlaps of peaks with pcHi-C data:
	o1 <- findOverlaps(peaks, oeGr, ignore.strand=T)
	o2 <- findOverlaps(peaks, baitGr, ignore.strand=T)
	geneOverlaps <- rbind(data.table(q=queryHits(o1), n=dM[subjectHits(o1),targetName]), data.table(q=queryHits(o2), n=dM[subjectHits(o2),targetName])) #data.table(q=queryHits(o1), n=dM[subjectHits(o1),targetName]) #
	geneOverlapsAgg <- geneOverlaps[,.(gnames=paste(unique(sort(n)),collapse=";")), by=q]
	peaksDt[geneOverlapsAgg$q, geneSymbols:=geneOverlapsAgg$gnames]
	peaksDt[, overlapByHiC:=!is.na(geneSymbols)]

	tAnnot <- lib$loadTranscriptAnnotation(genomeVer)
	o <- distanceToNearest(peaks, lib$dtToGr(tAnnot, "chrom", "tss", "tss"), ignore.strand=T)
	maxDist <- 50000
	o <- o[o@elementMetadata$distance<=maxDist,]
	peaksDt[queryHits(o), geneSymbols:=paste(geneSymbols, tAnnot[subjectHits(o), geneSymbol], sep=";")]
	peaksDt[queryHits(o), distToNextGene:=o@elementMetadata$distance]
	peaksDt[,geneSymbols:=sapply(geneSymbols, function(g) sort(unique(unlist(strsplit(gsub("^NA;", "", g),";"))) ))]

	peaksDt[,isPromoPeak:=F]
	peaksDt[unique(queryHits(o[o@elementMetadata$distance<=1000,])),isPromoPeak:=T]
	nPromoterPeaks <- peaksDt[,sum(isPromoPeak)]
	msgF("n(promoter)=%d [%.1f%%], n(distal)=%d", nPromoterPeaks, nPromoterPeaks/length(peaks)*100, length(peaks)-nPromoterPeaks)

	geneOverlaps <- rbindlist(list(
		hiC=geneOverlaps, 
		proximity=data.table(q=queryHits(o), n=tAnnot[subjectHits(o), geneSymbol])
	), idcol="linkType")

	write.table(peaksDt[,.(chrom,start,end,sprintf("%s:%s", rid, geneSymbols))], file=resultsDir("peaks.bed"), col.names=F, row.names=F, sep="\t", quote=F)

	lib$pdfPlot("hist_peak_len_filtered", 4, 4)
	hist(peaksDt[,end-start])
	dev.off()

	plotData <- melt(peaksDt, id.vars=c("rid","geneSymbols","chrom","start","end"))[value==1, .N, by=variable]
	plotData <- cbind(plotData, dA[as.character(plotData$variable),])
	p <- ggplot(plotData, aes(x=reorder(sample_name,-N), y=N, fill=gsub("LCH_","",sample_group))) + geom_bar(stat="identity") + defTheme(flipX=T) + scale_fill_manual(values=colorPalettes$transition_state)
	gg(p, "peak_num_filtered_merged", 4, 4, type="pdf", expand0=TRUE)

	list(peaks=peaks, peaksDt=peaksDt, geneOverlaps=geneOverlaps)
}, noload=T, recreate=F)

# count reads in peaks:
simpleCache("atac_counts", {
	simpleCache("atac_peaks", assignToVar="tmp")
	peaks <- tmp$peaks
	
	counts <- summarizeReads(peaks, dA$sample_name, dA$reads_file, useDuplicates=useDuplicates)
	counts[, rpm:=count / lib$dt2namedVec(counts[,sum(count),by=sampleName],"sampleName")[sampleName] * 1000000]
	counts
}, noload=T, recreate=F)
