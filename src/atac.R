#!/usr/bin/env $CODEBASE/Rscript
# Analysis of ATAC-seq data (chromatin accessibility).
# run("atac")

run("atac", "init")
run("atac", "load")
run("atac", "differential")
run("atac", "enrichment")
run("atac", "motifs")