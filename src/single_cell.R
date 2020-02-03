#!/usr/bin/env $CODEBASE/Rscript
# Analysis of 10x single-cell RNA-seq data (gene expression).
# run("single_cell")

doAllPlots <- T

run("single_cell", "init")
run("single_cell", "load")
run("single_cell", "entropy")
run("single_cell", "markers") 
run("single_cell", "enrichment")

