#!/usr/bin/env $CODEBASE/Rscript
# Master script executing the entire secondary analysis.
# 
# run()

# intialize general parameters:
run("common") 

# re-analyze legacy microarray data:
run("microarray")

# analysis of 10x single-cell RNA-seq data:
run("single_cell")

# analysis of ATAC-seq data:
run("atac")

# integrative analysis combining data types:
run("integrative")