#!/usr/bin/env $CODEBASE/Rscript
#
# Load in legacy microarray data.
#
# run("microarray", "load")


dt <- fread(lib$dataDir("array161108_ex_vivo_sorted_subsets.csv"))

dWide <- as.matrix(dt[,-(1:3)])
rownames(dWide) <- dt[,syms]
