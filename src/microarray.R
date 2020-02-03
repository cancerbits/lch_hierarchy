#!/usr/bin/env $CODEBASE/Rscript
# Re-analysis of old microarray data (from Hutter et al. 2012).
# run("microarray")

run("microarray", "init")
run("microarray", "load")
run("microarray", "differential")

