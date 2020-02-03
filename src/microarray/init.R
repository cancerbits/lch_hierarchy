#!/usr/bin/env $CODEBASE/Rscript
#
# Initialize microarray analysis (not much to do here).
#
# run("microarray", "init")


if(lib$projectName!="lch") run("common")
setCurrentAnalysis(paste0(cfg$ver,"/microarray"))
loadLibraries(c("limma"))


