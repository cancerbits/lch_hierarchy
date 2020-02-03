#!/usr/bin/env $CODEBASE/Rscript
# Integrative analysis combining data types.
#
# run("integrative")

doAllPlots <- F

run("common")

# make sure some variables from the ATAC-seq and single-cell analysis are in memory:
run("atac","init")
run("single_cell", c("init","load"))

# run the integrative part:
run("integrative", "init")
run("integrative", "external")
run("integrative", "cell_type_signatures")
run("integrative", "cor_expression")
run("integrative", "networks")

