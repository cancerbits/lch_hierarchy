# Epigenomics and Single-cell Sequencing Define a Developmental Hierarchy in Langerhans Cell Histiocytosis


# Abstract

Langerhans cell histiocytosis (LCH) is a rare neoplasm predominantly affecting children. It occupies a hybrid position between cancers and inflammatory diseases, and it provides an attractive model for studying cancer development. To explore the molecular mechanisms underlying the pathophysiology of LCH and its characteristic clinical heterogeneity, we investigated the transcriptomic and epigenomic diversity in primary LCH lesions. Using single-cell RNA sequencing, we identified multiple recurrent types of LCH cells within these biopsies, including putative LCH progenitor cells and several subsets of differentiated LCH cells. We confirmed the presence of proliferative LCH cells in all analysed biopsies using immunohistochemistry, and we defined an epigenomic and gene regulatory basis of the different LCH cell subsets by chromatin accessibility profiling. In summary, our single-cell analysis of LCH un-covered an unexpected degree of cellular, transcriptomic, and epigenomic heterogeneity among LCH cells, indicative of complex developmental hierarchies in LCH lesions.

This study sketches a molecular portrait of LCH lesions by combining single-cell transcriptomics with epigenome profiling. We uncovered extensive cellular heterogeneity, explained in part by an intrinsic developmental hierarchy of LCH cells. Our findings provide new insights and hypotheses for advancing LCH research and a starting point for personalising therapy.

# Website and Links

http://doi.org/10.1158/2159-8290.CD-19-0138

http://lch-hierarchy.computational-epigenetics.org

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133706

https://ega-archive.org/studies/EGAS00001003822

# Repository Structure

* `src` -- source code for pipelines and analysis (plus execution scripts in main folder). The code has been split into different components (`src/single_cell`, etc.) and there is a master script for most steps (`src/single_cell.R`, etc.).
* `metadata` -- sample annotation, additional inputs, and configurable parameters (e.g. colors)

Execution scripts:
* `start_pipeline.sh` -- set off basic processing pipeline for ATAC-seq data (using pyPiper/looper framework). These have been executed at CeMM via `slurm`.
* `run_all.sh` -- execute the complete scRNA-seq, ATAC-seq, and integrative analysis (`src/main.R`), plus re-analysis of legacy microarray data. This has been run at the CCRI via `nohup`.

Additional code archival (EGA, GEO, website) is available in `src/archival/`.

I used R 3.5.2 and Seurat 2.3.4 -- with a lot of other libraries.

```R
other attached packages:
 [1] corpcor_1.6.9               LandSCENT_0.99.2
 [3] SingleCellExperiment_1.4.1  SummarizedExperiment_1.12.0
 [5] GenomicRanges_1.34.0        GenomeInfoDb_1.18.2
 [7] Seurat_2.3.4                cowplot_0.9.4
 [9] org.Mm.eg.db_3.7.0          org.Hs.eg.db_3.7.0
[11] AnnotationDbi_1.44.0        DelayedMatrixStats_1.4.0
[13] DelayedArray_0.8.0          BiocParallel_1.16.6
[15] IRanges_2.16.0              S4Vectors_0.20.1
[17] matrixStats_0.54.0          monocle_2.10.1
[19] DDRTree_0.1.5               irlba_2.3.3
[21] VGAM_1.1-1                  Biobase_2.42.0
[23] BiocGenerics_0.28.0         Matrix_1.2-15
[25] tsne_0.1-3                  Rtsne_0.15
[27] doMC_1.3.5                  iterators_1.0.10
[29] foreach_1.4.4               fpc_2.1-11.1
[31] Cairo_1.5-10                viridis_0.5.1
[33] viridisLite_0.3.0           usethis_1.4.0
[35] devtools_2.0.1              BiocManager_1.30.4
[37] pheatmap_1.0.12             pryr_0.1.4
[39] dplyr_0.8.0.1               plyr_1.8.4
[41] ggrepel_0.8.0               tidyr_0.8.3
[43] digest_0.6.18               simpleCache_0.4.1
[45] fastcluster_1.1.25          data.table_1.12.0
[47] ggplot2_3.2.0               RColorBrewer_1.1-2

loaded via a namespace (and not attached):
  [1] prabclus_2.2-7         R.methodsS3_1.7.1      acepack_1.4.1
  [4] bit64_0.9-7            knitr_1.21             R.utils_2.8.0
  [7] rpart_4.1-13           RCurl_1.95-4.11        metap_1.1
 [10] snow_0.4-3             callr_3.1.1            RSQLite_2.1.1
 [13] RANN_2.6.1             combinat_0.0-8         proxy_0.4-22
 [16] bit_1.1-14             assertthat_0.2.0       isva_1.9
 [19] xfun_0.5               DEoptimR_1.0-8         caTools_1.17.1.1
 [22] igraph_1.2.4           DBI_1.0.0              htmlwidgets_1.3
 [25] sparsesvd_0.1-4        purrr_0.3.1            backports_1.1.3
 [28] trimcluster_0.1-2.1    gbRd_0.4-11            remotes_2.0.2
 [31] ROCR_1.0-7             withr_2.1.2            robustbase_0.93-3
 [34] checkmate_1.9.1        prettyunits_1.0.2      mclust_5.4.2
 [37] cluster_2.0.7-1        ape_5.2                segmented_0.5-3.0
 [40] lazyeval_0.2.1         crayon_1.3.4           hdf5r_1.2.0
 [43] edgeR_3.24.3           pkgconfig_2.0.2        slam_0.1-45
 [46] labeling_0.3           nlme_3.1-137           vipor_0.4.5
 [49] pkgload_1.0.2          nnet_7.3-12            rlang_0.4.0
 [52] diptest_0.75-7         dbscan_1.1-3           doSNOW_1.0.16
 [55] rprojroot_1.3-2        lmtest_0.9-36          Rhdf5lib_1.4.2
 [58] zoo_1.8-4              base64enc_0.1-3        beeswarm_0.2.3
 [61] ggridges_0.5.1         processx_3.2.1         png_0.1-7
 [64] bitops_1.0-6           R.oo_1.22.0            KernSmooth_2.23-15
 [67] blob_1.1.1             lars_1.2               qvalue_2.14.1
 [70] stringr_1.4.0          scales_1.0.0           memoise_1.1.0
 [73] magrittr_1.5           ica_1.0-2              gplots_3.0.1.1
 [76] bibtex_0.4.2           gdata_2.18.0           zlibbioc_1.28.0
 [79] compiler_3.5.2         HSMMSingleCell_1.2.0   lsei_1.2-0
 [82] clue_0.3-57            fitdistrplus_1.0-14    JADE_2.0-1
 [85] cli_1.0.1              dtw_1.20-1             XVector_0.22.0
 [88] pbapply_1.4-0          ps_1.3.0               htmlTable_1.13.1
 [91] Formula_1.2-3          MASS_7.3-51.1          tidyselect_0.2.5
 [94] stringi_1.3.1          densityClust_0.3       locfit_1.5-9.1
 [97] latticeExtra_0.6-28    grid_3.5.2             tools_3.5.2
[100] rstudioapi_0.9.0       foreign_0.8-71         gridExtra_2.3
[103] FNN_1.1.3              qlcMatrix_0.9.7        Rcpp_1.0.0
[106] SDMTools_1.1-221       httr_1.4.0             npsurv_0.4-0
[109] kernlab_0.9-27         Rdpack_0.10-1          colorspace_1.4-0
[112] fs_1.2.6               reticulate_1.11        scater_1.10.1
[115] flexmix_2.3-15         sessioninfo_1.1.1      jsonlite_1.6
[118] marray_1.58.0          dynamicTreeCut_1.63-1  modeltools_0.2-22
[121] R6_2.4.0               Hmisc_4.2-0            pillar_1.3.1
[124] htmltools_0.3.6        glue_1.3.0             BiocNeighbors_1.0.0
[127] class_7.3-15           codetools_0.2-16       pkgbuild_1.0.2
[130] mvtnorm_1.0-9          lattice_0.20-38        tibble_2.0.1
[133] mixtools_1.1.0         ggbeeswarm_0.6.0       gtools_3.8.1
[136] misc3d_0.8-4           survival_2.43-3        limma_3.38.3
[139] docopt_0.6.1           desc_1.2.0             fastICA_1.2-1
[142] munsell_0.5.0          rhdf5_2.26.2           GenomeInfoDbData_1.2.0
[145] plot3D_1.1.1           HDF5Array_1.10.1       reshape2_1.4.3
[148] gtable_0.2.0
```
