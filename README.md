Impact of time to deep-freezing on cancer tissue metabolomics - data set
================
Andreas Mock, National Center for Tumor Diseases (NCT) Heidelberg

Getting started
===============

This document requires R version 3.4.2 or higher. For more detailed information about the package requirements see the session information chapter at the end of this document.

Install package dependencies
============================

``` r
packages <- c("MultiAssayExperiment","MetaboDiff")
invisible(lapply(packages, library, character.only = TRUE))
```

Study set and time to deep-freezing
===================================

-   3 cortically located, bulk-resected *IDH1* wild-type glioblastomas
-   time to deep-freezing: 0-180 min
-   two samples per time point

| GBM 1     | GBM 2   | GBM 3     |
|-----------|---------|-----------|
| 0 - 0     | 0 - 0   | 0 - 0     |
| 5 - 5     | 5 - 5   | 5 - 5     |
| 10 - 10   | 10 - 10 | 10 - 10   |
| 30 - 30   | 30 - 30 | 30 - 30   |
| 60 - 60   | 60 - 60 | 60 - 60   |
| 180 - 180 |         | 180 - 180 |

Load data set
=============

The object `q` is a so called MultiAssayExperiment class[1] containing the metabolomic data, as well as all the associated metadata. Please refer to the [MetaboDiff package vignette](https://github.com/andreasmock/MetaboDiff) for more information about the rationale of using a MultiAssayExperiment class, as well as an introduction to basic usage.

``` r
load("quenching.RData")
```

We display the summary information of the object by:

``` r
q
```

    ## A MultiAssayExperiment object of 1 listed
    ##  experiment with a user-defined name and respective class. 
    ##  Containing an ExperimentList class object of length 1: 
    ##  [1] raw: SummarizedExperiment with 597 rows and 34 columns 
    ## Features: 
    ##  experiments() - obtain the ExperimentList instance 
    ##  colData() - the primary/phenotype DataFrame 
    ##  sampleMap() - the sample availability DataFrame 
    ##  `$`, `[`, `[[` - extract colData columns, subset, or experiment 
    ##  *Format() - convert into a long or wide DataFrame 
    ##  assays() - convert ExperimentList to a SimpleList of matrices

As we can see from the summary information, the object contains the raw relative metabolomic measurments as provided by Metabolon Inc. A total of 597 named metabolites were identified accross the 34 samples.

Imputation of missing values
============================

In contrast to other high-throughput technologies, missing values are common in quantitative metabolomic datasets. Imputation is performed by k-nearest neighbor imputation, which could be shown to minimize the effects on the normality and variance of the data as long as the number of missing data does not exceed 40%. The function `knn_impute` adds the slot “impute” to the MultiAssayExperiment object that contains the imputed relative metabolite measurements for all metabolites with raw measurements in more than 60% of cases. We recommend a cutoff of 40% (i.e. 0.4). However the cutoff might be changed according to the discretion of the user.

``` r
(q = knn_impute(q, 0.4))
```

    ## A MultiAssayExperiment object of 2 listed
    ##  experiments with user-defined names and respective classes. 
    ##  Containing an ExperimentList class object of length 2: 
    ##  [1] raw: SummarizedExperiment with 597 rows and 34 columns 
    ##  [2] imputed: SummarizedExperiment with 562 rows and 34 columns 
    ## Features: 
    ##  experiments() - obtain the ExperimentList instance 
    ##  colData() - the primary/phenotype DataFrame 
    ##  sampleMap() - the sample availability DataFrame 
    ##  `$`, `[`, `[[` - extract colData columns, subset, or experiment 
    ##  *Format() - convert into a long or wide DataFrame 
    ##  assays() - convert ExperimentList to a SimpleList of matrices

Normalization
=============

Variance stabilizing normalization (vsn) is used to ensure that the variance remains nearly constant over the measured spectrum.

``` r
(q = normalize_met(q))
```

    ## A MultiAssayExperiment object of 4 listed
    ##  experiments with user-defined names and respective classes. 
    ##  Containing an ExperimentList class object of length 4: 
    ##  [1] raw: SummarizedExperiment with 597 rows and 34 columns 
    ##  [2] imputed: SummarizedExperiment with 562 rows and 34 columns 
    ##  [3] norm: SummarizedExperiment with 597 rows and 34 columns 
    ##  [4] norm_imputed: SummarizedExperiment with 562 rows and 34 columns 
    ## Features: 
    ##  experiments() - obtain the ExperimentList instance 
    ##  colData() - the primary/phenotype DataFrame 
    ##  sampleMap() - the sample availability DataFrame 
    ##  `$`, `[`, `[[` - extract colData columns, subset, or experiment 
    ##  *Format() - convert into a long or wide DataFrame 
    ##  assays() - convert ExperimentList to a SimpleList of matrices

At this point the data processing is completed with the `MultiAssayExperiment` object containing 4 slots:

-   raw - raw relative metabolic measurements as provided by company or core facility
-   imputed - imputed relative metabolic measurements (k-nearest neighbor imputation)
-   norm - normalized relative metabolic measurements (vsn)
-   norm\_imputed - normalized and imputed relative metabolic measurements (vsn)

Exclude GBM3 @ 30 minutes from further analysis
===============================================

The two pieces of tissue GBM3 at 30 minutes represent two samples with a low tumor content, as well as a high fraction of necrosis and apoptotic cells and were hence excluded from further analysis.

``` r
q <- q[,!(q$tumor=="GBM3"&q$time==30)]
```

    ## harmonizing input:
    ##   removing 8 sampleMap rows with 'colname' not in colnames of experiments
    ##   removing 2 colData rownames not in sampleMap 'primary'

Working with the data object
============================

1.  Display sample metadata

``` r
colData(q)
```

    ## DataFrame with 32 rows and 4 columns
    ##      sample_id      time     piece    tumor
    ##       <factor> <integer> <integer> <factor>
    ## S72      1_0_1         0         1     GBM1
    ## S74      1_0_2         0         2     GBM1
    ## S75      1_5_1         5         1     GBM1
    ## S76      1_5_2         5         2     GBM1
    ## S78     1_10_1        10         1     GBM1
    ## ...        ...       ...       ...      ...
    ## S104    3_10_2        10         2     GBM3
    ## S107    3_60_1        60         1     GBM3
    ## S108    3_60_2        60         2     GBM3
    ## S109   3_180_1       180         1     GBM3
    ## S110   3_180_2       180         2     GBM3

1.  Normalized and imputed metabolic measurements

``` r
head(assay(q[["norm_imputed"]]))
```

    ##            S72      S74      S75      S76      S78      S79      S81
    ## 20675 21.91029 21.71348 21.65519 21.51916 21.68176 21.71669 21.68131
    ## 33971 24.80303 24.44446 24.25727 24.45003 24.87667 24.29809 24.72898
    ## 33972 23.35797 22.71748 22.29044 22.76144 22.77624 22.60618 23.25546
    ## 37752 21.40502 20.67796 20.66504 20.70528 20.95955 20.30374 20.37451
    ## 35186 24.20370 23.85083 24.49277 24.15916 24.16978 24.34175 24.43529
    ## 34214 21.31666 21.35067 21.75965 22.01847 22.39883 20.95779 21.84452
    ##            S82      S83      S85      S86      S87      S89      S90
    ## 20675 21.82947 21.65797 21.46906 21.47045 21.42752 24.29755 23.78600
    ## 33971 24.42530 26.08888 25.44599 25.57105 25.18796 23.99209 25.00149
    ## 33972 22.71081 24.46391 23.24237 23.37847 23.41343 22.40170 23.01346
    ## 37752 20.91805 20.20300 20.94806 21.20299 21.48732 18.96592 21.32424
    ## 35186 24.54709 24.05221 24.19894 24.13281 24.69871 24.46906 25.07940
    ## 34214 21.86775 20.95390 22.27382 21.90753 21.67450 21.90613 22.72705
    ##            S91      S92      S93      S94      S95      S96      S97
    ## 20675 24.25896 24.50416 23.87086 24.17165 23.94165 23.88324 23.87995
    ## 33971 23.91583 25.31405 25.18485 23.67085 24.73729 25.53193 24.52403
    ## 33972 22.71535 23.91949 23.39871 21.54296 22.90049 23.75941 23.15263
    ## 37752 20.56781 19.52185 20.63943 20.63474 21.03182 21.00767 21.11288
    ## 35186 24.84381 23.95223 24.65904 24.59431 24.02618 24.41599 24.85080
    ## 34214 21.55558 21.00962 22.87570 21.96245 22.44750 23.53190 22.53308
    ##            S98      S99     S100     S101     S102     S103     S104
    ## 20675 25.00204 23.89552 23.97743 23.83964 23.77682 23.90039 23.77349
    ## 33971 24.46351 24.68053 24.58138 24.40341 25.02742 24.67866 25.48721
    ## 33972 23.33977 23.54983 23.38348 23.08589 23.79771 23.18911 24.18006
    ## 37752 19.79200 19.75678 20.28635 19.85507 19.72430 19.95059 18.35137
    ## 35186 23.50485 24.74536 23.76994 23.90720 23.32903 24.23222 24.68284
    ## 34214 20.12588 22.05694 21.56192 21.55208 21.34292 21.69249 22.69806
    ##           S107     S108     S109     S110
    ## 20675 23.86506 24.10153 24.01145 23.76698
    ## 33971 24.72453 24.32189 25.13932 25.21019
    ## 33972 23.52851 23.28023 23.87641 23.72792
    ## 37752 20.25392 20.12406 20.61898 19.34000
    ## 35186 24.05495 23.92915 23.16864 22.41464
    ## 34214 22.25604 21.69776 21.44941 21.02457

1.  Metabolic annotation for normalized and imputed data

``` r
rowData(q[["norm_imputed"]])
```

    ## DataFrame with 562 rows and 11 columns
    ##     PATHWAY_SORTORDER                               BIOCHEMICAL
    ##             <integer>                               <character>
    ## 1                 515              1,5-anhydroglucitol (1,5-AG)
    ## 2                 636                10-heptadecenoate (17:1n7)
    ## 3                 648                 10-nonadecenoate (19:1n9)
    ## 4                 710                          13-HODE + 9-HODE
    ## 5                 919 1-arachidonoylglycerophosphoethanolamine*
    ## ...               ...                                       ...
    ## 558               839                                  pyruvate
    ## 559               388                                  spermine
    ## 560              3229                                 tartarate
    ## 561              2828                       thiamin diphosphate
    ## 562              2680             uridine 5'-triphosphate (UTP)
    ##              SUPER_PATHWAY
    ##                <character>
    ## 1             Carbohydrate
    ## 2                    Lipid
    ## 3                    Lipid
    ## 4                    Lipid
    ## 5                    Lipid
    ## ...                    ...
    ## 558           Carbohydrate
    ## 559             Amino Acid
    ## 560            Xenobiotics
    ## 561 Cofactors and Vitamins
    ## 562             Nucleotide
    ##                                              SUB_PATHWAY   COMP_ID
    ##                                              <character> <integer>
    ## 1   Glycolysis, Gluconeogenesis, and Pyruvate Metabolism     20675
    ## 2                                  Long Chain Fatty Acid     33971
    ## 3                                  Long Chain Fatty Acid     33972
    ## 4                                Fatty acid, monohydroxy     37752
    ## 5                                              Lysolipid     35186
    ## ...                                                  ...       ...
    ## 558 Glycolysis, Gluconeogenesis, and Pyruvate Metabolism     22250
    ## 559                                 Polyamine Metabolism       603
    ## 560                                 Food Component/Plant     15336
    ## 561                                  Thiamine Metabolism     35670
    ## 562             Pyrimidine Metabolism, Uracil containing     33448
    ##           PLATFORM          RI        MASS               CAS     PUBCHEM
    ##        <character> <character> <character>       <character> <character>
    ## 1            GC/MS      1788.7         217         154-58-5;            
    ## 2        LC/MS Neg        5558       267.3       29743-97-3;     5312435
    ## 3        LC/MS Neg        5775       295.4       73033-09-7;     5312513
    ## 4        LC/MS Neg        5247       295.2                              
    ## 5        LC/MS Neg        5731       500.3                              
    ## ...            ...         ...         ...               ...         ...
    ## 558    LC/MS Polar        2650   175.02481          127-17-3        1060
    ## 559 LC/MS Pos Late         758   203.22303           71-44-3        1103
    ## 560      LC/MS Neg       615.3   149.00916 87-69-4;6106-24-7      444305
    ## 561      LC/MS Neg      1195.1   423.02986          154-87-0        1132
    ## 562      LC/MS Neg         572   482.96126        19817-92-6        6133
    ##            KEGG
    ##     <character>
    ## 1        C07326
    ## 2              
    ## 3              
    ## 4              
    ## 5              
    ## ...         ...
    ## 558      C00022
    ## 559      C00750
    ## 560      C00898
    ## 561      C00068
    ## 562      C00075

Session information
===================

``` r
sessionInfo()
```

    ## R version 3.4.2 (2017-09-28)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: macOS  10.14.2
    ## 
    ## Matrix products: default
    ## BLAS: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8
    ## 
    ## attached base packages:
    ##  [1] grid      parallel  stats4    stats     graphics  grDevices utils    
    ##  [8] datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] MetaboDiff_0.9.2           forcats_0.2.0             
    ##  [3] stringr_1.3.1              dplyr_0.7.6               
    ##  [5] purrr_0.2.5                readr_1.1.1               
    ##  [7] tidyr_0.8.0                tibble_1.4.2              
    ##  [9] ggplot2_3.0.0              tidyverse_1.2.1           
    ## [11] ComplexHeatmap_1.17.1      SummarizedExperiment_1.8.1
    ## [13] DelayedArray_0.4.1         matrixStats_0.54.0        
    ## [15] Biobase_2.38.0             GenomicRanges_1.30.1      
    ## [17] GenomeInfoDb_1.14.0        IRanges_2.12.0            
    ## [19] S4Vectors_0.16.0           BiocGenerics_0.24.0       
    ## [21] MultiAssayExperiment_1.4.9
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] nlme_3.1-131           bitops_1.0-6           lubridate_1.7.2       
    ##  [4] RColorBrewer_1.1-2     httr_1.3.1             rprojroot_1.3-2       
    ##  [7] tools_3.4.2            backports_1.1.2        affyio_1.48.0         
    ## [10] R6_2.3.0               lazyeval_0.2.1         colorspace_1.3-2      
    ## [13] GetoptLong_0.1.6       withr_2.1.2            tidyselect_0.2.4      
    ## [16] mnormt_1.5-5           preprocessCore_1.40.0  compiler_3.4.2        
    ## [19] cli_1.0.0              rvest_0.3.2            xml2_1.2.0            
    ## [22] scales_0.5.0           psych_1.7.8            affy_1.56.0           
    ## [25] digest_0.6.18          foreign_0.8-69         rmarkdown_1.10        
    ## [28] XVector_0.18.0         pkgconfig_2.0.1        htmltools_0.3.6       
    ## [31] limma_3.34.8           rlang_0.3.0.1          GlobalOptions_0.0.12  
    ## [34] readxl_1.0.0           rstudioapi_0.7         impute_1.52.0         
    ## [37] BiocInstaller_1.28.0   shiny_1.1.0            shape_1.4.4           
    ## [40] bindr_0.1.1            jsonlite_1.5           RCurl_1.95-4.10       
    ## [43] magrittr_1.5           GenomeInfoDbData_1.0.0 Matrix_1.2-12         
    ## [46] Rcpp_0.12.19           munsell_0.4.3          vsn_3.46.0            
    ## [49] stringi_1.2.4          yaml_2.2.0             zlibbioc_1.24.0       
    ## [52] plyr_1.8.4             promises_1.0.1         shinydashboard_0.6.1  
    ## [55] crayon_1.3.4           lattice_0.20-35        haven_1.1.1           
    ## [58] circlize_0.4.3         hms_0.4.1              knitr_1.20            
    ## [61] pillar_1.1.0           rjson_0.2.15           reshape2_1.4.3        
    ## [64] glue_1.3.0             evaluate_0.11          modelr_0.1.1          
    ## [67] httpuv_1.4.5           cellranger_1.1.0       gtable_0.2.0          
    ## [70] assertthat_0.2.0       mime_0.6               xtable_1.8-3          
    ## [73] broom_0.4.3            later_0.7.5            bindrcpp_0.2.2

[1] Sig M (2017). MultiAssayExperiment: Software for the integration of multi-omics experiments in Bioconductor. R package version 1.2.1
