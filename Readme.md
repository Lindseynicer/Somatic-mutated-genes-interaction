---
title: "MAF_PRAD_Somatic_Interaction"
author: "Foong Lian Chee"
date: "4/18/2021"
output:
  html_document:
    keep_md: true
---
# Introduction

**Question:** How frequent is the co-occurrence of PTEN loss/mutate, TP53 loss/mutate, and the activation of "KRAS/MAPK/ERK" pathway genes in prostate cancer in the available dataset?

**Answer/output:**
See the last figure in this report, and also a csv file (maf_PRAD_mutect2.csv) is generated which showing all the somatic interaction of MAPK signaling pathway genes, together with PTEN and TP53 genes in TCGA PRAD dataset.

**Analysis/Work flow:**
The Cancer Genome Atlas Project (TCGA) has sequenced over 30 different cancers with sample size of each cancer type being over 200. Resulting data consisting of somatic variants are stored in the form of Mutation Annotation Format (MAF) and can be downloaded from National Cancer Institute Genomic Data Commons (NIH GDC) platform.

- A. Download MAF file
- B. Visualization of the MAF file
- C. Analysis of somatic interaction among the gene-of-interest

# A. Download MAF file

There are two ways to download MAF file:

## 1. via the GDC repository
<https://portal.gdc.cancer.gov/projects/TCGA-PRAD>
then go to "Files" > Repository > select "MAF" data format and "open" access > there are 4 types of MAF files to select (i.e. varscan, mutect, muse, and somaticsniper) 

## 2. through R "maftools" and "TCGAbiolinks" packages
"TCGAbiolinks" helps to retrieve the TCGA cancer database. Specifically, we use "GDCquery_Maf" (Mutation Annotation Format) to acquire the TCGA PRAD (shortform for prostate cancer) set.

In TCGA datasets, Variant calling is performed using five separate pipelines:MuSe, MuTect2, VarScan2, SomaticSniper, Pindel. According to GDC website, There is currently no scientific consensus on the best variant calling pipeline. Here, I use *Mutect2* pipelines because it contains relatively a richer mutation data compared to other pipelines.

* read here to find out the difference among these pipelines, https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/



```r
getwd()
```

```
## [1] "C:/Users/lianc/Documents/R_discovery"
```

```r
#setwd("R_discovery")

# refer to https://rpubs.com/connorH982/491637
## Search and download MAF data from cancer database 
library(maftools)
library(TCGAbiolinks)
library(dplyr)

# In GDC website, it mentions that MuTect2 pipeline employs a "Panel of Normals" to identify additional germline mutations.This panel is generated using TCGA blood normal genomes from thousands of individuals that were curated and confidently assessed to be cancer-free. This method allows for a higher level of confidence to be assigned to somatic variants that were called by the MuTect2 pipeline.

maf_mutect2 <- GDCquery_Maf("PRAD", pipelines = "mutect2") %>% subset(Sequencer == "Illumina HiSeq 2000") %>% read.maf
```

```
## -Validating
## -Silent variants: 9821 
## -Summarizing
## --Possible FLAGS among top ten genes:
##   TTN
##   MUC16
##   SYNE1
## -Processing clinical data
## --Missing clinical data
## -Finished in 18.2s elapsed (6.390s cpu)
```

```r
# read here to know what are the columns in a MAF file, https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/
```


```r
getGeneSummary(maf_mutect2)
```

```
##       Hugo_Symbol Frame_Shift_Del Frame_Shift_Ins In_Frame_Del In_Frame_Ins
##    1:        TP53               9               3            0            0
##    2:        SPOP               1               0            1            0
##    3:         TTN               2               5            0            0
##    4:       FOXA1               4               2           16            0
##    5:       KMT2D              11               2            2            0
##   ---                                                                      
## 9825:      ZSWIM5               0               0            0            0
## 9826:        ZW10               0               0            0            0
## 9827:      ZWILCH               0               0            0            0
## 9828:      ZYG11B               0               0            0            0
## 9829:        ZZZ3               0               0            0            0
##       Missense_Mutation Nonsense_Mutation Nonstop_Mutation Splice_Site
##    1:                39                 2                0           7
##    2:                54                 0                0           0
##    3:                62                 5                0           0
##    4:                10                 0                0           0
##    5:                14                 4                0           0
##   ---                                                                 
## 9825:                 1                 0                0           0
## 9826:                 1                 0                0           0
## 9827:                 1                 0                0           0
## 9828:                 0                 1                0           0
## 9829:                 0                 1                0           0
##       Translation_Start_Site total MutatedSamples AlteredSamples
##    1:                      0    60             57             57
##    2:                      0    56             55             55
##    3:                      0    74             52             52
##    4:                      0    32             32             32
##    5:                      0    33             27             27
##   ---                                                           
## 9825:                      0     1              1              1
## 9826:                      0     1              1              1
## 9827:                      0     1              1              1
## 9828:                      0     1              1              1
## 9829:                      0     1              1              1
```

```r
# This dataset shows some frame shift samples which is not found in the MuSe type data.
```

# B. Visualization of the MAF data
## Plotting MAF summary.
We can use plotmafSummary to plot the summary of the maf file, which displays number of variants in each sample as a stacked barplot and variant types as a boxplot summarized by Variant_Classification. We can add either mean or median line to the stacked barplot to display average/median number of variants across the cohort.

```r
# refer to https://bioconductor.statistik.tu-dortmund.de/packages/3.8/bioc/vignettes/maftools/inst/doc/maftools.html#7_visualization 
# plot the summary of the MAF file
plotmafSummary(maf = maf_mutect2 , rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
```

![](MAF_PRAD_Somatic_Interaction_20210418_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

## Oncoplots
Oncoplot is also known as waterfall plots.We draw oncoplots for top ten mutated genes and non-mutated samples are removed from the plot for better visualization.


```r
# Oncoplot
# Variants annotated as Multi_Hit are those genes which are mutated more than once in the same sample
oncoplot(maf = maf_mutect2 , top = 10, removeNonMutated = TRUE)
```

![](MAF_PRAD_Somatic_Interaction_20210418_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

## Transition and Transversions
titv function classifies SNPs into Transitions and Transversions.A boxplot showing overall distribution of six different conversions. The bottom is a stacked barplot showing fraction of conversions in each sample.


```r
titv = titv(maf = maf_mutect2 , plot = FALSE, useSyn = TRUE)
# plot titv summary
plotTiTv(res = titv)
```

![](MAF_PRAD_Somatic_Interaction_20210418_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

## Compare mutation load against TCGA PRAD cohorts
The default is to compare with all 30++ cancer cohorts. But here, I just wanna check if my data look similar to the TCGA PRAD dataset.


```r
# Compare mutation load against TCGA cohorts
# prad.mutload = tcgaCompare(maf = maf_mutect2 , cohortName = 'Mydata-PRAD')
prad.mutload = tcgaCompare(maf = maf_mutect2 , cohortName = 'Mydata-PRAD', tcga_cohorts = "PRAD", logscale = TRUE, capture_size = 50)
```

```
## Warning in FUN(X[[i]], ...): Removed 1 samples with zero mutations.
```

![](MAF_PRAD_Somatic_Interaction_20210418_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

## Plotting VAF (Variant Allele Frequencies)
This plot helps toestimate clonal status of top mutated genes (clonal genes usually have mean allele frequency around ~50% assuming pure sample).


```r
# Plotting VAF (Variant Allele Frequencies)
plotVaf(maf = maf_mutect2 )
```

![](MAF_PRAD_Somatic_Interaction_20210418_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

## Oncogenic Signaling Pathways
Now, have a glimpse on the pathways. OncogenicPathways function checks for enrichment of known Oncogenic Signaling Pathways in TCGA PRAD. In this figure, tumor suppressor genes are in red, and oncogenes are in blue font.


```r
OncogenicPathways(maf = maf_mutect2 )
```

![](MAF_PRAD_Somatic_Interaction_20210418_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

```
##        Pathway  N n_affected_genes fraction_affected Mutated_samples
##  1:       TP53  6                3         0.5000000              74
##  2:       NRF2  3                3         1.0000000               9
##  3:        MYC 13                6         0.4615385              10
##  4:   TGF-Beta  7                6         0.8571429              13
##  5: Cell_Cycle 15                8         0.5333333              13
##  6:       PI3K 29               22         0.7586207              49
##  7:      Hippo 38               26         0.6842105              73
##  8:      NOTCH 71               40         0.5633803              63
##  9:        WNT 68               48         0.7058824              71
## 10:    RTK-RAS 85               62         0.7294118              80
##     Fraction_mutated_samples
##  1:               0.14949495
##  2:               0.01818182
##  3:               0.02020202
##  4:               0.02626263
##  5:               0.02626263
##  6:               0.09898990
##  7:               0.14747475
##  8:               0.12727273
##  9:               0.14343434
## 10:               0.16161616
```

```r
PlotOncogenicPathways(maf = maf_mutect2 , pathways = "RTK-RAS")
```

![](MAF_PRAD_Somatic_Interaction_20210418_files/figure-html/unnamed-chunk-8-2.png)<!-- -->

```r
# Tumor suppressor genes are in red, and oncogenes are in blue font.
```

# C. Analysis of somatic interaction among the gene-of-interestAnalysis
## Somatic Interactions
This "somaticInteractions" package performs Pair-wise Fisher's Exact test to detect mutually exclusive or co-occuring events. I listed the genes in KEGG MAPK Signaling pathway from GSEA platform, and then added PTEN in the list manually (because it is not in the MAPK pathway), the file name is geneset_KEGG_MAPK.txt.


```r
# Genes to investigate
## https://www.gsea-msigdb.org/gsea/msigdb/cards/KEGG_MAPK_SIGNALING_PATHWAY
## added manually PTEN in the txt file
geneset <- read.delim("geneset_KEGG_MAPK.txt")
geneset <- geneset[-1,]

#subset MAF file that only contains geneset
maf_mutect2_subset <- subsetMaf(maf = maf_mutect2, genes = geneset)

# exclusive/co-occurance event analysis on top 25 mutated genes.
# this is with default parameters
# somaticInteractions(maf = maf_mutect2, top = 25, pvalue = c(0.05, 0.01))

# Now, to save a csv file, I make countType as "all" which include cooccurrence and mutually exclusive results
maf_PRAD_mutect2 <- somaticInteractions(maf = maf_mutect2_subset, top = 25, countStats = "sig", countType = "all", fontSize = 0.6, pvalue = c(0.05, 0.01))
```

![](MAF_PRAD_Somatic_Interaction_20210418_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

```r
# output as csv file
write.csv(maf_PRAD_mutect2, 
          file = "maf_PRAD_mutect2.csv", 
          quote = TRUE, 
          row.names = FALSE)
```


```r
sessionInfo()
```

```
## R version 4.0.4 (2021-02-15)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 19041)
## 
## Matrix products: default
## 
## locale:
## [1] LC_COLLATE=English_Malaysia.1252  LC_CTYPE=English_Malaysia.1252   
## [3] LC_MONETARY=English_Malaysia.1252 LC_NUMERIC=C                     
## [5] LC_TIME=English_Malaysia.1252    
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] dplyr_1.0.5         TCGAbiolinks_2.18.0 maftools_2.6.05    
## 
## loaded via a namespace (and not attached):
##  [1] MatrixGenerics_1.2.1        Biobase_2.50.0             
##  [3] httr_1.4.2                  tidyr_1.1.3                
##  [5] sass_0.3.1                  bit64_4.0.5                
##  [7] jsonlite_1.7.2              splines_4.0.4              
##  [9] R.utils_2.10.1              bslib_0.2.4                
## [11] assertthat_0.2.1            askpass_1.1                
## [13] highr_0.8                   BiocFileCache_1.14.0       
## [15] stats4_4.0.4                blob_1.2.1                 
## [17] GenomeInfoDbData_1.2.4      yaml_2.2.1                 
## [19] progress_1.2.2              pillar_1.6.0               
## [21] RSQLite_2.2.4               lattice_0.20-41            
## [23] downloader_0.4              glue_1.4.2                 
## [25] digest_0.6.27               GenomicRanges_1.42.0       
## [27] RColorBrewer_1.1-2          XVector_0.30.0             
## [29] rvest_1.0.0                 colorspace_2.0-0           
## [31] plyr_1.8.6                  htmltools_0.5.1.1          
## [33] Matrix_1.2-18               R.oo_1.24.0                
## [35] XML_3.99-0.6                pkgconfig_2.0.3            
## [37] biomaRt_2.46.3              zlibbioc_1.36.0            
## [39] purrr_0.3.4                 scales_1.1.1               
## [41] openssl_1.4.3               tibble_3.1.0               
## [43] generics_0.1.0              TCGAbiolinksGUI.data_1.10.0
## [45] IRanges_2.24.1              ggplot2_3.3.3              
## [47] ellipsis_0.3.1              cachem_1.0.4               
## [49] SummarizedExperiment_1.20.0 BiocGenerics_0.36.0        
## [51] survival_3.2-7              magrittr_2.0.1             
## [53] crayon_1.4.1                memoise_2.0.0              
## [55] evaluate_0.14               R.methodsS3_1.8.1          
## [57] fansi_0.4.2                 xml2_1.3.2                 
## [59] prettyunits_1.1.1           tools_4.0.4                
## [61] data.table_1.14.0           hms_1.0.0                  
## [63] lifecycle_1.0.0             matrixStats_0.58.0         
## [65] stringr_1.4.0               S4Vectors_0.28.1           
## [67] munsell_0.5.0               DelayedArray_0.16.3        
## [69] AnnotationDbi_1.52.0        compiler_4.0.4             
## [71] jquerylib_0.1.3             GenomeInfoDb_1.26.4        
## [73] rlang_0.4.10                grid_4.0.4                 
## [75] RCurl_1.98-1.3              rappdirs_0.3.3             
## [77] bitops_1.0-6                rmarkdown_2.7              
## [79] gtable_0.3.0                curl_4.3                   
## [81] DBI_1.1.1                   R6_2.5.0                   
## [83] knitr_1.31                  fastmap_1.1.0              
## [85] bit_4.0.4                   utf8_1.1.4                 
## [87] readr_1.4.0                 stringi_1.5.3              
## [89] parallel_4.0.4              Rcpp_1.0.6                 
## [91] vctrs_0.3.6                 dbplyr_2.1.0               
## [93] tidyselect_1.1.0            xfun_0.21
```

