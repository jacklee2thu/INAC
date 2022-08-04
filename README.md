# INAC (Integrated analysis toolkit for dissecting whole genome wide features of cell-free DNA)
## *Jie Li & Xun Lan*

Patterns in whole genome wide features of cell-free DNA (cfDNA) in human plasma provide a noninvasive diagnostic approach for cancer detection. Here, we provide integrated analysis toolkit for whole genome wide features of cfDNA (INAC) to explore **fragment size, fragment size ratio, copy number variance, TSS relative coverage, Promoter fragmentation entropy**. a collected independent dataset (50 patients with gastric cancer and 50 healthy controls) with ~10x sequence depth were estimated by INAC. INAC also explores the relationships between these cfDNA features and other genomic or transcriptomic features including chromatin accessibility, gene expression levels and tumor derived programs, and reveals the aberrant patterns in patients with cancer.
we also recommend reading the [publication]() to get a better idea of what INAC can do.

## INAC workflow
![](https://github.com/jacklee2thu/INAC/blob/main/image/workflow.jpg)


## INAC function
Taking into consideration of large BAM file size of cfDNA sequencing files and long running process, We used sample case file to show the input and output of INAC function. the prepared [downsample BAM file](https://github.com/jacklee2thu/INAC/blob/main/Datasets/downsample_cfDNA.bam) and [index file](https://github.com/jacklee2thu/INAC/blob/main/Datasets/downsample_cfDNA.bam.bai) were supported to users to reproduct the INAC functions. Part of INAC functions is designed to be run on Unix-based operating systems such as macOS and linux.

|Function name|Type|Input files|main out|
|:--:|:--:|:--:|:--:|
|[INAC_initial_QC](#inac_initial_qc)|QC|BAM|the mapped reads count, coverage, mean depth, mean MAPQ in 24 chromatin and mitochondria|
|[INAC_QC](#inac_qc)|QC|BAM|the fraction of cfDNA fragment size on 30-80, 80-150, 150-220, 220-1000, 1000-longer (bp)|
|[INAC_FR](#inac_fr)|Feature|BAM|the counts and fraction of short and long cfDNA fragments on whole genome wide|
|[INAC_FR_visibility](#inac_fr_visibility)|visibility|feature matrix|the cfDNA fragments ratio on whole genome wide|
|[INAC_CNV](#inac_cnv)|Feature|BAM|the number of copy number variance of cfDNA on whole genome wide|
|[INAC_TSS_NDR](#inac_tss_ndr)|Feature|BAM|the relative coverage of NDR around TSS locations|
|[INAC_TSS_2K](#inac_tss_2k)|Feature|BAM|the relative coverage of 2K region around TSS locations|
|[INAC_PFE](#inac_pfe)|Feature|BAM|the PFE values around TSS locations|
|[INAC_ML](#inac_ml)|visibility|feature matrix|model performance of each feature|


### [INAC_initial_QC]()
##### Description
‘INAC_initial_QC’ could indicate cfDNA genome mapped reads count, coverage, mean depth, mean MAPQ in 24 chromatin and mitochondria.

##### Usage
INAC_initial_QC (sample_name,input_dir,output_dir,samtools)  
##### Arguments
|arguments|meaning|
|:--:|:--:|
|sample_name|the name of BAM file|
|input_dir|the dir path of BAM file|
|output_dir|the output file path of BAM file|  
|samtools|the path of samtools|

##### output file
1. [sample_name_coverage.txt](https://github.com/jacklee2thu/INAC/blob/main/Datasets/LJ_10_12_he_coverage.txt) the text file contains the cfDNA genome mapped reads count, coverage, mean depth, mean MAPQ in all chromatins.  
2. [Initial_QC.pdf](https://github.com/jacklee2thu/INAC/blob/main/Datasets/Initial_QC.pdf) the four barplots contain the cfDNA genome mapped reads count, coverage, mean depth, mean MAPQ in 24 chromatin and mitochondria.

##### Examples
```
Rscript INAC_initial_QC.R sample_name input_dir output_dir samtools
```

### [INAC_QC](https://github.com/jacklee2thu/INAC/blob/main/Rscripts/INAC_QC.R)
##### Description
‘INAC_QC’ could indicate the fraction of cfDNA fragment size on 30-80, 80-150, 150-220, 220-1000, 1000-longer (bp).

##### Usage
INAC_QC (sample_name,input_dir,output_dir,samtools,consensusBlacklist)
##### Arguments
|arguments|meaning|
|:--:|:--:|
|sample_name|the name of BAM file|
|input_dir|the dir path of BAM file|
|output_dir|the output file path of BAM file|  
|samtools|the path of samtools|
|[consensusBlacklist](https://github.com/jacklee2thu/INAC/blob/main/materials/consensusBlacklist.txt)|consensus Blacklist in the materials|

##### output file
1. [cfDNA_density.pdf](https://github.com/jacklee2thu/INAC/blob/main/Datasets/cfDNA_density.pdf) the cfDNA fragment size density plot om 0-500 bp.  
2. [INAC_fraction_QC_barplot.pdf](https://github.com/jacklee2thu/INAC/blob/main/Datasets/INAC_fraction_QC_barplot.pdf) the fraction of cfDNA fragment size on 30-80, 80-150, 150-220, 220-1000, 1000-longer (bp).

##### Examples
```
Rscript INAC_QC.R sample_name input_dir output_dir samtools consensusBlacklist
```

### [INAC_FR](https://github.com/jacklee2thu/INAC/blob/main/Rscripts/INAC_FR.R)
##### Description
‘INAC_FR’ could indicate the counts and fraction of short and long cfDNA fragments on whole genome wide.

##### Usage
INAC_FR (sample_name, input_dir, output_dir, samtools, consensusBlacklist, bin)
##### Arguments
|arguments|meaning|
|:--:|:--:|
|sample_name|the name of BAM file|
|input_dir|the dir path of BAM file|
|output_dir|the output file path of fragment ratio table|  
|samtools|the path of samtools|
|[consensusBlacklist](https://github.com/jacklee2thu/INAC/blob/main/materials/consensusBlacklist.txt)|consensus Blacklist in the materials|
|[bin](https://github.com/jacklee2thu/INAC/blob/main/materials/bin.txt)|5 mb bins bed file in the materials|

##### output file
1. [summary_table.txt](https://github.com/jacklee2thu/INAC/blob/main/Datasets/summary_table.txt) a sample value of short (<150 bp), long (>150 bp) cfDNA fragment, ratio and ratio.center in each 5 mb on whole genome wide.  
2. sample_name_table.Rdata summary_table.txt transfored into Rdata form.  


##### Examples
```
Rscript INAC_FR.R sample_name input_dir output_dir samtools consensusBlacklist bin
```

### [INAC_FR_visibility](https://github.com/jacklee2thu/INAC/blob/main/Rscripts/INAC_FR_visibility.R)
##### Description
‘INAC_FR_visibility’ could integrate the fragment ratio table and show the different pattern of fragment ratio between cancer patients and healthy controls.

##### Usage
INAC_FR_visibility (input_dir, output_dir, bin)
##### Arguments
|arguments|meaning|
|:--:|:--:|
|input_dir|the dir path of fragment ratio table|
|output_dir|the output file path of integrated fragment ratio table and plots|  
|[bin](https://github.com/jacklee2thu/INAC/blob/main/materials/bin.txt)|5 mb bins bed file in the materials|

##### output file
1. [all_table.Rdata](https://github.com/jacklee2thu/INAC/blob/main/Datasets/all_table.Rdata) a sample value of short (<150 bp), long (>150 bp) cfDNA fragment, ratio and other factors in each 5 mb on whole genome wide.  
2. [INAC_FR_visibility.pdf](https://github.com/jacklee2thu/INAC/blob/main/Datasets/INAC_FR_visibility.pdf) integrated fragment ratio plots.  

##### Examples
```
Rscript INAC_FR_visibility input_dir output_dir bin
```


### [INAC_CNV](https://github.com/jacklee2thu/INAC/blob/main/Rscripts/INAC_CNV.R)
##### Description
‘INAC_CNV’ could indicate the number of copy number variance of cfDNA on whole genome wide.

##### Usage
INAC_CNV (sample_name, input_dir, output_dir, samtools, consensusBlacklist, bin_gc, healthy_standard_copy)
##### Arguments
|arguments|meaning|
|:--:|:--:|
|sample_name|the name of BAM file|
|input_dir|the dir path of BAM file|
|output_dir|the output file path of log2 transformed copy number variance and plots|
|samtools|the path of samtools|
|[consensusBlacklist](https://github.com/jacklee2thu/INAC/blob/main/materials/consensusBlacklist.txt)|consensus Blacklist in the materials|
|[bin_gc](https://github.com/jacklee2thu/INAC/blob/main/materials/bin_gc.txt)|50 kb bins bed file in the materials|
|[healthy_standard_copy](https://github.com/jacklee2thu/INAC/blob/main/materials/healthy_standard_copy.Rdata)|golden standard file of healthy controls in the materials|

##### output file
1. [gastric_copy_number_table.Rdata](https://github.com/jacklee2thu/INAC/blob/main/Datasets/gastric_copy_number_table.Rdata) a sample value of log2 transformed copy number variance in each 50 kb on whole genome wide.  
2. [gastric_cancer_copy_num_density.pdf](https://github.com/jacklee2thu/INAC/blob/main/Datasets/gastric_cancer_copy_num_density.pdf) the density plot of log2 transformed copy number variance.  

##### Examples
```
Rscript INAC_CNV sample_name input_dir output_dir samtools consensusBlacklist bin_gc healthy_standard_copy
```


### [INAC_TSS_NDR](https://github.com/jacklee2thu/INAC/blob/main/Rscripts/INAC_TSS_NDR.R)
##### Description
‘INAC_TSS_NDR’ could indicate the relative coverage of NDR around TSS locations.

##### Usage
INAC_TSS_NDR (sample_name, input_dir, output_dir, samtools, tss_bed_up3000_1000, tss_bed_down1000_3000, tss_bed_up1000_down1000, tss_bed_up150_down50, tss_pro_table)
##### Arguments
|arguments|meaning|
|:--:|:--:|
|sample_name|the name of BAM file|
|input_dir|the dir path of BAM file|
|output_dir|the output file path of the relative coverage of NDR file|
|samtools|the path of samtools|
|[tss_bed_up3000_1000](https://github.com/jacklee2thu/INAC/blob/main/materials/tss_bed_up3000_1000.bed)|the bed file contains upstream 3000 bp to downstream 1000 bp relative to TSS locations|
|[tss_bed_down1000_3000](https://github.com/jacklee2thu/INAC/blob/main/materials/tss_bed_down1000_3000.bed)|the bed file contains upstream 1000 bp to downstream 3000 bp relative to TSS locations|
|[tss_bed_up1000_down1000](https://github.com/jacklee2thu/INAC/blob/main/materials/tss_bed_up1000_down1000.bed)|the bed file contains upstream 1000 bp to downstream 1000 bp relative to TSS locations|
|[tss_bed_up150_down50](https://github.com/jacklee2thu/INAC/blob/main/materials/tss_bed_up150_down50.bed)|the bed file contains upstream 150 bp to downstream 50 bp relative to TSS locations|
|[tss_pro_table](https://github.com/jacklee2thu/INAC/blob/main/materials/tss_pro_table.Rdata)|TSS locations information|

##### output file
1. [sample_name_tss_NDR.Rdata](https://github.com/jacklee2thu/INAC/blob/main/Datasets/sample_name_tss_NDR.Rdata) a sample value of the relative coverage of NDR around TSS locations on whole genome wide.  


##### Examples
```
Rscript INAC_TSS_NDR sample_name input_dir output_dir samtools tss_bed_up3000_1000 tss_bed_down1000_3000 tss_bed_up1000_down1000 tss_bed_up150_down50 tss_pro_table
```


### [INAC_TSS_2K](https://github.com/jacklee2thu/INAC/blob/main/Rscripts/INAC_TSS_2K.R)
##### Description
‘INAC_TSS_2K’ could indicate the relative coverage of 2K region around TSS locations.

##### Usage
INAC_TSS_2K (sample_name, input_dir, output_dir, samtools, tss_bed_up3000_1000, tss_bed_down1000_3000, tss_bed_up1000_down1000, tss_bed_up150_down50, tss_pro_table)
##### Arguments
|arguments|meaning|
|:--:|:--:|
|sample_name|the name of BAM file|
|input_dir|the dir path of BAM file|
|output_dir|the output file path of the relative coverage of 2K region file|
|samtools|the path of samtools|
|[tss_bed_up3000_1000](https://github.com/jacklee2thu/INAC/blob/main/materials/tss_bed_up3000_1000.bed)|the bed file contains upstream 3000 bp to downstream 1000 bp relative to TSS locations|
|[tss_bed_down1000_3000](https://github.com/jacklee2thu/INAC/blob/main/materials/tss_bed_down1000_3000.bed)|the bed file contains upstream 1000 bp to downstream 3000 bp relative to TSS locations|
|[tss_bed_up1000_down1000](https://github.com/jacklee2thu/INAC/blob/main/materials/tss_bed_up1000_down1000.bed)|the bed file contains upstream 1000 bp to downstream 1000 bp relative to TSS locations|
|[tss_bed_up150_down50](https://github.com/jacklee2thu/INAC/blob/main/materials/tss_bed_up150_down50.bed)|the bed file contains upstream 150 bp to downstream 50 bp relative to TSS locations|
|[tss_pro_table](https://github.com/jacklee2thu/INAC/blob/main/materials/tss_pro_table.Rdata)|TSS locations information|

##### output file
1. [sample_name_tss_2K.Rdata](https://github.com/jacklee2thu/INAC/blob/main/Datasets/sample_name_tss_2K.Rdata) a sample value of the relative coverage of 2K region around TSS locations on whole genome wide.  


##### Examples
```
Rscript INAC_TSS_2K sample_name input_dir output_dir samtools tss_bed_up3000_1000 tss_bed_down1000_3000 tss_bed_up1000_down1000 tss_bed_up150_down50 tss_pro_table
```


### [INAC_PFE](https://github.com/jacklee2thu/INAC/blob/main/Rscripts/INAC_PFE.R)
##### Description
‘INAC_PFE’ could indicate the PFE values around TSS locations.

##### Usage
INAC_PFE (sample_name, input_dir, output_dir, samtools, bedtools, consensusBlacklist, tss_bed_up1000_down1000, tss_bed_up1000_up750, tss_bed_down750_down1000, tss_pro_table, report_NBT)
##### Arguments
|arguments|meaning|
|:--:|:--:|
|sample_name|the name of BAM file|
|input_dir|the dir path of BAM file|
|output_dir|the output file path of PFE file|
|samtools|the path of samtools|
|bedtools|the path of bedtools|
|[consensusBlacklist](https://github.com/jacklee2thu/INAC/blob/main/materials/consensusBlacklist.txt)|consensus Blacklist in the materials|
|[tss_bed_up1000_down1000](https://github.com/jacklee2thu/INAC/blob/main/materials/tss_bed_up1000_down1000.bed)|the bed file contains upstream 1000 bp to downstream 1000 bp relative to TSS locations|
|[tss_bed_up1000_up750](https://github.com/jacklee2thu/INAC/blob/main/materials/tss_bed_up1000_up750.bed)|the bed file contains upstream 1000 bp to upstream 750 bp relative to TSS locations|
|[tss_bed_down750_down1000](https://github.com/jacklee2thu/INAC/blob/main/materials/tss_bed_down750_down1000.bed)|the bed file contains downstream 750 bp to downstream 1000 bp relative to TSS locations|
|[tss_pro_table](https://github.com/jacklee2thu/INAC/blob/main/materials/tss_pro_table.Rdata)|TSS locations information|
|[report_NBT](https://github.com/jacklee2thu/INAC/blob/main/materials/all.tss.genes.canonical.ensembl75.txt)| NBT reported negative control genes|

##### output file
1. [gastric_cancer_adjust_final_ent.Rdata](https://github.com/jacklee2thu/INAC/blob/main/Datasets/gastric_cancer_adjust_final_ent.Rdata) a sample value of PFE around TSS locations on whole genome wide.  


##### Examples
```
Rscript INAC_PFE sample_name input_dir output_dir samtools bedtools consensusBlacklist tss_bed_up1000_down1000 tss_bed_up1000_up750 tss_bed_down750_down1000 tss_pro_table report_NBT
```

### [INAC_ML](https://github.com/jacklee2thu/INAC/blob/main/Rscripts/INAC_ML.R)
##### Description
‘INAC_ML’ could indicate the model performance of each feature.


##### Usage
INAC_ML (sample_name, input_dir, output_dir, cancer_sample, healthy_sample, seed, divide_fraction, choose_method)
##### Arguments
|arguments|meaning|
|:--:|:--:|
|sample_name|the file of feature matrix, [example data](https://github.com/jacklee2thu/INAC/blob/main/Features/fragment_ratio_matrix.Rdata) is showed. [published data]() containg 50 gastric cancer patients and 50 healthy controls with each feature is here|
|input_dir|the dir path of feature matrix|
|output_dir|the output file path of model performance file|
|cancer_sample|the list of cancer sample names|
|seed|the seed used to divide dataset into training and test datasets|
|divide_fraction|the fraction used to divide dataset into training and test datasets, such as 1/2, 2/3, 3/4|
|choose_method|machine learning method, such as 'gbm','rf'. which could be found in [here](http://topepo.github.io/caret/train-models-by-tag.html)|


##### output file
1. [train_ROC.pdf](https://github.com/jacklee2thu/INAC/blob/main/Datasets/fragment_size_train_ROC.pdf) the training dataset ROC plot contains AUC value.  
2. [test_ROC.pdf](https://github.com/jacklee2thu/INAC/blob/main/Datasets/fragment_size_test_ROC.pdf) the test dataset ROC plot contains AUC value. 
3. [confusion_matrix_test_data.Rdata](https://github.com/jacklee2thu/INAC/blob/main/Datasets/confusion_matrix_fragment_size.Rdata) the test dataset confusion matrix. 


##### Examples
```
Rscript INAC_ML sample_name input_dir output_dir cancer_sample healthy_sample seed divide_fraction choose_method
```


### sessionInfo
```
R version 3.6.0 (2019-04-26)
Platform: x86_64-redhat-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS/LAPACK: /usr/lib64/R/lib/libRblas.so

locale:
 [1] LC_CTYPE=zh_CN.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=zh_CN.UTF-8        LC_COLLATE=zh_CN.UTF-8    
 [5] LC_MONETARY=zh_CN.UTF-8    LC_MESSAGES=zh_CN.UTF-8   
 [7] LC_PAPER=zh_CN.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=zh_CN.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] ROCR_1.0-11                            
 [2] gbm_2.1.8                              
 [3] pROC_1.17.0.1                          
 [4] caret_6.0-88                           
 [5] lattice_0.20-44                        
 [6] forcats_0.5.1                          
 [7] stringr_1.4.0                          
 [8] dplyr_1.0.7                            
 [9] purrr_0.3.4                            
[10] readr_2.0.1                            
[11] tidyr_1.1.3                            
[12] tibble_3.1.6                           
[13] tidyverse_1.3.1                        
[14] gtools_3.9.2                           
[15] limma_3.40.6                           
[16] clusterProfiler_3.12.0                 
[17] biovizBase_1.32.0                      
[18] BSgenome.Hsapiens.UCSC.hg19_1.4.0      
[19] BSgenome_1.52.0                        
[20] Biostrings_2.52.0                      
[21] XVector_0.24.0                         
[22] Homo.sapiens_1.3.1                     
[23] TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2
[24] org.Hs.eg.db_3.8.2                     
[25] GO.db_3.8.2                            
[26] OrganismDbi_1.26.0                     
[27] GenomicFeatures_1.36.4                 
[28] AnnotationDbi_1.46.1                   
[29] Biobase_2.44.0                         
[30] rtracklayer_1.44.4                     
[31] GenomicRanges_1.36.1                   
[32] GenomeInfoDb_1.20.0                    
[33] IRanges_2.18.3                         
[34] S4Vectors_0.22.1                       
[35] BiocGenerics_0.30.0                    
[36] pheatmap_1.0.12                        
[37] digest_0.6.27                          
[38] ggrepel_0.9.1                          
[39] ggpubr_0.4.0                           
[40] ggplot2_3.3.5                          

loaded via a namespace (and not attached):
  [1] utf8_1.2.2                  tidyselect_1.1.1           
  [3] RSQLite_2.2.8               htmlwidgets_1.5.3          
  [5] grid_3.6.0                  BiocParallel_1.18.1        
  [7] munsell_0.5.0               codetools_0.2-18           
  [9] withr_2.4.2                 colorspace_2.0-2           
 [11] GOSemSim_2.10.0             knitr_1.33                 
 [13] rstudioapi_0.13             ggsignif_0.6.2             
 [15] DOSE_3.10.2                 urltools_1.7.3             
 [17] GenomeInfoDbData_1.2.1      polyclip_1.10-0            
 [19] bit64_4.0.5                 farver_2.1.0               
 [21] vctrs_0.3.8                 generics_0.1.0             
 [23] ipred_0.9-11                xfun_0.25                  
 [25] R6_2.5.1                    graphlayouts_0.7.1         
 [27] AnnotationFilter_1.8.0      bitops_1.0-7               
 [29] cachem_1.0.6                fgsea_1.10.1               
 [31] gridGraphics_0.5-1          DelayedArray_0.10.0        
 [33] assertthat_0.2.1            scales_1.1.1               
 [35] ggraph_2.0.5                nnet_7.3-16                
 [37] enrichplot_1.4.0            gtable_0.3.0               
 [39] tidygraph_1.2.0             ensembldb_2.8.1            
 [41] timeDate_3043.102           rlang_0.4.11               
 [43] splines_3.6.0               rstatix_0.7.0              
 [45] lazyeval_0.2.2              ModelMetrics_1.2.2.2       
 [47] dichromat_2.0-0             broom_0.7.9                
 [49] europepmc_0.4               checkmate_2.0.0            
 [51] BiocManager_1.30.16         reshape2_1.4.4             
 [53] abind_1.4-5                 modelr_0.1.8               
 [55] backports_1.2.1             qvalue_2.16.0              
 [57] Hmisc_4.5-0                 RBGL_1.60.0                
 [59] lava_1.6.9                  tools_3.6.0                
 [61] ggplotify_0.0.9             ellipsis_0.3.2             
 [63] RColorBrewer_1.1-2          ggridges_0.5.3             
 [65] Rcpp_1.0.7                  plyr_1.8.6                 
 [67] base64enc_0.1-3             progress_1.2.2             
 [69] zlibbioc_1.30.0             RCurl_1.98-1.4             
 [71] prettyunits_1.1.1           rpart_4.1-15               
 [73] viridis_0.6.1               cowplot_1.1.1              
 [75] SummarizedExperiment_1.14.1 haven_2.4.3                
 [77] cluster_2.1.2               fs_1.5.0                   
 [79] magrittr_2.0.1              data.table_1.14.0          
 [81] DO.db_2.9                   openxlsx_4.2.4             
 [83] reprex_2.0.1                triebeard_0.3.0            
 [85] ProtGenerics_1.16.0         matrixStats_0.60.1         
 [87] hms_1.1.0                   XML_3.99-0.3               
 [89] rio_0.5.27                  jpeg_0.1-9                 
 [91] readxl_1.3.1                gridExtra_2.3              
 [93] compiler_3.6.0              biomaRt_2.40.5             
 [95] crayon_1.4.1                htmltools_0.5.1.1          
 [97] tzdb_0.1.2                  Formula_1.2-4              
 [99] lubridate_1.7.10            DBI_1.1.1                  
[101] tweenr_1.0.2                dbplyr_2.1.1               
[103] MASS_7.3-54                 Matrix_1.3-4               
[105] car_3.0-11                  cli_3.0.1                  
[107] gower_0.2.2                 igraph_1.2.6               
[109] pkgconfig_2.0.3             rvcheck_0.1.8              
[111] GenomicAlignments_1.20.1    foreign_0.8-71             
[113] recipes_0.1.16              foreach_1.5.1              
[115] xml2_1.3.2                  prodlim_2019.11.13         
[117] rvest_1.0.1                 yulab.utils_0.0.2          
[119] VariantAnnotation_1.30.1    graph_1.62.0               
[121] cellranger_1.1.0            fastmatch_1.1-3            
[123] htmlTable_2.2.1             curl_4.3.2                 
[125] Rsamtools_2.0.3             nlme_3.1-152               
[127] lifecycle_1.0.0             jsonlite_1.7.2             
[129] carData_3.0-4               viridisLite_0.4.0          
[131] fansi_0.5.0                 pillar_1.6.2               
[133] fastmap_1.1.0               httr_1.4.2                 
[135] survival_3.2-12             glue_1.4.2                 
[137] zip_2.2.0                   UpSetR_1.4.0               
[139] iterators_1.0.13            png_0.1-7                  
[141] bit_4.0.4                   class_7.3-19               
[143] ggforce_0.3.3               stringi_1.7.3              
[145] blob_1.2.2                  latticeExtra_0.6-29        
[147] memoise_2.0.0
```


## INAC citation

Jie Li, Xun Lan.2022-7. Integrated analysis toolkit for dissecting whole genome wide features of cell-free DNA

## INAC license

RECIPIENT acknowledges that the Program is a research tool still in the development stage and that it is being supplied as is, without any accompanying services, support, or improvements.Any accompanying information, materials, or manuals, free of charge for non-commercial use only and bound by the above [license](https://github.com/jacklee2thu/INAC/blob/main/License%20and%20Terms%20of%20Use.txt) agreement.

