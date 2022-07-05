# INAC (Integrated analysis toolkit for dissecting whole genome wide features of cell-free DNA)
## *Jie Li & Xun Lan*

## INAC function
Taking into consideration of large BAM file size of cfDNA sequencing files and long running process, We used a downsample example cfDNA BAM file to show the input and output of INAC function. Of course, the final processed output files of each feature also be supported in the [Datasets](https://github.com/jacklee2thu/INAC/tree/main/Datasets). (Part of INAC functions is designed to be run on Unix-based operating systems such as macOS and linux)

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
Rscript INAC_initial_QC.R sample_name,input_dir,output_dir,samtools
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
Rscript INAC_QC.R sample_name,input_dir,output_dir,samtools,consensusBlacklist
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


##### Examples
```
Rscript INAC_FR.R sample_name, input_dir, output_dir, samtools, consensusBlacklist, bin
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
Rscript INAC_FR_visibility input_dir, output_dir, bin
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
Rscript INAC_CNV sample_name, input_dir, output_dir, samtools, consensusBlacklist, bin_gc, healthy_standard_copy
```



