# INAC (Integrated analysis toolkit for dissecting whole genome wide features of cell-free DNA)
## *Jie Li & Xun Lan*

## INAC function
Taking into consideration of large BAM file size of cfDNA sequencing files and long running process, We used sample case file to show the input and output of INAC function. Part of INAC functions is designed to be run on Unix-based operating systems such as macOS and linux.

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
Rscript INAC_TSS_NDR sample_name, input_dir, output_dir, samtools, tss_bed_up3000_1000, tss_bed_down1000_3000, tss_bed_up1000_down1000, tss_bed_up150_down50, tss_pro_table
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
Rscript INAC_TSS_2K sample_name, input_dir, output_dir, samtools, tss_bed_up3000_1000, tss_bed_down1000_3000, tss_bed_up1000_down1000, tss_bed_up150_down50, tss_pro_table
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
Rscript INAC_PFE sample_name, input_dir, output_dir, samtools, bedtools, consensusBlacklist, tss_bed_up1000_down1000, tss_bed_up1000_up750, tss_bed_down750_down1000, tss_pro_table, report_NBT
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
Rscript INAC_ML sample_name, input_dir, output_dir, cancer_sample, healthy_sample, seed, divide_fraction, choose_method
```


### sessionInfo





