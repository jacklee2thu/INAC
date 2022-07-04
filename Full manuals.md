# INAC (Integrated analysis toolkit for dissecting whole genome wide features of cell-free DNA)
## *Jie Li & Xun Lan*

## INAC function
Taking into consideration of large BAM file size of cfDNA sequencing files and long running process, We used a downsample example cfDNA BAM file to show the input and output of INAC function. Of course, the final processed output files of each feature also be supported in the [Datasets](https://github.com/jacklee2thu/INAC/tree/main/Datasets). (Part of INAC functions is designed to be run on Unix-based operating systems such as macOS and linux)

|Function name|Type|Input files|main out|
|:--:|:--:|:--:|:--:|
|[INAC_QC](#inac_initial_qc)|QC|BAM|the fraction of cfDNA fragment size on 30-80, 80-150, 150-220, 220-1000, 1000-longer (bp)|
|INAC_FR|Feature|BAM|the counts and fraction of short and long cfDNA fragments on whole genome wide|
|INAC_CNV|Feature|BAM|the number of copy number variance of cfDNA on whole genome wide|
|INAC_TSS_NDR|Feature|BAM|the relative coverage of NDR around TSS locations|
|INAC_TSS_2K|Feature|BAM|the relative coverage of 2K region around TSS locations|
|INAC_PFE|Feature|BAM|the PFE values around TSS locations|
|INAC_ML|machine learning|feature matrix|model performance of each feature|

### [INAC_initial_QC]()
##### Description
‘INAC_initial_QC’ could indicate cfDNA human genome mapped reads count, coverage, mean depth, mean MAPQ in 24 chromatin and mitochondria, and the fraction of cfDNA fragment size on 30-80, 80-150, 150-220, 220-1000, 1000-longer (bp).

Usage
addGCBias(object, ...)

output
(Ranged)SummarizedExperiment object with new column in row metadata with the gc content of the peak in question

Examples
data(example_counts, package = "chromVAR")
# show example on small part of data 
subset_counts <- example_counts[1:500,]
library(BSgenome.Hsapiens.UCSC.hg19)
example_counts <- addGCBias(subset_counts, 
                              genome = BSgenome.Hsapiens.UCSC.hg19)
