# INAC (Integrated analysis toolkit for dissecting whole genome wide features of cell-free DNA)
## *Jie Li & Xun Lan*

Patterns in whole genome wide features of cell-free DNA (cfDNA) in human plasma provide a noninvasive diagnostic approach for cancer detection. Here, we provide integrated analysis toolkit for whole genome wide features of cfDNA (INAC) to explore **fragment size, fragment size ratio, copy number variance, TSS relative coverage, Promoter fragmentation entropy**. a collected independent dataset (50 patients with gastric cancer and 50 healthy controls) with ~10x sequence depth were estimated by INAC. INAC also explores the relationships between these cfDNA features and other genomic or transcriptomic features including chromatin accessibility, gene expression levels and tumor derived programs, and reveals the aberrant patterns in patients with cancer.
we also recommend reading the [publication]() to get a better idea of what INAC can do.

## INAC workflow
![](https://github.com/jacklee2thu/INAC/blob/main/image/workflow.jpg)
## INAC functions

|Function name|Type|Input files|main out|
|:--:|:--:|:--:|:--:|
|INAC_initial_QC|QC|BAM|the mapped reads count, coverage, mean depth, mean MAPQ in 24 chromatin and mitochondria|
|INAC_QC|QC|BAM|the fraction of cfDNA fragment size on 30-80, 80-150, 150-220, 220-1000, 1000-longer (bp)|
|INAC_FR|Feature|BAM|the counts and fraction of short and long cfDNA fragments on whole genome wide|
|INAC_FR_visibility|visibility|feature matrix|the cfDNA fragments ratio on whole genome wide|
|INAC_CNV|Feature|BAM|the number of copy number variance of cfDNA on whole genome wide|
|INAC_TSS_NDR|Feature|BAM|the relative coverage of NDR around TSS locations|
|INAC_TSS_2K|Feature|BAM|the relative coverage of 2K region around TSS locations|
|INAC_PFE|Feature|BAM|the PFE values around TSS locations|
|INAC_ML|visibility|feature matrix|model performance of each feature|
the details and usage of INAC functions are shown in the [full manuals](https://github.com/jacklee2thu/INAC/blob/main/Full%20manuals.md).

## INAC citation

Jie Li, Xun Lan.2022-7. Integrated analysis toolkit for dissecting whole genome wide features of cell-free DNA

## INAC license

RECIPIENT acknowledges that the Program is a research tool still in the development stage and that it is being supplied as is, without any accompanying services, support, or improvements.Any accompanying information, materials, or manuals, free of charge for non-commercial use only and bound by the above [license](https://github.com/jacklee2thu/INAC/blob/main/License%20and%20Terms%20of%20Use.txt) agreement.

