# INAC (Integrated analysis toolkit for dissecting whole genome wide features of cell-free DNA)
## *Jie Li & Xun Lan*

|Function name|Type|Input files|main out|
|:--:|:--:|:--:|:--:|
|INAC_QC|QC|BAM|the fraction of cfDNA fragment size on 30-80, 80-150, 150-220, 220-1000, 1000-longer (bp)|


addGCBias
Description
Computes GC content for peaks

Usage
addGCBias(object, ...)

## S4 method for signature 'RangedSummarizedExperiment'
addGCBias(object,
  genome = GenomeInfoDb::genome(object))

## S4 method for signature 'SummarizedExperiment'
addGCBias(object, peaks,
  genome = GenomeInfoDb::genome(peaks))
Arguments
object	
(Ranged)SummarizedExperiment

...	
additional arguments

genome	
BSgenome object, by defualt hg19

peaks	
GenomicRanges with peaks, needed if object is SummarizedExperiment and not RangedSummarizedExperiment

Value
(Ranged)SummarizedExperiment object with new column in row metadata with the gc content of the peak in question

Methods (by class)
RangedSummarizedExperiment: method for RangedSummarizedExperiment

SummarizedExperiment: method for SummarizedExperiment

Examples

data(example_counts, package = "chromVAR")
# show example on small part of data 
subset_counts <- example_counts[1:500,]
library(BSgenome.Hsapiens.UCSC.hg19)
example_counts <- addGCBias(subset_counts, 
                              genome = BSgenome.Hsapiens.UCSC.hg19)
