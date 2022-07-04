
###### INAC initial quality control
argv <- commandArgs(TRUE)
sample_name<-as.character(argv[1])
input_dir<-as.character(argv[2])
output_dir<-as.character(argv[3])
samtools<-as.character(argv[4])

options(stringsAsFactors=F)
setwd(output_dir)
system(paste0(samtools,' coverage ',input_dir,'/',sample_name,'.bam',' > ',sample_name,'_coverage.txt'),intern=TRUE, wait=TRUE)

LJ_10_12_he_coverage<-read.table(file=paste0(sample_name,'_coverage.txt'),sep="\t",header=F)
colnames(LJ_10_12_he_coverage)<-c('rname','startpos','endpos','numreads','covbases','coverage','meandepth','meanbaseq','meanmapq')
library(ggplot2)
library(ggpubr)
library(ggrepel)

pdf('Initial_QC.pdf',width=15,height=10)
LJ_10_12_he_coverage<-LJ_10_12_he_coverage[1:25,]
ggbarplot(LJ_10_12_he_coverage, x = "rname", y = "numreads", 
          add = c("mean_se"),sort.by.groups=T,
          color = "red2",fill = "red2", palette = "npg",
position = position_dodge(0.8))
ggbarplot(LJ_10_12_he_coverage, x = "rname", y = "meanmapq", 
          add = c("mean_se"),sort.by.groups=T,
          color = "red2",fill = "red2", palette = "npg",
position = position_dodge(0.8))
ggbarplot(LJ_10_12_he_coverage, x = "rname", y = "coverage", 
          add = c("mean_se"),sort.by.groups=T,
          color = "red2",fill = "red2", palette = "npg",
position = position_dodge(0.8))
ggbarplot(LJ_10_12_he_coverage, x = "rname", y = "meandepth", 
          add = c("mean_se"),sort.by.groups=T,
          color = "red2",fill = "red2", palette = "npg",
position = position_dodge(0.8))
dev.off()

