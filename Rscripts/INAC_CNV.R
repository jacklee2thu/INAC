
###### INAC CNV
argv <- commandArgs(TRUE)
sample_name<-as.character(argv[1])
input_dir<-as.character(argv[2])
output_dir<-as.character(argv[3])
samtools<-as.character(argv[4])
consensusBlacklist<-as.character(argv[5])
bin_gc<-as.character(argv[6])
healthy_standard_copy<-as.character(argv[7])

options(stringsAsFactors=F)
setwd(output_dir)
system(paste0(samtools,' view -q 30 ',input_dir,'/',sample_name,'.bam|cut -f 3,4,9',' > ',sample_name,'_length.txt'),intern=TRUE, wait=TRUE)

cfDNA_length<-read.table(file=paste0(sample_name,'_length.txt'),sep="\t",header=F)
blacklist<-read.table(file=consensusBlacklist,sep="\t",header=F)
colnames(blacklist)<-c("chr","start","end","lable","num","dot")
colnames(cfDNA_length)<-c("chr","start","length")

aoto_chr<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7",
"chr8","chr9","chr10","chr11","chr12","chr13","chr14",
"chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")

new_cfDNA_length<-cfDNA_length[cfDNA_length$length>0,]
new_cfDNA_length<-new_cfDNA_length[new_cfDNA_length$chr%in%aoto_chr,]
#save(new_cfDNA_length,file="new_cfDNA_length.Rdata")

chr_exclude_cfDNA<-list()
for(i in 1:length(aoto_chr)){
chr_cfDNA<-new_cfDNA_length[new_cfDNA_length$chr%in%aoto_chr[i],]
chr_blacklist<-blacklist[blacklist$chr%in%aoto_chr[i],]
exclude_index<-c()
for(j in 1:dim(chr_blacklist)[1]){
exclude_index_temp<-which(chr_cfDNA$start>chr_blacklist[j,2]&chr_cfDNA$start<chr_blacklist[j,3])
exclude_index<-c(exclude_index,exclude_index_temp)
}
chr_exclude_cfDNA[[i]]<-chr_cfDNA[-exclude_index,]
}
names(chr_exclude_cfDNA)<-aoto_chr

all_exclude_cfDNA<-data.frame()
for(k in 1:length(chr_exclude_cfDNA)){
all_exclude_cfDNA<-rbind(all_exclude_cfDNA,chr_exclude_cfDNA[[k]])
}

#######统计长度
bin_gc<-read.table(file=bin_gc,sep="\t",header=T)##50kb


proper_cfDNA<-all_exclude_cfDNA[which(all_exclude_cfDNA$length>100&all_exclude_cfDNA$length<220),]

summary_copy<-data.frame()
for(i in 1:length(aoto_chr)){
proper_chr_cfDNA<-proper_cfDNA[proper_cfDNA$chr%in%aoto_chr[i],]
chr_copy<-bin_gc[bin_gc$chromosome%in%aoto_chr[i],]
temp_info_copy<-data.frame()
summary_chr_cfDNA<-data.frame()

for(k in 1:dim(chr_copy)[1]){
new_copy<-proper_chr_cfDNA[which(proper_chr_cfDNA$start>chr_copy$start[k]&proper_chr_cfDNA$start<chr_copy$end[k]),]
short_copy<-length(which(new_copy$length>100&new_copy$length<=150))
long_copy<-length(which(new_copy$length>150&new_copy$length<=220))
new_info_copy<-c(chr_copy[k,],short_copy,long_copy,short_copy+long_copy,short_copy/long_copy)
names(new_info_copy)<-c(colnames(chr_copy),"short","long","frag","ratio")
temp_info_copy<-rbind(temp_info_copy,new_info_copy)
}
summary_copy<-rbind(summary_copy,temp_info_copy)
}


##50Kb为一个bin统计量copy_num

gc.correct <- function(coverage, bias) {
    i <- seq(min(bias, na.rm=TRUE), max(bias, na.rm=TRUE), by = 0.001)
    coverage.trend <- loess(coverage ~ bias)
    coverage.model <- loess(predict(coverage.trend, i) ~ i)
    coverage.pred <- predict(coverage.model, bias)
    coverage.corrected <- coverage - coverage.pred + median(coverage)
}

##gastric
gastric_table<-data.frame()
for(j in 1:length(sample_name)){
load(sample_name[j])
summary_copy$ratio.center<-scale(summary_copy$ratio,scale = F)
temp_table<-data.frame(rep(last_sample[j],51120),summary_copy,bin_gc$gc)
temp_table$short.corrected=gc.correct(temp_table$short, temp_table$bin_gc.gc)
temp_table$short.corrected_fre<-temp_table$short.corrected/sum(na.omit(temp_table$short.corrected))
temp_table$long.corrected=gc.correct(temp_table$long, temp_table$bin_gc.gc)
temp_table$long.corrected_fre<-temp_table$long.corrected/sum(na.omit(temp_table$long.corrected))
temp_table$nfrags.corrected=gc.correct(temp_table$short+temp_table$long, temp_table$bin_gc.gc)
temp_table$nfrags.corrected_fre<-temp_table$nfrags.corrected/sum(na.omit(temp_table$nfrags.corrected))
#temp_table$ratio.corrected=gc.correct(temp_table$ratio, temp_table$bin_gc.gc)

gastric_table<-rbind(gastric_table,temp_table)
}
colnames(gastric_table)[1]<-"sample_name"

#
load(healthy_standard_copy)
gastric_sample<-names(table(gastric_table$sample_name))
gastric_copy_number<-list()
for(i in 1:length(gastric_sample)){
gastric_copy_number[[i]]<-log2(gastric_table[gastric_table$sample_name%in%gastric_sample[i],]$nfrags.corrected_fre/healthy_standard_copy)
}
names(gastric_copy_number)<-gastric_sample

gastric_table$copy_num<-unlist(gastric_copy_number)
save(gastric_table,file="gastric_copy_number_table.Rdata")

library(ggplot2)
pdf('gastric_cancer_copy_num_density.pdf',width=50,height=10)
p <- ggplot(gastric_table[gastric_table$sample_name%in%sample_name,], aes(x=bin, y=copy_num))+ylim(-1,1)
p + stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white")+scale_fill_gradient(low="yellow", high="red")+
theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
dev.off()

