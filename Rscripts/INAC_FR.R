
###### INAC fragment size ratio
argv <- commandArgs(TRUE)
sample_name<-as.character(argv[1])
input_dir<-as.character(argv[2])
output_dir<-as.character(argv[3])
samtools<-as.character(argv[4])
consensusBlacklist<-as.character(argv[5])
bin<-as.character(argv[6])


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
bin_dataframe<-read.table(file=bin,sep="\t",header=T)


proper_cfDNA<-all_exclude_cfDNA[which(all_exclude_cfDNA$length>100&all_exclude_cfDNA$length<220),]

summary_table<-data.frame()
for(i in 1:length(aoto_chr)){
proper_chr_cfDNA<-proper_cfDNA[proper_cfDNA$chr%in%aoto_chr[i],]
chr_bin<-bin_dataframe[bin_dataframe$chromosome%in%aoto_chr[i],]
chr_copy<-bin_copy[bin_copy$chromosome%in%aoto_chr[i],]
temp_info_chr<-data.frame()###每条染色体上的summary
temp_info_copy<-data.frame()
summary_chr_cfDNA<-data.frame()
for(j in 1:dim(chr_bin)[1]){
new_bin<-proper_chr_cfDNA[which(proper_chr_cfDNA$start>chr_bin$start[j]&proper_chr_cfDNA$start<chr_bin$end[j]),]
short_num<-length(which(new_bin$length>100&new_bin$length<=150))
long_num<-length(which(new_bin$length>150&new_bin$length<=220))
new_info_chr<-c(chr_bin[j,],short_num,long_num,short_num+long_num,short_num/long_num)
names(new_info_chr)<-c(colnames(chr_bin),"short","long","frag","ratio")
temp_info_chr<-rbind(temp_info_chr,new_info_chr)
summary_chr_cfDNA<-rbind(new_bin,summary_chr_cfDNA)
}
summary_table<-rbind(summary_table,temp_info_chr)
}
summary_table$ratio.center<-scale(summary_table$ratio,scale=F)

write.table(summary_table,file='summary_table.txt',sep='\t',row.names = F,col.names = T,quote = F)


