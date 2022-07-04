
###### INAC quality control
argv <- commandArgs(TRUE)
sample_name<-as.character(argv[1])
input_dir<-as.character(argv[2])
output_dir<-as.character(argv[3])
samtools<-as.character(argv[4])
consensusBlacklist<-as.character(argv[5])


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
cfDNA_frame<-data.frame()
cfDNA_frame<-all_exclude_cfDNA
L_length<-dim(all_exclude_cfDNA)[1]
library(ggplot2)
pdf("cfDNA_density.pdf")
ggplot(cfDNA_frame,aes(x=length,colour='red'))+geom_density()+xlim(0,500)+theme_bw()+
theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+labs(title = "cfDNA_size_density")
dev.off()
######百分比
proper_cfDNA_30_80<-all_exclude_cfDNA[which(all_exclude_cfDNA$length>30&all_exclude_cfDNA$length<80),]
proper_cfDNA_80_150<-all_exclude_cfDNA[which(all_exclude_cfDNA$length>80&all_exclude_cfDNA$length<150),]
proper_cfDNA_150_220<-all_exclude_cfDNA[which(all_exclude_cfDNA$length>150&all_exclude_cfDNA$length<220),]
proper_cfDNA_220_1000<-all_exclude_cfDNA[which(all_exclude_cfDNA$length>220&all_exclude_cfDNA$length<1000),]
proper_cfDNA_1000_longer<-all_exclude_cfDNA[which(all_exclude_cfDNA$length>1000),]

dfm<-data.frame()
temp_length<-L_length
temp_dfm<-data.frame(c(dim(proper_cfDNA_30_80)[1]/temp_length,dim(proper_cfDNA_80_150)[1]/temp_length,
dim(proper_cfDNA_150_220)[1]/temp_length,dim(proper_cfDNA_220_1000)[1]/temp_length,dim(proper_cfDNA_1000_longer)[1]/temp_length),
c('proper_cfDNA_30_80','proper_cfDNA_80_150',
'proper_cfDNA_150_220','proper_cfDNA_220_1000','proper_cfDNA_1000_longer'))
colnames(temp_dfm)<-c('reads fraction','group')
dfm<-temp_dfm
dfm[,2]<-factor(dfm[,2],levels=c('proper_cfDNA_30_80','proper_cfDNA_80_150',
'proper_cfDNA_150_220','proper_cfDNA_220_1000','proper_cfDNA_1000_longer'))


library(ggpubr)
library(ggrepel)
pdf('INAC_fraction_QC_barplot.pdf',width=12,height=10)
ggbarplot(dfm, x = "group", y = "reads fraction", 
          add = c("mean_se"),sort.by.groups=T,
          color = "group",fill = "group", palette = "npg",
position = position_dodge(0.8))
dev.off()