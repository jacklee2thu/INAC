
###### INAC PFE
argv <- commandArgs(TRUE)
sample_name<-as.character(argv[1])
input_dir<-as.character(argv[2])
output_dir<-as.character(argv[3])
samtools<-as.character(argv[4])
bedtools<-as.character(argv[4])
consensusBlacklist<-as.character(argv[5])
tss_pro_bed_up1000_down1000<-as.character(argv[5])
tss_pro_bed_up1000_up750<-as.character(argv[6])
tss_pro_bed_down750_down1000<-as.character(argv[7])
tss_pro_table<-as.character(argv[9])
report_NBT<-as.character(argv[9])

options(stringsAsFactors=F)
##跑出每个人要转录的TSS

setwd(input_dir)
if(is.na(match(paste0(sample_name,".bam.bai"),dir()))){
system(paste0(samtools,' index ',sample_name,".bam"))
}

system(paste0(samtools,' view -q 30 ',input_dir,'/',sample_name,'.bam|cut -f 3,4,9',' > ',sample_name,'_length.txt'),intern=TRUE, wait=TRUE)


cfDNA_length<-read.table(file=paste0(sample_name,"_length.txt"),sep="\t",header=F)
blacklist<-read.table(file=consensusBlacklist,sep="\t",header=F)
colnames(blacklist)<-c("chr","start","end","lable","num","dot")
colnames(cfDNA_length)<-c("chr","start","length")

aoto_chr<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7",
"chr8","chr9","chr10","chr11","chr12","chr13","chr14",
"chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")

new_cfDNA_length<-cfDNA_length[cfDNA_length$length>0,]
new_cfDNA_length<-new_cfDNA_length[new_cfDNA_length$chr%in%aoto_chr,]
##排除黑箱子区域
chr_exclude_cfDNA<-list()
for(i in 1:length(aoto_chr)){
chr_cfDNA<-new_cfDNA_length[new_cfDNA_length$chr%in%aoto_chr[i],]
chr_blacklist<-blacklist[blacklist$chr%in%aoto_chr[i],]
exclude_index<-c()###将黑箱子区域索引
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

########计算样本区间内的count数,只要有交集就纳入统计
all_exclude_cfDNA<-all_exclude_cfDNA[which(all_exclude_cfDNA$length<=300 & all_exclude_cfDNA$length>=100),]
tss_pro_bed_up1000_down1000<-read.table(file=tss_pro_bed_up1000_down1000,header=F,sep='\t')
colnames(tss_pro_bed_up1000_down1000)<-c('chr','start','end')
tss_pro_bed_up1000_down1000$combine_name<-paste0(tss_pro_bed_up1000_down1000$chr,'_',tss_pro_bed_up1000_down1000$start,'_',tss_pro_bed_up1000_down1000$end)
tss_pro_bed_up1000_up750<-read.table(file=tss_pro_bed_up1000_up750,header=F,sep='\t')
colnames(tss_pro_bed_up1000_up750)<-c('chr','start','end')
tss_pro_bed_up1000_up750$combine_name<-paste0(tss_pro_bed_up1000_up750$chr,'_',tss_pro_bed_up1000_up750$start,'_',tss_pro_bed_up1000_up750$end)
tss_pro_bed_down750_down1000<-read.table(file=tss_pro_bed_down750_down1000,header=F,sep='\t')
colnames(tss_pro_bed_down750_down1000)<-c('chr','start','end')
tss_pro_bed_down750_down1000$combine_name<-paste0(tss_pro_bed_down750_down1000$chr,'_',tss_pro_bed_down750_down1000$start,'_',tss_pro_bed_down750_down1000$end)
new_length_bed<-data.frame(all_exclude_cfDNA$chr,all_exclude_cfDNA$start,
c(all_exclude_cfDNA$start+all_exclude_cfDNA$length),all_exclude_cfDNA$length)
colnames(new_length_bed)<-c('chr','start','end','length')
write.table(new_length_bed,file=paste0(sample_name,"_new_length_bed.bed"),sep="\t",quote = F,col.names =F,row.names=F)

system(paste0('mkdir -p 777 intersect_bed'),intern=TRUE, wait=TRUE)
system(paste0('mkdir -p 777 length_count'),intern=TRUE, wait=TRUE)
system(paste0(bedtools,
' intersect -a ',tss_pro_bed_up1000_down1000,' -b ',sample_name,'_new_length_bed.bed',' -wao > intersect_bed/',sample_name,'_intersect_up1000_down1000.txt'),
intern=TRUE, wait=TRUE)
system(paste0(bedtools,
' intersect -a ',tss_pro_bed_up1000_up750,' -b ',sample_name,'_new_length_bed.bed',' -wao > intersect_bed/',sample_name,'_intersect_up1000_up750.txt'),
intern=TRUE, wait=TRUE)
system(paste0(bedtools,
' intersect -a ',tss_pro_bed_down750_down1000,' -b ',sample_name,'_new_length_bed.bed',' -wao > intersect_bed/',sample_name,'_intersect_down750_down1000.txt'),
intern=TRUE, wait=TRUE)
system(paste0('rm ',sample_name,"_new_length_bed.bed"),
intern=TRUE, wait=TRUE)
intersect_up1000_down1000<-read.table(file=paste0('intersect_bed/',sample_name,'_intersect_up1000_down1000.txt'),sep='\t',header=F)
colnames(intersect_up1000_down1000)<-c('intra_chr','intra_start','intra_end','tar_chr','tar_start','tar_end','tar_length','tar_overlap')
intersect_up1000_down1000$combine_name<-paste0(intersect_up1000_down1000$intra_chr,'_',intersect_up1000_down1000$intra_start,'_',intersect_up1000_down1000$intra_end)
length(table(intersect_up1000_down1000$combine_name))
intersect_up1000_up750<-read.table(file=paste0('intersect_bed/',sample_name,'_intersect_up1000_up750.txt'),sep='\t',header=F)
colnames(intersect_up1000_up750)<-c('intra_chr','intra_start','intra_end','tar_chr','tar_start','tar_end','tar_length','tar_overlap')
intersect_up1000_up750$combine_name<-paste0(intersect_up1000_up750$intra_chr,'_',intersect_up1000_up750$intra_start,'_',intersect_up1000_up750$intra_end)
length(table(intersect_up1000_up750$combine_name))
intersect_down750_down1000<-read.table(file=paste0('intersect_bed/',sample_name,'_intersect_down750_down1000.txt'),sep='\t',header=F)
colnames(intersect_down750_down1000)<-c('intra_chr','intra_start','intra_end','tar_chr','tar_start','tar_end','tar_length','tar_overlap')
intersect_down750_down1000$combine_name<-paste0(intersect_down750_down1000$intra_chr,'_',intersect_down750_down1000$intra_start,'_',intersect_down750_down1000$intra_end)
length(table(intersect_down750_down1000$combine_name))
save(intersect_up1000_down1000,intersect_up1000_up750,intersect_down750_down1000,file=paste0('intersect_bed/',sample_name,'_intersect_data.Rdata'))

#####建立100-300的长度矩阵，行为tss name，列为长度区间。
load(paste0('intersect_bed/',sample_name,'_intersect_data.Rdata'))
count_frame_up1000_down1000<-list()
tss_name<-names(table(intersect_up1000_down1000$combine_name))
for(i in 1:length(tss_name)){
count_frame_up1000_down1000[[i]]<-intersect_up1000_down1000[intersect_up1000_down1000$combine_name%in%tss_name[i],]
}
names(count_frame_up1000_down1000)<-tss_name
count_frame_up1000_up750<-list()
tss_name<-names(table(intersect_up1000_up750$combine_name))
for(i in 1:length(tss_name)){
count_frame_up1000_up750[[i]]<-intersect_up1000_up750[intersect_up1000_up750$combine_name%in%tss_name[i],]
}
names(count_frame_up1000_up750)<-tss_name
count_frame_down750_down1000<-list()
tss_name<-names(table(intersect_down750_down1000$combine_name))
for(i in 1:length(tss_name)){
count_frame_down750_down1000[[i]]<-intersect_down750_down1000[intersect_down750_down1000$combine_name%in%tss_name[i],]
print(c(i,'count_frame'))
}
names(count_frame_down750_down1000)<-tss_name

length_frame_up1000_down1000<-data.frame()
for(i in 1:length(count_frame_up1000_down1000)){
if(count_frame_up1000_down1000[[i]]$tar_length=='.'){
temp_count<-rep(0,300)
names(temp_count)<-c(1:300)
length_frame_up1000_down1000<-rbind(length_frame_up1000_down1000,t(temp_count))
}else{
temp_table<-table(count_frame_up1000_down1000[[i]]$tar_length)
temp_count<-rep(0,300)
names(temp_count)<-c(1:300)
temp_count[names(temp_table)]<-temp_table[names(temp_table)]
length_frame_up1000_down1000<-rbind(length_frame_up1000_down1000,t(temp_count))
}
}
rownames(length_frame_up1000_down1000)<-names(count_frame_up1000_down1000)
length_frame_up1000_up750<-data.frame()
for(i in 1:length(count_frame_up1000_up750)){
if(count_frame_up1000_up750[[i]]$tar_length=='.'){
temp_count<-rep(0,300)
names(temp_count)<-c(1:300)
length_frame_up1000_up750<-rbind(length_frame_up1000_up750,t(temp_count))
}else{
temp_table<-table(count_frame_up1000_up750[[i]]$tar_length)
temp_count<-rep(0,300)
names(temp_count)<-c(1:300)
temp_count[names(temp_table)]<-temp_table[names(temp_table)]
length_frame_up1000_up750<-rbind(length_frame_up1000_up750,t(temp_count))
}
}
rownames(length_frame_up1000_up750)<-names(count_frame_up1000_up750)
length_frame_down750_down1000<-data.frame()
for(i in 1:length(count_frame_down750_down1000)){
if(count_frame_down750_down1000[[i]]$tar_length=='.'){
temp_count<-rep(0,300)
names(temp_count)<-c(1:300)
length_frame_down750_down1000<-rbind(length_frame_down750_down1000,t(temp_count))
}else{
temp_table<-table(count_frame_down750_down1000[[i]]$tar_length)
temp_count<-rep(0,300)
names(temp_count)<-c(1:300)
temp_count[names(temp_table)]<-temp_table[names(temp_table)]
length_frame_down750_down1000<-rbind(length_frame_down750_down1000,t(temp_count))
}
}
rownames(length_frame_down750_down1000)<-names(count_frame_down750_down1000)
save(length_frame_up1000_down1000,length_frame_up1000_up750,length_frame_down750_down1000,file=paste0('length_count/',sample_name,'_length_frame.Rdata'))

#######计算样本熵，the longest distance from the center of the TSS −1 kbp to −750bp (upstream) and +750bp to +1 kbp (downstream)
calculate_ent <- function(x,range){x=x/sum(x);ent=-sum(x*log2(x));return(ent)}
raw_ent<-c()
for(i in 1:dim(length_frame_up1000_down1000)[1]){
raw_gene<-length_frame_up1000_down1000[i,]
temp_raw_ent<-calculate_ent(raw_gene)
raw_ent<-c(final_ent,temp_raw_ent)
print(i)
}
names(final_ent)<-rownames(length_frame_up1000_down1000)

#####背景矫正
load(paste0('length_count/',sample_name,'_length_frame.Rdata'))
length_frame_up1000_up750<-length_frame_up1000_up750[,100:300]
length_frame_down750_down1000<-length_frame_down750_down1000[,100:300]
length_frame_up1000_down1000<-length_frame_up1000_down1000[,100:300]
length_sample<-apply(rbind(length_frame_up1000_up750,length_frame_down750_down1000),2,sum)
sample_ent<-calculate_ent(length_sample)

###计算5个背景基因集的熵,每个基因集的基因个数为n
n=20
report_NBT<-read.table(file='all.tss.genes.canonical.ensembl75.txt',sep='\t',header=T)##negativeControl gene
load(tss_pro_table)
NC_tss<-tss_pro_table[tss_pro_table$gene_name%in%intersect(report_NBT[report_NBT$Category=='negativeControl',]$Gene.Symbol,tss_pro_table$gene_name),]
NC_tss$combine_name<-paste0(NC_tss$chromosome,'_',NC_tss$tss-1000,'_',NC_tss$tss+1000)
set.seed(2022)
index0<-sample(1:dim(NC_tss)[1],n,replace =T)
index1<-sample(1:dim(NC_tss)[1],n,replace =T)
index2<-sample(1:dim(NC_tss)[1],n,replace =T)
index3<-sample(1:dim(NC_tss)[1],n,replace =T)
index4<-sample(1:dim(NC_tss)[1],n,replace =T)
length_negcon0<-apply(length_frame_up1000_down1000[NC_tss$combine_name[index0],],2,sum) + 1e-5
length_negcon1<-apply(length_frame_up1000_down1000[NC_tss$combine_name[index1],],2,sum) + 1e-5
length_negcon2<-apply(length_frame_up1000_down1000[NC_tss$combine_name[index2],],2,sum) + 1e-5
length_negcon3<-apply(length_frame_up1000_down1000[NC_tss$combine_name[index3],],2,sum) + 1e-5
length_negcon4<-apply(length_frame_up1000_down1000[NC_tss$combine_name[index4],],2,sum) + 1e-5
baseents0 = calculate_ent(length_negcon0)
baseents1 = calculate_ent(length_negcon1)
baseents2 = calculate_ent(length_negcon2)
baseents3 = calculate_ent(length_negcon3)
baseents4 = calculate_ent(length_negcon4)
baseents = c(baseents0,baseents1,baseents2,baseents3,baseents4)

#####计算熵
library(gtools)
calculate_ent_bayesian <- function(counts,distfull, alpha0,n=250, baseents = NULL){
  countsnew = counts
  distfullnew = distfull
  distfullnew = distfullnew/sum(distfullnew)
  updatedalpha = 1+(distfullnew*alpha0+countsnew)
  if (sum(updatedalpha)<10)
  {
    updatedalpha = updatedalpha/sum(updatedalpha)*10
  }
  N = sum(counts)
  mydirichlet <- rdirichlet(n=n, alpha=as.numeric(updatedalpha)) + 1e-5
  mynewent = digamma(sum(updatedalpha)+1) - sum(updatedalpha/sum(updatedalpha)*digamma(updatedalpha+1))
  newnet = mynewent/log(2)
  base_ent = calculate_ent(distfullnew)
  gammanew <- function(arrayin){o=pgamma(arrayin,shape=.5,rate=1);return(o)}
  if (is.null(baseents))
  {
    allents = (apply(mydirichlet,1,function(x){x=x+1e-5;
    ent = calculate_ent(x/sum(x));return((ent/base_ent))}))- 1
    bayesian_sig = mean(gammanew(allents))
  }else{
    bayesian_sig = 0
    dirichletEnts = (apply(mydirichlet,1,function(x){x=x+1e-5;ent = calculate_ent(x/sum(x));return((ent))}))
    for (base_ent in baseents)
    {
      allents = dirichletEnts/base_ent - 1
      bayesian_sig = bayesian_sig+mean(gammanew(allents)) 
    }
    bayesian_sig = bayesian_sig/length(baseents)
  }
  return(bayesian_sig)
}
final_ent<-c()
for(i in 1:dim(length_frame_up1000_down1000)[1]){
temp_gene<-length_frame_up1000_down1000[i,]
temp_ent<-calculate_ent_bayesian(temp_gene,length_sample,20,n=200,baseents)
final_ent<-c(final_ent,temp_ent)
print(i)
}
names(final_ent)<-rownames(length_frame_up1000_down1000)
setwd(outputdir)
save(final_ent,file=paste0(sample_name,'_adjust_final_ent.Rdata'))



