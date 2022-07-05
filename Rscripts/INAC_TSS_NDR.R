
###### INAC TSS NDR
argv <- commandArgs(TRUE)
sample_name<-as.character(argv[1])
input_dir<-as.character(argv[2])
output_dir<-as.character(argv[3])
samtools<-as.character(argv[4])
tss_bed_up3000_1000<-as.character(argv[6])
tss_bed_down1000_3000<-as.character(argv[7])
tss_bed_up1000_down1000<-as.character(argv[8])
tss_bed_up150_down50<-as.character(argv[9])
tss_pro_table<-as.character(argv[10])


options(stringsAsFactors=F)
##跑出每个人要转录的TSS

setwd(input_dir)
if(is.na(match(paste0(sample_name,".bam.bai"),dir()))){
system(paste0(samtools,' index ',sample_name,".bam"))
}

tss_bed_up3000_1000<-read.table(file=tss_bed_up3000_1000,sep="\t",header=F)
tss_bed_down1000_3000<-read.table(file=tss_bed_down1000_3000,sep="\t",header=F)
tss_bed_up1000_down1000<-read.table(file=tss_bed_up1000_down1000,sep="\t",header=F)
tss_bed_up150_down50<-read.table(file=tss_bed_up150_down50,sep="\t",header=F)

patient_tss1000<-list()
patient_tss150<-list()
for(i in 1:dim(tss_bed_up3000_1000)[1]){
##上游3000-1000和下游1000-3000用于做金标准
temp_file<-system(paste0(samtools,
" depth -q 30"," -r ",tss_bed_up3000_1000[i,1],":",tss_bed_up3000_1000[i,2],"-",tss_bed_up3000_1000[i,3]," ",sample_name,".bam"),intern=TRUE, wait=TRUE)

if(length(temp_file)==0){
temp_up3000_1000<-NA
}else{
temp_up3000_1000<-data.frame()
for(j in 1:length(temp_file)){
temp_frame<-unlist(strsplit(temp_file[j],'\t'))
temp_up3000_1000<-rbind(temp_up3000_1000,temp_frame)
}
colnames(temp_up3000_1000)<-c('chr','site','coverage')
temp_up3000_1000[,3]<-as.numeric(temp_up3000_1000[,3])
}

temp_file<-system(paste0(samtools,
" depth -q 30"," -r ",tss_bed_down1000_3000[i,1],":",tss_bed_down1000_3000[i,2],"-",tss_bed_down1000_3000[i,3]," ",sample_name,".bam"),intern=TRUE, wait=TRUE)
if(length(temp_file)==0){
temp_down1000_3000<-NA
}else{
temp_down1000_3000<-data.frame()
for(j in 1:length(temp_file)){
temp_frame<-unlist(strsplit(temp_file[j],'\t'))
temp_down1000_3000<-rbind(temp_down1000_3000,temp_frame)
}
colnames(temp_down1000_3000)<-c('chr','site','coverage')
temp_down1000_3000[,3]<-as.numeric(temp_down1000_3000[,3])
}

if(unique(is.na(temp_down1000_3000))&unique(is.na(temp_up3000_1000))){
normal_tss_value<-NA
}
if(unique(is.na(temp_up3000_1000))&unique(!is.na(temp_down1000_3000))){
normal_tss_value<-mean(temp_down1000_3000[,3])
}
if(unique(is.na(temp_down1000_3000))&unique(!is.na(temp_up3000_1000))){
normal_tss_value<-mean(temp_up3000_1000[,3])
}

if(unique(!is.na(temp_down1000_3000))&unique(!is.na(temp_up3000_1000))){
normal_tss<-rbind(temp_up3000_1000,temp_down1000_3000)
normal_tss_value<-mean(normal_tss[,3])
}


##上游1000-下游1000
temp_file<-system(paste0(samtools,
" depth -q 30"," -r ",tss_bed_up1000_down1000[i,1],":",tss_bed_up1000_down1000[i,2],"-",tss_bed_up1000_down1000[i,3]," ",sample_name,".bam"),intern=TRUE, wait=TRUE)
if(length(temp_file)==0){
temp_up1000_down1000<-NA
}else{
temp_up1000_down1000<-data.frame()
for(j in 1:length(temp_file)){
temp_frame<-unlist(strsplit(temp_file[j],'\t'))
temp_up1000_down1000<-rbind(temp_up1000_down1000,temp_frame)
}
temp_up1000_down1000[,3]<-as.numeric(temp_up1000_down1000[,3])
temp_up1000_down1000[,4]<-temp_up1000_down1000[,3]/normal_tss_value
colnames(temp_up1000_down1000)<-c('chr','site','coverage','mean_coverage')
}

##上游150-下游50
temp_file<-system(paste0(samtools,
" depth -q 30"," -r ",tss_bed_up150_down50[i,1],":",tss_bed_up150_down50[i,2],"-",tss_bed_up150_down50[i,3]," ",sample_name,".bam"),intern=TRUE, wait=TRUE)
if(length(temp_file)==0){
temp_up150_down50<-NA
}else{
temp_up150_down50<-data.frame()
for(j in 1:length(temp_file)){
temp_frame<-unlist(strsplit(temp_file[j],'\t'))
temp_up150_down50<-rbind(temp_up150_down50,temp_frame)
}
temp_up150_down50[,3]<-as.numeric(temp_up150_down50[,3])
temp_up150_down50[,4]<-temp_up150_down50[,3]/normal_tss_value
colnames(temp_up150_down50)<-c('chr','site','coverage','mean_coverage')
}

patient_tss1000[[i]]<-temp_up1000_down1000
patient_tss150[[i]]<-temp_up150_down50
print(i)
}
setwd(output_dir)

###
load(tss_pro_table)
all_sample_tss<-list()
for(j in 1:length(sample_name)){
load(paste0(sample_name,"_alltss.Rdata"))
tss_1000<-c()
tss_150<-c()

for(i in 1:length(patient_tss1000)){
if(!is.na(patient_tss1000[[i]])){
tss_1000[i]<-mean(patient_tss1000[[i]][,4])
}else{tss_1000[i]<-NA}
if(!is.na(patient_tss150[[i]])){
tss_150[i]<-mean(patient_tss150[[i]][,4])
}else{tss_150[i]<-NA}
}
tss_relative_coverage<-data.frame(tss_150,tss_1000)
colnames(tss_relative_coverage)<-c("tss_150","tss_1000")
all_sample_tss[[j]]<-tss_relative_coverage
}
names(all_sample_tss)<-sample_name


###########
temp_tss_NDR<-all_sample_tss[[1]][,1]
save(temp_tss_NDR,file=paste0(sample_name,'_tss_NDR.Rdata'))



