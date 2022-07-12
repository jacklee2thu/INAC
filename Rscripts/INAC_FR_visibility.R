
###### INAC fragment size ratio
argv <- commandArgs(TRUE)
input_dir<-as.character(argv[1])
output_dir<-as.character(argv[2])
bin<-as.character(argv[3])


options(stringsAsFactors=F)
setwd(input_dir)

##5Mb为一个bin统计量
options(stringsAsFactors=F)
bin_gc<-read.table(file=bin,sep="\t",header=T)
gc.correct <- function(coverage, bias) {
    i <- seq(min(bias, na.rm=TRUE), max(bias, na.rm=TRUE), by = 0.001)####bias的范围超出了i的范围，导致出现NA值
    coverage.trend <- loess(coverage ~ bias)
    coverage.model <- loess(predict(coverage.trend, i) ~ i)
    coverage.pred <- predict(coverage.model, bias)
    coverage.corrected <- coverage - coverage.pred + median(coverage)
}

##integrate summary_table 
sample_name<-dir()
last_sample<-c()
for(i in 1:length(sample_name)){
temp_name<-strsplit(sample_name,".R")
last_sample<-c(last_sample,temp_name[[i]][1])
}
gastric_table<-data.frame()
for(j in 1:length(sample_name)){
load(sample_name[j])
summary_table$ratio.center<-scale(summary_table$ratio,scale = F)
temp_table<-data.frame(rep(last_sample[j],504),summary_table,bin_gc$gc)
temp_table$short.corrected=gc.correct(temp_table$short, temp_table$bin_gc.gc)
temp_table$short.corrected_fre<-temp_table$short.corrected/sum(na.omit(temp_table$short.corrected))
temp_table$long.corrected=gc.correct(temp_table$long, temp_table$bin_gc.gc)
temp_table$long.corrected_fre<-temp_table$long.corrected/sum(na.omit(temp_table$long.corrected))
temp_table$nfrags.corrected=gc.correct(temp_table$short+temp_table$long, temp_table$bin_gc.gc)
temp_table$nfrags.corrected_fre<-temp_table$nfrags.corrected/sum(na.omit(temp_table$nfrags.corrected))
temp_table$ratio.corrected=gc.correct(temp_table$ratio, temp_table$bin_gc.gc)

gastric_table<-rbind(gastric_table,temp_table)
}
colnames(gastric_table)[1]<-"sample_name"
setwd(output_dir)
save(gastric_table,file="all_table.Rdata")

###
load('all_table.Rdata')
all_summary_table<-gastric_table
all_summary_table<-cbind(c(rep("gastric",dim(gastric_table)[1])),all_summary_table)
colnames(all_summary_table)[1]<-"type"
all_summary_table<-as.data.frame(all_summary_table)

##ratio&ratio.center
library(ggplot2)
pdf('INAC_FR_visibility.pdf',width=70,height=10)
gg <- ggplot(all_summary_table, aes(x=bin, y=ratio, group=sample_name,color=type)) +
scale_color_manual(breaks = c("gastric","healthy"),values=c("#D45252","#6EA6D9"))+geom_line(size=0.001)
gg+theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
axis.title.x =element_text(size=20), axis.title.y=element_text(size=20),axis.text.x=element_text(size=20),
axis.text.y=element_text(size=20),legend.text=element_text(size=20))+labs(title = "fragment ratio")


gg <- ggplot(gastric_table, aes(x=bin, y=ratio.center, group=sample_name)) + geom_line( size=0.001,color="#D45252")
gg+theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
axis.title.x =element_text(size=20), axis.title.y=element_text(size=20),axis.text.x=element_text(size=20),
axis.text.y=element_text(size=20),legend.text=element_text(size=20))+labs(title = "fragment ratio.center of gastric cancer patients")

dev.off()

