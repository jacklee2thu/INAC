
###### INAC ML
argv <- commandArgs(TRUE)
sample_name<-as.character(argv[1])
input_dir<-as.character(argv[2])
output_dir<-as.character(argv[3])
cancer_sample<-as.character(argv[4])
healthy_sample<-as.character(argv[5])
seed<-as.character(argv[6])
divide_fraction<-as.character(argv[7])
choose_method<-as.character(argv[8])

######输入数据
options(stringsAsFactors=F)
setwd(input_dir)
load(sample_name)
gastric_sample<-cancer_sample
xiehe_sample<-healthy_sample


######gbm模型

library(tidyverse)
library(caret)
library(pROC)
library(gbm)

xiehe_index<-match(xiehe_sample,rownames(features.sl))
gastric_index<-match(gastric_sample,rownames(features.sl))
type_name<-rep(0,dim(features.sl)[1])
type_name[xiehe_index]<-'Healthy'
type_name[gastric_index]<-'Cancer'
type_name<-factor(type_name,levels=c('Cancer','Healthy'))

###训练集和测试集
set.seed(seed)
inTrain = createDataPartition(type_name, divide_fraction, list = FALSE)
trainx = features.sl[inTrain,]
testx = features.sl[-inTrain,]
trainy = type_name[inTrain]
testy = type_name[-inTrain]

####若变量过多则进行特征选择

subsets = c(20,50,100)####特征个数选择
ctrl= rfeControl(functions = rfFuncs,method = 'cv',verbose = FALSE,returnResamp = 'final')
Profile = rfe(trainx, trainy, sizes = subsets, rfeControl = ctrl)


##########
if(length(Profile$optVariables)>=101){
trainx=trainx[,Profile$optVariables[1:100]]
testx=testx[,Profile$optVariables[1:100]]
}else{
trainx=trainx[,Profile$optVariables]
testx=testx[,Profile$optVariables]
}


######模型训练
features<-trainx

features$type<-trainy

ctrl <- trainControl(method = "repeatedcv",
                     number = 10,
                     repeats = 10,
                     verboseIter = FALSE,
                     savePredictions=TRUE,
                     classProbs=TRUE,
#                      preProcOptions=list(thres = 0.90),
                     summaryFunction = twoClassSummary)
					 
gbmGrid = expand.grid(interaction.depth = 3,n.trees = 150,shrinkage = 0.1,n.minobsinnode = 10)
#gbmFit1 = train(features,method = 'gbm',trControl = ctrl,tuneGrid = gbmGrid,verbose = FALSE)
	
					 
					 
set.seed(1234)
model_gbm <- caret::train(type ~ .,
                               data = features,
                               method = choose_method,
							   tuneGrid=gbmGrid,
							   preProcess = c("corr", "nzv"),
                         trControl = ctrl)


pred.tbl <- model_gbm$pred %>%
group_by(rowIndex) %>% dplyr::summarize(obs=obs[1], Cancer=mean(Cancer))
pred.tbl$sample <- rownames(features)
combine_info<-cbind(rownames(features),features)
colnames(combine_info)[1]<-'sample'
pred.tbl <- inner_join(pred.tbl, combine_info)

## 90% specificity
cutoff <- (pred.tbl %>% filter(type=="Healthy") %>%
           arrange(desc(Cancer)))$Cancer[15]
## 95% specificity cutoff to be used in tissue prediction.
cutoff95 <- (pred.tbl %>% filter(type=="Healthy") %>%
             arrange(desc(Cancer)))$Cancer[9]

pred.tbl <- pred.tbl %>%
    mutate(detected90 = ifelse(Cancer > cutoff, "Cancer", "Healthy"),
       detected95 = ifelse(Cancer > cutoff95, "Cancer", "Healthy"))
	   

###
library(ROCR)
pred.tbl$lable=ifelse(pred.tbl$obs=='Cancer',yes=1,0)
pred1 = prediction(pred.tbl$Cancer,pred.tbl$lable)
perf1 = performance(pred1, measure='tpr', x.measure='fpr')
auc_ROCR <-performance(pred1,measure ="auc")
auc_ROCR@y.values[[1]]
train_standard<-cbind(pred.tbl,rownames(trainx))[,c(1,2,7)]

library(ggplot2)
df <- data.frame(Curve=as.factor(rep(c(1), each=length(perf1@x.values[[1]]))), 
                 FalsePositive=c(perf1@x.values[[1]]),
                 TruePositive=c(perf1@y.values[[1]]))
setwd(output_dir)
pdf('train_ROC.pdf')
ggplot(df, aes(x=FalsePositive, y=TruePositive, color=Curve)) +ggtitle(paste0("ROC Curve  AUC=", auc_ROCR@y.values[[1]]))+ geom_line()+theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position='none')
dev.off()

######测试集

models<-list(model_gbm)
predValues = extractPrediction(models,testX = testx, testY = testy)
testValues = subset(predValues, dataType == 'Test')
probValues = extractProb(models,testX = testx, testY = testy)
testProbs = subset(probValues, dataType == 'Test')
pred.tbl = subset(testValues, model == 'gbm')
confusion_matrix<-confusionMatrix(pred.tbl$pred, pred.tbl$obs)
save(confusion_matrix,file='confusion_matrix_test_data.Rdata')

library(ROCR)
pred.tbl$lable=ifelse(pred.tbl$obs=='Cancer',yes=1,0)
pred1 = prediction(testProbs$Cancer,pred.tbl$lable)
perf1 = performance(pred1, measure='tpr', x.measure='fpr')
auc_ROCR <-performance(pred1,measure ="auc")

library(ggplot2)
df <- data.frame(Curve=as.factor(rep(c(1), each=length(perf1@x.values[[1]]))), 
                 FalsePositive=c(perf1@x.values[[1]]),
                 TruePositive=c(perf1@y.values[[1]]))
pdf('test_ROC.pdf')
ggplot(df, aes(x=FalsePositive, y=TruePositive, color=Curve))+ggtitle(paste0("ROC Curve  AUC=", auc_ROCR@y.values[[1]])) + geom_line()+theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position='none')
dev.off()

