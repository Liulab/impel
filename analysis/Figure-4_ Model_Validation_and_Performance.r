#####################################################################################################
###Figure 4A Kaplan-Meier curves were plotted to analyze the correlations between three groups      # 
###and DFS in GSE39279.                                                                             #
###Figure 4B Kaplan-Meier curves were plotted to analyze the correlations between two groups        #
###and DFS in GSE39279.                                                                             #
###Figure 4C Kaplan-Meier curves were plotted to analyze the correlations between three groups      #
###and DFS in TCGA.                                                                                 #
###Figure 4D Kaplan-Meier curves were plotted to analyze the correlations between two groups        #
###and DFS in TCGA.                                                                                 #
###Figure 4E-F T The AUC of model with eight methylations was higher than those constructed         #
###by less features in GSE39279 and TCGA.                                                           #
###Figure 4G ROC curves were generated to predict the recurrence of lung cancer in validation data. #
#####################################################################################################

############
##Figure 4A#
############
options(stringsAsFactors=F)
options(scipen=200)
library(glmnet)
library(pROC)
library(survival)
library("survminer")
load("/data1/dengyl/project/S1/chennan/new3/GSE39279_mt100.RData")
load("/data1/dengyl/project/S1/chennan/new3/labelGSE39279.RData")
load("/data1/dengyl/project/S1/chennan/new3/cliniGSE39279.RData")
load("cvfit.lasso.RData")
tcgaMtF <- GSE39279_mt100[rownames(model)[-1],]
lasso_C <- predict(cvfit.lasso, newx = t(tcgaMtF),s = "lambda.min",type = "response")
lasso_C  <- lasso_C [,1]
names(lasso_C) <- colnames(tcgaMtF)
cliniGSE39279 <- cliniGSE39279[match(names(lasso_C),cliniGSE39279[,"Sample_geo_accession"]),]
label <- rep(0,nrow(cliniGSE39279))
label[cliniGSE39279$"recurrence"==" yes"] <- 1
y <-Surv(as.numeric(cliniGSE39279$"time"),label)
Clusters <-rep(2,length(lasso_C))
Clusters[lasso_C>=0.6] <- 3
Clusters[lasso_C<=0.4] <- 1
Clusters <- factor(as.character(Clusters),levels=c("1","2","3"))
cliniT <- data.frame(recurrenceTime=as.numeric(cliniGSE39279$"time"),
label=label,Clusters=Clusters,stringsAsFactors=F)
Sur_Curve<-survfit(Surv(recurrenceTime,label)~Clusters,data = cliniT)
Diff_Survival<-survdiff(Surv(as.numeric(cliniGSE39279$"time"),label)~Clusters)
P_Value<-pchisq(Diff_Survival$chisq,length(table(Clusters))-1,lower.tail=F)#0.003807353
ggsurvplot(Sur_Curve,data=cliniT,pvalue=TRUE,surv.median.line="hv",
risk.table = T,break.x.by=2,ylab="Disease-Free Survival")

############
##Figure 4B#
############
options(stringsAsFactors=F)
options(scipen=200)
library(glmnet)
library(pROC)
library(RColorBrewer)
library(ROCR)
library(survival)
library("survminer")
library(survcomp) 
load("/data1/dengyl/project/S1/chennan/new3/GSE39279_mt100.RData")
load("/data1/dengyl/project/S1/chennan/new3/labelGSE39279.RData")
load("/data1/dengyl/project/S1/chennan/new3/cliniGSE39279.RData")
load("cvfit.lasso.RData")
tcgaMtF <- GSE39279_mt100[rownames(model)[-1],]
lasso_C <- predict(cvfit.lasso, newx = t(tcgaMtF),s = "lambda.min",type = "response")
lasso_C  <- lasso_C [,1]
lasso_C <- lasso_C[(lasso_C>=0.6)|(lasso_C<=0.4)]
names(lasso_C) <- colnames(tcgaMtF)
cliniGSE39279 <- cliniGSE39279[match(names(lasso_C),cliniGSE39279[,"Sample_geo_accession"]),]
label <- rep(0,nrow(cliniGSE39279))
label[cliniGSE39279$"recurrence"==" yes"] <- 1
y <-Surv(as.numeric(cliniGSE39279$"time"),label)
Clusters <-rep(1,length(lasso_C))
Clusters[lasso_C>=0.6] <- 2
Clusters <- factor(as.character(Clusters),levels=c("1","2"))
cliniT <- data.frame(recurrenceTime=as.numeric(cliniGSE39279$"time"),
label=label,Clusters=Clusters,stringsAsFactors=F)
Sur_Curve<-survfit(Surv(recurrenceTime,label)~Clusters,data = cliniT)
Diff_Survival<-survdiff(Surv(as.numeric(cliniGSE39279$"time"),label)~Clusters)
P_Value<-pchisq(Diff_Survival$chisq,length(table(Clusters))-1,lower.tail=F)#0.004
ggsurvplot(Sur_Curve,data=cliniT,pvalue=TRUE,surv.median.line="hv",
risk.table = T,break.x.by=2,ylab="Disease-Free Survival")


############
##Figure 4C#
############
options(stringsAsFactors=F)
options(scipen=200)
library(GenomicRanges)
library(glmnet)
library(pROC)
library(sva)
library(survival)
library("survminer")
load("/data1/dengyl/project/S1/chennan/new3/labeltcga.RData")
load("/data1/dengyl/project/S1/chennan/new3/tcga_mt100.RData")
load("/data1/dengyl/project/S1/chennan/new3/cliniTCGA.RData")
load("cvfit.lasso.RData")
xiao_patient <- read.table(file="/data1/dengyl/project/S1/chennan/new2/nationwidechildrens.org_clinical_patient_luad.txt",
sep="\t",quote = "",stringsAsFactors=F,header=T)
xiao_drug <- read.table(file="/data1/dengyl/project/S1/chennan/new2/nationwidechildrens.org_clinical_drug_luad.txt",
sep="\t",quote = "",stringsAsFactors=F,header=T)
tcgaMtF <- tcga_mt100[rownames(model)[-1],]
lasso_C <- predict(cvfit.lasso, newx = t(tcgaMtF),s = "lambda.min",type = "response")
lasso_C  <- lasso_C [,1]
names(lasso_C) <- colnames(tcgaMtF)
cliniTCGA <- cliniTCGA[match(names(lasso_C),cliniTCGA[,1]),]
label <- rep(0,nrow(cliniTCGA))
label[(cliniTCGA$"recurrence"%in%c("Distant Metastasis","Locoregional Recurrence"))] <- 1
y <-Surv(as.numeric(cliniTCGA$"recurrenceTime"),label)
Clusters <-rep(2,length(lasso_C))
Clusters[lasso_C>=0.6] <- 3
Clusters[lasso_C<=0.4] <- 1
Clusters <- factor(as.character(Clusters),levels=c("1","2","3"))
xiao_patient[,1] <- substr(xiao_patient[,1],9,12)
smoke <- rep(1,length(lasso_C))
names(smoke) <- names(lasso_C)
smoke[names(smoke)%in%xiao_patient[xiao_patient[,"tobacco_smoking_history_indicator"]=="Lifelong Non-smoker",1]] <- 0
xiao_drug[,1] <- substr(xiao_drug[,1],9,12)
drug <- rep(0,length(lasso_C))
names(drug) <- names(lasso_C)
drug[names(drug)%in%xiao_drug[xiao_drug[,"pharmaceutical_therapy_type"]=="Chemotherapy",1]] <- 1
cliniT <- data.frame(recurrenceTime=as.numeric(cliniTCGA$"recurrenceTime"),
label=label,Clusters=Clusters,stringsAsFactors=F)
Sur_Curve<-survfit(Surv(recurrenceTime,label)~Clusters,data = cliniT)
print(Sur_Curve)
summary(Sur_Curve)
Diff_Survival<-survdiff(Surv(as.numeric(cliniTCGA$"recurrenceTime"),label)~Clusters)
P_Value<-pchisq(Diff_Survival$chisq,length(table(Clusters))-1,lower.tail=F)#0.032
ggsurvplot(Sur_Curve,data=cliniT,pvalue=TRUE,surv.median.line="hv",
risk.table = T,break.x.by=2,ylab="Disease-Free Survival")


############
##Figure 4D#
############
options(stringsAsFactors=F)
options(scipen=200)
library(glmnet)
library(pROC)
library(RColorBrewer)
library(ROCR)
library(survival)
library("survminer")
library(survcomp) 
load("/data1/dengyl/project/S1/chennan/new3/labeltcga.RData")
load("/data1/dengyl/project/S1/chennan/new3/tcga_mt100.RData")
load("/data1/dengyl/project/S1/chennan/new3/cliniTCGA.RData")
load("cvfit.lasso.RData")
tcgaMtF <- tcga_mt100[rownames(model)[-1],]
lasso_C <- predict(cvfit.lasso, newx = t(tcgaMtF),s = "lambda.min",type = "response")
lasso_C  <- lasso_C [,1]
names(lasso_C) <- colnames(tcgaMtF)
lasso_C <- lasso_C[(lasso_C>=0.6)|(lasso_C<=0.4)]
cliniTCGA <- cliniTCGA[match(names(lasso_C),cliniTCGA[,1]),]
label <- rep(0,nrow(cliniTCGA))
label[(cliniTCGA$"recurrence"%in%c("Distant Metastasis","Locoregional Recurrence"))] <- 1
y <-Surv(as.numeric(cliniTCGA$"recurrenceTime"),label)
Clusters <-rep(1,length(lasso_C))
Clusters[lasso_C>=0.6] <- 2
Clusters <- factor(as.character(Clusters),levels=c("1","2"))
cliniT <- data.frame(recurrenceTime=as.numeric(cliniTCGA$"recurrenceTime"),
label=label,Clusters=Clusters,stringsAsFactors=F)
Sur_Curve<-survfit(Surv(recurrenceTime,label)~Clusters,data = cliniT)
Diff_Survival<-survdiff(Surv(as.numeric(cliniTCGA$"recurrenceTime"),label)~Clusters)
P_Value<-pchisq(Diff_Survival$chisq,length(table(Clusters))-1,lower.tail=F)#0.009990579
ggsurvplot(Sur_Curve,data=cliniT,pvalue=TRUE,surv.median.line="hv",
risk.table = T,break.x.by=2,ylab="Disease-Free Survival") 

############
##Figure 4E#
############
#subsample
options(stringsAsFactors=F)
options(scipen=200)
library(glmnet)
library(pROC)
library(ROCR)
load("/data1/dengyl/project/S1/chennan/new3/labelGSE39279.RData")
load("/data1/dengyl/project/S1/chennan/new3/GSE39279_mt100.RData")
load("cvfit.lasso.RData")
load("/data1/dengyl/project/S1/chennan/new3/cliniGSE39279.RData")
regionAll <- rownames(model)[-1]
pos1 <- rep(c(1,0),128)
pos2 <- rep(c(1,1,0,0),64)
pos3 <- rep(c(1,1,1,1,0,0,0,0),32)
pos4 <- rep(c(1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0),16)
pos5 <- rep(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),8)
pos6 <- rep(c(rep(1,32),rep(0,32)),4)
pos7 <- rep(c(rep(1,64),rep(0,64)),2)
pos8 <- rep(c(1,0),c(128,128))
para_mt <- cbind(pos1,cbind(pos2,cbind(pos3,cbind(pos4,cbind(pos5,cbind(pos6,cbind(pos7,pos8)))))))
para_mt <- para_mt[apply(para_mt,1,sum)>1,]
per1000 <- lapply(1:1000,function(y){
	res_mt <-t(sapply(seq(nrow(para_mt)),function(x){ 
		tmp <- tran_fun2(mt_tr=m100_mt[regionAll[as.logical(para_mt[x,])],],
			label_tr=label,
			mt_our=our_mt100,
			label_our=label_our,
			mt_tcga=tcga_mt100,
			label_tcga=labeltcga,
			mt_gse=GSE39279_mt100,
			labe_gse=labelGSE39279,
			seed=y)
		return(tmp)	
	}))
	print(y)
	return(res_mt)
})
model_gse <- sapply(per1000,function(x) x[1,3])
p_gse <- sapply(1:246,function(y){
	model_our_sub1 <- sapply(per1000,function(x) x[y+1,3])
	return(t.test(model_our_sub1,model_gse,alternative="less")$"p.value")
})

#plot
options(scipen=200)
library(glmnet)
library(pROC)
library(ROCR)
library(RColorBrewer)
load("/data1/dengyl/project/S1/chennan/new3/labelGSE39279.RData")
load("/data1/dengyl/project/S1/chennan/new3/GSE39279_mt100.RData")
load("cvfit.lasso.RData")
load("per1000.RData")
model_gse <- sapply(per1000,function(x) x[1,3])
p_gse <- sapply(1:246,function(y){
	model_our_sub1 <- sapply(per1000,function(x) x[y+1,3])
	return(t.test(model_our_sub1,model_gse,alternative="less")$"p.value")
})
p_gse_list <- lapply(1:246,function(y){
	model_our_sub1 <- sapply(per1000,function(x) x[y+1,3])
	return(model_our_sub1)
}) 
pos1 <- rep(c(1,0),128)
pos2 <- rep(c(1,1,0,0),64)
pos3 <- rep(c(1,1,1,1,0,0,0,0),32)
pos4 <- rep(c(1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0),16)
pos5 <- rep(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),8)
pos6 <- rep(c(rep(1,32),rep(0,32)),4)
pos7 <- rep(c(rep(1,64),rep(0,64)),2)
pos8 <- rep(c(1,0),c(128,128))
para_mt <- cbind(pos1,cbind(pos2,cbind(pos3,cbind(pos4,cbind(pos5,cbind(pos6,cbind(pos7,pos8)))))))
para_mt <- para_mt[apply(para_mt,1,sum)>1,]
colo <- brewer.pal(8, "Reds")[2:7]
p_our_index <- p_gse_list[order(sapply(p_gse_list,median),decreasing=T)]
boxplot(p_our_index,lwd=0.1,
col=colo[apply(para_mt,1,sum)[order(sapply(p_gse_list,median),decreasing=T)]],
main="AUC in GSE39279",outline = F,ylim=c(0.54,0.72))
abline(h=0.71,col=brewer.pal(8, "Reds")[8],lwd=2)

############
##Figure 4F#
############
options(stringsAsFactors=F)
options(scipen=200)
library(glmnet)
library(pROC)
library(ROCR)
load("/data1/dengyl/project/S1/chennan/new3/labeltcga.RData")
load("cvfit.lasso.RData")
load("/data1/dengyl/project/S1/chennan/new3/cliniTCGA.RData")
regionAll <- rownames(model)[-1]
pos1 <- rep(c(1,0),128)
pos2 <- rep(c(1,1,0,0),64)
pos3 <- rep(c(1,1,1,1,0,0,0,0),32)
pos4 <- rep(c(1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0),16)
pos5 <- rep(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),8)
pos6 <- rep(c(rep(1,32),rep(0,32)),4)
pos7 <- rep(c(rep(1,64),rep(0,64)),2)
pos8 <- rep(c(1,0),c(128,128))
para_mt <- cbind(pos1,cbind(pos2,cbind(pos3,cbind(pos4,cbind(pos5,cbind(pos6,cbind(pos7,pos8)))))))
para_mt <- para_mt[apply(para_mt,1,sum)>1,]
per1000 <- lapply(1:1000,function(y){
	res_mt <-t(sapply(seq(nrow(para_mt)),function(x){ 
		tmp <- tran_fun2(mt_tr=m100_mt[regionAll[as.logical(para_mt[x,])],],
			label_tr=label,
			mt_our=our_mt100,
			label_our=label_our,
			mt_tcga=tcga_mt100,
			label_tcga=labeltcga,
			mt_gse=GSE39279_mt100,
			labe_gse=labelGSE39279,
			seed=y)
		return(tmp)	
	}))
	print(y)
	return(res_mt)
})

model_tcga <- sapply(per1000,function(x) x[1,2])
p_tcga <- sapply(1:246,function(y){
	model_our_sub1 <- sapply(per1000,function(x) x[y+1,2])
	return(t.test(model_our_sub1,model_tcga,alternative="less")$"p.value")
})
#plot
options(stringsAsFactors=F)
options(scipen=200)
library(glmnet)
library(pROC)
library(ROCR)
library(RColorBrewer)
load("/data1/dengyl/project/S1/chennan/new3/tcga_mt100.RData")
load("cvfit.lasso.RData")
load("/data1/dengyl/project/S1/chennan/new3/cliniTCGA.RData")
load("per1000.RData")
model_tcga <- sapply(per1000,function(x) x[1,2])
p_tcga <- sapply(1:246,function(y){
	model_our_sub1 <- sapply(per1000,function(x) x[y+1,2])
	return(t.test(model_our_sub1,model_tcga,alternative="less")$"p.value")
})
p_tcga_list <- lapply(1:246,function(y){
	model_our_sub1 <- sapply(per1000,function(x) x[y+1,2])
	return(model_our_sub1)
}) 
pos1 <- rep(c(1,0),128)
pos2 <- rep(c(1,1,0,0),64)
pos3 <- rep(c(1,1,1,1,0,0,0,0),32)
pos4 <- rep(c(1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0),16)
pos5 <- rep(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),8)
pos6 <- rep(c(rep(1,32),rep(0,32)),4)
pos7 <- rep(c(rep(1,64),rep(0,64)),2)
pos8 <- rep(c(1,0),c(128,128))
para_mt <- cbind(pos1,cbind(pos2,cbind(pos3,cbind(pos4,cbind(pos5,cbind(pos6,cbind(pos7,pos8)))))))
para_mt <- para_mt[apply(para_mt,1,sum)>1,]
colo <- brewer.pal(8, "Reds")[2:7]
p_our_index <- p_tcga_list[order(sapply(p_tcga_list,median),decreasing=T)]
boxplot(p_our_index,lwd=0.1,
col=colo[apply(para_mt,1,sum)[order(sapply(p_tcga_list,median),decreasing=T)]],
main="AUC in TCGA",outline = F,ylim=c(0.5,0.67))
abline(h=0.66,col=brewer.pal(8, "Reds")[8],lwd=2)

############
##Figure 4G#
############
options(stringsAsFactors=F)
options(scipen=200)
library(GenomicRanges)
library(glmnet)
library(gplots)
library(pROC)
library(RColorBrewer)
library(ROCR)
colo=brewer.pal(4, "Set1")[1:4]
load("cvfit.lasso.RData")
load("/data1/dengyl/project/S1/chennan/new2/labeltcga.RData")
load("/data1/dengyl/project/S1/chennan/new2/tcga_mt100.RData")
load("/data1/dengyl/project/S1/chennan/new2/cliniTCGA.RData")
load("/data1/dengyl/project/S1/chennan/new2/GSE39279_mt100.RData")
load("/data1/dengyl/project/S1/chennan/new2/labelGSE39279.RData")
load("/data1/dengyl/project/S1/chennan/new2/cliniGSE39279.RData")
tcgaMtF <- tcga_mt100[rownames(model)[-1],]
lasso_C1 <- predict(cvfit.lasso, newx = t(tcgaMtF),s = "lambda.min",type = "response")
lasso_C1  <- lasso_C1[,1]
names(lasso_C1) <- colnames(tcgaMtF)
lasso_C1 <- lasso_C1[(lasso_C1>=0.6)|(lasso_C1<0.4)]
cliniTCGA <- cliniTCGA[match(names(lasso_C1),cliniTCGA[,1]),]
juan1 <- read.table(file="/data1/syt/project/00.Lung_meth/02.bismark/07.450k/endversion/stage/juan.tcga.txt",
sep="\t",quote = "",stringsAsFactors=F,header=T)
juan_score1 <- apply(juan1[,1:5],1,function(x) as.numeric(sum(unlist(x))>1)+1)
names(juan_score1) <- substr(juan1[,6],9,12)
juan_score1 <- juan_score1[names(lasso_C1)]
label1 <- rep(0,nrow(cliniTCGA))
label1[(cliniTCGA$"recurrence"%in%c("Distant Metastasis","Locoregional Recurrence"))] <- 1
tcgaMtF <- GSE39279_mt100[rownames(model)[-1],]
lasso_C2 <- predict(cvfit.lasso, newx = t(tcgaMtF),s = "lambda.min",type = "response")
lasso_C2  <- lasso_C2 [,1]
names(lasso_C2) <- colnames(tcgaMtF)
lasso_C2 <- lasso_C2[(lasso_C2>=0.6)|(lasso_C2<0.4)]
cliniGSE39279 <- cliniGSE39279[match(names(lasso_C2),cliniGSE39279[,"Sample_geo_accession"]),]
juan2 <- read.table(file="/data1/syt/project/00.Lung_meth/02.bismark/07.450k/endversion/stage/juan.gse.txt",
sep="\t",quote = "",stringsAsFactors=F,header=T)
score2 <- apply(juan2[,1:5],1,function(x) as.numeric(sum(unlist(x))>1)+1)
names(score2) <- juan2[,6]
score2 <- score2[names(lasso_C2)]
names(labelGSE39279) <- colnames(GSE39279_mt100)
labelGSE39279 <- labelGSE39279[names(lasso_C2)]
label2 <- rep(0,nrow(cliniGSE39279))
label2[cliniGSE39279$"recurrence"==" yes"] <- 1
pred_list1 <- prediction( c(lasso_C1,lasso_C2), c(label1,label2) )
roc1 <- roc(c(label1,label2),c(lasso_C1,lasso_C2), percent=TRUE,ci=TRUE, boot.n=100, ci.alpha=0.9)
pred_list2 <- prediction( c(as.numeric(cliniTCGA$"stage"=="Stage II"),
as.numeric(cliniGSE39279$"Stage"==" II")), c(label1,label2) )
roc2 <- roc(c(label1,label2),c(as.numeric(cliniTCGA$"stage"=="Stage II"),
as.numeric(cliniGSE39279$"Stage"==" II")), percent=TRUE,ci=TRUE, boot.n=100, ci.alpha=0.9)
pred_list3 <- prediction( c(juan_score1,score2), c(label1,label2) )
roc3 <- roc(c(label1,label2),c(juan_score1,score2),
percent=TRUE,ci=TRUE, boot.n=100, ci.alpha=0.9)
#plot
perf <- performance( pred_list1, "tpr", "fpr" )
plot( perf,col=colo[1],main="validation_method_compare_auc",lwd=2)
auc2 <- performance( pred_list1, "auc")@"y.values"[[1]]
text(0.5,0.5,label=paste0("model auc=",substr(as.character(auc2),1,5)),col = colo[1])
perf <- performance( pred_list2, "tpr", "fpr" )
plot(perf, add = TRUE,col=colo[2],lwd=2)
auc2 <- performance( pred_list2, "auc")@"y.values"[[1]]
text(0.5,0.5-0.05*1,label=paste0("Stage auc=",substr(as.character(auc2),1,5)),col = colo[2])
perf <- performance( pred_list3, "tpr", "fpr" )
plot(perf, add = TRUE,col=colo[3],lwd=2)
auc2 <- performance( pred_list3, "auc")@"y.values"[[1]]
text(0.5,0.5-0.05*2,label=paste0("Sandoval auc=",substr(as.character(auc2),1,5)),col = colo[3])
