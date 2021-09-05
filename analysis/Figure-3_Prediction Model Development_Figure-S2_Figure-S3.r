#########################################################################################
###Figure 3A The constructing DNA methylation-based recurrence model.                   #
###Figure 3B  Unsupervised hierarchical clustering and heatmap of our cohort based on   # 
###the methylation patterns of the methylation eight locus.                             # 
###in functional regions.                                                               #
###Figure 3C Kaplan-Meier curves were plotted to analyze the correlations               # 
###between three groups and DFS in training data.                                       # 
###were enriched in lung/stem cell specific enhancers                                   #
###Figure 3D Kaplan-Meier curves were plotted to analyze the correlations between       # 
###two groups and DFS in training data.                                                 #
###Figure 3E  The AUC of model with eight methylations was higher than those            #
###constructed by less features.                                                        #
###Figure 3F ROC curves were generated to predict the recurrence of lung cancer in      #
###training data.                                                                       #
###Figure S2 The bimodal distribution of model score in the training data.              #
###Figure S3 The nomogram of the methylation model.                                     #
#########################################################################################

############
##Figure 3A#
############
#The codes for annotation of lung-specific or stem cell-specific enhancer regions
#were included in codes in Figure 2B

#gene annotation
options(stringsAsFactors=F)
options(scipen=200)
load("detaBeta_our.RData")
load("auc.RData")
ano <- read.table(file="/data1/dengyl/download/gsea/db/c5.bp.v6.1.symbols.gmt",sep="\n",quote = "",stringsAsFactors=F,header=F)
gene_list <- strsplit(ano[,1],"\t")
names(gene_list) <- sapply(gene_list,function(y) y[1])
metastasis_associated <- c("GO_NEGATIVE_REGULATION_OF_CELL_ADHESION",
"GO_CELLULAR_LIPID_METABOLIC_PROCESS",
"GO_REGULATION_OF_LYMPHOCYTE_APOPTOTIC_PROCESS")
gene_anno <- unique(unlist(gene_list[c(metastasis_associated)]))
anno_enhancer <- rep(0,length(names_enhancer))
names(anno_enhancer) <- meth_screen_all_enhancerd_exp[,"gene"]
anno_enhancer[names_enhancer%in%gene_anno] <- 1
load("meth_screen_all_enhancerd_exp.RData")
gtfl <- read.table(file="/data1/ref/hg38/gencode.v28.annotation.gtf",sep="\t",quote = "",stringsAsFactors=F,header=F)
genegtf <- gtfl[gtfl[,3]=="gene",]
genegtf_v9 <- strsplit(genegtf[,9],"\"| |;")
pos_g <- sapply(genegtf_v9,function(x) which(x=="gene_id"))
gene_id <- sapply(seq(length(genegtf_v9)),function(x) genegtf_v9[[x]][pos_g[x]+2])
pos_t <- sapply(genegtf_v9,function(x) which(x=="gene_name"))
gene_name <- sapply(seq(length(genegtf_v9)),function(x) genegtf_v9[[x]][pos_t[x]+2])
id_name <- data.frame(gene_id=gene_id,gene_name=gene_name,stringsAsFactors=F)
id_name[,"gene_id"] <- substr(id_name[,"gene_id"],1,15)
names_enhancer <- id_name[match(meth_screen_all_enhancerd_exp[,"gene"],id_name[,"gene_id"]),"gene_name"]
enhancer_candidate <- data.frame(region=meth_screen_all_enhancerd_exp[,"region"],
chr=meth_screen_all_enhancerd_exp[,2],start=meth_screen_all_enhancerd_exp[,3],
end=meth_screen_all_enhancerd_exp[,4],
element="enhancer",source=meth_screen_all_enhancerd_exp[,"V5"],
cellline=meth_screen_all_enhancerd_exp[,"V4"],
gene_id=meth_screen_all_enhancerd_exp[,"gene"],
gene_name=names_enhancer,
exp_coeff=meth_screen_all_enhancerd_exp[,"estiamte"],exp_p=meth_screen_all_enhancerd_exp[,"pvalue"],
detlaBeta=detaBeta_our[match(meth_screen_all_enhancerd_exp[,"region"],detaBeta_our[,1]),"deltaBeta"],
deltaBeta_P=detaBeta_our[match(meth_screen_all_enhancerd_exp[,"region"],detaBeta_our[,1]),"pvalue"],
anno=anno_enhancer,stringsAsFactors=F)
feature_candidate <- enhancer_candidate
feature_candidate <- feature_candidate[feature_candidate[,"anno"]==1,]
rownames(feature_candidate) <- NULL

#model training
options(stringsAsFactors=F)
options(scipen=200)
library(glmnet)
library(gplots)
library(pROC)
library(ROCR)
load("/data1/dengyl/project/S1/chennan/new3/GSE39279_mt100.RData")
load("/data1/dengyl/project/S1/chennan/new3/m100_mt.RData")
#training label
label <- rep(0,ncol(m100_mt))
label[substr(colnames(m100_mt),1,nchar(colnames(m100_mt))-3)=="tumor_relapse"] <- 1
our_mt100 <- m100_mt[,sapply(strsplit(colnames(m100_mt),"_"),function(x) x[1])=="tumor"]
label_our <- rep(0,ncol(our_mt100))
label_our[substr(colnames(our_mt100),1,nchar(colnames(our_mt100))-3)=="tumor_relapse"] <- 1
set.seed(1)
cvfit.lasso=cv.glmnet(t(m100_mt[enhancer_candidate[,1],]), label, family = "binomial",alpha=0)
model <- coef(cvfit.lasso, s = "lambda.min")
feature <- paste(rownames(model)[model[,1]!=0][-1],collapse=";")

############
##Figure 3B#
############
options(stringsAsFactors=F)
options(scipen=200)
library(gplots)
library(RColorBrewer)
library(glmnet)
library(pheatmap)
load("cvfit.lasso.RData")
label1 <- sapply(strsplit(colnames(m100_mt),"_"),function(x) paste(x[-length(x)],collapse="_"))
feature <- rownames(model)[-1]
mt100F <- m100_mt[feature,order(label1)]
label1 <- sort(label1)
colo <- rev(c(brewer.pal(3, "Set1"),brewer.pal(7, "Greens")[3]))
names(colo) <- sort(unique(label1))
annotation_col = data.frame(clss = factor(label1))
rownames(annotation_col) = colnames(mt100F)
pheatmap(mt100F,annotation_col = annotation_col,
clustering_distance_col = "correlation",
clustering_method = 'complete',
cellwidth = 2, cellheight = 12,
color = colorRampPalette(c("navy","white", "firebrick3"))(50),scale = "row")

##############
##Figure 3C-D#
##############
#Figure 3C
options(stringsAsFactors=F)
options(scipen=200)
library(glmnet)
library(pROC)
library(RColorBrewer)
library(ROCR)
library(survival)
library("survminer")
library(survcomp) 
setwd("/data1/dengyl/project/S1/chennan/new3/enumer_lung_enhancer")
load("/data1/dengyl/project/S1/chennan/new3/m100_mt.RData")
load("cvfit.lasso.RData")
load("/data1/dengyl/project/S1/chennan/new3/clini_ori.RData")
our_mt100 <- m100_mt[rownames(model)[-1],sapply(strsplit(colnames(m100_mt),"_"),function(x) x[1])=="tumor"]
label_our <- rep(0,ncol(our_mt100))
label_our[substr(colnames(our_mt100),1,nchar(colnames(our_mt100))-3)=="tumor_relapse"] <- 1
lasso_C <- predict(cvfit.lasso, newx = t(our_mt100),s = "lambda.min",type = "response")
lasso_C  <- lasso_C[,1]
names(lasso_C) <- colnames(our_mt100)
clini_ori <- clini_ori[match(names(lasso_C),clini_ori[,"label"]),]
nrow(clini_ori)
label <- rep(0,nrow(clini_ori))
label[(clini_ori$"Recurrence"=="Yes")] <- 1
y <-Surv(as.numeric(clini_ori$"DFS"),label)
Clusters <-rep(2,length(lasso_C))
Clusters[lasso_C>=0.6] <- 3
Clusters[lasso_C<=0.4] <- 1
Clusters <- factor(as.character(Clusters),levels=c("1","2","3"))
cliniT <- data.frame(recurrenceTime=as.numeric(clini_ori$"DFS")/12,
label=label,Clusters=Clusters,stringsAsFactors=F)
Sur_Curve<-survfit(Surv(recurrenceTime,label)~Clusters,data = cliniT)
Diff_Survival<-survdiff(Surv(as.numeric(clini_ori$"DFS"),label)~Clusters)
P_Value<-pchisq(Diff_Survival$chisq,length(table(Clusters))-1,lower.tail=F)#0.000005392187
ggsurvplot(Sur_Curve,data=cliniT,pvalue=TRUE,surv.median.line="hv",
risk.table = T,break.x.by=2,ylab="Disease-Free Survival")

#Figure 3D
options(stringsAsFactors=F)
options(scipen=200)
library(glmnet)
library(pROC)
library(RColorBrewer)
library(ROCR)
library(survival)
library("survminer")
library(survcomp) 
setwd("/data1/dengyl/project/S1/chennan/new3/enumer_lung_enhancer")
load("/data1/dengyl/project/S1/chennan/new3/m100_mt.RData")
load("cvfit.lasso.RData")
load("/data1/dengyl/project/S1/chennan/new3/clini_ori.RData")
our_mt100 <- m100_mt[rownames(model)[-1],sapply(strsplit(colnames(m100_mt),"_"),function(x) x[1])=="tumor"]
label_our <- rep(0,ncol(our_mt100))
label_our[substr(colnames(our_mt100),1,nchar(colnames(our_mt100))-3)=="tumor_relapse"] <- 1
lasso_C <- predict(cvfit.lasso, newx = t(our_mt100),s = "lambda.min",type = "response")
lasso_C  <- lasso_C[,1]
names(lasso_C) <- colnames(our_mt100)
auc_our_continus_all_sample <- pROC::roc(label_our, lasso_C,direction= "<" )$"auc"[[1]]##
lasso_C <- lasso_C[(lasso_C>=0.6)|(lasso_C<=0.4)]
clini_ori <- clini_ori[match(names(lasso_C),clini_ori[,"label"]),]
names(label_our) <- colnames(our_mt100)
label_our <- label_our[names(lasso_C)]
tt2 <- rep(0,length(lasso_C))
tt2[lasso_C>=0.6] <- 1
label <- rep(0,nrow(clini_ori))
label[(clini_ori$"Recurrence"=="Yes")] <- 1
y <-Surv(as.numeric(clini_ori$"DFS"),label)
Clusters <-rep(1,length(lasso_C))
Clusters[lasso_C>=0.6] <- 2
Clusters <- factor(as.character(Clusters),levels=c("1","2"))
cliniT <- data.frame(recurrenceTime=as.numeric(clini_ori$"DFS")/12,
label=label,Clusters=Clusters,stringsAsFactors=F)
Sur_Curve<-survfit(Surv(recurrenceTime,label)~Clusters,data = cliniT)
Diff_Survival<-survdiff(Surv(as.numeric(clini_ori$"DFS"),label)~Clusters)
P_Value<-pchisq(Diff_Survival$chisq,length(table(Clusters))-1,lower.tail=F)#0.0000009714635
ggsurvplot(Sur_Curve,data=cliniT,pvalue=TRUE,surv.median.line="hv",
risk.table = T,break.x.by=2,ylab="Disease-Free Survival")

############
##Figure 3E#
############
#subsample
options(stringsAsFactors=F)
options(scipen=200)
library(glmnet)
library(pROC)
library(ROCR)
load("/data1/dengyl/project/S1/chennan/new3/m100_mt.RData")
load("cvfit.lasso.RData")
label <- rep(0,ncol(m100_mt))
label[substr(colnames(m100_mt),1,nchar(colnames(m100_mt))-3)=="tumor_relapse"] <- 1
our_mt100 <- m100_mt[,sapply(strsplit(colnames(m100_mt),"_"),function(x) x[1])=="tumor"]
label_our <- rep(0,ncol(our_mt100))
label_our[substr(colnames(our_mt100),1,nchar(colnames(our_mt100))-3)=="tumor_relapse"] <- 1
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
model_our <- sapply(per1000,function(x) x[1,1])
p_our <- sapply(1:246,function(y){
	model_our_sub1 <- sapply(per1000,function(x) x[y+1,1])
	return(t.test(model_our_sub1,model_our,alternative="less")$"p.value")
}) 
#for Figure 3E
options(stringsAsFactors=F)
options(scipen=200)
library(glmnet)
library(pROC)
library(ROCR)
library(RColorBrewer)
load("/data1/dengyl/project/S1/chennan/new3/m100_mt.RData")
load("cvfit.lasso.RData")
load("per1000.RData")
label <- rep(0,ncol(m100_mt))
label[substr(colnames(m100_mt),1,nchar(colnames(m100_mt))-3)=="tumor_relapse"] <- 1
our_mt100 <- m100_mt[,sapply(strsplit(colnames(m100_mt),"_"),function(x) x[1])=="tumor"]
label_our <- rep(0,ncol(our_mt100))
label_our[substr(colnames(our_mt100),1,nchar(colnames(our_mt100))-3)=="tumor_relapse"] <- 1
model_our <- sapply(per1000,function(x) x[1,1])
p_our <- sapply(1:246,function(y){
	model_our_sub1 <- sapply(per1000,function(x) x[y+1,1])
	return(t.test(model_our_sub1,model_our,alternative="less")$"p.value")
}) 
p_our_list <- lapply(1:246,function(y){
	model_our_sub1 <- sapply(per1000,function(x) x[y+1,1])
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
p_our_index <- p_our_list[order(sapply(p_our_list,median),decreasing=T)]
boxplot(p_our_index,lwd=0.1,
col=colo[apply(para_mt,1,sum)[order(sapply(p_our_list,median),decreasing=T)]],
main="AUC in our Data",ylim=c(0.6,0.9),outline = F)
abline(h=0.89,col=brewer.pal(8, "Reds")[8],lwd=2)

############
##Figure 3F#
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
load("/data1/dengyl/project/S1/chennan/new2/Tmt.RData")
load("cvfit.lasso.RData")
load("/data1/dengyl/project/S1/chennan/new3/m100_mt.RData")
load("/data1/dengyl/project/S1/chennan/new3/clini_ori.RData")
jama_model <- read.table(file="/data1/dengyl/project/S1/chennan/new2/jama_model.txt",sep="\t",quote = "",stringsAsFactors=F,header=T)
our_mt100 <- m100_mt[rownames(model)[-1],sapply(strsplit(colnames(m100_mt),"_"),function(x) x[1])=="tumor"]
label_our <- rep(0,ncol(our_mt100))
label_our[substr(colnames(our_mt100),1,nchar(colnames(our_mt100))-3)=="tumor_relapse"] <- 1
lasso_C <- predict(cvfit.lasso, newx = t(our_mt100),s = "lambda.min",type = "response")
lasso_C  <- lasso_C[,1]
names(lasso_C) <- colnames(our_mt100)
label_our <- label_our[(lasso_C>=0.6)|(lasso_C<0.4)]
lasso_C <- lasso_C[(lasso_C>=0.6)|(lasso_C<0.4)]
clini_ori <- clini_ori[match(names(lasso_C),clini_ori[,"label"]),]
rownames(Tmt) <- sapply(strsplit(rownames(Tmt),".",fixed=T),function(x) x[1])
Tmt <- Tmt[,names(lasso_C)]
score <- apply(Tmt,2,function(x) {
	names(x) <- rownames(Tmt)
	res <- sum(as.numeric(x[jama_model[,2]] < x[jama_model[,4]])*jama_model[,5])
	return(res)
}) 
pred_list1 <- prediction( lasso_C, label_our )
roc1 <- roc(label_our,lasso_C, percent=TRUE,ci=TRUE, boot.n=100, ci.alpha=0.9)
pred_list2 <- prediction( score, label_our )
roc2 <- roc(label_our,score, percent=TRUE,ci=TRUE, boot.n=100, ci.alpha=0.9)
pred_list3 <- prediction( as.numeric(clini_ori$"Pleural.Invasion"=="Yes"), label_our )
roc3 <- roc(label_our,as.numeric(clini_ori$"Pleural.Invasion"=="Yes"),
percent=TRUE,ci=TRUE, boot.n=100, ci.alpha=0.9)
pred_list4 <- prediction( as.numeric(clini_ori$"Stage"=="IIA"), label_our )
perf <- performance( pred_list1, "tpr", "fpr" )
plot( perf,col=colo[1],main="our_method_compare_auc",lwd=2)
auc2 <- performance( pred_list1, "auc")@"y.values"[[1]]
text(0.5,0.5,label=paste0("model auc=",substr(as.character(auc2),1,5)),col = colo[1])
##jama
perf <- performance( pred_list2, "tpr", "fpr" )
plot(perf, add = TRUE,col=colo[2],lwd=2)
auc2 <- performance( pred_list2, "auc")@"y.values"[[1]]
text(0.5,0.5-0.05*1,label=paste0("Li auc=",substr(as.character(auc2),1,5)),col = colo[2])
##Pleural.Invasion
perf <- performance( pred_list3, "tpr", "fpr" )
plot(perf, add = TRUE,col=colo[3],lwd=2)
auc2 <- performance( pred_list3, "auc")@"y.values"[[1]]
text(0.5,0.5-0.05*2,label=paste0("pleural invasion auc=",substr(as.character(auc2),1,5)),col = colo[3])
##stage
perf <- performance( pred_list4, "tpr", "fpr" )
plot(perf, add = TRUE,col=colo[4],lwd=2)
auc2 <- performance( pred_list2, "auc")@"y.values"[[1]]
text(0.5,0.5-0.05*3,label=paste0("stage auc=",substr(as.character(auc2),1,5)),col = colo[4])


############
##Figure S2#
############
options(stringsAsFactors=F)
options(scipen=200)
library(GenomicRanges)
library(glmnet)
library(gplots)
library(pROC)
library(RColorBrewer)
library(ROCR)
load("/data1/dengyl/project/S1/chennan/new2/m100_mt.RData")
load("cvfit.lasso.RData")
label <- rep(0,ncol(m100_mt))
label[substr(colnames(m100_mt),1,nchar(colnames(m100_mt))-3)=="tumor_relapse"] <- 1
m100_mt <- m100_mt[rownames(model)[-1],]
lasso_C <- predict(cvfit.lasso, newx = t(m100_mt),s = "lambda.min",type = "response")
lassoC <- lasso_C[,1]
names(lassoC) <- sapply(strsplit(rownames(lasso_C),"_"),function(x) as.numeric(x[length(x)]))
plot(density(lassoC),main="distribution of model score",xlab="",axes="n")
abline(v=0.6,col="black",lty="longdash")
abline(v=0.4,col="black",lty="longdash")
axis(2)
axis(1,c(0,0.2,0.4,0.6,0.8,1),as.character(c(0,0.2,0.4,0.6,0.8,1)))
box()

############
##Figure S3#
############
options(stringsAsFactors=F)
options(scipen=200)
library(glmnet)
load("cvfit.lasso.RData")

#formula
sigmod_F2 <- function(x)
{
	return(exp(x)/(1+exp(x)))
}

#function
nomogram_F2 <- function(cvfit)
{
	library(glmnet)
	eff <- coef(cvfit,s = "lambda.min")

	##2.2.text
	scl <- ceiling(max(abs(eff))*2)/2
	ncov <- sum(eff!=0)
	plot(c(-2,5),c(0,ncov+3),type="n",axes = FALSE,xlab="",ylab="")
	text(-2, ncov+2, labels = "Points",pos=4)
	lines(c(0,5),c(ncov+2,ncov+2))
	for(i in seq(from=0,to=5,by=0.1))
	{
		lines(c(i,i),c(ncov+2,ncov+2.1))
	}
	for(i in seq(from=0,to=5,by=0.5))
	{
		lines(c(i,i),c(ncov+2,ncov+2.25))
		text(i, ncov+2.6, labels = i*scl*20,cex=0.6)
	}

	eff1 <- eff[-1]
	names(eff1) <- rownames(eff)[-1]
	eff1 <- eff1[eff1!=0]

	for(j in seq(length(eff1)))
	{
		lines(c(0,5*abs(eff1[j])/scl),c(ncov+2-j,ncov+2-j))
		text(-2, ncov+2-j, labels = names(eff1)[j],pos=4)
		seq1 <- seq(from=0,to=5*abs(eff1[j])/scl,by=1*abs(eff1[j])/scl)
		if(eff1[j]>0)
		{
			for(i in 1:6)
			{
				lines(c(seq1[i],seq1[i]),c(ncov+2-j,ncov+1.8-j))
				text(seq1[i], ncov+1.5-j, labels = (i-1)*20,cex=0.4)
			}
		}else{
			for(i in 1:6)
			{
				lines(c(seq1[i],seq1[i]),c(ncov+2-j,ncov+1.8-j))
				text(seq1[i], ncov+1.5-j, labels = (6-i)*20,cex=0.4)
			}
		}
	}

	Rsigmod_F <- function(x)
	{
		return(-log((1-x)/x))
	}
	lines(c(0,5),c(2,2))
	text(-2, 2, labels = "Total Points",pos=4)
	seq2 <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
	seq1 <- seq(from=0,to=sum(abs(eff1)),by=sum(abs(eff1))/5)
	lines(c(0,5),c(2,2))
	for(i in seq(from=0,to=5,by=0.1))
	{
		lines(c(i,i),c(2,1.9))
	}
	for(i in seq(from=0,to=5,by=1))
	{
		lines(c(i,i),c(2,1.75))
		text(i, 1.6, labels = round(seq1[i+1]*100),cex=0.6)
	}
	
	seq3 <- (Rsigmod_F(seq2)-eff[1]+sum(abs(eff1[eff1<0])))*5/sum(abs(eff1))

	lines(c(seq3[1],seq3[length(seq3)]),c(1,1))
	text(-2, 1, labels = "Risk",pos=4)
	for(i in seq(from=0,to=8,by=1))
	{
		lines(c(seq3[i+1],seq3[i+1]),c(1,0.75))
		text(seq3[i+1], 0.6, labels = seq2[i+1],cex=0.6)
	}
}

#plot
nomogram_F2(cvfit.lasso)


