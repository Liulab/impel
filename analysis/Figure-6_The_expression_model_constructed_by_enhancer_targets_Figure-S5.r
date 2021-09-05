##########################################################################################################
###Figure 6A The locus zoom for enhancer targets of chr16:1079894-1080188 with link (top),               #
###ENCODE cCREs (middle) and DNase signals (bottom).                                                     #
###Figure 6B The locus zoom for enhancer targets of chr11:9752161- 9816165 with link (top),              #
###ENCODE cCREs (middle) and DNase signals (bottom).                                                     #
###Figure 6C Kaplan-Meier curves were plotted to analyze the correlations between expression model and   #
###DFS in training data.                                                                                 #
###Figure 6D ROC curves were generated to predict the recurrence of lung cancer in training data.        #
###Figure 6E Kaplan-Meier curves were plotted to analyze the correlations between expression model       #
###and DFS in TCGA.                                                                                      #
###Figure 6F ROC curves were generated to predict the recurrence of lung cancer in TCGA.                 #
###Figure 6G Violin plot for the mean expression of target genes in epithelial cells from nLung, tLung,  #
###mLN and mBrain.                                                                                       #
###Figure 6H UMAP plot of cancer cells color-coded by tissue origins in patient 19.                      # 
###Figure 6I Diffusion map showed the difference between c1 and other classes.                           #
###Figure 6J Diffusion map, where c1 and mBrain were excluded, showed that c2 was close to mLN at        #
###diffusion componeDC2.                                                                                 #
###Figure 6K Diffusion map, where c1 and mLN were excluded, showed that c2 was close to mBrain at DC3    #
###Figure 6L The distribution of cells that showed positive signals of enhancer targets, where           #
###positive signals mean that expression of any target genes.                                            #
###Figure 6J Volin plot for the mean expression of target genes in c2 and other cancer cells from tLung. #
###Figure S5A Overview of tissue origins in the single cell RNA sequencing. nLung, epithelial cells      #
###from normal lung tissues                                                                              #
###Figure S5B UMAP plot of epithelial cells color-coded by the mean expression of target genes.          #
###Figure S5C UMAP plot of cancer cells color-coded by tissue origins in patient 19.                     #
###Figure S5D The cancer cells from tLung were clustered into 5 classes. The numbers of cells for each   # 
###classes were labeled in parentheses.                                                                  #
###Figure S5E The cumulative distribution function of DC2 for c2, c3, c4 and mLN, respectively.          # 
###Figure S5F The significantly enriched pathways for DC2 by fgsea.                                      #
###Figure S5G The cumulative distribution function of DC3 for c2, c3, c4 and mBrain, respectively.       #
###Figure S5H The significantly enriched pathways for DC3 by fgsea.                                      #
###Figure S5I Volin plot for the mean expression of target genes in c1, c2, c3 and c4.                   #
##########################################################################################################

############
##Figure 6A#
############
#visulation by UCSC


############
##Figure 6B#
############
#visulation by UCSC

############
##Figure 6C#
############
#The enhancer target genes were obtained though EpiMap(10.1038/s41586-020-03145-z) with 
#the RNA sequencing data from our cohort
#deal with epimap datasets
#bash
gunzip -c links_by_group.adipose.tsv.gz > links_by_group.adipose.tsv
gunzip -c links_by_group.lung.tsv.gz > links_by_group.lung.tsv
gunzip -c links_by_group.blood_t_cell.tsv.gz > links_by_group.blood_t_cell.tsv
gunzip -c links_by_group.bone.tsv.gz > links_by_group.bone.tsv
gunzip -c links_by_group.brain.tsv.gz > links_by_group.brain.tsv
gunzip -c links_by_group.digestive.tsv.gz > links_by_group.digestive.tsv
gunzip -c links_by_group.endocrine.tsv.gz > links_by_group.endocrine.tsv
gunzip -c links_by_group.endothelial.tsv.gz > links_by_group.endothelial.tsv
gunzip -c links_by_group.epithelial.tsv.gz > links_by_group.epithelial.tsv
gunzip -c links_by_group.esc.tsv.gz > links_by_group.esc.tsv
gunzip -c links_by_group.es_deriv.tsv.gz > links_by_group.es_deriv.tsv
gunzip -c links_by_group.eye.tsv.gz > links_by_group.eye.tsv
gunzip -c links_by_group.heart.tsv.gz > links_by_group.heart.tsv
gunzip -c links_by_group.hsc_b_cell.tsv.gz > links_by_group.hsc_b_cell.tsv
gunzip -c links_by_group.ipsc.tsv.gz > links_by_group.ipsc.tsv
gunzip -c links_by_group.kidney.tsv.gz > links_by_group.kidney.tsv
gunzip -c links_by_group.liver.tsv.gz > links_by_group.liver.tsv
gunzip -c links_by_group.lymphoblastoid.tsv.gz > links_by_group.lymphoblastoid.tsv
gunzip -c links_by_group.mesench.tsv.gz > links_by_group.mesench.tsv
gunzip -c links_by_group.muscle.tsv.gz > links_by_group.muscle.tsv
gunzip -c links_by_group.myosat.tsv.gz > links_by_group.myosat.tsv
gunzip -c links_by_group.neurosph.tsv.gz > links_by_group.neurosph.tsv
gunzip -c links_by_group.pancreas.tsv.gz > links_by_group.pancreas.tsv
gunzip -c links_by_group.placenta_eem.tsv.gz > links_by_group.placenta_eem.tsv
gunzip -c links_by_group.pns.tsv.gz > links_by_group.pns.tsv
gunzip -c links_by_group.reproductive.tsv.gz > links_by_group.reproductive.tsv
gunzip -c links_by_group.smmuscle.tsv.gz > links_by_group.smmuscle.tsv
gunzip -c links_by_group.spleen.tsv.gz > links_by_group.spleen.tsv
gunzip -c links_by_group.stromal.tsv.gz > links_by_group.stromal.tsv
gunzip -c links_by_group.thymus.tsv.gz > links_by_group.thymus.tsv
gunzip -c links_by_group.urinary.tsv.gz > links_by_group.urinary.tsv
cat *.tsv > alle.tsv
#R target gene
options(stringsAsFactors=F)
options(scipen=200)
library(glmnet)
library(gplots)
library(pROC)
library(ROCR)
library(survival)
library(GenomicRanges)
load("cvfit.lasso.RData")
load("m100_mt.RData")
load("Tmt.RData")
our_mt100 <- m100_mt[,sapply(strsplit(colnames(m100_mt),"_"),function(x) x[1])=="tumor"]
label_our <- rep(0,ncol(our_mt100))
label_our[substr(colnames(our_mt100),1,nchar(colnames(our_mt100))-3)=="tumor_relapse"] <- 1
Tmt2 <- Tmt[!duplicated(t1),]
rownames(Tmt2) <- sapply(strsplit(rownames(Tmt2),"_"),function(x) x[2])
Tmt2 <- log2(Tmt2+1)
Tmt2 <- Tmt2[,sapply(strsplit(colnames(Tmt2),"_"),function(x) x[1])=="tumor"]
label <- rep(0,ncol(Tmt2))
label[substr(colnames(Tmt2),1,nchar(colnames(Tmt2))-3)=="tumor_relapse"] <- 1
house<-c('RPL38', 'UBA52', 'RPL4', 'RPS29', 'SLC25A3', 'CLTC', 'RPL37', 'PSMA1', 'RPL8', 'PPP2CA','TXNL1', 'MMADHC', 'PSMC1', 'RPL13A', 'MRFAP1')
enhancer_target <- read.table(file="/NAS/dyl/project/chennan/5.expression/EpiMap_enhancer/alle.tsv",sep="\t",quote = "",stringsAsFactors=F,header=T) 
t2 <- strsplit(rownames(Tmt),"_",fixed=T)
gene_id <- matrix(unlist(strsplit(rownames(Tmt)[sapply(t2,length)==2],"_",fixed=T)),ncol=2,byrow=T)
gene_id[,1] <- substr(gene_id[,1],1,15)
enhancer_target <- enhancer_target[enhancer_target[,1]%in%gene_id[,1],]
Gregion <- GRanges(seqnames = Rle(enhancer_target[,5]),
ranges = IRanges(as.numeric(enhancer_target[,6]),end = as.numeric(enhancer_target[,7])))
model <- coef(cvfit.lasso, s = "lambda.min")
er <- rownames(model)[-1]
ermt <- matrix(unlist(strsplit(er,".",fixed=T)),ncol=3,byrow=T)
Eregion <- GRanges(seqnames = Rle(ermt[,1]),
ranges = IRanges(as.numeric(ermt[,2])-1,end = as.numeric(ermt[,3]))) 
mtch1 <- findOverlaps(Gregion,Eregion,maxgap=20000L)
pair1 <- as.matrix(mtch1)
enhancer_target_pair <- cbind(enhancer_target[pair1[,1],],er[pair1[,2]],stringsAsFactors=F)
target <- unique(enhancer_target_pair[,1])
target_gene <- gene_id[match(target,gene_id[,1]),2]
deltaBeta <- apply(Tmt2[target_gene,],1,function(x) mean(as.numeric(x)[label==1]/apply(Tmt2[house,label==1],2,mean))-mean(as.numeric(x)[label==0]/apply(Tmt2[house,label==0],2,mean)))
pvalue <- apply(Tmt2[target_gene,],1,function(x) {
	t1 <- as.numeric(x)[label==1]/apply(Tmt2[house,label==1],2,mean)
	t2 <- as.numeric(x)[label==0]/apply(Tmt2[house,label==0],2,mean)
	t3 <- wilcox.test(t1,t2)$"p.value"
	return(t3)
})
pvalue[is.na(pvalue)] <- 1
auc_our <- c()
for(i in target_gene){
	roc1 <- pROC::roc(label, unlist(Tmt2[i,]-apply(Tmt2[house,],2,mean)),direction= "<" ) 
	auc_our <- c(auc_our,roc1$"auc")
}
names(auc_our) <- target_gene
gene_list <- target_gene[(auc_our>0.5)&(deltaBeta>0)&(pvalue<0.05)]

#develop expression model
options(stringsAsFactors=F)
options(scipen=200)
library(glmnet)
library(gplots)
library(pROC)
library(ROCR)
library(survival)
load("Tmt.RData")
t1 <- sapply(strsplit(rownames(Tmt),"_"),function(x) x[2])
Tmt2 <- Tmt[!duplicated(t1),]
rownames(Tmt2) <- sapply(strsplit(rownames(Tmt2),"_"),function(x) x[2])
Tmt2 <- log2(Tmt2+1)
Tmt2 <- Tmt2[,sapply(strsplit(colnames(Tmt2),"_"),function(x) x[1])=="tumor"]
label <- rep(0,ncol(Tmt2))
label[substr(colnames(Tmt2),1,nchar(colnames(Tmt2))-3)=="tumor_relapse"] <- 1
house<-c('RPL38', 'UBA52', 'RPL4', 'RPS29', 'SLC25A3', 'CLTC', 'RPL37', 'PSMA1', 'RPL8', 'PPP2CA','TXNL1', 'MMADHC', 'PSMC1', 'RPL13A', 'MRFAP1')
mt_tr= t(apply(Tmt2[gene_list,],1,function(x) x/apply(Tmt2[house,],2,mean))
label_tr=label
alpha=1
set.seed(1)
cvfit.lasso=cv.glmnet(t(mt_tr), label_tr, family = "binomial",alpha=alpha)
model <- coef(cvfit.lasso, s = "lambda.min")

#plot Figure 6C
options(stringsAsFactors=F)
options(scipen=200)
library(glmnet)
library(pROC)
library(RColorBrewer)
library(ROCR)
library(survival)
library("survminer")
library(survcomp) 
load("Tmt.RData")
load("expression_model.RData")
load("clini_ori.RData")
house<-c('RPL38', 'UBA52', 'RPL4', 'RPS29', 'SLC25A3', 'CLTC', 'RPL37', 'PSMA1', 'RPL8', 'PPP2CA','TXNL1', 'MMADHC', 'PSMC1', 'RPL13A', 'MRFAP1')
t1 <- sapply(strsplit(rownames(Tmt),"_"),function(x) x[2])
Tmt2 <- Tmt[!duplicated(t1),]
rownames(Tmt2) <- sapply(strsplit(rownames(Tmt2),"_"),function(x) x[2])
Tmt2 <- log2(Tmt2+1)
Tmt2 <- Tmt2[,sapply(strsplit(colnames(Tmt2),"_"),function(x) x[1])=="tumor"]
mt_tr=t(apply(Tmt2[gene_list,],1,function(x) as.numeric(x)-apply(Tmt2[house,],2,mean)))
lasso_C <- predict(cvfit.lasso, newx = t(mt_tr), s = "lambda.min",type = "response")
clini_ori <- clini_ori[match(colnames(mt_tr),clini_ori[,"label"]),]
label <- rep(0,nrow(clini_ori))
label[(clini_ori$"Recurrence"=="Yes")] <- 1
y <-Surv(as.numeric(clini_ori$"DFS"),label)
Clusters <-rep(2,length(lasso_C))
Clusters[lasso_C<=0.5] <- 1
Clusters <- factor(as.character(Clusters),levels=c("1","2"))
cliniT <- data.frame(recurrenceTime=as.numeric(clini_ori$"DFS")/12,
label=label,Clusters=Clusters,stringsAsFactors=F)
Sur_Curve<-survfit(Surv(recurrenceTime,label)~Clusters,data = cliniT)
Diff_Survival<-survdiff(Surv(as.numeric(clini_ori$"DFS"),label)~Clusters)
P_Value<-pchisq(Diff_Survival$chisq,length(table(Clusters))-1,lower.tail=F)#0.0001701026
ggsurvplot(Sur_Curve,data=cliniT,pvalue=TRUE,surv.median.line="hv",
risk.table = T,break.x.by=2,ylab="Disease-Free Survival")

############
##Figure 6D#
############
options(stringsAsFactors=F)
options(scipen=200)
library(glmnet)
library(gplots)
library(pROC)
library(RColorBrewer)
library(ROCR)
colo=brewer.pal(4, "Set1")[1:4]
#expression matrix
load("Tmt.RData")
load("expression_model.RData")
expcvfit.lasso <- cvfit.lasso
#methylation model
load("cvfit.lasso.RData")
#methylation matrix
load("m100_mt.RData")
load("clini_ori.RData")

model <- coef(cvfit.lasso, s = "lambda.min")
our_mt100 <- m100_mt[rownames(model)[-1],sapply(strsplit(colnames(m100_mt),"_"),function(x) x[1])=="tumor"]
label_our <- rep(0,ncol(our_mt100))
label_our[substr(colnames(our_mt100),1,nchar(colnames(our_mt100))-3)=="tumor_relapse"] <- 1
house<-c('RPL38', 'UBA52', 'RPL4', 'RPS29', 'SLC25A3', 'CLTC', 'RPL37', 'PSMA1', 'RPL8', 'PPP2CA','TXNL1', 'MMADHC', 'PSMC1', 'RPL13A', 'MRFAP1')
lasso_C1 <- predict(cvfit.lasso, newx = t(our_mt100),s = "lambda.min",type = "response")
lasso_C  <- lasso_C1[,1]
names(lasso_C) <- colnames(our_mt100)
label_our <- label_our[(lasso_C>=0.6)|(lasso_C<0.4)]
lasso_C <- lasso_C[(lasso_C>=0.6)|(lasso_C<0.4)]
clini_ori <- clini_ori[match(names(lasso_C),clini_ori[,"label"]),]
t1 <- sapply(strsplit(rownames(Tmt),"_"),function(x) x[2])
Tmt2 <- Tmt[!duplicated(t1),]
rownames(Tmt2) <- sapply(strsplit(rownames(Tmt2),"_"),function(x) x[2])
Tmt2 <- log2(Tmt2+1)
Tmt2 <- Tmt2[,sapply(strsplit(colnames(Tmt2),"_"),function(x) x[1])=="tumor"]
label <- rep(0,ncol(Tmt2))
label[substr(colnames(Tmt2),1,nchar(colnames(Tmt2))-3)=="tumor_relapse"] <- 1
mt_tr=t(apply(Tmt2[gene_list[[78]],],1,function(x) as.numeric(x)-apply(Tmt2[house,],2,mean)))
lasso_Cexp <- predict(expcvfit.lasso, newx = t(mt_tr), s = "lambda.min",type = "response")
lasso_Cexp <- lasso_Cexp[(lasso_C1[,1]>=0.6)|(lasso_C1[,1]<0.4)]
roc1 <- roc(label_our,lasso_C, percent=TRUE,ci=TRUE, boot.n=100, ci.alpha=0.9)
pred_list1 <- prediction( lasso_C, label_our )
roc2 <- roc(label_our,lasso_Cexp, percent=TRUE,ci=TRUE, boot.n=100, ci.alpha=0.9)
pred_list2 <- prediction( lasso_Cexp, label_our )
perf <- performance( pred_list1, "tpr", "fpr" )
plot( perf,col=colo[1],main="our dataset",lwd=2)
auc2 <- performance( pred_list1, "auc")@"y.values"[[1]]
text(0.5,0.5,label=paste0("methylation auc=",substr(as.character(auc2),1,5)),col = colo[1])
perf <- performance( pred_list2, "tpr", "fpr" )
plot(perf, add = TRUE,col=colo[2],lwd=2)
auc2 <- performance( pred_list2, "auc")@"y.values"[[1]]
text(0.5,0.5-0.05*1,label=paste0("expression auc=",substr(as.character(auc2),1,5)),col = colo[2])

	
############
##Figure 6E#
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

load("expression_model.RData")
load("FPKMmt_geneName_log2.RData")
load("cliniTCGA.RData")
house<-c('RPL38', 'UBA52', 'RPL4', 'RPS29', 'SLC25A3', 'CLTC', 'RPL37', 'PSMA1', 'RPL8', 'PPP2CA','TXNL1', 'MMADHC', 'PSMC1', 'RPL13A', 'MRFAP1')
inter <- intersect(colnames(FPKMmt),cliniTCGA[,1])##122
FPKMmt <- FPKMmt[,inter]
cliniTCGA <- cliniTCGA[match(inter,cliniTCGA[,1]),]
labeltcga <- rep(0,nrow(cliniTCGA))
labeltcga[(cliniTCGA$"recurrence"%in%c("Distant Metastasis","Locoregional Recurrence"))] <- 1

mt_tcga=t(apply(FPKMmt[gene_list[[78]],],1,function(x) as.numeric(x)-apply(FPKMmt[house,],2,mean)))

lasso_C <- predict(cvfit.lasso, newx = t(mt_tcga),s = "lambda.min",type = "response")
lasso_C  <- lasso_C [,1]
names(lasso_C) <- colnames(mt_tcga)
cliniTCGA <- cliniTCGA[match(names(lasso_C),cliniTCGA[,1]),]
label <- rep(0,nrow(cliniTCGA))
label[(cliniTCGA$"recurrence"%in%c("Distant Metastasis","Locoregional Recurrence"))] <- 1
y <-Surv(as.numeric(cliniTCGA$"recurrenceTime"),label)
Clusters <-rep(2,length(lasso_C))
Clusters[lasso_C<=0.5] <- 1
Clusters <- factor(as.character(Clusters),levels=c("1","2"))
cliniT <- data.frame(recurrenceTime=as.numeric(cliniTCGA$"recurrenceTime"),
label=label,Clusters=Clusters,stringsAsFactors=F)
Sur_Curve<-survfit(Surv(recurrenceTime,label)~Clusters,data = cliniT)
Diff_Survival<-survdiff(Surv(as.numeric(cliniTCGA$"recurrenceTime"),label)~Clusters)
P_Value<-pchisq(Diff_Survival$chisq,length(table(Clusters))-1,lower.tail=F)#0.03400936
ggsurvplot(Sur_Curve,data=cliniT,pvalue=TRUE,surv.median.line="hv",
risk.table = T,break.x.by=2,ylab="Disease-Free Survival")

############
##Figure 6F#
############
options(stringsAsFactors=F)
options(scipen=200)
library(glmnet)
library(gplots)
library(pROC)
library(RColorBrewer)
library(ROCR)
colo=brewer.pal(4, "Set1")[1:4]
load("expression_model.RData")
expcvfit.lasso <- cvfit.lasso
load("cvfit.lasso.RData")
load("labeltcga.RData")
load("tcga_mt100.RData")
load("cliniTCGA.RData")
load("FPKMmt_geneName_log2.RData")

inter <- intersect(colnames(tcga_mt100),colnames(FPKMmt))
tcga_mt100 <- tcga_mt100[,inter]
FPKMmt <- FPKMmt[,inter]
house<-c('RPL38', 'UBA52', 'RPL4', 'RPS29', 'SLC25A3', 'CLTC', 'RPL37', 'PSMA1', 'RPL8', 'PPP2CA','TXNL1', 'MMADHC', 'PSMC1', 'RPL13A', 'MRFAP1')
model <- coef(cvfit.lasso, s = "lambda.min")
tcga_mt100 <- tcga_mt100[rownames(model)[-1],]
lasso_C1 <- predict(cvfit.lasso, newx = t(tcga_mt100),s = "lambda.min",type = "response")
lasso_C  <- lasso_C1[,1]
names(lasso_C) <- colnames(tcga_mt100)
lasso_C <- lasso_C[(lasso_C>=0.6)|(lasso_C<0.4)]
cliniTCGA <- cliniTCGA[match(names(lasso_C),cliniTCGA[,1]),]
label <- rep(0,nrow(cliniTCGA))
label[(cliniTCGA$"recurrence"%in%c("Distant Metastasis","Locoregional Recurrence"))] <- 1
mt_tcga <- t(apply(FPKMmt[gene_list[[57]],],1,function(x) as.numeric(x)-apply(FPKMmt[house,],2,mean)))
lasso_Cexp <- predict(expcvfit.lasso, newx = t(mt_tcga), s = "lambda.min",type = "response")
lasso_Cexp <- lasso_Cexp[(lasso_C1[,1]>=0.6)|(lasso_C1[,1]<0.4)]
roc1 <- roc(label,lasso_C, percent=TRUE,ci=TRUE, boot.n=100, ci.alpha=0.9)#61.6%
tt2 <- rep(0,length(lasso_C))
tt2[lasso_C>=0.6] <- 1
pred_list1 <- prediction( lasso_C, label )
pred_list2 <- prediction( lasso_Cexp, label )
perf <- performance( pred_list1, "tpr", "fpr" )
plot( perf,col=colo[1],main="TCGA dataset",lwd=2)
auc2 <- performance( pred_list1, "auc")@"y.values"[[1]]
text(0.5,0.5,label=paste0("methylation auc=",substr(as.character(auc2),1,5)),col = colo[1])
perf <- performance( pred_list2, "tpr", "fpr" )
plot(perf, add = TRUE,col=colo[2],lwd=2)
auc2 <- performance( pred_list2, "auc")@"y.values"[[1]]
text(0.5,0.5-0.05*1,label=paste0("expression auc=",substr(as.character(auc2),1,5)),col = colo[2])


##################
##Figure 6G,S5A-B#
##################
#process of scRNAseq
#raw counts and phenotype information could obtain from GEO (GSE131907)
library(Seurat)
library(dplyr)
library(patchwork)
ref.data <- read.table(file = "/NAS/dyl/project/singlecell/GGO/3.seraut/3.integrate/GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt",sep="\t",quote = "",stringsAsFactors=F,header=T,row.names=1) 
pheno <- read.table(file = "/NAS/dyl/project/singlecell/GGO/3.seraut/3.integrate/GSE131907/GSE131907_Lung_Cancer_cell_annotation.txt",sep="\t",quote = "",stringsAsFactors=F,header=T)
smpsel <- pheno[pheno[,"Cell_type"]%in%c("Epithelial cells"),1]
ref.data <- ref.data[,smpsel]
ref <- CreateSeuratObject(counts = ref.data, project = "NC", min.cells = 3, min.features = 200)
ref[["percent.mt"]] <- PercentageFeatureSet(ref, pattern = "^MT-")
ref <- subset(ref, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
ref <-subset(ref,features=setdiff(rownames(ref@assays[[1]]),
grep("MT-",rownames(ref@assays[[1]]), value = T))) 
ref <- NormalizeData(ref, verbose = FALSE)
expr_mt <- GetAssayData(ref,assay = "RNA")
load("FPKMmt_geneName_log2.RData")
load("Tmt.RData")
t1 <- sapply(strsplit(rownames(Tmt),"_"),function(x) x[2])
Tmt2 <- Tmt[!duplicated(t1),]
rownames(Tmt2) <- sapply(strsplit(rownames(Tmt2),"_"),function(x) x[2])
intergene <- intersect(rownames(Tmt2),intersect(rownames(FPKMmt),rownames(expr_mt)))
expr_mt <- expr_mt[intergene,]
labelscRNAseq <- pheno[match(colnames(expr_mt),pheno[,1]),4]

#plot Figure 6G
library(Seurat)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(glmnet)
library(ggplot2)
colo=brewer.pal(6, "Set1")[1:6]
ref.data <- read.table(file = "/NAS/dyl/project/singlecell/GGO/3.seraut/3.integrate/GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt",sep="\t",quote = "",stringsAsFactors=F,header=T,row.names=1) 
pheno <- read.table(file = "/NAS/dyl/project/singlecell/GGO/3.seraut/3.integrate/GSE131907/GSE131907_Lung_Cancer_cell_annotation.txt",sep="\t",quote = "",stringsAsFactors=F,header=T)
smpsel <- pheno[pheno[,"Cell_type"]%in%c("Epithelial cells"),1]
ref.data <- ref.data[,smpsel]
ref <- CreateSeuratObject(counts = ref.data, project = "NC", min.cells = 3, min.features = 200)
ref[["percent.mt"]] <- PercentageFeatureSet(ref, pattern = "^MT-")
ref <- subset(ref, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
ref <-subset(ref,features=setdiff(rownames(ref@assays[[1]]),
grep("MT-",rownames(ref@assays[[1]]), value = T))) 
ref <- NormalizeData(ref, verbose = FALSE)
ref <- ScaleData(ref)
ref <- FindVariableFeatures(ref)
ref <- RunPCA(ref)
print(ElbowPlot(ref,ndims=50))
ref <- RunUMAP(ref, dims = 1:50, verbose = FALSE)
ref <- RunTSNE(ref, dims = 1:50, nthreads = 4, max_iter = 2000)
ref$cluster<- pheno[match(rownames(ref@meta.data),pheno[,1]),"Sample_Origin"]
load("expr_mt.RData")
expr_mt <- as.matrix(expr_mt)
load("model/expression_model.RData")
model <- coef(cvfit.lasso, s = "lambda.min")
lasso_C <- apply(expr_mt[rownames(model)[model[,1]!=0][-1],],2,mean)
tt1 <- as.numeric(lasso_C)
names(tt1) <- colnames(expr_mt)
ref$signature<- tt1[rownames(ref@meta.data)]
refsub <- subset(ref,cells=rownames(ref@meta.data)[ref@meta.data[,"cluster"]%in%c("nLung","tLung","mLN","mBrain")])
score_df <- data.frame(score=refsub@meta.data[,"signature"],
label=factor(refsub@meta.data[,"cluster"],levels=c("nLung","tLung","mLN","mBrain")),stringsAsFactors=F)
colo2 <- brewer.pal(6, "Set1")[c(3,4,2,1)]

#plot Figure 6G
dp <- ggplot(score_df, aes(x=label, y=score, fill=label)) + 
  geom_violin()+scale_fill_manual(values=colo2)+geom_boxplot(width=0.1,fill="white")+
  labs(title="target signature",y = "mean expression") + theme_classic()
dp

#plot Figure S5A
DimPlot(ref, group.by="cluster",reduction = "umap",cols=colo,label = FALSE) 

#plot Figure S5B
refsub <- subset(ref,cells=rownames(ref@meta.data)[ref@meta.data[,"cluster"]%in%c("nLung","tLung","mLN",
ref$signature<- tt1[rownames(ref@meta.data)]
FeaturePlot(refsub, features = "signature",min.cutoff="q10") + 
    RotatedAxis()


##################
##Figure 6H,S5C-D#
##################
#cluster of patient 19
library(Seurat)
library(dplyr)
library(patchwork)
ref.data <- readRDS(file = "/NAS/dyl/project/singlecell/GGO/3.seraut/3.integrate/GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.rds")
pheno <- read.table(file = "/NAS/dyl/project/singlecell/GGO/3.seraut/3.integrate/GSE131907/GSE131907_Lung_Cancer_cell_annotation.txt",sep="\t",quote = "",stringsAsFactors=F,header=T)
smpsel19 <- pheno[(pheno[,"Cell_type"]%in%c("Epithelial cells"))&(pheno[,"Sample"]%in%c("LUNG_N19",
"LUNG_T19","EBUS_19","NS_19")),1]
smp19.data <- ref.data[,smpsel19]
pheno_sub <- pheno[match(smpsel19,pheno[,1]),] 
sample_list <- c("LUNG_N19","LUNG_T19","EBUS_19","NS_19")
data_list <- lapply(seq(length(sample_list)),function(x){ 
	smp19.data_sub <- smp19.data[,pheno_sub[,"Sample"]==sample_list[x]]
	pbmc <- CreateSeuratObject(counts = smp19.data_sub, project = sample_list[x], min.cells = 3, min.features = 200)
	pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
	pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
	pbmc <-subset(pbmc,features=setdiff(rownames(pbmc@assays[[1]]),
	grep("MT-",rownames(pbmc@assays[[1]]), value = T))) 
	pbmc <- NormalizeData(pbmc, verbose = FALSE)
    pbmc <- FindVariableFeatures(pbmc, selection.method = "vst",  verbose = FALSE)
	return(pbmc)
})
names(data_list) <- sample_list
anchors <- FindIntegrationAnchors(object.list = data_list,  dims = 1:50,k.filter = 20)
smp19.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)
DefaultAssay(smp19.integrated) <- "integrated"
smp19.integrated <- ScaleData(smp19.integrated, verbose = FALSE)
smp19.integrated <- RunPCA(smp19.integrated,npcs=50,verbose = FALSE)
print(ElbowPlot(smp19.integrated,ndims=50))
smp19.integrated <- RunUMAP(smp19.integrated, dims = 1:40, verbose = FALSE)
smp19.integrated <- RunTSNE(smp19.integrated, dims = 1:40, nthreads = 4, max_iter = 2000)
smp19.integrated <- FindNeighbors(smp19.integrated, dims = 1:40, verbose = FALSE)
smp19.integrated <- FindClusters(smp19.integrated, verbose = FALSE,resolution=0.8)
smp19.integrated$cluster<- pheno[match(rownames(smp19.integrated@meta.data),pheno[,1]),"Sample_Origin"]

#plot Figure S5C
DimPlot(smp19.integrated, group.by="cluster",reduction = "umap",label = T) 

#copy number alteration of LUNG_T19 by copykat
library(Seurat)
library(copykat)
ref.data <- readRDS(file = "/NAS/dyl/project/singlecell/GGO/3.seraut/3.integrate/GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.rds") 
pheno <- read.table(file = "/NAS/dyl/project/singlecell/GGO/3.seraut/3.integrate/GSE131907/GSE131907_Lung_Cancer_cell_annotation.txt",sep="\t",quote = "",stringsAsFactors=F,header=T)
source("/NAS/dyl/project/singlecell/GGO/13.copyKAT/1.study/copykat1.R")
smpsel19 <- pheno[(pheno[,"Sample"]%in%c("LUNG_T19")),1]
smp19.data <- ref.data[,smpsel19]
pheno_sub <- pheno[match(smpsel19,pheno[,1]),]
exp.rawdata <- as.matrix(smp19.data)
copykat.smp19t <- copykat1(rawmat=exp.rawdata, id.type="S", ngene.chr=5, win.size=15, KS.cut=0.05, sam.name="LUNG_T19", distance="euclidean", 
norm.cell.names=colnames(exp.rawdata)[pheno_sub[,"Cell_type"]!="Epithelial cells"],n.cores=32)
saveRDS(copykat.smp19t, file = "copykat.smp19t.rds")

#copy number alteration of EBUS_19 by copykat
library(Seurat)
library(copykat)
ref.data <- readRDS(file = "/NAS/dyl/project/singlecell/GGO/3.seraut/3.integrate/GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.rds")
pheno <- read.table(file = "/NAS/dyl/project/singlecell/GGO/3.seraut/3.integrate/GSE131907/GSE131907_Lung_Cancer_cell_annotation.txt",sep="\t",quote = "",stringsAsFactors=F,header=T)
source("/NAS/dyl/project/singlecell/GGO/13.copyKAT/1.study/copykat1.R")
smpsel19 <- pheno[(pheno[,"Sample"]%in%c("EBUS_19")),1]
smp19.data <- ref.data[,smpsel19]
pheno_sub <- pheno[match(smpsel19,pheno[,1]),]
exp.rawdata <- as.matrix(smp19.data)
copykat.smp19mln <- copykat1(rawmat=exp.rawdata, id.type="S", ngene.chr=5, win.size=15, KS.cut=0.05, sam.name="EBUS_19", distance="euclidean", 
norm.cell.names=colnames(exp.rawdata)[pheno_sub[,"Cell_type"]!="Epithelial cells"],n.cores=32)
saveRDS(copykat.smp19mln, file = "copykat.smp19mln.rds")

#copy number alteration of NS_19 by copykat
library(Seurat)
library(copykat)
ref.data <- readRDS(file = "/NAS/dyl/project/singlecell/GGO/3.seraut/3.integrate/GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.rds") 
pheno <- read.table(file = "/NAS/dyl/project/singlecell/GGO/3.seraut/3.integrate/GSE131907/GSE131907_Lung_Cancer_cell_annotation.txt",sep="\t",quote = "",stringsAsFactors=F,header=T)
source("/NAS/dyl/project/singlecell/GGO/13.copyKAT/1.study/copykat1.R")
smpsel19 <- pheno[(pheno[,"Sample"]%in%c("NS_19")),1]
smp19.data <- ref.data[,smpsel19]
pheno_sub <- pheno[match(smpsel19,pheno[,1]),]
exp.rawdata <- as.matrix(smp19.data)
copykat.smp19brain <- copykat1(rawmat=exp.rawdata, id.type="S", ngene.chr=5, win.size=15, KS.cut=0.05, sam.name="NS_19", distance="euclidean", 
norm.cell.names=colnames(exp.rawdata)[pheno_sub[,"Cell_type"]!="Epithelial cells"],n.cores=32)
saveRDS(copykat.smp19brain, file = "copykat.smp19brain.rds")

#cluster of cancer cell from LUNG_T19
library(Seurat)
library(dplyr)
library(patchwork)
library(RColorBrewer)
colo=c(brewer.pal(6, "Set1")[1:4],brewer.pal(6, "Set2")[1:6])

smp19 <- readRDS(file = "/NAS/dyl/project/chennan/5.expression/sc19_cluster/smp19.integrated.rds")
copykat.smp19t <- readRDS(file = "/NAS/dyl/project/chennan/5.expression/sc19_cluster/copykat.smp19t.rds")
copykat.smp19mln <- readRDS(file = "/NAS/dyl/project/chennan/5.expression/sc19_cluster/copykat.smp19mln.rds")
copykat.smp19brain <- readRDS(file = "/NAS/dyl/project/chennan/5.expression/sc19_cluster/copykat.smp19brain.rds")
pheno <- read.table(file = "/NAS/dyl/project/singlecell/GGO/3.seraut/3.integrate/GSE131907/GSE131907_Lung_Cancer_cell_annotation.txt",sep="\t",quote = "",stringsAsFactors=F,header=T)
smp19$sample<- pheno[match(rownames(smp19@meta.data),pheno[,1]),"Sample"]

sample_list <- list()
sample_list[[1]] <- copykat.smp19t
sample_list[[2]] <- copykat.smp19mln
sample_list[[3]] <- copykat.smp19brain
names(sample_list) <- c("LUNG_T19","EBUS_19","NS_19")

anno_label_list <- c() 
for(smp in names(sample_list)) 
{
	pbmc <- sample_list[smp][[1]]
	tumor.cells <- pbmc["prediction"][[1]][which(pbmc["prediction"][[1]][,"copykat.pred"]=="aneuploid"),"cell.names"]
	anno_label <- rownames(smp19@meta.data)[as.character(smp19@meta.data[,"sample"])==smp]
	intersmp <- intersect(anno_label,tumor.cells)
	names(intersmp) <- NULL
	anno_label_list <- c(anno_label_list,intersmp)
}

label <- rep("normal",nrow(smp19@meta.data))
names(label) <- rownames(smp19@meta.data)
label[names(label)%in%anno_label_list] <- "tumor"

smp19$cnv <- label

smp19cancer <-subset(smp19,cells=rownames(smp19@meta.data)[(as.character(smp19@meta.data[,"cnv"])=="tumor")]) 

smp19cancer <- RunPCA(smp19cancer, features = VariableFeatures(object = smp19cancer))
smp19cancer <- RunUMAP(smp19cancer, dims = 1:40, verbose = FALSE)
smp19cancer <- RunTSNE(smp19cancer, dims = 1:40, nthreads = 4, max_iter = 2000)
smp19cancer <- FindNeighbors(smp19cancer, dims = 1:40, verbose = FALSE)
smp19cancer <- FindClusters(smp19cancer, verbose = FALSE,resolution=3)

smp19pri <-subset(smp19cancer,cells=rownames(smp19cancer@meta.data)[(as.character(smp19cancer@meta.data[,"cluster"])=="tLung")]) 
smp19pri <- RunPCA(smp19pri, features = VariableFeatures(object = smp19pri))
smp19pri <- RunUMAP(smp19pri, dims = 1:40, verbose = FALSE)
smp19pri <- RunTSNE(smp19pri, dims = 1:40, nthreads = 4, max_iter = 2000)
smp19pri <- FindNeighbors(smp19pri, dims = 1:40, verbose = FALSE)
smp19pri <- FindClusters(smp19pri, verbose = FALSE,resolution=1.2)
smp19pri$cl <- paste0("c",as.numeric(as.character(smp19pri@meta.data[,"seurat_clusters"]))+1)

#plot Figure S5D
DimPlot(smp19pri,group.by="cl", reduction = "umap",cols=colo[c(1:4,6)],label = T) 

#plot Figure 6H
clss <- as.character(smp19pri@meta.data[,"cl"])
names(clss) <- rownames(smp19pri@meta.data)
cluster <- as.character(smp19cancer@meta.data[,"cluster"])
names(cluster) <- rownames(smp19cancer@meta.data)
cluster[names(clss)] <- clss
smp19cancer$"cluster" <- cluster
DimPlot(smp19cancer, reduction = "umap",group.by="cluster",
cells=rownames(smp19cancer@meta.data)[cluster!="c5"],label = T,cols=colo[c(5:8,1,2)]) 

############
##Figure 6I#
############
library(destiny)
library(ggplot2)
library(RColorBrewer)
colo=c(brewer.pal(6, "Set1")[1:4],brewer.pal(6, "Set2")[1:6])[c(5:8,1,2)]
load("smp19cancer.integrated_pri_ln_brain.RData")
expMt <- as.matrix(smp19cancer.integrated[,!cluster%in%c("c5")])
cluster <- cluster[!cluster%in%c("c5")]
dif <- DiffusionMap(t(expMt),n_pcs=40)
dpt <- DPT(dif)
names(colo) <- unique(cluster)
colov <- colo[cluster]
qplot(DC1, DC2, data = dif, colour = cluster) 

##################
##Figure 6J,S5E-F#
##################
setwd("/NAS/dyl/project/chennan/5.expression/sc19_cluster/plot")
library(destiny)
library(ggplot2)
library(RColorBrewer)
library(sROC)
library(fgsea)
colo=c(brewer.pal(6, "Set1")[1:4],brewer.pal(6, "Set2")[1:6])[c(5:8,1,2)]
load("smp19cancer.integrated_pri_ln_brain.RData")
expMt <- as.matrix(smp19cancer.integrated[,!cluster%in%c("c5","c1","mBrain")])
cluster <- cluster[!cluster%in%c("c5","c1","mBrain")]
dif <- DiffusionMap(t(expMt),n_pcs=40)
dpt <- DPT(dif)
names(colo) <- unique(cluster)
colov <- colo[cluster]
#plot Figure 6J
qplot(DC2, DC3, data = dif, colour = cluster) 

#plot Figure S5E
dc_exp <- as.data.frame(dif)
dc2 <- dc_exp[,2]
dc2_stage <- split(dc2,f=factor(cluster))
plot(kCDF(dc2_stage[[1]]),col=brewer.pal(6, "Set2")[2],CI=F,xlim=c(min(dc_exp[,2]),max(dc_exp[,2])),
main="cumulative distribution function",xlab="DC2")
par(new=T)
plot(kCDF(dc2_stage[[2]]),col=brewer.pal(6, "Set2")[3],CI=F,xlim=c(min(dc_exp[,2]),max(dc_exp[,2])),
axes = F,main="",xlab="",ylab="")
par(new=T)
plot(kCDF(dc2_stage[[3]]),col=brewer.pal(6, "Set2")[4],CI=F,xlim=c(min(dc_exp[,2]),max(dc_exp[,2])),axes = F,main="",xlab="",ylab="")
par(new=T)
plot(kCDF(dc2_stage[[4]]),col=brewer.pal(6, "Set2")[2],CI=F,xlim=c(min(dc_exp[,2]),max(dc_exp[,2])),axes = F,
main="",xlab="",ylab="")
legend("bottomright",names(dc2_stage),col=c(brewer.pal(6, "Set2")[2:4],brewer.pal(6, "Set2")[2]))

#plot Figure S5F
gene_cor <- sapply(rownames(expMt),function(x){
	return(cor(unlist(expMt[x,]),dc_exp[,2]))
})
names(gene_cor) <- rownames(expMt)
gene_cor_sel <- gene_cor[order(abs(gene_cor),decreasing=T)][1:200]
gmt.file_kegg <- "/NAS/dyl/source/gsea/db/c2.cp.kegg.v6.1.symbols.gmt"
pathways_go <- gmtPathways(gmt.file_go)
fgseaRes_list_go <- as.data.frame(fgsea(pathways = pathways_go, stats = gene_cor_sel,minSize=3,maxSize=500,nperm=1000))
fgseaRes_list_go <- fgseaRes_list_go[order(abs(fgseaRes_list_go[,"NES"]),
decreasing=T),]
topPathways <- fgseaRes_list_go[(fgseaRes_list_go[,"pval"]<0.05),1]
seltopPathways <- c("GO_OXIDATION_REDUCTION_PROCESS","GO_CELL_CYCLE_PROCESS","GO_INNATE_IMMUNE_RESPONSE","GO_MICROTUBULE_BASED_PROCESS","GO_CELLULAR_GLUCOSE_HOMEOSTASIS")
plotGseaTable(pathways_go[seltopPathways], gene_cor_sel, fgseaRes_list_go, 
              gseaParam = 0.5)


##################
##Figure 6K,S5G-H#
##################
library(destiny)
library(ggplot2)
library(RColorBrewer)
library(sROC)
library(fgsea)
colo=c(brewer.pal(6, "Set1")[1:4],brewer.pal(6, "Set2")[1:6])[c(5:8,1,2)]
load("smp19cancer.integrated_pri_ln_brain.RData")
expMt <- as.matrix(smp19cancer.integrated[,!cluster%in%c("c5","c1","mLN")])
cluster <- cluster[!cluster%in%c("c5","c1","mLN")]
dif <- DiffusionMap(t(expMt),n_pcs=40)
dpt <- DPT(dif)
names(colo) <- unique(cluster)
colov <- colo[cluster]
#plot Figure 6K
qplot(DC1, DC3, data = dif, colour = cluster) 

#plot Figure S5G
dc_exp <- as.data.frame(dif)
dc3 <- dc_exp[,3]
dc3_stage <- split(dc3,f=factor(cluster))
plot(kCDF(dc3_stage[[1]]),col=1,CI=F,xlim=c(min(dc_exp[,3]),max(dc_exp[,3])),
main="cumulative distribution function",xlab="DC3")
par(new=T)
plot(kCDF(dc3_stage[[2]]),col=2,CI=F,xlim=c(min(dc_exp[,3]),max(dc_exp[,3])),
axes = F,main="",xlab="",ylab="")
par(new=T)
plot(kCDF(dc3_stage[[3]]),col=3,CI=F,xlim=c(min(dc_exp[,3]),max(dc_exp[,3])),axes = F,main="",xlab="",ylab="")
par(new=T)
plot(kCDF(dc3_stage[[4]]),col=4,CI=F,xlim=c(min(dc_exp[,3]),max(dc_exp[,3])),axes = F,
main="",xlab="",ylab="")
legend("bottomright",names(dc3_stage),col=c(brewer.pal(6, "Set2")[2:4],brewer.pal(6, "Set2")[2]))

#plot Figure S5H
gene_cor <- sapply(rownames(expMt),function(x){
	return(cor(unlist(expMt[x,]),dc_exp[,3]))
})
names(gene_cor) <- rownames(expMt)
gene_cor_sel <- gene_cor[order(abs(gene_cor),decreasing=T)][1:200]


gmt.file_go <- "/NAS/dyl/source/gsea/db/c5.bp.v6.1.symbols.gmt"
gmt.file_kegg <- "/NAS/dyl/source/gsea/db/c2.cp.kegg.v6.1.symbols.gmt"
##2.3.runGSEA
pathways_go <- gmtPathways(gmt.file_go)
fgseaRes_list_go <- as.data.frame(fgsea(pathways = pathways_go, stats = gene_cor_sel,minSize=3,maxSize=500,nperm=1000))
fgseaRes_list_go <- fgseaRes_list_go[order(abs(fgseaRes_list_go[,"NES"]),
decreasing=T),]
topPathways <- fgseaRes_list_go[(fgseaRes_list_go[,"pval"]<0.05),1]
seltopPathways <- c("GO_INNATE_IMMUNE_RESPONSE","GO_CELL_MOTILITY",
"GO_POSITIVE_REGULATION_OF_CELL_PROLIFERATION","GO_NEGATIVE_REGULATION_OF_CELL_CELL_ADHESION")
pdf("plotGseaTable_mBrain_DC3.pdf")
plotGseaTable(pathways_go[seltopPathways], gene_cor_sel, fgseaRes_list_go, 
              gseaParam = 0.5)

############
##Figure 6L#
############
library(Seurat)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(ggplot2)
colo=brewer.pal(6, "Set2")[1:6]
smp19 <- readRDS(file = "/NAS/dyl/project/chennan/5.expression/sc19_cluster/smp19.integrated.rds")
copykat.smp19t <- readRDS(file = "/NAS/dyl/project/chennan/5.expression/sc19_cluster/copykat.smp19t.rds")
copykat.smp19mln <- readRDS(file = "/NAS/dyl/project/chennan/5.expression/sc19_cluster/copykat.smp19mln.rds")
copykat.smp19brain <- readRDS(file = "/NAS/dyl/project/chennan/5.expression/sc19_cluster/copykat.smp19brain.rds")
pheno <- read.table(file = "/NAS/dyl/project/singlecell/GGO/3.seraut/3.integrate/GSE131907/GSE131907_Lung_Cancer_cell_annotation.txt",sep="\t",quote = "",stringsAsFactors=F,header=T)
smp19$sample<- pheno[match(rownames(smp19@meta.data),pheno[,1]),"Sample"]
sample_list <- list()
sample_list[[1]] <- copykat.smp19t
sample_list[[2]] <- copykat.smp19mln
sample_list[[3]] <- copykat.smp19brain
names(sample_list) <- c("LUNG_T19","EBUS_19","NS_19")
anno_label_list <- c() 
for(smp in names(sample_list)) 
{
	pbmc <- sample_list[smp][[1]]
	tumor.cells <- pbmc["prediction"][[1]][which(pbmc["prediction"][[1]][,"copykat.pred"]=="aneuploid"),"cell.names"]
	anno_label <- rownames(smp19@meta.data)[as.character(smp19@meta.data[,"sample"])==smp]
	intersmp <- intersect(anno_label,tumor.cells)
	names(intersmp) <- NULL
	anno_label_list <- c(anno_label_list,intersmp)
}
label <- rep("normal",nrow(smp19@meta.data))
names(label) <- rownames(smp19@meta.data)
label[names(label)%in%anno_label_list] <- "tumor"
smp19$cnv <- label
smp19cancer <-subset(smp19,cells=rownames(smp19@meta.data)[(as.character(smp19@meta.data[,"cnv"])=="tumor")]) 
smp19cancer <- RunPCA(smp19cancer, features = VariableFeatures(object = smp19cancer))
smp19cancer <- RunUMAP(smp19cancer, dims = 1:40, verbose = FALSE)
smp19cancer <- RunTSNE(smp19cancer, dims = 1:40, nthreads = 4, max_iter = 2000)
smp19cancer <- FindNeighbors(smp19cancer, dims = 1:40, verbose = FALSE)
smp19cancer <- FindClusters(smp19cancer, verbose = FALSE,resolution=3)
smp19pri <-subset(smp19cancer,cells=rownames(smp19cancer@meta.data)[(as.character(smp19cancer@meta.data[,"cluster"])=="tLung")]) 
smp19pri <- RunPCA(smp19pri, features = VariableFeatures(object = smp19pri))
smp19pri <- RunUMAP(smp19pri, dims = 1:40, verbose = FALSE)
smp19pri <- RunTSNE(smp19pri, dims = 1:40, nthreads = 4, max_iter = 2000)
smp19pri <- FindNeighbors(smp19pri, dims = 1:40, verbose = FALSE)
smp19pri <- FindClusters(smp19pri, verbose = FALSE,resolution=1.2)
smp19pri$cl <- paste0("c",as.numeric(as.character(smp19pri@meta.data[,"seurat_clusters"]))+1)
score_df2 <- data.frame(score=smp19pri@meta.data[,"signature"],
label=as.character(smp19pri@meta.data[,"cl"]),stringsAsFactors=F)
score_df2 <- score_df2[score_df2[,2]!="c5",]
score_df2[,2] <- factor(score_df2[,2],levels=c("c2","c1","c3","c4"))
score_df3 <- score_df2
score_df3[,2] <- as.character(score_df3[,2])
score_df3[score_df3[,2]%in%c("c1","c3","c4"),2] <- "c1+c3+c4"
score_df3[,2] <- factor(score_df3[,2],levels=c("c2","c1+c3+c4"))
pie(table(score_df3[as.character(score_df3[,2])=="c2",1]>0))
pie(table(score_df3[as.character(score_df3[,2])=="c1+c3+c4",1]>0))

################
##Figure 6M,S5I#
################
library(Seurat)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(ggplot2)
colo=brewer.pal(6, "Set2")[1:6]
smp19 <- readRDS(file = "/NAS/dyl/project/chennan/5.expression/sc19_cluster/smp19.integrated.rds")
copykat.smp19t <- readRDS(file = "/NAS/dyl/project/chennan/5.expression/sc19_cluster/copykat.smp19t.rds")
copykat.smp19mln <- readRDS(file = "/NAS/dyl/project/chennan/5.expression/sc19_cluster/copykat.smp19mln.rds")
copykat.smp19brain <- readRDS(file = "/NAS/dyl/project/chennan/5.expression/sc19_cluster/copykat.smp19brain.rds")
pheno <- read.table(file = "/NAS/dyl/project/singlecell/GGO/3.seraut/3.integrate/GSE131907/GSE131907_Lung_Cancer_cell_annotation.txt",sep="\t",quote = "",stringsAsFactors=F,header=T)
smp19$sample<- pheno[match(rownames(smp19@meta.data),pheno[,1]),"Sample"]
sample_list <- list()
sample_list[[1]] <- copykat.smp19t
sample_list[[2]] <- copykat.smp19mln
sample_list[[3]] <- copykat.smp19brain
names(sample_list) <- c("LUNG_T19","EBUS_19","NS_19")
anno_label_list <- c() 
for(smp in names(sample_list)) 
{
	pbmc <- sample_list[smp][[1]]
	tumor.cells <- pbmc["prediction"][[1]][which(pbmc["prediction"][[1]][,"copykat.pred"]=="aneuploid"),"cell.names"]
	anno_label <- rownames(smp19@meta.data)[as.character(smp19@meta.data[,"sample"])==smp]
	intersmp <- intersect(anno_label,tumor.cells)
	names(intersmp) <- NULL
	anno_label_list <- c(anno_label_list,intersmp)
}
label <- rep("normal",nrow(smp19@meta.data))
names(label) <- rownames(smp19@meta.data)
label[names(label)%in%anno_label_list] <- "tumor"
smp19$cnv <- label
smp19cancer <-subset(smp19,cells=rownames(smp19@meta.data)[(as.character(smp19@meta.data[,"cnv"])=="tumor")]) 
smp19cancer <- RunPCA(smp19cancer, features = VariableFeatures(object = smp19cancer))
smp19cancer <- RunUMAP(smp19cancer, dims = 1:40, verbose = FALSE)
smp19cancer <- RunTSNE(smp19cancer, dims = 1:40, nthreads = 4, max_iter = 2000)
smp19cancer <- FindNeighbors(smp19cancer, dims = 1:40, verbose = FALSE)
smp19cancer <- FindClusters(smp19cancer, verbose = FALSE,resolution=3)
smp19pri <-subset(smp19cancer,cells=rownames(smp19cancer@meta.data)[(as.character(smp19cancer@meta.data[,"cluster"])=="tLung")]) 
smp19pri <- RunPCA(smp19pri, features = VariableFeatures(object = smp19pri))
smp19pri <- RunUMAP(smp19pri, dims = 1:40, verbose = FALSE)
smp19pri <- RunTSNE(smp19pri, dims = 1:40, nthreads = 4, max_iter = 2000)
smp19pri <- FindNeighbors(smp19pri, dims = 1:40, verbose = FALSE)
smp19pri <- FindClusters(smp19pri, verbose = FALSE,resolution=1.2)
smp19pri$cl <- paste0("c",as.numeric(as.character(smp19pri@meta.data[,"seurat_clusters"]))+1)

#plot Figure S5I
score_df2 <- data.frame(score=smp19pri@meta.data[,"signature"],
label=as.character(smp19pri@meta.data[,"cl"]),stringsAsFactors=F)
score_df2 <- score_df2[score_df2[,2]!="c5",]
score_df2[,2] <- factor(score_df2[,2],levels=c("c2","c1","c3","c4"))
dp <- ggplot(score_df2, aes(x=label, y=score, fill=label)) + 
  geom_violin()+scale_fill_manual(values=colo[c(2,1,3,4)])+geom_boxplot(width=0.1,fill="white")+
  labs(title=levels(score_df2[,2]),y = "expression") + theme_classic()
dp

#plot Figure 6M
score_df3 <- score_df2
score_df3[,2] <- as.character(score_df3[,2])
score_df3[score_df3[,2]%in%c("c1","c3","c4"),2] <- "c1+c3+c4"
score_df3[,2] <- factor(score_df3[,2],levels=c("c2","c1+c3+c4"))
dp <- ggplot(score_df3, aes(x=label, y=score, fill=label)) + 
  geom_violin()+scale_fill_manual(values=brewer.pal(8, "Set2")[c(2,8)])+geom_boxplot(width=0.1,fill="white")+
  labs(title=levels(score_df3[,2]),y = "expression") + theme_classic()
dp
