#################################################################################################
###Figure 5A The most significantly enriched Gene Ontology (GO) terms for the genes that showed #
###significant correlation with model score.                                                    # 
###Figure 5B The genes in the most significantly enriched GO terms.                             #
###Figure 5C The locus zoom for a methylation locus (chr17:76561412).                           #
###Figure 5D-E The expression of CYGB was negative correlated with the methylation of           #
###chr17:76561412 in our data and TCGA.                                                         #
#################################################################################################

############
##Figure 5A#
############
#co-expression network
options(stringsAsFactors=F)
options(scipen=200)
library(GenomicRanges)
library(glmnet)
library(gplots)
library(pROC)
library(RColorBrewer)
library(ROCR)
load("/data1/dengyl/project/S1/chennan/new2/Tmt.RData")
load("cvfit.lasso.RData")
load("/data1/dengyl/project/S1/chennan/new3/m100_mt.RData")
rownames(Tmt) <- sapply(strsplit(rownames(Tmt),".",fixed=T),function(x) x[1])
m100_mt <- m100_mt[rownames(model)[-1],]
lasso_C <- predict(cvfit.lasso, newx = t(m100_mt),s = "lambda.min",type = "response")
lasso_C  <- lasso_C[,1]
cor_p <- t(sapply(seq(nrow(Tmt)),function(x){
	t1 <- cor.test(lasso_C,unlist(Tmt[x,]))
	return(c(t1$"estimate",t1$"p.value"))
}))
colnames(cor_p) <- c("estiamte","pvalue")
cor_p <- cbind(rownames(Tmt),cor_p)
padj <- p.adjust(cor_p[,"pvalue"],method="fdr")
cor_p <- cbind(cor_p,padj)
cor_p <- cor_p[cor_p[,"padj"]<0.05,]

#enrichment was performed by DAVID

#plot
library(plotrix)
options(stringsAsFactors=F)
enrich <- read.table(file="DAVID.txt",
sep="\t",quote = "",stringsAsFactors=F,header=T)
enrich <- enrich[order(enrich$"Benjamini"),]
GOname <- enrich[,"Term"]
GOnameMt <- matrix(unlist(strsplit(GOname,"~",fixed=T)),byrow=T,ncol=2)
GOnm <- GOnameMt[,2]
fc <- enrich[,"Fold.Enrichment"]
names(fc) <- GOnm
barplot(rev(fc),xlab="Index",xlim=c(0,4),las=2,
  ylab="Fold Enrichment",main="Enrich GO Terms",col=rep("#bebebe",length(fc)),horiz = T)
p_adjust <- enrich[,"Benjamini"]
log10p <- -log10(p_adjust)
names(log10p) <- GOnm
barplot(rev(log10p),  axes=F,  ylab="", xlab="adjusted P value",las=2,col="#7e8000",xlim=c(0,8),horiz = T)
box()
axis(1,at=c(0,  2,  4,6,8),labels=c("1","0.01","0.0001","10-6","10-8"))


############
##Figure 5B#
############
#prepare for cytoscape
options(stringsAsFactors=F)
options(scipen=200)
library(gplots)
library(RColorBrewer)
cor_p <- read.table(file="cor_p.txt",
sep="\t",quote = "",stringsAsFactors=F,header=T)
DAVID <- read.table(file="DAVID.txt",
sep="\t",stringsAsFactors=F,T)
gtfl <- read.table(file="/data1/ref/hg38/gencode.v28.annotation.gtf",sep="\t",quote = "",stringsAsFactors=F,header=F)
genegtf <- gtfl[gtfl[,3]=="gene",]
genegtf_v9 <- strsplit(genegtf[,9],"\"| |;")
pos_g <- sapply(genegtf_v9,function(x) which(x=="gene_id"))
gene_id <- sapply(seq(length(genegtf_v9)),function(x) genegtf_v9[[x]][pos_g[x]+2])
pos_t <- sapply(genegtf_v9,function(x) which(x=="gene_name"))
gene_name <- sapply(seq(length(genegtf_v9)),function(x) genegtf_v9[[x]][pos_t[x]+2])
id_name <- data.frame(gene_id=gene_id,gene_name=gene_name,stringsAsFactors=F)
id_name[,"gene_id"] <- substr(id_name[,"gene_id"],1,15)
cor_p <- cor_p[order(abs(cor_p[,2]),decreasing=T),]
gene_list <- strsplit(DAVID[,"Genes"],", ")
names(gene_list) <- DAVID[,1]
gene_listsub <- lapply(gene_list,function(x) 
{
	t1 <- x[x%in%cor_p[1:200,1]]
	t2 <- intersect(t1,id_name[,1])
	if(length(t2)>0)
	{
		return(id_name[match(t2,id_name[,1]),2])
	}
})
names(gene_listsub) <- names(gene_list)
gene_listsub <- gene_listsub[sapply(gene_listsub,length)>5]
gene2 <- unique(unlist(gene_listsub))
gattri <- sapply(gene2,function(x) names(gene_listsub)[sapply(gene_listsub,function(y) x%in%y)][1])
geneS_df <- data.frame(gene=gene2,att=gattri)
write.table(geneS_df, file="node.attri.txt", sep="\t" , quote=FALSE, row.names=F, col.names=T )
CorS <- cor_p[match(id_name[match(gene2,id_name[,2]),1],cor_p[,1]),2]
score <- rep(1,length(CorS))
score[CorS<0] <- 0
CorS_df <- data.frame(lasso="model",gene=gene2,score=score)
write.table(CorS_df, file="edge.attri.txt", sep="\t" , quote=FALSE, row.names=F, col.names=T )

#visulasition by cytoscape

############
##Figure 5C#
############
#visulation in UCSC

##############
##Figure 5D-E#
##############
options(stringsAsFactors=F)
options(scipen=200)
library(glmnet)
library(GenomicRanges)
library(RColorBrewer)
setwd("/data1/dengyl/project/S1/chennan/new3/enumer_lung_enhancer")
load("/data1/dengyl/project/S1/chennan/new3/m100_mt.RData")
load("cvfit.lasso.RData")
load("/data1/dengyl/project/S1/chennan/new2/Tmt.RData")
load("/data1/dengyl/project/S1/chennan/new2/tcga_mt100.RData")
load("/data1/dengyl/download/TCGA/LUAD/mRNAexpression/FPKMmt_geneName_log2.RData")

gtfl <- read.table(file="/data1/ref/hg38/gencode.v28.annotation.gtf",sep="\t",quote = "",stringsAsFactors=F,header=F)
genegtf <- gtfl[gtfl[,3]=="gene",]
genegtf_v9 <- strsplit(genegtf[,9],"\"| |;")
pos_g <- sapply(genegtf_v9,function(x) which(x=="gene_id"))
gene_id <- sapply(seq(length(genegtf_v9)),function(x) genegtf_v9[[x]][pos_g[x]+2])
pos_t <- sapply(genegtf_v9,function(x) which(x=="gene_type"))
gene_type <- sapply(seq(length(genegtf_v9)),function(x) genegtf_v9[[x]][pos_t[x]+2])
pos_n <- sapply(genegtf_v9,function(x) which(x=="gene_name"))
gene_name <- sapply(seq(length(genegtf_v9)),function(x) genegtf_v9[[x]][pos_n[x]+2])
gtf <- data.frame(gene=gene_id,gene_type=gene_type,gene_name=gene_name,
chr=genegtf[,1],start=genegtf[,4],end=genegtf[,5],strand=genegtf[,7],stringsAsFactors=F)
width <- gtf[,"end"]-gtf[,"start"]+1
gtf <- gtf[width>200,]
gtf_p <- gtf[gtf[,"strand"]=="+",]
enhancerd_p <- data.frame(gene=gtf_p[,"gene"],gene_type=gtf_p[,"gene_type"],gene_name=gtf_p[,"gene_name"],
chr=gtf_p[,"chr"],start=gtf_p[,"start"]-100000,end=gtf_p[,"end"]+100000,strand=gtf_p[,"strand"],stringsAsFactors=F)
gtf_n <- gtf[gtf[,"strand"]=="-",]
enhancerd_n <- data.frame(gene=gtf_n[,"gene"],gene_type=gtf_n[,"gene_type"],gene_name=gtf_n[,"gene_name"],
chr=gtf_n[,"chr"],start=gtf_n[,"start"]-100000,end=gtf_n[,"end"]+100000,strand=gtf_n[,"strand"],stringsAsFactors=F)
enhancerd <- rbind(enhancerd_p,enhancerd_n,stringsAsFactors=F)
enhancerd[,"gene"] <- substr(enhancerd[,"gene"],1,15)
enhancerd[enhancerd[,"start"] < 1,"start"] <- 1
region <- rownames(model)[-1]
region_bed <- matrix(unlist(strsplit(region,".",fixed=T)),byrow=T,ncol=3)
Gregion <- GRanges(seqnames = Rle(region_bed[,1]),
ranges = IRanges(as.numeric(region_bed[,2])-1,end = as.numeric(region_bed[,3]))) 
Gscreen_all_enhancerd <- GRanges(seqnames = Rle(enhancerd[,"chr"]),
ranges = IRanges(as.numeric(enhancerd[,"start"]),end = as.numeric(enhancerd[,"end"]))) 
mtch1 <- findOverlaps(Gregion,Gscreen_all_enhancerd)
pair1 <- as.matrix(mtch1)
meth_screen_all_enhancerd <- cbind(region[pair1[,1]],enhancerd[pair1[,2],],stringsAsFactors=F)
colnames(meth_screen_all_enhancerd)[1] <- "region"
rownames(Tmt) <- sapply(strsplit(rownames(Tmt),".",fixed=T),function(x) x[1])
meth_screen_all_enhancerd_exp <- meth_screen_all_enhancerd[meth_screen_all_enhancerd[,"gene"]%in%rownames(Tmt),]
gtfl <- read.table(file="/data1/ref/hg38/gencode.v28.annotation.gtf",sep="\t",quote = "",stringsAsFactors=F,header=F)
genegtf <- gtfl[gtfl[,3]=="gene",]
genegtf_v9 <- strsplit(genegtf[,9],"\"| |;")
pos_g <- sapply(genegtf_v9,function(x) which(x=="gene_id"))
gene_id <- sapply(seq(length(genegtf_v9)),function(x) genegtf_v9[[x]][pos_g[x]+2])
pos_t <- sapply(genegtf_v9,function(x) which(x=="gene_name"))
gene_name <- sapply(seq(length(genegtf_v9)),function(x) genegtf_v9[[x]][pos_t[x]+2])
id_name <- data.frame(gene_id=gene_id,gene_name=gene_name,stringsAsFactors=F)
id_name[,"gene_id"] <- substr(id_name[,"gene_id"],1,15)
cor_p <- t(sapply(seq(nrow(meth_screen_all_enhancerd_exp)),function(x){
	methp <- meth_screen_all_enhancerd_exp[x,1]
	genep <- meth_screen_all_enhancerd_exp[x,"gene"]
	t1 <- cor.test(unlist(m100_mt[methp,]),unlist(Tmt[genep,]), method =  "spearman",alternative="less")
	return(c(t1$"estimate",t1$"p.value"))
}))
colnames(cor_p) <- c("estiamte","pvalue")
meth_screen_all_enhancerd_exp <- cbind(meth_screen_all_enhancerd_exp,cor_p)
meth_screen_all_enhancerd_exp <- meth_screen_all_enhancerd_exp[(meth_screen_all_enhancerd_exp[,"pvalue"]<0.05),]
inter <- intersect(colnames(tcga_mt100),colnames(FPKMmt))
tcga_mt100 <- tcga_mt100[,inter]
FPKMmt <- FPKMmt[,inter]
#our data
x <- m100_mt["chr17.76561412.76561412",]
y <- log2(unlist(Tmt["ENSG00000161544",])+1)
plot(x,y,main="chr17.76561412.76561412_CYGB",pch=19,xlab="methylation (beta)",ylab="expression log2(TPM+1)")
abline(lm(formula = y ~ x),col="red")
Cor <- cor.test(x,y)
text(0.5,1,label=paste0("cor=",Cor$"estimate",",P=",Cor$"p.value"))
#TCGA data
x <- tcga_mt100["chr17.76561412.76561412",]
y <- FPKMmt["CYGB",]
plot(x,y,main="chr17.76561412.76561412_CYGB",pch=19,xlab="methylation (beta)",ylab="expression log2(TPM+1)")
abline(lm(formula = y ~ x),col="red")
Cor <- cor.test(x,y)
text(0.5,1,label=paste0("cor=",Cor$"estimate",",P=",Cor$"p.value"))
