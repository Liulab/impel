#########################################################################################
###Figure 2A The volcano plot for differentially methylated locus.                      #
###Figure 2B The The distribution of differentially methylated locus                    # 
###in functional regions.                                                               #
###Figure 2C The differentially methylated locus within enhancers                       # 
###were enriched in lung/stem cell specific enhancers                                   #
###Figure 2D The methylation locus in enhancer regions shower more discriminative power #
#########################################################################################

############
##Figure 2A#
############
options(stringsAsFactors=F)
options(scipen=200)
library(ggplot2)
library(RColorBrewer)

#matrix of beta value
matrix_dir = "m100_mt.RData"
load(matrix_dir)

#group label
label <- rep(0,ncol(m100_mt))
label[substr(colnames(m100_mt),1,nchar(colnames(m100_mt))-3)=="tumor_relapse"] <- 1

#The differentially methylated loci
deltaBeta <- apply(m100_mt,1,function(x) mean(as.numeric(x)[label==1])-mean(as.numeric(x)[label==0]))
pvalue <- apply(m100_mt,1,function(x) {
	t1 <- as.numeric(x)[label==1]
	t2 <- as.numeric(x)[label==0]
	t3 <- wilcox.test(t1,t2,alternative="greater")$"p.value"
	t4 <- wilcox.test(t1,t2,alternative="less")$"p.value"
	return(min(c(t3,t4)))
})
pvalue[is.na(pvalue)] <- 1
detaBeta_our <- data.frame(pos=rownames(m100_mt),deltaBeta=deltaBeta,pvalue=pvalue,stringsAsFactors=F)

#plot Figure 2A
fc_p <- detaBeta_our[,"deltaBeta"]
q_e <- -log10(detaBeta_our[,"pvalue"])
colo <- rep("grey",length(fc_p))
colo[(abs(fc_p)>0.1)&(q_e>2)] <- "red"
plot(fc_p,q_e,xlab="delta Beta",ylab="-log10(p value)",col=colo,pch=19)
abline(v=0.1,lty=2)
abline(v=-0.1,lty=2)
abline(h=2,lty=2)


##############
##Figure 2B-C#
##############
#The enhancer annotation of enhancerAtlas was downloaded from http://www.enhanceratlas.org/
#The enhancer annotation of dbSUPER was downloaded from http://asntech.org/dbsuper/
#The enhancer annotation of SCREEN was downloaded from http://screen.encodeproject.org

#The annotation of promoter
options(stringsAsFactors=F)
options(scipen=200)
library(GenomicRanges)
gtfl <- read.table(file="/data1/ref/hg38/gencode.v28.annotation.gtf",sep="\t",quote = "",stringsAsFactors=F,header=F)
genegtf <- gtfl[gtfl[,3]=="gene",]
genegtf_v9 <- strsplit(genegtf[,9],"\"| |;")
pos_g <- sapply(genegtf_v9,function(x) which(x=="gene_id"))
gene_id <- sapply(seq(length(genegtf_v9)),function(x) genegtf_v9[[x]][pos_g[x]+2])
pos_t <- sapply(genegtf_v9,function(x) which(x=="gene_type"))
gene_type <- sapply(seq(length(genegtf_v9)),function(x) genegtf_v9[[x]][pos_t[x]+2])
gtf <- data.frame(gene=gene_id,gene_type=gene_type,
chr=genegtf[,1],start=genegtf[,4],end=genegtf[,5],strand=genegtf[,7],stringsAsFactors=F)
width <- gtf[,"end"]-gtf[,"start"]+1
gtf <- gtf[width>200,]
#consider strand information
gtf_p <- gtf[gtf[,"strand"]=="+",]
promoter_p <- data.frame(gene=gtf_p[,"gene"],gene_type=gtf_p[,"gene_type"],
chr=gtf_p[,"chr"],start=gtf_p[,"start"]-2000,end=gtf_p[,"start"]+2000,strand=gtf_p[,"strand"],stringsAsFactors=F)
gtf_n <- gtf[gtf[,"strand"]=="-",]
promoter_n <- data.frame(gene=gtf_n[,"gene"],gene_type=gtf_n[,"gene_type"],
chr=gtf_n[,"chr"],start=gtf_n[,"end"]-2000,end=gtf_n[,"end"]+2000,strand=gtf_n[,"strand"],stringsAsFactors=F)
promoter <- rbind(promoter_p,promoter_n,stringsAsFactors=F)
promoter[,"gene"] <- substr(promoter[,"gene"],1,15)
screen_all <- read.table(file="/data1/dengyl/project/S1/chennan/new/screen/GRCh38-ccREs.PLS.bed",sep="\t",quote = "",stringsAsFactors=F,header=F)
Gscreen_all <- GRanges(seqnames = Rle(screen_all[,1]),
ranges = IRanges(as.numeric(screen_all[,2]),end = as.numeric(screen_all[,3]))) 
Gpromoter <- GRanges(seqnames = Rle(promoter[,"chr"]),
ranges = IRanges(as.numeric(promoter[,"start"]),end = as.numeric(promoter[,"end"]))) 
mtch1 <- findOverlaps(Gscreen_all,Gpromoter)
pair1 <- as.matrix(mtch1)
screen_all_promoter <- cbind(screen_all[pair1[,1],],promoter[pair1[,2],],stringsAsFactors=F)
#matrix of beta value
matrix_dir = "m100_mt.RData"
load(matrix_dir)
m100_mt <- m100_mt[detaBeta_our[(detaBeta_our[,"pvalue"]<0.01)&(abs(detaBeta_our[,"deltaBeta"])>0.1),1],]
region <- rownames(m100_mt)
region_bed <- matrix(unlist(strsplit(region,".",fixed=T)),byrow=T,ncol=3)
Gregion <- GRanges(seqnames = Rle(region_bed[,1]),
ranges = IRanges(as.numeric(region_bed[,2])-1,end = as.numeric(region_bed[,3]))) 
Gscreen_all_promoter <- GRanges(seqnames = Rle(screen_all_promoter[,1]),
ranges = IRanges(as.numeric(screen_all_promoter[,2]),end = as.numeric(screen_all_promoter[,3]))) 
mtch1 <- findOverlaps(Gregion,Gscreen_all_promoter)
pair1 <- as.matrix(mtch1)
meth_screen_all_promoter <- cbind(region[pair1[,1]],screen_all_promoter[pair1[,2],],stringsAsFactors=F)
colnames(meth_screen_all_promoter)[1] <- "region"
meth_promoter <- unique(meth_screen_all_promoter[,1])

#enhancer regions from dbSUPER
options(stringsAsFactors=F)
options(scipen=200)
setwd("/data1/dengyl/project/S1/chennan/new/dbSUPER/all_hg19_bed/")
fl <- dir()
fl <- grep("bed",fl,value=T)
celltype <- sub(".bed","",fl)
dbSUPER_list <- lapply(celltype,function(x) {
	t1 <- read.table(file=paste0("/data1/dengyl/project/S1/chennan/new/dbSUPER/all_hg19_bed/",x,".bed"),sep="\t",quote = "",stringsAsFactors=F,header=F)
	res <- data.frame(chr=t1[,1],start=t1[,2],end=t1[,3],
	cellline=gsub(" ","_",x),source="dbSUPER",stringsAsFactors=F)
	return(res)
})
dbSUPER <- do.call(rbind,dbSUPER_list)
setwd("/data1/dengyl/project/S1/chennan/new3")
write.table(dbSUPER, file="dbSUPER_enhancer_hg19.bed", sep="\t" , quote=FALSE, row.names=F, col.names=F )
#exit R, and run liftOver
/data1/luohao/software/vep86/liftOver/liftOver dbSUPER_enhancer_hg19.bed \
/data1/ref/hg19/hg19ToHg38.over.chain \
dbSUPER_enhancer_hg38.bed unmap.bed

#lung-specific or stem cell-specific enhancer regions from dbSUPER 
#For dbSUPER,convert hg19 to hg 38
options(stringsAsFactors=F)
options(scipen=200)
celltype <- c("H1","IMR90","Lung","NHLF","H2171","GLC16","NCI-H82","NCI-H69")
dbSUPER_list <- lapply(celltype,function(x) {
	t1 <- read.table(file=paste0("/data1/dengyl/project/S1/chennan/new/dbSUPER/all_hg19_bed/",x,".bed"),sep="\t",quote = "",stringsAsFactors=F,header=F)
	res <- data.frame(chr=t1[,1],start=t1[,2],end=t1[,3],cellline=x,source="dbSUPER",stringsAsFactors=F)
	return(res)
})
dbSUPER_lung <- do.call(rbind,dbSUPER_list)
write.table(dbSUPER_lung, file="dbSUPER_lung_hg19.bed", sep="\t" , quote=FALSE, row.names=F, col.names=F )
#exit R, and run liftOver
/data1/luohao/software/vep86/liftOver/liftOver dbSUPER_lung_hg19.bed \
/data1/ref/hg19/hg19ToHg38.over.chain \
dbSUPER_lung_hg38.bed unmap.bed

#annotation of lung-specific or stem cell-specific enhancer regions from dbSUPER 
options(stringsAsFactors=F)
options(scipen=200)
library(GenomicRanges)
screen_all_enhancerd_dbSUPER <- read.table(file="/data1/dengyl/project/S1/chennan/new/dbSUPER_hg38.bed",sep="\t",quote = "",stringsAsFactors=F,header=F) 
load("m100_mt.RData")
region <- rownames(m100_mt)
region_bed <- matrix(unlist(strsplit(region,".",fixed=T)),byrow=T,ncol=3)
Gregion <- GRanges(seqnames = Rle(region_bed[,1]),
ranges = IRanges(as.numeric(region_bed[,2])-1,end = as.numeric(region_bed[,3]))) 
Gscreen_all_enhancerd <- GRanges(seqnames = Rle(screen_all_enhancerd_dbSUPER[,1]),
ranges = IRanges(as.numeric(screen_all_enhancerd_dbSUPER[,2]),
end = as.numeric(screen_all_enhancerd_dbSUPER[,3]))) 
mtch1 <- findOverlaps(Gregion,Gscreen_all_enhancerd)
pair1 <- as.matrix(mtch1)
meth_screen_all_enhancerd_dbSUPER <- cbind(region[pair1[,1]],screen_all_enhancerd_dbSUPER[pair1[,2],],stringsAsFactors=F)
colnames(meth_screen_all_enhancerd_dbSUPER)[1] <- "region"
dbSUPER <- unique(meth_screen_all_enhancerd_dbSUPER[,1])

#enhancer regions from enhancerAtlas
options(stringsAsFactors=F)
options(scipen=200)
load("/data1/dengyl/project/S1/chennan/new/enhancerAtlas_allspecies.RData")
celltype <- colnames(HS)[-1]
pos <- apply(HS[,celltype],1,sum)>0
bed_info <- HS[pos,1]
bed_mt <- matrix(unlist(strsplit(as.character(bed_info),":|-")),byrow=T,ncol=3)
tissue <- apply(HS[pos,celltype],1,function(x) paste0(celltype[x!=0],collapse=","))
bed <- data.frame(chr=bed_mt[,1],start=as.numeric(bed_mt[,2]),end=as.numeric(bed_mt[,3]),
tissue=tissue,source="enhancerAtlas")
write.table(bed, file="enhancerAtlas_enhancer_hg19.bed", sep="\t" , quote=FALSE, row.names=F, col.names=F )
#exit R, and run liftover
/data1/luohao/software/vep86/liftOver/liftOver enhancerAtlas_enhancer_hg19.bed \
/data1/ref/hg19/hg19ToHg38.over.chain \
enhancerAtlas_enhancer_hg38.bed unmap.bed

#lung-specific or stem cell-specific enhancer regions from enhancerAtlas
#For enhancerAtlas, convert hg19 to hg 38
options(stringsAsFactors=F)
options(scipen=200)
load("/data1/dengyl/project/S1/chennan/new/enhancerAtlas_allspecies.RData")
celltype <- c("CiPSC","IMR90","ESC","iPSC_FX52","NHLF","NHBE","HSC","iPSC",
"Fetal_lung","NCC_H9","FiPSC","H2171","Bronchia_epithelial","H9",
"HUES64","H1299","H128","PC.9","H1703","NCI.H524","EiPSC","Lung","A549","H1")
pos <- apply(HS[,celltype],1,sum)>0
bed_info <- HS[pos,1]
bed_mt <- matrix(unlist(strsplit(as.character(bed_info),":|-")),byrow=T,ncol=3)
tissue <- apply(HS[pos,celltype],1,function(x) paste0(celltype[x!=0],collapse=","))
bed <- data.frame(chr=bed_mt[,1],start=as.numeric(bed_mt[,2]),end=as.numeric(bed_mt[,3]),
tissue=tissue,source="enhancerAtlas")
write.table(bed, file="enhancerAtlas_lung_hg19.bed", sep="\t" , quote=FALSE, row.names=F, col.names=F )
#exit R, and run liftOver
/data1/luohao/software/vep86/liftOver/liftOver enhancerAtlas_lung_hg19.bed \
/data1/ref/hg19/hg19ToHg38.over.chain \
enhancerAtlas_lung_hg38.bed unmap.bed

#the annotation of lung-specific or stem cell-specific enhancer regions from enhancerAtlas
options(stringsAsFactors=F)
options(scipen=200)
library(GenomicRanges)
load("m100_mt.RData")
screen_all <- read.table(file="/data1/dengyl/project/S1/chennan/new/enhancerAtlas_hg38.bed",sep="\t",quote = "",stringsAsFactors=F,header=F) 
screen_all_enhancerd_enhancerAtlas <- screen_all
region <- rownames(m100_mt)
region_bed <- matrix(unlist(strsplit(region,".",fixed=T)),byrow=T,ncol=3)
Gregion <- GRanges(seqnames = Rle(region_bed[,1]),
ranges = IRanges(as.numeric(region_bed[,2])-1,end = as.numeric(region_bed[,3]))) 
Gscreen_all_enhancerd <- GRanges(seqnames = Rle(screen_all_enhancerd_enhancerAtlas[,1]),
ranges = IRanges(as.numeric(screen_all_enhancerd_enhancerAtlas[,2]),end = as.numeric(screen_all_enhancerd_enhancerAtlas[,3]))) 
mtch1 <- findOverlaps(Gregion,Gscreen_all_enhancerd)
pair1 <- as.matrix(mtch1)
meth_screen_all_enhancerd_enhancerAtlas <- cbind(region[pair1[,1]],screen_all_enhancerd_enhancerAtlas[pair1[,2],],stringsAsFactors=F)
colnames(meth_screen_all_enhancerd_enhancerAtlas)[1] <- "region"
enhancerAtlas <- unique(meth_screen_all_enhancerd_enhancerAtlas[,1])

#enhancer regions from screen
options(stringsAsFactors=F)
options(scipen=200)
setwd("/data1/dengyl/project/S1/chennan/new/screenAll")
fl <- dir()
flname <- grep("bed",fl,value=T)
celltype <- sub(".bed","",fl)
celltype <- setdiff(celltype,c("download.log","download.txt"))
screen_list <- lapply(celltype,function(x) {
	t1 <- read.table(file=paste0("/data1/dengyl/project/S1/chennan/new/screenAll/",x,".bed"),sep="\t",quote = "",stringsAsFactors=F,header=F)
	t1 <- t1[t1[,10]%in%c("pELS","dELS","High-H3K27ac","High-H3K27ac,High-CTCF",
	"dELS,CTCF-bound","pELS,CTCF-bound"),]
	print(x)
	if(nrow(t1)>0)
	{
		res <- data.frame(chr=t1[,1],start=t1[,2],end=t1[,3],cellline=x,source="screen",stringsAsFactors=F)
		return(res)
	}
})
screen_list <- screen_list[sapply(screen_list,length)>0]
setwd("/data1/dengyl/project/S1/chennan/new3")
screen_enhancer <- do.call(rbind,screen_list)
save(screen_enhancer,file="screen_enhancer.RData")

#lung-specific or stem cell-specific enhancer regions from screen
options(stringsAsFactors=F)
options(scipen=200)
library(GenomicRanges)
screen_all1 <- read.table(file="/data1/dengyl/project/S1/chennan/new/screen/GRCh38-ccREs.pELS.bed",sep="\t",quote = "",stringsAsFactors=F,header=F)
screen_all2 <- read.table(file="/data1/dengyl/project/S1/chennan/new/screen/GRCh38-ccREs.dELS.bed",sep="\t",quote = "",stringsAsFactors=F,header=F)
screen_all_enhancerp <- rbind(screen_all1,screen_all2,stringsAsFactors=F)
load("m100_mt.RData")
region <- rownames(m100_mt)
region_bed <- matrix(unlist(strsplit(region,".",fixed=T)),byrow=T,ncol=3)
Gregion <- GRanges(seqnames = Rle(region_bed[,1]),
ranges = IRanges(as.numeric(region_bed[,2])-1,end = as.numeric(region_bed[,3]))) 
Gscreen_all_enhancerp <- GRanges(seqnames = Rle(screen_all_enhancerp[,1]),
ranges = IRanges(as.numeric(screen_all_enhancerp[,2]),end = as.numeric(screen_all_enhancerp[,3]))) 
mtch1 <- findOverlaps(Gregion,Gscreen_all_enhancerp)
pair1 <- as.matrix(mtch1)
meth_screen_all_enhancerp <- cbind(region[pair1[,1]],screen_all_enhancerp[pair1[,2],],stringsAsFactors=F)
colnames(meth_screen_all_enhancerp)[1] <- "region"
meth_screen <- unique(meth_screen_all_enhancerp[,1])

#plot for figure 2B-C
options(stringsAsFactors=F)
options(scipen=200)
library(GenomicRanges)
library(Vennerable)
load("screen_enhancer.RData")
load("m100_mt.RData")
load("detaBeta_our.RData")
load("enhancer_lung.RData")
load("meth_promoter.RData")

dbSUPER_lung <- read.table(file="dbSUPER_enhancer_hg38.bed",sep="\t",quote = "",stringsAsFactors=F,header=F) 
colnames(screen_enhancer) <- colnames(dbSUPER_lung)
enhancerAtlas_lung <- read.table(file="enhancerAtlas_enhancer_hg38.bed",sep="\t",quote = "",stringsAsFactors=F,header=F) 
screen_all_enhancerd <- rbind(screen_enhancer,rbind(dbSUPER_lung,enhancerAtlas_lung))
region <- rownames(m100_mt)
region_bed <- matrix(unlist(strsplit(region,".",fixed=T)),byrow=T,ncol=3)
Gregion <- GRanges(seqnames = Rle(region_bed[,1]),
ranges = IRanges(as.numeric(region_bed[,2])-1,end = as.numeric(region_bed[,3]))) 
Gscreen_all_enhancerd <- GRanges(seqnames = Rle(screen_all_enhancerd[,1]),
ranges = IRanges(as.numeric(screen_all_enhancerd[,2]),end = as.numeric(screen_all_enhancerd[,3]))) 
mtch1 <- findOverlaps(Gregion,Gscreen_all_enhancerd)
pair1 <- as.matrix(mtch1)
meth_screen_all_enhancerd <- cbind(region[pair1[,1]],screen_all_enhancerd[pair1[,2],],stringsAsFactors=F)
colnames(meth_screen_all_enhancerd)[1] <- "region"
meth_screen_all_enhancerd_exp <- meth_screen_all_enhancerd
colnames(meth_screen_all_enhancerd_exp) <- c("region","chromosome","start","end","sample_name","resource")
enhancer_all <- unique(meth_screen_all_enhancerd_exp[,1])
meth_diff <- detaBeta_our[(detaBeta_our[,"pvalue"]<0.01)&(abs(detaBeta_our[,"deltaBeta"])>0.1),1]
enhancer_all_diff = intersect(meth_diff,enhancer_all)
anno <- rep("others",length(meth_diff))
anno[meth_diff%in%c(enhancer_all)] <- "enhancer_in_other_tissue"
anno[meth_diff%in%c(enhancer_lung)] <- "enhancer_in_lung_or_stem_cell"
anno[meth_diff%in%meth_promoter] <- "promoter"
anno[meth_diff%in%intersect(meth_promoter,enhancer_all)] <- "enhancer&promoter"
#Figure 2B
pie(table(anno))
#figure 2C
t1 <- rep(0,nrow(detaBeta_our))
t1[detaBeta_our[,1]%in%intersect(enhancer_all,meth_diff)] <- 1
t2 <- rep(0,nrow(detaBeta_our))
t2[detaBeta_our[,1]%in%enhancer_lung] <- 1
x<- which(t1==1)    
y<- which(t2==1)  
dataD <-Venn(list("diff_enhancer"=x,"lung_or_stem_cell_enhancer"=y))    
plot(dataD,doWeight=T)

############
##Figure 2D#
############
#calculation of AUC
options(scipen=200)
library(pROC)
library(ggplot2)
library(RColorBrewer)
auc_our <- c()
for(i in seq(nrow(m100_mt))){
	roc1 <- pROC::roc(label, m100_mt[i,],direction= "<" ) 
	auc_our <- c(auc_our,roc1$"auc")
}
names(auc_our) <- rownames(m100_mt)

#plot figure 2D
auc_p <- auc_our[meth_promoter]
auc_e <- auc_our[enhancer_lung]
auc_ea <- auc_our[setdiff(enhancer_all,enhancer_lung)]
auc_p[auc_p<0.5] <- 1-auc_p[auc_p<0.5]
auc_e[auc_e<0.5] <- 1-auc_e[auc_e<0.5]
auc_ea[auc_ea<0.5] <- 1-auc_ea[auc_ea<0.5]
auc <- c(auc_p,auc_e,auc_ea)
label <- c(rep("promoter",length(auc_p)),rep("lung_or_stem_cell_enhancer",length(auc_e)),
rep("other_enhancer",length(auc_ea)))
MHI_df <- data.frame(MHI=auc,label=factor(label))
colo <- brewer.pal(3, "Set1")
dp <- ggplot(MHI_df, aes(x=label, y=MHI, fill=label)) + 
  geom_violin()+scale_fill_manual(values=colo)+geom_boxplot(width=0.1,fill="white")+
  labs(title="auc_promoter_vs_enhancer",y = "auc") + theme_classic()
dp
