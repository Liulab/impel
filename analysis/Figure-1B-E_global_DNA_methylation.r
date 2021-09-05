###############################################################################################
###Figure 1B The distribution of methylation beta values across the genome                    #
###Figure 1C The methylation heterogeneity index(MHI) value                                   #
###Figure 1D The global methylation                                                           #
###Figure 1E The average methylation of Solo-WCGW CpGs in  partial methylation domains (PMDs) #
###############################################################################################

##############
##Figure 1B-C#
##############
library(ggplot2)
library(RColorBrewer)
options(stringsAsFactors=F)
options(scipen=200)

#matrix of beta value
matrix_dir = "/data1/dengyl/project/S1/chennan/new/data/our_mt.imputed.RData"
load(matrix_dir)

#group label
label = substr(colnames(our_mt.imputed),1,nchar(colnames(our_mt.imputed))-3)
ex1 <- as.numeric(our_mt.imputed[,label==sort(unique(label),decreasing=T)[1]])
ex2 <- as.numeric(our_mt.imputed[,label==sort(unique(label),decreasing=T)[2]])
ex3 <- as.numeric(our_mt.imputed[,label==sort(unique(label),decreasing=T)[3]])
ex4 <- as.numeric(our_mt.imputed[,label==sort(unique(label),decreasing=T)[4]])

#plot Figure 1B
plot(density(ex1,adjust=2),col=colo[1],xlab="Methylation beta-value",main="",xlim=c(0.2,0.8),ylim=c(0.15,0.5))
lines(density(ex2,adjust=2),col=colo[2],xlim=c(0.2,0.8),ylim=c(0.15,0.5))
lines(density(ex3,adjust=2),col=colo[3],xlim=c(0.2,0.8),ylim=c(0.15,0.5))
lines(density(ex4,adjust=2),col=colo[4],xlim=c(0.2,0.8),ylim=c(0.15,0.5))

#threshold for Figure 1C
thUp=0.85
thDown=0.25

#MHI value
MHI <- apply(our_mt.imputed,2,function(x) sum((x<thUp)&(x>thDown))/length(x))

#plot Figure 1C
MHI_list <- split(MHI,f=factor(label))
MHI_df <- data.frame(MHI=MHI,label=factor(label,levels=sort(unique(label),decreasing=T)))
colo <- c(brewer.pal(3, "Set1"),brewer.pal(7, "Greens")[3])
dp <- ggplot(MHI_df, aes(x=label, y=MHI, fill=label)) + 
  geom_violin()+scale_fill_manual(values=colo)+geom_boxplot(width=0.1,fill="white")+
  labs(title="MHI_reccurence",y = "Methylation Heterogeneity Index") + theme_classic()
dp


############
##Figure 1D#
############
options(stringsAsFactors=F)
options(scipen=200)
library(GenomicRanges)
library(ggplot2)
library(RColorBrewer)

#matrix of beta value
matrix_dir = "/data1/dengyl/project/S1/chennan/new/data/our_mt.imputed.RData"
load(matrix_dir)

#group label
label = substr(colnames(our_mt.imputed),1,nchar(colnames(our_mt.imputed))-3)

#The LINE-1 family annotation was downloaded from the RepeatMasker track of the UCSC genome browser. 
rmsk <- read.table(file="/data1/ref/hg38/UCSC/rmsk.txt",
sep="\t",quote = "",stringsAsFactors=F,header=F)
rmsk_line <- rmsk[rmsk[,12]=="LINE",]

#overlap between methylation locus and LINE-1 family annotation
our_bed <- matrix(unlist(strsplit(rownames(our_mt.imputed),".",fixed=T)),byrow=T,ncol=2)
Gour <- GRanges(seqnames = Rle(our_bed[,1]),
ranges = IRanges(as.numeric(our_bed[,2]),end = as.numeric(our_bed[,2])))
Grmsk <- GRanges(seqnames = Rle(rmsk_line[,6]),
ranges = IRanges(as.numeric(rmsk_line[,7]),end = as.numeric(rmsk_line[,8])))
mtch1 <- findOverlaps(Gour, Grmsk)
pair1 <- as.matrix(mtch1)

#global methylation
gmeth <- apply(our_mt.imputed[unique(pair1[,1]),],2,mean)
gmeth_list <- split(gmeth,f=factor(label))

#plot Figure 1D
gmeth_df <- data.frame(gmeth=gmeth,label=factor(label,levels=sort(unique(label),decreasing=T)))
colo <- c(brewer.pal(3, "Set1"),brewer.pal(7, "Greens")[3])
dp <- ggplot(gmeth_df, aes(x=label, y=gmeth, fill=label)) + 
  geom_violin()+scale_fill_manual(values=colo)+geom_boxplot(width=0.1,fill="white")+
  labs(title="global_methylation_reccurence",y = "") + theme_classic()
dp


############
##Figure 1E#
############
options(stringsAsFactors=F)
options(scipen=200)
library(GenomicRanges)

#matrix of beta value
matrix_dir = "/data1/dengyl/project/S1/chennan/new/data/our_mt.imputed.RData"
load(matrix_dir)

#group label
label = substr(colnames(our_mt.imputed),1,nchar(colnames(our_mt.imputed))-3)

#The common PMDs were dowmload from https://zwdzwd.github.io/pmd.
rmsk <- read.table(file="/data1/dengyl/project/S1/chennan/platform/sd.binmean.bed",
sep="\t",quote = "",stringsAsFactors=F,header=F)
rmsk_line <- rmsk[rmsk[,6]=="commonPMD",]

#overlap between methylation locus and common PMDs
our_bed <- matrix(unlist(strsplit(rownames(our_mt.imputed),".",fixed=T)),byrow=T,ncol=2)
Gour <- GRanges(seqnames = Rle(our_bed[,1]),
ranges = IRanges(as.numeric(our_bed[,2]),end = as.numeric(our_bed[,2])))
Grmsk <- GRanges(seqnames = Rle(rmsk_line[,1]),
ranges = IRanges(as.numeric(rmsk_line[,2]),end = as.numeric(rmsk_line[,3])))
mtch1 <- findOverlaps(Gour, Grmsk)
pair1 <- as.matrix(mtch1)
our_bed_in_pmd <- our_bed[unique(pair1[,1]),]

#obtain Solo-WCGW CpGs
our_bed_in_pmd_flank35 <- data.frame(chr=c(our_bed_in_pmd[,1],our_bed_in_pmd[,1]),
start=c(as.numeric(our_bed_in_pmd[,2])-36,as.numeric(our_bed_in_pmd[,2])-36),
end=c(as.numeric(our_bed_in_pmd[,2])+35,as.numeric(our_bed_in_pmd[,2])+35),
name=c(paste(our_bed_in_pmd[,1],our_bed_in_pmd[,2],"+",sep="."),
paste(our_bed_in_pmd[,1],our_bed_in_pmd[,2],"+",sep=".")),
value=1,strand=c(rep("+",nrow(our_bed_in_pmd)),rep("-",nrow(our_bed_in_pmd))),stringsAsFactors=F)
write.table(our_bed_in_pmd_flank35 , file="our_bed_in_pmd_flank35.bed", sep="\t" , quote=FALSE, row.names=F, col.names=F )

#exit R and run bedtools
bedtools getfasta -fi /data1/ref/hg38/GATK_bundle/Homo_sapiens_assembly38.fasta \
-bed our_bed_in_pmd_flank35.bed -s -name > our_bed_in_pmd_flank35.reduce.fasta

#R again
library(ggplot2)
library(RColorBrewer)
options(stringsAsFactors=F)
options(scipen=200)

#matrix of beta value
matrix_dir = "/data1/dengyl/project/S1/chennan/new/data/our_mt.imputed.RData"
load(matrix_dir)

#obtain Solo-WCGW CpGs
fl <- read.table(file="our_bed_in_pmd_flank35.reduce.fasta",sep="\n",quote = "",stringsAsFactors=F,header=F)[,1]
pos1 <- grep("^>",fl)
poss <- pos1+1
pose <- c(pos1[-1]-1,length(fl))
sequen <- lapply(seq(pos1),function(x) unlist(strsplit(fl[poss[x]:pose[x]],"")))
posc <- sapply(sequen,function(x) x[36]=="C")
possolo <- sapply(sequen,function(x){
	sum(sapply(1:70,function(y) paste0(x[y],x[y+1]))=="CG")==1
})
fasta_name <- grep("^>",fl,value=T)
solo_WGCW <- fasta_name[posc&possolo]
solo_WGCW <- sub(">","",solo_WGCW )
solo_WGCW  <- sapply(strsplit(solo_WGCW,":"),function(x) x[1])
solo_WGCW <- sapply(strsplit(solo_WGCW,".",fixed=T),function(x) paste(x[1],x[2],sep="."))

#PMD
out_mt.solo <- our_mt.imputed[solo_WGCW,]
pmd <- apply(out_mt.solo,2,mean)
pmd_list <- split(pmd,f=factor(label))

#plot Figure 1E
pmd_df <- data.frame(pmd=pmd,label=factor(label,levels=sort(unique(label),decreasing=T)))
colo <- c(brewer.pal(3, "Set1"),brewer.pal(7, "Greens")[3])
dp <- ggplot(pmd_df, aes(x=label, y=pmd, fill=label)) + 
  geom_violin()+scale_fill_manual(values=colo)+geom_boxplot(width=0.1,fill="white")+
  labs(title="partial methylation domains",y = "") + theme_classic()
dp
