#####################################################
###Figure S1 Comparison of MHI at various parameters#
#####################################################

options(stringsAsFactors=F)
options(scipen=200)
library(gplots)
library(RColorBrewer)

#matrix of beta value
matrix_dir = "/data1/dengyl/project/S1/chennan/new/data/our_mt.imputed.RData"
load(matrix_dir)

#group label
label = substr(colnames(our_mt.imputed),1,nchar(colnames(our_mt.imputed))-3)

#grid search
thUp_list <- seq(from=0.7,to=0.9,by=0.01)
thDown_list <- seq(from=0.1,to=0.3,by=0.01)
Combn <- data.frame(up=rep(thUp_list,length(thDown_list)),down=rep(thDown_list,rep(length(thDown_list),length(thDown_list))))

#Comparison of MHI at various parameters
p_34 <- c()
p_421 <- c()
p_321 <- c()
for(i in seq(nrow(Combn)))
{
thUp=Combn[i,1]
thDown=Combn[i,2]
MHI <- apply(our_mt.imputed,2,function(x) sum((x<thUp)&(x>thDown))/length(x))
MHI_list <- split(MHI,f=factor(label))
p_34 <- c(p_34,wilcox.test(MHI_list[[4]],MHI_list[[3]],alternative="greater")$"p.value")
p_421 <- c(p_421,wilcox.test(MHI_list[[4]],c(MHI_list[[1]],MHI_list[[2]]),alternative="greater")$"p.value")
p_321 <- c(p_321,wilcox.test(MHI_list[[3]],c(MHI_list[[1]],MHI_list[[2]]),alternative="greater")$"p.value")
}
names(p_34) <- paste(Combn[,1],Combn[,2],sep="_")
names(p_421) <- paste(Combn[,1],Combn[,2],sep="_")
names(p_321) <- paste(Combn[,1],Combn[,2],sep="_")

#plot Figure S1A
p_34mt <- matrix(p_34,21,21)
rownames(p_34mt) <- seq(from=0.7,to=0.9,by=0.01) 
colnames(p_34mt) <- seq(from=0.1,to=0.3,by=0.01) 
heatmap.2(p_34mt,Rowv = F,Colv=F,dendrogram = "none",scale = "none",trace="none",
           breaks=c(0.0001,0.001,0.01,0.05,0.2),
           col=c(brewer.pal(4, "Reds")[4:2],"grey"))

#plot Figure S1B		   
p_421mt <- matrix(p_421,21,21)
rownames(p_421mt) <- seq(from=0.7,to=0.9,by=0.01) 
colnames(p_421mt) <- seq(from=0.1,to=0.3,by=0.01) 
heatmap.2(p_421mt,Rowv = F,Colv=F,dendrogram = "none",scale = "none",trace="none",
           col=rev(colorRampPalette(brewer.pal(5, "Reds")[2:4])(20)))

#plot Figure S1C 
p_321mt <- matrix(p_321,21,21)
rownames(p_321mt) <- seq(from=0.7,to=0.9,by=0.01) 
colnames(p_321mt) <- seq(from=0.1,to=0.3,by=0.01) 
heatmap.2(p_321mt,Rowv = F,Colv=F,dendrogram = "none",scale = "none",trace="none",
           col=rev(colorRampPalette(brewer.pal(5, "Reds")[2:4])(20)))
