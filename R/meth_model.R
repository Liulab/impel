meth_model <-
function(methylation.matrix,methylation_model)
{
if(sum(!c(rownames(methylation_model)[methylation_model[,1]!=0][-1])%in%rownames(methylation.matrix))>0)
{
"Features are missing!"
}else{
mt_tr=methylation.matrix[rownames(methylation_model)[methylation_model[,1]!=0][-1],]
lasso_C2 <- apply(mt_tr,2,function(x) {
tt1 <- 1/(exp(-1*(sum(x*(methylation_model[-1,1]))+methylation_model[1,1]))+1)
return(tt1)
})
risk.class <- rep("intermediate",ncol(methylation.matrix))
risk.class[lasso_C2>=0.6] <- "high risk"
risk.class[lasso_C2<=0.4] <- "low risk"
res_df <- data.frame(sample.name=colnames(methylation.matrix),
risk.score=lasso_C2,risk.class=risk.class,stringsAsFactors=F)
rownames(res_df) <- NULL
return(res_df)
}
}
