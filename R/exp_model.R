exp_model <-
function(expression.matrix,expression_model,house.keeping.genes)
{
if(sum(!c(house.keeping.genes,rownames(expression_model)[expression_model[,1]!=0][-1])%in%rownames(expression.matrix))>0)
{
"Feature genes are missing!"
}else{
mt_tr=t(apply(expression.matrix[rownames(expression_model)[expression_model[,1]!=0][-1],],1,function(x) as.numeric(x)-apply(expression.matrix[house.keeping.genes,],2,mean)))
lasso_C2 <- apply(mt_tr,2,function(x) {
tt1 <- 1/(exp(-1*(sum(x*(expression_model[-1,1]))+expression_model[1,1]))+1)
return(tt1)
})
risk.class <- rep("intermediate",ncol(expression.matrix))
risk.class[lasso_C2>=0.5] <- "high risk"
risk.class[lasso_C2<=0.5] <- "low risk"
res_df <- data.frame(sample.name=colnames(expression.matrix),
risk.score=lasso_C2,risk.class=risk.class,stringsAsFactors=F)
rownames(res_df) <- NULL
return(res_df)
}
}
