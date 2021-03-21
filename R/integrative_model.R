integrative_model <-
function(methylation.matrix,methylation_model,expression.matrix,expression_model,house.keeping.genes)
{
res_meth <- meth_model(methylation.matrix=methylation.matrix,methylation_model=methylation_model)
res_exp <- exp_model(expression.matrix=expression.matrix,expression_model=expression_model,house.keeping.genes=house.keeping.genes)
res_df <- data.frame(sample.name=res_meth[,1],
meth.score=res_meth[,2],exp_score=res_exp[,2],
integrative_score=0.5*(res_meth[,2]+res_meth[,2]),stringsAsFactors=F)
return(res_df)
}
