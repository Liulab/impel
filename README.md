# impel
Integrative Methylation model to Predict recurrence of Early Lung cancer. A model constructed by DNA methylation of enhancer region, a model constructed by RNA expression of target genes, and a model integrating above two models. These three models could be used to predict the recurrence of early lung cancer.




#
### Data Availability

The data analyzed in this study are available from the National Genomics Data Center (accession numbers: HRA000208).




#
### install
```bash
library(devtools)
install_github("Liulab/impel")
```



#
### Quick start

#### meth_model 

The DNA methylationn model constructed by enhancer region to predict the recurrence of early lung cancer
```bash
data(methylation_matrix)
data(methylation_model)
res = meth_model(methylation.matrix=methylation_matrix,methylation_model=methylation_model)
```



#### exp_model

The expression model constructed by enhacer target to predict the recurrence of early lung cancer
```bash
data(expression.matrix)
data(expression_model)
data(house.keeping.genes)
res =  exp_model(expression.matrix=expression.matrix,expression_model=expression_model,house.keeping.genes=house.keeping.genes)
```



#### integrative_model

A model integrating  DNA methylation of enhancer regionand RNA expression of target genes, to predict the recurrence of early lung cancer.

```bash
data(methylation_matrix)
data(methylation_model)
data(expression.matrix)
data(expression_model)
data(house.keeping.genes)
res =  integrative_model(methylation.matrix=methylation_matrix,methylation_model=methylation_model,
expression.matrix=expression.matrix,expression_model=expression_model,
house.keeping.genes=house.keeping.genes)
```

