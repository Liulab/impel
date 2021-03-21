# impel
Integrative Methylation model to Predict recurrence of Early Lung cancer. A model constructed by DNA methylation of enhancer region, a model constructed by RNA expression of target genes, and a model integrating above two models. These three models could be used to predict the recurrence of early lung cancer.

# install
library(devtools)

install_github("Liulab/impel")

# functions
meth_model 

The DNA methylationn model constructed by enhancer region to predict the recurrence of early lung cancer


exp_model

The expression model constructed by enhacer target to predict the recurrence of early lung cancer


integrative_model

A model integrating  DNA methylation of enhancer regionand RNA expression of target genes, to predict the recurrence of early lung cancer.
