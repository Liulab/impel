\name{integrative_model}
\alias{integrative_model}
\title{
Calculating integrative model
}
\description{
A model integrating  DNA methylation of enhancer regionand RNA expression of target genes, to predict the recurrence of early lung cancer.
}
\usage{
integrative_model(methylation.matrix, methylation_model, expression.matrix, expression_model, house.keeping.genes)
}
\arguments{
  \item{methylation.matrix}{
A DNA methylation matrix contains locus (in rows) and sample (in columns). 
The annotation of methylation locus was based on the hg38 reference genome.
}
  \item{methylation_model}{
The coeffecients of feature methylation locus.
}
  \item{expression.matrix}{
A gene expression matrix contains each gene (in rows) and sample (in columns). 
The expression matrix should be log-normalized.
}
  \item{expression_model}{
The coeffecients of feature genes.
}
  \item{house.keeping.genes}{
The house keeping genes used to normalized gene expression.
}
}
\examples{
data(methylation_matrix)
data(methylation_model)
data(expression.matrix)
data(expression_model)
data(house.keeping.genes)
res <-  integrative_model(methylation.matrix=methylation_matrix,methylation_model=methylation_model,
expression.matrix=expression.matrix,expression_model=expression_model,
house.keeping.genes=house.keeping.genes)
}
