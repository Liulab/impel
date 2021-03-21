\name{meth_model}
\alias{meth_model}
\title{
Calculation of DNA methylation model constructed by enhancer region
}
\description{
The DNA methylationn model constructed by enhancer region to predict the recurrence of early lung cancer
}
\usage{
meth_model(methylation.matrix, methylation_model)
}
\arguments{
  \item{methylation.matrix}{
A DNA methylation matrix contains locus (in rows) and sample (in columns). 
The annotation of methylation locus was based on the hg38 reference genome.
}
  \item{methylation_model}{
The coeffecients of feature methylation locus.
}
}
\examples{
data(methylation_matrix)
data(methylation_model)
res <- meth_model(methylation.matrix=methylation_matrix,methylation_model=methylation_model)
}
