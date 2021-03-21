\name{exp_model}
\alias{exp_model}
\title{
Calculation of expression model constructed by enhacer target
}
\description{
The expression model constructed by enhacer target to predict the recurrence of early lung cancer
}
\usage{
exp_model(expression.matrix, expression_model, house.keeping.genes)
}
\arguments{
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
data(expression.matrix)
data(expression_model)
data(house.keeping.genes)
res <-  exp_model(expression.matrix=expression.matrix,expression_model=expression_model,house.keeping.genes=house.keeping.genes)
}
