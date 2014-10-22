\name{ImportQcols}
\alias{ImportQcols}
\title{Imports Qcols from STRUCTURE output files}
\usage{
ImportQcols(input.data.folder, compress = 0.99)
}
\arguments{
  \item{input.data.folder}{Path to input data folder}

  \item{compress}{By default any Q-colums with a higher
  pearson-rank correlation of 0.98 are pooled}
}
\value{
A list with a matrix of all concatenated Q-columns
(accessed thrugh $Qcols) and the a vector that gives K's
for each Q-column (accessed through $K)
}
\description{
From a path to a data folder containing STRUCTURE output
data, Q-columns are concatenated into a single matrix that
can be used by \code{\link{Qcols2KDAG}}
}
\details{
More info to be added
}

