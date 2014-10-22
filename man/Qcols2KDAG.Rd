\name{Qcols2KDAG}
\alias{Qcols2KDAG}
\title{Parses STRUCTURE output to a K-directed acyclig graph (K-DAG)}
\usage{
Qcols2KDAG(Qcols, K)
}
\arguments{
  \item{Qcols}{A matrix of concatenated Qcols from
  replicate STRUCTURE runs across many K values}

  \item{K}{A vector of K's for each Qcolumn in
  \code{Qcols}}
}
\value{
returns a list with two objects, \code{KDAG} which contains
an edgelist containing the informatino of the K-directed
acyclig graph and \code{toMerge}, which contains a list of
nodes from \code{KDAG} and which Qcolumns they represent.
}
\description{
Takes a matrix of concatenated Qcolumns from mutltiple
STRUCTURE runs and transforms these into a K-DAG
}
\details{
More info to come
}
\author{
Petri Kemppainen, \email{petrikemppainen2@gmail.com} and
Stuart Baird, \email{stuartj.e.baird@gmail.com}
}
\keyword{Qcols2DAG}

