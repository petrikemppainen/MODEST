% Generated by roxygen2 (4.0.2): do not edit by hand
\name{plotKDAG}
\alias{plotKDAG}
\title{Plots K-directed acyclic graph (KDAG)}
\usage{
plotKDAG(modest, Qdist, KDAG, min.vertex.size = 4, max.vertex.size = 20,
  intensities = c(1, 0.6, 0.3), plot.PCA = NULL, hue.saturation = c(1, 2))
}
\arguments{
\item{modest}{A file produced by \code{importData} that contains parsed data from replicate STRUCTURE runs}

\item{Qdist}{A file produced by \code{getQdist} that contains all pairwise distances between the Q-columns in file \code{modest}}

\item{KDAG}{A file produced by \code{Qdist2KDAG} that contains an edge list describing a KDAG.}

\item{min.vertex.size}{To come}

\item{max.vertex.size}{To come}

\item{intensities}{To come}

\item{plot.PCA}{To come}

\item{hue.saturation}{To come}
}
\value{
a list with the following objects. \code{vertices}: a list with each vertex in the DAG as an object containing index of each Q-column that have been merged into that vertex. \code{K_vertices}: the K for each vertex in file \code{vertices}. \code{color}: a vector of color codes that were used for each vertex in the KDAG plot.
}
\description{
Plots KDAG based on information produced by \code{\link{importData}}, \code{\link{getQdist}} and \code{\link{Qdist2KDAG}} and exports vertex information for function \code{\link{plotHist}}.
}
\details{
To come
}
\examples{
data(MODEST)
KDAG <- Qdist2KDAG(Qdist_MCMV, modest_MCMV$K, merge.by="topology")
plotKDAG(modest_MCMV, Qdist_MCMV, KDAG) # default, more options available
}
\author{
Petri Kemppainen \email{petrikemppainen2@gmail.com}
}
\seealso{
\code{\link{importData}}, \code{\link{getQdist}} and \code{\link{Qdist2KDAG}}
}
\keyword{STRUCTURE,}
\keyword{getQdist}
\keyword{importData,}
\keyword{plotDAG,}
\keyword{plotKDAG,}

