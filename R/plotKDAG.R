#' Plots KDAG
#' 
#' Plots output from \code{\link{Qcols2KDAG}} 
#' 
#' More info to come
#' 
#' @param KDAG 0utput \code{KDAG} from \code{\link{Qcols2KDAG}}
#' @param ... Additional graphical parameters to be passed on the the plot.igraph functions
#' @return NULL
#' @export
#' @keywords plotDAG
#' @author Petri Kemppainen, \email{petrikemppainen2@@gmail.com} and Stuart Baird, \email{stuartj.e.baird@@gmail.com}

plotKDAG <- function(KDAG, ...){
  DC <- graph.edgelist(apply(KDAG, 2, as.character))
  
  V(DC)$shape <- "circle"
  strsplit(V(DC)$name, ",")
  V(DC)$shape[sapply(strsplit(V(DC)$name, ","), function(x) all(x<0))] <- "square"
  length(V(DC))
  V(DC)$label.cex <- 0.75
  
  V(DC)$color <- "yellow"
  V(DC)$color[sapply(strsplit(V(DC)$name, ","), function(x) all(x<0))] <- "white"
  V(DC)$size <- 5
  V(DC)$size[sapply(strsplit(V(DC)$name, ","), function(x) all(x<0))] <- 1
  E(DC)$color <- "black"
  E(DC)$lty <- 1
  
  lay2 <- layout.sugiyama(DC, attributes="all", maxiter = 1000)
  plot(DC, layout=lay2$layout, vertex.label.color="black", edge.arrow.size=0.1, edge.label.family="Helvetica", edge.label.color="black", ...)
}
