#' Parses STRUCTURE output to a K-directed acyclig graph (K-DAG)
#' 
#' Takes a matrix of concatenated Qcolumns from mutltiple STRUCTURE runs and transforms these into a K-DAG
#' 
#' More info to come
#' 
#' @param Qcols A matrix of concatenated Qcols from replicate STRUCTURE runs across many K values
#' @param K A vector of K's for each Qcolumn in \code{Qcols}
#' @return returns a list with two objects, \code{KDAG} which contains an edgelist containing the informatino of the K-directed acyclig graph and \code{toMerge}, which contains a list of nodes from \code{KDAG} and which Qcolumns they represent.
#' @export
#' @keywords Qcols2DAG
#' @author Petri Kemppainen, \email{petrikemppainen2@@gmail.com} and Stuart Baird, \email{stuartj.e.baird@@gmail.com}

Qcols2KDAG <- function(Qcols, K){
  mergeDAG(Qcols2KDAG_inner(Qcols, K), K)
} 

Qcols2KDAG_inner <- function(Qcols, K){ # the imputs for this function are a matrix of Qcols and a vector of K for each Qcolumn
  ### prepare slink
  colnames(Qcols) <- c(1:ncol(Qcols))
  correlations <- cor(Qcols, method="pearson") # calculate correlations
  slink <- as.phylo(hclust(as.dist(1-correlations), method="single")) # get slink from distance matrix and parse it to a 'phylo' object
  ### prepare some useful info 
  Ntips <- length(slink$tip.label) # number of tips
  Nclust <- slink$Nnode # number of clusters
  el <- slink$edge # shorter name for the edge list which contains the slink tree structure
  
  ### prepare info to be used by the parser
  Cl <- lapply((Ntips+1):max(el[,1]), function(x) extract.clade(slink, x)) # a list of all possible sub-cluster from slink using function 'extract.clade' from 'ape' package
  Cl <- lapply(Cl, function(x) as.numeric(x$tip.label)) # list with tip labels (Qindex) for each cluster 
  Cl <- c(1:Ntips, Cl) # now add all tips as well 
  Cl_K <- lapply(Cl, function(x) K[x]) # List with K:s for each cluster/tip
  Cl_minK <- sapply(Cl_K, min) # minimum K for each cluster/tip
  Cl_minK_Qcol <- lapply(1:length(Cl), function(x) sort(Cl[[x]][Cl_K[[x]] == Cl_minK[[x]]])) # Q cols that have minK for each cluster; these will be vertices in the DAG
  DAG <- matrix(NA, 0, 2) # file to append to
  Eindex <- 0
  index <- NULL # gives the names for each cluster for which the parser was used
  #x <- 360
  ### The parser
  parser <- function(x){ 
    ## prepare info
    children <- el[,2][el[,1]==x] # get children
    children <- children[order(Cl_minK[children])] # order children such that delta K is never negative
    # get Ks for children
    child1_K <- unlist(Cl_minK[children][1])
    child2_K <- unlist(Cl_minK[children][2])
    # Do not create edges when K for children is the same, this is the only exception (!)
    if(child1_K!=child2_K){
      index <<- c(index, x)
      # get grandchildren and vertex names
      grandChildren1 <- unlist(Cl_minK_Qcol[children[1]])
      grandChildren2 <- unlist(Cl_minK_Qcol[children[2]])
      child1_names <- grandChildren1
      child2_names <- grandChildren2
      #two nested functions that create all edges
      parseInner <- function(y){
        ## prepare some info
        # get Ks for children
        createEdges <- function(x) {
          createEdgesInner <- function(y){
            from <- child1_names[y]
            to <- child2_names[x]      
            if(diff(c(child1_K, child2_K))==1){ # if delta K ==1
              DAG <- rbind(DAG, c(from, to))
            }else{ # if delta K >1 create necessary empty vertices
              Eindex <<- Eindex-1
              DAG <- rbind(DAG, c(Eindex, to))
              if(diff(c(child1_K, child2_K))>2){
                for(i in 1:(diff(c(child1_K, child2_K))-2)){
                  DAG <- rbind(DAG, c(Eindex-1, Eindex))
                  Eindex <<- Eindex-1                
                }
              }
              DAG <- rbind(DAG, c(from, Eindex))
            }
          }
          DAG <- rbind(DAG, do.call('rbind', lapply(1:length(child1_names), createEdgesInner)))  
          return(DAG)
        }
        DAG <- do.call('rbind', lapply(1:length(child2_names), createEdges))
        return(DAG)
      }
      DAG <- parseInner()  
    }
  }
  ### now parse
  DAGlist <- lapply((Ntips+1):(Ntips+Nclust), parser) # use parser for all clusters
  names(DAGlist) <- index # give each object relevant names (equals to node numbers for which the parser was used), used for debugging
  DAG <- do.call('rbind', DAGlist) # append 
  DAG <- DAG[!duplicated(DAG),]
  ##create root
  names <- unlist(Cl[which(K==2)])
  DAG <- rbind(DAG, cbind(0, names)) # add root
  return(DAG)
}
##################
mergeDAG <- function(DAG, K){
  #start with root
  out <- list()
  out[[1]] <- list(0)
  names(out[[1]][[1]]) <- "root"
  currentK <- 1
  MDAG <- matrix(NA, 0, 2)
  while((currentK)<max(K)){
    parents <- out[[currentK]]
    currentK <- currentK+1
    temp <- list()
    x <- 0
    #i <- 11
    #collapse nodes with comparable topologies and build merged DAG
    for(i in 1:length(parents)){
      if(any(DAG[,1] %in% unlist(parents[[i]]))){
        x <- x+1
        from <- paste(sort(unlist(parents[[i]])), collapse=",")
        children <- unique(DAG[,2][DAG[,1] %in% unlist(parents[[i]])])
        subDAGs <- lapply(children, function(x) extractDAG(DAG, x))
        subDAGsTop <- lapply(subDAGs, function(x) collapseSinglesEl(x))
        trees <- subDAGsTop
        which.not.tip <- which(subDAGsTop != "tip")
        trees[which.not.tip] <- lapply(subDAGsTop[which.not.tip], DAG2Phylo)
        d.tree = matrix( nrow=length(trees), ncol=length(trees), dimnames=list(children, children) )
        #decision tree for what two trees have comparable topologies
        for(h in 1:length(trees)){
          for(k in h:length(trees)){
            if(any(c(any(trees[[k]]=="tip"), any(trees[[h]]=="tip")))){
              if(all(trees[[k]]=="tip") && any(trees[[h]]=="tip")){
                d.tree[k,h] <- 0 
              }else{
                d.tree[k,h] <- 1 
              }
            }else{
              if(trees[[k]]$Nnode+trees[[h]]$Nnode==2){
                if(length(trees[[k]]$tip.label)==length(trees[[h]]$tip.label)){
                  d.tree[k,h] <- 0  
                }else{
                  d.tree[k,h] <- 1
                }
              }
              if(trees[[k]]$Nnode+trees[[h]]$Nnode==3) d.tree[k,h] <- 1  
              if(trees[[k]]$Nnode+trees[[h]]$Nnode>3){
                if(any(c(trees[[k]]$Nnode,trees[[h]]$Nnode)==1)){
                  d.tree[k,h] <- 1
                }else{
                  d.tree[k,h] <- dist.topo(trees[[k]], trees[[h]])  
                }
              } 
            }
          }
        }
        #find clusters of the same topology, here I use networks
        d.tree[d.tree!=0] <- -1
        d.tree[d.tree==0] <- 1
        g <- graph.adjacency(d.tree, mode="lower", diag=FALSE, weighted=T)
        g <- delete.edges(g, which(E(g)$weight==-1))
        temp2 <- lapply(decompose.graph(g), function(x) V(x)$name)
        #
        for(j in 1:length(temp2)){
          if(any(unlist(temp2[j]) %in% unlist(temp))){
            shared <- unlist(temp2[j])[unlist(temp2[j]) %in% unlist(temp)]
            which.shared <- which(unlist(lapply(temp, function(x) any(x %in% shared))))
            new.name <- paste(sort(unique(c(unlist(temp[which.shared]), unlist(temp2[j])))), collapse=",")
            for(y in which.shared){
              MDAG[,2][MDAG[,2] == paste(unlist(temp[y]), collapse=",")] <- new.name
              temp[[y]] <- sort(unique(c(unlist(temp[which.shared]), unlist(temp2[j]))))
            }
            MDAG <- rbind(MDAG, c(from, new.name))  
            if(j==length(temp2)) x <- x-1
          }else{
            columns <- sort(unlist(temp2[j]))
            names(columns) <- NULL
            temp[[x]] <- columns
            to <- paste(columns, collapse=",")
            MDAG <- rbind(MDAG, c(from, to))
            if(j<length(temp2)) x <- x+1
          }
        }
        #find which columns should be merged and build merged DAG
      }
    }
    duplicated(temp)
    out[[currentK]] <- temp
  }
  out <- list(out, MDAG)
  names(out) <- c("toMerge", "KDAG")
  out
}
############## Functions
extractDAG <- function(DAG, node){
  el <- DAG
  node.first <- node
  vertices <- NULL
  while(!identical(el[,2][el[,1] %in% node], numeric(0))){
    vertices <- c(vertices, el[,2][el[,1] %in% node])
    node <- el[,2][el[,1] %in% node]
  } 
  vertices <- unique(c(node.first, vertices))
  subDAG <- el[el[,1] %in% vertices & el[,2] %in% vertices,]
  subDAG  
}
###
DAG2Phylo <- function(DAG){
  el <- DAG
  el.new <- el
  tips <- unique(el[,2][!el[,2] %in% el[,1]])
  nodes <- unique(as.vector(el)[!as.vector(el) %in% tips])
  nTips <- length(tips)
  nNodes <- length(nodes)
  el.new[,2] <- c((nTips+1):(nNodes+nTips))[match(el[,2], nodes)]
  el.new[,1] <- c((nTips+1):(nNodes+nTips))[match(el[,1], nodes)]
  el.new[,2][is.na(el.new[,2])] <- na.omit(c(1:nTips)[match(el[,2], tips)])
  tree <- list()
  tree$edge <- apply(el.new, 2, as.numeric)
  tree$tip.label <- tips
  tree$Nnode <- nNodes
  class(tree) <- "phylo"
  tree
}
#####
collapseSinglesEl <- function(el){
  
  if(is.vector(el)){return("tip")
  }else{
    temp <- table(el[,1])
    singles <- as.numeric(names(temp)[which(temp==1)])
    
    if(length(temp)==length(singles)){return("tip")
    }else{
      for(i in singles){
        el[,2][el[,2] == i] <- el[,2][el[,1] == i]
        el <- el[el[,1] != i,]
      }
      el <- el[!duplicated(el),]
      if(is.vector(el)){
        return("tip")
      }else{
        return(el)  
      } 
    }
  }
}
############
