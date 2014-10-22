#' Imports Qcols from STRUCTURE output files
#' 
#' From a path to a data folder containing STRUCTURE output data, Q-columns are concatenated into a single matrix that can be used by \code{\link{Qcols2KDAG}}
#' 
#' More info to be added
#' 
#' @param input.data.folder Path to input data folder
#' @param compress By default any Q-colums with a higher pearson-rank correlation of 0.98 are pooled
#' @return A list with a matrix of all concatenated Q-columns (accessed thrugh $Qcols) and the a vector that gives K's for each Q-column (accessed through $K)
#' @export

ImportQcols <- function(input.data.folder, compress=0.99){
  raw <- Q.extract(input.data.folder)
  fColnames <- lapply(raw[[1]], function(x) colnames(x))
  extractIndex <- lapply(fColnames, function(y) which(sapply(strsplit(y, "", fixed=T), function(x) x[1])=="C"))
  temp <- lapply(1:length(raw[[1]]), function(x) raw[[1]][[x]][,extractIndex[[x]]])
  Qcols <- do.call('cbind', temp)
  K <- rep(sapply(extractIndex, length), sapply(extractIndex, length))
  Qcols <- Qcols[,K!=1]
  K <- K[K!=1]
  
  if(compress<1){
    colnames(Qcols) <- c(1:ncol(Qcols))
    correlations <- cor(Qcols, method="pearson") # calculate correlations
    g <- graph.adjacency(correlations, mode="upper", diag=FALSE, weighted=TRUE)
    g <- delete.edges(g, which(E(g)$weight<compress))
    temp <- lapply(decompose.graph(g), function(x) V(x)$name)
    temp <- temp[sapply(temp, function(x) length(x)>1)]
    remove <- lapply(temp, function(x) sample(x, length(x)-1))
    Qcols <- Qcols[,-as.numeric(unlist(remove))]
    K <- K[-as.numeric(unlist(remove))]
  }
  out <- list(Qcols, K)
  names(out) <- c("Qcols", "K")
  out
}

Q.extract=function(input.data.folder,Kcrit=0.55)
  {
  
  # check format of the input.data.folder for further use, as it needs a / at the end.
  if(read.fwf(textConnection(input.data.folder),widths=c(nchar(input.data.folder)-1,nchar(input.data.folder)))[2] != "/")
  {
    input.data.folder=paste(input.data.folder,"/",sep="") 
  }
  # end check format
  
  input.files=paste(input.data.folder,list.files(input.data.folder),sep="")
  
  output.list=list()
  subset.list=list()
  list.teller=1
  info.data=data.frame(filenumber=numeric(0),columnnumber=numeric(0),subsetnumber=numeric(0),subsetsize=numeric(0),K=numeric(0))
  # i=1;ii=1
  for(i in 1:length(input.files))
  {
    
    inp=suppressWarnings(readLines(input.files[i]))
    
    N=as.numeric(gsub(" |individuals","",inp[grep("individuals",inp)[which(grep("individuals",inp)>grep("Run parameters",inp))][1]]))
    
    rstart=grep("(%Miss)",inp)+1
    rstop=grep("(%Miss)",inp)+N
    
    inp2=gsub(",|\\(|\\)"," ",inp[rstart:rstop])
    dat=read.table(text=inp2)
    #   head(dat)
    #   which(dat[1,]==":")
    #  if(grepl("Label",inp[rstart-1])) {idx.id=2;idx.first.cluster.column=6} else {idx.id=1;idx.first.cluster.column=5}    # check if there are individual labels, which changes the columns that need to be extracted
    
    if(grepl("Label",inp[rstart-1])) {idx.id=2} else {idx.id=1}    # check if there are individual labels, which changes the columns that need to be extracted
    idx.first.cluster.column=which(dat[1,]==":")+1
    
    dat.export=dat[,c(idx.id,idx.first.cluster.column:ncol(dat))]
    
    colnames(dat.export)=c("Label",paste(c(rep("C",(ncol(dat.export)-1)/3),rep(c("LCI","UCI"),(ncol(dat.export)-1)/3)),c(1:((ncol(dat.export)-1)/3),rep(1:((ncol(dat.export)-1)/3),each=2)),sep=""));colnames(dat.export)
    dat.export$Label=as.character(dat.export$Label)
    output.list[[i]]=dat.export
    for(ii in 1:((ncol(dat.export)-1)/3))
    {
      if(length(which(dat.export[ii+1]>Kcrit))>0)
      {
      subset.list[[list.teller]]=dat.export[which(dat.export[ii+1]>Kcrit),1]
      info.data[list.teller,]=c(i,ii,list.teller,length(subset.list[[list.teller]]),(ncol(dat.export)-1)/3)
      }
      if(length(which(dat.export[ii+1]>Kcrit))==0)
      {
        subset.list[[list.teller]]=character(0)
        info.data[list.teller,]=c(i,ii,list.teller,length(subset.list[[list.teller]]),(ncol(dat.export)-1)/3)
      }
      list.teller=list.teller+1
    }
  
    print(paste(i," of ",length(input.files)))
  } 
  print(paste("       ",N," individuals",sep=""))
  print(paste("       ",length(subset.list)," confident subsets",sep=""))
  
  
  info.data=info.data[order(info.data[,4],decreasing=T),]
  out=list(output.list,subset.list,info.data)
  
  return(out)
}
