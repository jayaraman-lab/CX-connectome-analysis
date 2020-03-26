connectivityMatrix <- function(connObj,
                               slctROIs,
                               allToAll=FALSE,
                               from="name.from",
                               to="name.to",
                               value="weightRelative",
                               ref=c("inputs","outputs")){UseMethod("connectivityMatrix")}

connectivityMatrix.data.frame <- function(connObj,
                                          slctROIs,
                                          allToAll=FALSE,
                                          from="name.from",
                                          to="name.to",
                                          value="weightRelative",
                                          ref=c("inputs","outputs")){
  #' Builds a connectivity matrix from a connection object
  #' @param connObj : a connection table
  #' @param slctROIs : which ROIs to consider
  #' @param allToAll : whether to build a square matrix or just a from -> to matrix
  #' @param from : which field to use as a "source" (default "name.from")
  #' @param to : which field to use as a "target" (default "name.to")
  #' @param value : which field to use to fill the matrix (default "weightRelative")
  #' @param ref : which channel will be used as the "reference" (to be the columns of the output). The
  #' other channel gets .roi affixed to their names in case several ROIs are considered
  #' 
  
  ref <- match.arg(ref)
  connObj <- filter(connObj,roi %in% slctROIs)
  if (any(is.na(c(connObj[[to]],connObj[[from]])))){
    warning("Some to/from entries are NA, using retype.na function.")
    connObj <- retype.na(connObj)
  }
  
  connObj[[to]] <- as.character(connObj[[to]])
  connObj[[from]] <- as.character(connObj[[from]])
  if (nrow(distinct_at(connObj,c(from,to,"roi"))) != nrow(connObj)){
    stop("Multiple entries for some of the from/to combinations. You need to either 
         use different from/to or summarize your data.frame beforehand.")}
  
  if (allToAll){
    bare <- unique(c(connObj[[from]],connObj[[to]]))
    if (length(slctROIs)==1){
      ins <- bare
      outs <- bare
    }else{
      if (ref=="inputs"){
        ins <- unique(c(paste(connObj[[from]],connObj[["roi"]],sep="."),c(paste(connObj[[to]],connObj[["roi"]],sep="."))))
        outs <- bare
      }else{
        ins <- bare
        outs <- unique(c(paste(connObj[[from]],connObj[["roi"]],sep="."),c(paste(connObj[[to]],connObj[["roi"]],sep="."))))
      }
    }
  }else{
    ins <- unique(connObj[[from]])
    outs <- unique(connObj[[to]])
  }
  
  if (length(slctROIs)>1){
    if (ref=="inputs"){
      ins <- unique(paste(connObj[[from]],connObj[["roi"]],sep="."))
      outMat <- matrix(0,nrow=length(ins),ncol=length(outs),dimnames=list("Inputs"=ins,"Outputs"=outs))
      
      for (l in 1:nrow(connObj)){
        outMat[paste0(connObj[[from]][l],".",connObj[["roi"]][l]),connObj[[to]][l]] <- connObj[[value]][l]
      }
    }else{
      outs <- unique(paste(connObj[[to]],connObj[["roi"]],sep="."))
      
      outMat <- matrix(0,nrow=length(ins),ncol=length(outs),dimnames=list("Inputs"=ins,"Outputs"=outs))
      
      for (l in 1:nrow(connObj)){
        outMat[connObj[[from]][l],paste0(connObj[[to]][l],".",connObj[["roi"]][l])] <- connObj[[value]][l]
      }
      outMat <- t(outMat)
    }
  }else{
    for (l in 1:nrow(connObj)){
      outMat <- matrix(0,nrow=length(ins),ncol=length(outs),dimnames=list("Inputs"=ins,"Outputs"=outs))
      outMat[connObj[[from]][l],connObj[[to]][l]] <- connObj[[value]][l]
    }
    if (ref=="outputs") outMat <-  t(outMat)
  }
  
  
  outMat
}

hClust_connMat <- function(connMat,rowSelect=rownames(connMat)){
  connMatS <- connMat[rownames(connMat) %in% rowSelect,]
  connCor <- cor(connMatS,method="spearman")
  connCor[is.na(connCor)] <- 1
  d <- as.dist(1-connCor)
  hc <- hclust(d)
  hc
}