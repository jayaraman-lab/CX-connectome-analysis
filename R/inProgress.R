connectivityMatrix <- function(connObj,
                               slctROI,
                               allToAll=FALSE,
                               from="name.from",
                               to="name.to",
                               value="weightRelative"){UseMethod("connectivityMatrix")}

connectivityMatrix.data.frame <- function(connObj,
                                          slctROI,
                                          allToAll=FALSE,
                                          from="name.from",
                                          to="name.to",
                                          value="weightRelative"){
  #' Builds a connectivity matrix from a connection object
  #' @param connObj : a connection table
  #' @param slctROI : which ROI to consider
  #' @param allToAll : whether to build a square matrix or just a from -> to matrix
  #' @param from : which field to use as a "source" (default "name.from")
  #' @param to : which field to use as a "target" (default "name.to")
  #' @param value : which field to use to fill the matrix (default "weightRelative")
  #' 
  
  connObj <- filter(connObj,roi==slctROI)
  if (any(is.na(c(connObj[[to]],connObj[[from]])))){
    warning("Some to/from entries are NA, using retype.na function.")
    connObj <- retype.na(connObj)
  }
  
  connObj[[to]] <- as.character(connObj[[to]])
  connObj[[from]] <- as.character(connObj[[from]])
  if (nrow(distinct_at(connObj,c(from,to))) != nrow(connObj)){
    stop("Multiple entries for some of the from/to combinations. You need to either 
         use different from/to or summarize your data.frame beforehand.")}
  
  if (allToAll){
    ins <- unique(c(connObj[[from]],connObj[[to]]))
    outs <- ins
  }else{
    ins <- unique(connObj[[from]])
    outs <- unique(connObj[[to]])
  }
  
  
  outMat <- matrix(0,nrow=length(ins),ncol=length(outs),dimnames=list("Inputs"=ins,"Outputs"=outs))
  
  for (l in 1:nrow(connObj)){
    outMat[connObj[[from]][l],connObj[[to]][l]] <- connObj[[value]][l]
  }
  outMat
}