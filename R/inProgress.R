connectivityMatrix <- function(connObj,slctROI,allToAll=FALSE,from="name.from",to="name.to",value="weightRelative"){UseMethod("connectivityMatrix")}

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
  
  
  connObj <- retype.na(connObj) %>% filter(roi==slctROI)
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