library(spatstat)

nncrossCross <- function(synapseSet,by="type"){
  tt <- unique(synapseSet[[by]])
  distTable <- dplyr::bind_rows(lapply(tt,function(tRef){
    matRef <-sapply(tt,function(tTo){
     k <- ifelse(tRef==tTo,2,1)
     nncrossSynapses(synapseSet[synapseSet[[by]]==tRef,],synapseSet[synapseSet[[by]]==tTo,],k=k) 
    })
    colnames(matRef) <- tt
    mutate(as.data.frame(matRef),type=tRef)
  }))
  pivot_longer(distTable,-type,names_to="to",values_to="distance")
}

nncrossSynapses <- function(sA,sB,...){
  nncross(ppx(sA[,c("x","y","z")]),ppx(sB[,c("x","y","z")]),what="dist",...)
}