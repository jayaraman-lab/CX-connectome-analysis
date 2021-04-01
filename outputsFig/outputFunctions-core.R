library(igraph)
library(Matrix)


#' Collect all synapses from neurons in `neuronsTable` in a set of ROIs `ROIs` of a given polarity
#' 
collectSynapses <- function(neuronsTable,ROIs,polarity=c("post","pre")){
  polarity <- match.arg(polarity)
  if(is.null(ROIs)) {ROIs <- list("1"=NULL)}
  syns <- dplyr::bind_rows(pbapply::pblapply(ROIs,function(r){
    synSet <- neuprint_get_synapses(neuronsTable$bodyid,roi=r,chunk=5,progress=FALSE) 
    if (!is.null(synSet)){
      if(is.null(r)){r <- "all"}
      synSet <- filter(mutate(synSet,roi=r),prepost==(ifelse(polarity=="pre",0,1)))
    }                
    return(synSet)                
  }))
  syns <- mutate(syns,
                 type=neuronsTable$type[match(bodyid,neuronsTable$bodyid)],
                 databaseType=neuronsTable$databaseType[match(bodyid,neuronsTable$bodyid)],
                 supertype = neuronsTable$supertype2[match(bodyid,neuronsTable$bodyid)]
  )
  
  syns
}

#' A custom retyping that splits DNb01 in 2 groups (personal communication rom Shigehiro Namiki that they are indeed
#' two different types)
customRetyping <- function(connections,postfix=c("raw", "to", "from")){
  # DNb01
  connections <- redefineTypeByBodyId(connections,sets=list(1566597156,1655997973),nameModifiers=c("_1_R","_2_R"),postfix=postfix,redefinePartners=T)
  # PLP042_a divided in two groups
  connections <- redefineTypeByBodyId(connections,sets=list(951018323,1224400894),nameModifiers=c("_2_R","_2_R"),postfix=postfix,redefinePartners=T)
  connections <- cxRetyping(connections,postfix=postfix)
}

#' This function takes an adjacency matrix and multiplies it recursively to convergence, restricting it to a given set 
#' of sources or targets
#' @param connmat An adjacency matrix (as returned by igraph::as_adjacency_matrix for example)
#' @param endpoint Indexes of the nodes to be considered as sources or targets
#' @param polarity Whether to run the algorithm from a set of sources (for example the CX output neurons) or a set of 
#' targets (for example the known neurons in the network). If running from source, the adjacency matrix should be in 
#' relative weights (or a metric <=1 when summing all the inputs to a given neuron), whereas when running from targets, 
#' it should be in a metric <=1 when summing all the outputs of a given neuron (we used a relative weights normalized)
#' @param eps The convergence threshold on the norm of the matrix
#' @param maxIt The maximum number of iterations
#' @return A list of n matrices, corresponding to the connectivity matrices for pathways of length 1:n
#' @details If you run this function with a lot of iterations, this could create a very large object. 
#' @seealso 'endpointConnections' does the same but only returns one matrix, the sum over the pathways of all lenghts
endpointConnections_raw <- function(connmat,endpoint=1,polarity=c("to_targets","from_sources"),eps=10^-3,maxIt=100){
  polarity <- match.arg(polarity)
  if(polarity=="to_targets") connmat <- t(connmat)
  
  simple <- vector(mode="list")
  
  simple[[1]] <- connmat[endpoint,]
  i <- 1
  while(i<maxIt & norm(simple[[i]])>eps){
    multip <- simple[[i]]
    multip[,endpoint] <- 0
    simple[[i+1]] <- multip %*% connmat
    simple[[i+1]] <- simple[[i+1]][endpoint,]
    i <- i+1  
  }
  
  if(polarity=="to_targets"){ 
    for(i in seq(length(simple))){
      simple[[i]] <- t(simple[[i]])
    }}
  return(simple)
}


#' This function takes an adjacency matrix and multiplies it recursively to convergence, restricting it to a given set 
#' of sources or targets
#' #' @seealso 'endpointConnections_raw'
endpointConnections <- function(connmat,endpoint=1,polarity=c("to_targets","from_sources"),eps=10^-3,maxIt=100){
  polarity <- match.arg(polarity)
  if(polarity=="to_targets") connmat <- t(connmat)
  
  simple <- connmat[endpoint,]
  res <- simple
  i <- 1
  while(i<maxIt & norm(simple)>eps){
    print(norm(simple))
    multip <- simple
    multip[,endpoint] <- 0
    simple <- multip %*% connmat
    simple <- simple[endpoint,]
    res <- res+simple
    i <- i+1  
  }
  res <- res[endpoint,]
  if(polarity=="to_targets") res <- t(res)
  
  return(res)
}

#' Supertyping customized for the outputs section.
customOutputSupertype <- function(typesTable,postfix = c("raw", "to", "from"),extraTypes=""){
  postfix <- match.arg(postfix)
  if (postfix == "raw") postfix <- "" else postfix <- paste0(".",postfix)
  
  postSp <- paste0("supertype",postfix)
  postDT <- paste0("databaseType",postfix)
  postSp3 <- paste0("supertype3",postfix)
  postT <- paste0("type",postfix)
  
  mutate(typesTable,
         !!paste0("customSupertype",postfix) := case_when(
           !!(sym(postSp)) %in% c("PFL","PFR") ~ as.character(!!(sym(postDT))),
           grepl("ExR[2-6]",!!(sym(postSp))) ~ "ExR2-6",
           !!(sym(postSp3))=="EB Columnar" ~ "EB Columnar",
           grepl("FB[6-7]",!!(sym(postSp))) ~ "FB6-7",
           grepl("FB[8-9]",!!(sym(postSp))) ~ "FB8-9",
           !!(sym(postT)) %in% extraTypes ~ as.character(!!(sym(postSp3))),
           TRUE ~ as.character(!!(sym(postSp)))))
  
}



simulate_roi_contra <- function(roiSummary,roiSet){
  simulated <- mutate(roiSummary,type=lrInvert(type),
                      roi=gsub("(R)","(L)",roi,fixed=T),
                      side=ifelse(side=="Right","Left","Central")) %>% filter(grepl("_L$|_L[1-9]$|_L[1-9]/[1-9]$|_L[1-9]_C[1-9]$|_L[1-9]_C[1-9]_irreg$|_L_C[1-9]_irreg$|_L_small$|_R$|_R[1-9]$|_R[1-9]/[1-9]$|_R[1-9]_C[1-9]$|_R[1-9]_C[1-9]_irreg$|_R_C[1-9]_irreg$|_R_small$",type))
  rbind(roiSummary,simulated)
}


