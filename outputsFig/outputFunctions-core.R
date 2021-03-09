library(igraph)
library(Matrix)


## A custom retyping that splits DNb01 in 2 groups
customRetyping <- function(connections,postfix=c("raw", "to", "from")){
  # DNb01
  connections <- redefineTypeByBodyId(connections,sets=list(1566597156,1655997973),nameModifiers=c("_1_R","_2_R"),postfix=postfix,redefinePartners=T)
  # PLP042_a divided in two groups
  connections <- redefineTypeByBodyId(connections,sets=list(951018323,1224400894),nameModifiers=c("_2_R","_2_R"),postfix=postfix,redefinePartners=T)
  connections <- cxRetyping(connections,postfix=postfix)
}

## Supertyping customized
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

endpointConnections_raw <- function(connmat,endpoint=1,polarity=c("downstream","upstream"),eps=10^-3,maxIt=100){
  polarity <- match.arg(polarity)
  if(polarity=="downstream") connmat <- t(connmat)
  
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
  
  if(polarity=="downstream"){ 
    for(i in seq(length(simple))){
      simple[[i]] <- t(simple[[i]])
  }}
  return(simple)
}

endpointConnections <- function(connmat,endpoint=1,polarity=c("downstream","upstream"),eps=10^-3,maxIt=100){
  polarity <- match.arg(polarity)
  if(polarity=="downstream") connmat <- t(connmat)
  
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
  if(polarity=="downstream") res <- t(res)
  
  return(res)
}

collectSynapses <- function(neuronsTable,ROIs,polarity=c("post","pre")){
  polarity <- match.arg(polarity)
  if(is.null(ROIs)) {ROIs <- list("1"=NULL)}
  syns <- dplyr::bind_rows(pbapply::pblapply(ROIs,function(r){
    synSet <- neuprint_get_synapses(neuronsTable$bodyid,roi=r,chunk=5,progress=FALSE) 
    if (!is.null(synSet)){
      if(is.null(r)){r <- "all"}
      synSet <- filter(mutate(synSet,roi=r),prepost==(ifelse(polarity=="post",0,1)))
    }                
    return(synSet)                
  }))
  syns <- mutate(syns,
                 type=neuronsTable$type[match(bodyid,neuronsTable$bodyid)],
                 databaseType=neuronsTable$databaseType[match(bodyid,neuronsTable$bodyid)],
                 supertype = neuronsTable$supertype2[match(bodyid,neuronsTable$bodyid)],
                 cluster=neuronsTable$cluster[match(bodyid,neuronsTable$bodyid)]
  )
  
  syns
}

