

collectSynapses <- function(neuronsTable,ROIs,polarity=c("post","pre")){
  polarity <- match.arg(polarity)
  syns <- dplyr::bind_rows(pbapply::pblapply(ROIs,function(r){
    synSet <- neuprint_get_synapses(neuronsTable$bodyid,roi=r,chunk=5,progress=FALSE) 
    if (!is.null(synSet)){
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

