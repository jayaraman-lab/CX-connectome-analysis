library(parallel)

getNeuronsInRoiTable <- function(slctROI,minTypePercentage=0.5) {
  #' Returns a table of all instances of neurons of types with non zero pre/post counts in slctROI.
  #' @param slctROI : the ROI to look for
  #' @return : a table of metadata for all neurons in the ROI, with extra columns \code{ROI_pre}
  #'  and \code{ROI_post}, the counts in the queried ROI.
  
  roi_Innervate <- neuprint_bodies_in_ROI(slctROI) %>%
    mutate(originalInstance = TRUE)
  metaRoi <- neuprint_get_meta(roi_Innervate) %>% drop_na(type)
 
  ## Get all instances of the types touching the ROI
  all_neurons <- getTypesTable(unique(metaRoi$type))
  
  ## Join to the full table of type members (and fill with zero where the extra instances do not
  ## innervate)
  roi_Innervate <- left_join(all_neurons,roi_Innervate,by=c("bodyid","pre","post","voxels"))
  roi_Innervate <- roi_Innervate %>% select(-c(voxels,cellBodyFiber)) %>%
    replace_na(list(ROI_pre = 0,ROI_post = 0,originalInstance=FALSE)) %>%
    mutate(databaseType = type) ## Convenience column for when types are changed

  roi_Innervate <-roi_Innervate %>% group_by(type) %>%
                   mutate(typePercentage = sum(originalInstance)/n()) %>%
                   ungroup() %>%
                   filter(typePercentage > minTypePercentage)
  
  return(roi_Innervate)
}

getTypesInRoiTable <- function(ROI,lateralize=FALSE,big=TRUE){
  neuronTable <- getNeuronsInRoiTable(ROI,minTypePercentage=ifelse(lateralize,0.25,0.5)) ## Remove types if less than 
  ## 25% of the instances touch (l/R)
  typesUnfiltered <- unique(neuronTable$type)
  
  if (big){
    roiConnections <- pblapply(typesUnfiltered,buildInputsOutputsByType,fixed=TRUE,cl=5)
    roiConnections <- do.call(bind_InoutLists,roiConnections)
  }else{
    roiConnections <- buildInputsOutputsByType(typesUnfiltered,fixed=TRUE)
  }
  
  if (lateralize == TRUE){
    roiConnections <- lateralizeInputOutputList(roiConnections)
  }
  
  roiConnections
}

typesInROI <- function(roiConnections,ROI){
  typesUnfiltered <- roiConnections$names$type
  inputs <- roiConnections$inputs %>% filter((roi == ROI) & (type.to %in% typesUnfiltered) &
                                              (type.from %in% typesUnfiltered))
  outputs <- roiConnections$outputs %>% filter((roi == ROI) & (type.to %in% typesUnfiltered) &
                                                (type.from %in% typesUnfiltered))
  
  roiTypes <- data.frame(type = unique(c(inputs$type.to,outputs$type.from))) %>%
    mutate(databaseType = roiConnections$names$databaseType[match(type,roiConnections$names$type)])
  
  return(roiTypes)
}