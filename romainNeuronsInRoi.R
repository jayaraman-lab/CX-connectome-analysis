library(parallel)

getNeuronsInRoiTable <- function(slctROI,minTypePercentage=0.5) {
  #' Returns a table of all instances of neurons of types with non zero pre/post counts in slctROI.
  #' @param slctROI : the ROI to look for
  #' @param minTypePercentage : the minimum proportion of the instances of a type that should be innervating the ROI for 
  #' it to be considered
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

getTypesInRoiTable <- function(ROI,lateralize=FALSE,big=TRUE,clN=5){
  #' Returns a neuronBag object of all the neurons in the ROI
  #' @param ROI : the ROI to consider
  #' @param lateralize : should the neuron types be divided in left/right (default FALSE)
  #' @param big : if TRUE, run through a pblapply call
  #' @param clN : if big is TRUE, this is the number of cores to use.
  #' 
  neuronTable <- getNeuronsInRoiTable(ROI,minTypePercentage=ifelse(lateralize,0.25,0.5)) ## Remove types if less than 
  ## 25% of the instances touch (l/R)
  typesUnfiltered <- unique(neuronTable$type)
  
  roiConnections <- buildInputsOutputsByType(neuronTable,fixed=TRUE,big=big,nc=clN)

  
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

getRoiTree <- function(){
  #' Build a more useful roi table, with things ordered in a semi useful way
  roiH <- neuprint_ROI_hierarchy() %>% mutate_all(as.character)
  roiT <- data.frame(level1 = roiH$roi[roiH$parent == "hemibrain"],stringsAsFactors = F) %>% filter(!(level1 %in% c("hemibrain","AOT(R)","GC","GF(R)","mALT(R)","POC")))
  roiT <- left_join(roiT,roiH,by=c("level1"="parent")) %>% rename(level2 = roi) %>% mutate(level2 = ifelse(is.na(level2),level1,level2))
  roiT <- left_join(roiT,roiH,by=c("level2"="parent")) %>% rename(level3 = roi) %>% mutate(level3 = ifelse(is.na(level3),level2,level3))
  roiT <- left_join(roiT,roiH,by=c("level3"="parent")) %>% rename(level4 = roi) %>% mutate(level4 = ifelse(is.na(level4),level3,level4))
  roiT <- roiT %>% mutate(side4 = "Central",
                          side2 = "Central")
  roiT$side4[grepl("(L",roiT$level4,fixed=T)] <- "Left"
  roiT$side4[grepl("(R",roiT$level4,fixed=T)] <- "Right"
  roiT$side2[grepl("(L",roiT$level2,fixed=T)] <- "Left"
  roiT$side2[grepl("(R",roiT$level2,fixed=T)] <- "Right"
  
  roiT$side4 <- factor(roiT$side4,levels=(c("Right","Central","Left")))
  roiT$side2 <- factor(roiT$side2,levels=(c("Right","Central","Left")))
  roiT$level1 <- factor(roiT$level1,levels= c("OL(R)","AL(R)","MB(+ACA)(R)","LH(R)","PENP","GNG","VLNP(R)","SNP(R)","VMNP","INP","LX(R)","CX","LX(L)","SNP(L)","MB(L)","AL(L)"))
  
  roiT <- arrange(roiT,side2,level1)
  #fineOrder <- c(which(roiT$side2!="Left"),rev(which(roiT$side2 == "Left")))
  roiT <- roiT %>% mutate_at(c("level2","level3","level4"),function(a) factor(a,levels=unique(a)))
  #roiT <- arrange(roiT,level4)
  roiT
}

selectRoiSet <- function(roiTree,default_level=2,exceptions=NULL,exceptionLevelMatch = 2){
  if (!is.null(exceptions)){
    normalRois <- roiTree %>% filter(!((!!as.name(paste0("level",exceptionLevelMatch))) %in% names(exceptions))) %>%
       mutate(customRois = (!!as.name(paste0("level",default_level))))
    
    exceptionsRois <- roiTree %>% filter(((!!as.name(paste0("level",exceptionLevelMatch))) %in% names(exceptions)))
    exceptionsRois$customRois <- unlist(sapply(names(exceptions),function(n){
      roiTree[[paste0("level",exceptions[[n]])]][roiTree[[paste0("level",exceptionLevelMatch)]] == n] 
    }))
    rois <- bind_rows(normalRois,exceptionsRois)
  }else{
    rois <- rois %>% mutate(customRois = (!!(paste0("level",default_level))))
  }
  
  rois <- rois %>% arrange(level4) %>% 
    mutate(customRois = factor(customRois,levels=unique(customRois)))
  
  return(distinct(rois))
}
 

