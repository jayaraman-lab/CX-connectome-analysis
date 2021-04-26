library(paletteer)
library(alphahull)
source("InputOutputByTypeUtils.R")

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
  roi_Innervate <- roi_Innervate %>% select(-c(voxels)) %>% #,cellBodyFiber
    replace_na(list(ROI_pre = 0,ROI_post = 0,originalInstance=FALSE)) %>%
    mutate(databaseType = type) ## Convenience column for when types are changed

  roi_Innervate <-roi_Innervate %>% group_by(type) %>%
                   mutate(typePercentage = sum(originalInstance)/n()) %>%
                   ungroup() %>%
                   filter(typePercentage > minTypePercentage)
  
  return(roi_Innervate)
}

getTypesInRoiTable <- function(ROI,lateralize=FALSE,...){
  #' Returns a neuronBag object of all the neurons in the ROI
  #' @param ROI : the ROI to consider
  #' @param lateralize : should the neuron types be divided in left/right (default FALSE)
  #' @param big : if TRUE, run through a pblapply call
  #' @param clN : if big is TRUE, this is the number of cores to use.
  #' 
  neuronTable <- getNeuronsInRoiTable(ROI,minTypePercentage=ifelse(lateralize,0.25,0.5)) ## Remove types if less than 
  ## 25% of the instances touch (l/R)
  typesUnfiltered <- unique(neuronTable$type)
  
  roiConnections <- buildInputsOutputsByType(neuronTable,fixed=TRUE,...)

  
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
  roiT <- data.frame(level1 = roiH$roi[roiH$parent == "hemibrain"],stringsAsFactors = F) %>% filter(!(level1 %in% c("hemibrain","AOT(R)","GC","GF(R)","mALT(R)","POC","mALT(L)")))
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
  roiT$level0 <- delateralize(roiT$level1)
  #fineOrder <- c(which(roiT$side2!="Left"),rev(which(roiT$side2 == "Left")))
  roiT <- roiT %>% mutate_at(c("level2","level3","level4"),function(a) factor(a,levels=unique(a)))
  #roiT <- arrange(roiT,level4)
  roiT
}

delateralize <- function(roiName){
  gsub("(L)","",gsub("(R)","",roiName,fixed=TRUE),fixed=TRUE)
}

selectRoiSet <- function(roiTree,default_level=2,exceptions=NULL,exceptionLevelMatch = default_level){
  if (!is.null(exceptions)){
    levelEx <- paste0("level",exceptionLevelMatch) 
    normalRois <- roiTree %>% filter(!((!!as.name(levelEx)) %in% names(exceptions))) %>%
       mutate(roi = (!!as.name(paste0("level",default_level))))
    
    exceptionsRois <- roiTree %>% filter(((!!as.name(levelEx)) %in% names(exceptions))) 
    
    roisEx <- as.character(exceptionsRois[[levelEx]])
    customLev <- sapply(roisEx,function(r) paste0("level",exceptions[[r]]))
    exceptionsRois$roi <- sapply(1:length(customLev),function(i) as.character(exceptionsRois[[customLev[i]]][i]))
    
    rois <- bind_rows(normalRois,exceptionsRois)
  }else{
    rois <- roiTree %>% mutate(roi = (!!(as.name(paste0("level",default_level)))))
  }
  
  rois <- rois %>% arrange(side2,level1) %>% 
    mutate(roi = factor(roi,levels=unique(roi)))
  
  return(distinct(rois))
}
 
roisPalette <- function(favoriteRegion="CX",my_palette=paletteer_d("Polychrome::palette36")){
  rois <- getRoiTree()
  roiL <- unique(delateralize(c(as.character(rois$level1),as.character(rois$level2[rois$level1==favoriteRegion]))))
  pal <- my_palette[1:length(roiL)]
  names(pal) <- roiL
  pal
}

roiOutline <- function(roiMesh,axis=c("x","y"),alpha=100,roiName){UseMethod("roiOutline")}

roiOutline.mesh3d <- function(roiMesh,alpha=100,roiName =deparse(substitute(roiMesh))){
  roiPts <-  data.frame(dotprops(roiMesh)$points)
  names(roiPts) <- c("x","y","z")
  roiHullxy <- ahull(x=roiPts$x,y=roiPts$y,alpha=alpha)
  roiHullxz <- ahull(x=roiPts$x,y=roiPts$z,alpha=alpha)
  
  roiOutxy <- data.frame(roiHullxy$arcs) %>% mutate(x=c1,y=c2,proj="xy",roi=roiName) %>% select(x,y,proj,roi)
  roiOutxy <-  bind_rows(roiOutxy,roiOutxy[1,])
  
  roiOutxz <- data.frame(roiHullxz$arcs) %>% mutate(x=c1,y=c2,proj="xz",roi=roiName) %>% select(x,y,proj,roi)
  roiOutxz <-  bind_rows(roiOutxz,roiOutxz[1,])
  
  bind_rows(roiOutxy,roiOutxz)
}

roiOutline.character <- function(roi,alpha=100){
  roiMesh <- neuprint_ROI_mesh(roi)
  roiOutline(roiMesh,alpha=alpha,roiName=roi)
}


getCompRoiInfo <- function(bodyidsIn,
                           bodyidsOut,
                           roi,
                           newDataset,
                           oldDataset){
  roiTableIn <- getCompRoiInfoInternal(bodyidsIn,roi,polarity="inputs",newDataset=newDataset,oldDataset=oldDataset)
  roiTableOut <- getCompRoiInfoInternal(bodyidsOut,roi,polarity="outputs",newDataset=newDataset,oldDataset=oldDataset)
  
  bind_rows(roiTableIn,roiTableOut) 
}

getCompRoiInfoInternal <- function(bodyids,ROI,polarity=c("inputs","outputs"),newDataset,oldDataset){
  polarity <- match.arg(polarity)
  newRefTable <- neuprint_get_meta(bodyids,dataset=newDataset)
  
  newRoiTable <- getRoiInfo(bodyids,dataset=newDataset) %>% filter(roi == ROI) %>% mutate(version="new")
  oldRoiTable <- getRoiInfo(bodyids,dataset=oldDataset) %>% filter(roi == ROI) %>% mutate(version="old")
  
  roiTable <- bind_rows(newRoiTable,oldRoiTable) 
  if (polarity == "inputs")
    roiTable <- pivot_wider(roiTable,names_from = version,values_from=upstream,id_cols = bodyid)
  else
    roiTable <- pivot_wider(roiTable,names_from = version,values_from=downstream,id_cols = bodyid)
  
  mutate(roiTable,databaseType=newRefTable$type[match(bodyid,newRefTable$bodyid)]) %>% supertype() %>% mutate(side=polarity)
}