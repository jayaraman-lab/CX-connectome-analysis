source("neuprintQueryUtils.R")
source("R/supertypeUtils.R")
library(pbapply)
library(ggnewscale)
#library(ggiraph)
library(paletteer)

## Define an S3 class, neuronBag, to hold the connectivity information of a bunch of neurons
neuronBag <- function(outputs,inputs,names,outputs_raw,inputs_raw,outputsTableRef){
  #' neuronBag objects contain information about the connectivity of a bunch of neurons
  #' @field names : a table of metadata associated with the neurons in the bag -- with an extra column, 'databaseType'
  #' which keeps the type name used in the database
  #' @field outputs : a type to type connectivity table of the outputs of the set of neurons
  #' @field inputs : a type to type connectivity table of the inputs of the set of neurons
  #' @field outputs_raw : a neuron to neuron connection table of the outputs of the set of neurons
  #' @field inputs_raw : a neuron to neuron connection table of the inputs of the set of neurons
  #' @field outputsTableRef : a table holding all instances of all the output types and their associated metadata
  res <- list(outputs=outputs,
              inputs=inputs,
              names=names,
              outputs_raw=outputs_raw,
              inputs_raw=inputs_raw,
              outputsTableRef=outputsTableRef)
  attr(res,"class") <- "neuronBag"
  return(res)
}

is.neuronBag <- function(x) inherits(x,"neuronBag")

buildInputsOutputsByType <- function(typeQuery,fixed=FALSE,by.roi=TRUE,...){
  #' Builds a neuronBag object either from a vector of query strings or a metadata data.frame.
  #' @param typeQuery : either a vector of queries (similar to neuprint_search queries) or a 
  #' metadata data.frame for a set of neurons
  #' @param fixed : if typeQuery is a query string, is it fixed?
  #' @param by.roi : return results by ROI or just global weights?
  #' @param verbose : Inform about progress if TRUE
  #' @param ... : to be passed to getConnectionTable
  UseMethod("buildInputsOutputsByType")}

buildInputsOutputsByType.character <- function(typeQuery,fixed=FALSE,by.roi=TRUE,verbose=FALSE,...){
  TypeNames <- distinct(bind_rows(lapply(typeQuery,neuprint_search,field="type",fixed=fixed))) %>%
                  mutate(databaseType = type)
  buildInputsOutputsByType(TypeNames,fixed=FALSE,by.roi=by.roi,verbose=verbose,...)
}
  
buildInputsOutputsByType.data.frame <- function(typeQuery,fixed=FALSE,selfRef=FALSE,by.roi=TRUE,verbose=FALSE,...){
  #'
  #'@param selfRef : Should the input data.frame be used as the type reference (use if you already renamed
  #'neurons/types in that data frame)
  #'
  
  if (verbose) message("Calculate raw outputs")
  outputsR <- getConnectionTable(typeQuery,synapseType = "POST",by.roi=by.roi,verbose=verbose,...)
  
  if (verbose) message("Calculate raw inputs")
  inputsR <- getConnectionTable(typeQuery,synapseType = "PRE",by.roi=by.roi,verbose=verbose,...)
  
  if (verbose) message("Calculate type to type outputs")
  if (length(outputsR)==0){OUTByTypes <- NULL
                           outputsTableRef <- NULL
                           unknowns <- NULL
                         }else{
    OUTByTypes <- getTypeToTypeTable(outputsR)
    outputsR <- retype.na(outputsR)
    outputsTableRef <- getTypesTable(unique(outputsR$type.to))
    unknowns <- retype.na_meta(neuprint_get_meta(unique(outputsR$to[!(outputsR$to %in% outputsTableRef$bodyid)])))
                         }
  
  if (verbose) message("Calculate type to type inputs")
  if (nrow(inputsR)==0){INByTypes <- NULL}else{
    if (selfRef){
      INByTypes <- getTypeToTypeTable(inputsR,typesTable = typeQuery)
    }else{
      INByTypes <- getTypeToTypeTable(inputsR)
    }
    inputsR <- retype.na(inputsR)}
  
  return(neuronBag(outputs = OUTByTypes,
              inputs = INByTypes,
              names = typeQuery,
              outputs_raw = outputsR,
              inputs_raw = inputsR,
              outputsTableRef = bind_rows(outputsTableRef,unknowns)
              ))
}

redefineTypeByNameInList <- function(IOList,
                                     typeList,
                                     pattern,
                                     newPostFixes,
                                     perl=FALSE){

  oldOutputs <- IOList$outputs
  oldInputs <- IOList$inputs
  
  for (t in typeList){
    for (col in c("from","to")){
      for (df in c("inputs_raw","outputs_raw")){
        IOList[[df]] = redefineTypeByName(IOList[[df]],
                                         type=t,
                                         pattern=pattern,
                                         newPostFixes=newPostFixes,
                                         type_col=paste0("type.",col),
                                         name_col=paste0("name.",col),
                                         perl=perl)
                  
      }
    }
    IOList$names = redefineTypeByName(IOList$names,
                                      type=t,
                                      pattern=pattern,
                                      newPostFixes=newPostFixes,
                                      type_col="type",
                                      name_col="name",
                                      perl=perl)
    
    IOList$outputsTableRef <- redefineTypeByName(IOList$outputsTableRef,
                                                 type=t,
                                                 pattern=pattern,
                                                 newPostFixes=newPostFixes,
                                                 type_col="type",
                                                 name_col="name",
                                                 perl=perl)
  }
  
 
  ## In case recursive modifs have been made
  
  IOList$outputs <- getTypeToTypeTable(IOList$outputs_raw,typesTable = IOList$outputsTableRef,oldTable = oldOutputs)
  IOList$inputs <- getTypeToTypeTable(IOList$inputs_raw,typesTable = IOList$names,oldTable = oldInputs)
  
  ## Exception: we want to keep connections that were lost through division of input types.
  

  return(IOList)
}

lateralizeInputOutputList <- function(inputOutputList,typeList=NULL){
  
  outputsLat <- lrSplit(inputOutputList$outputs_raw,nameCol = "name.from",typeCol = "type.from",typeList=typeList)
  outputsLat <- lrSplit(outputsLat,typeList=typeList,nameCol = "name.to",typeCol = "type.to")
  
  inputsLat <- lrSplit(inputOutputList$inputs_raw,nameCol = "name.from",typeCol = "type.from",typeList=typeList)
  inputsLat <- lrSplit(inputsLat,typeList=typeList,nameCol = "name.to",typeCol = "type.to")
 
  outputsRef <- lrSplit(inputOutputList$outputsTableRef,nameCol="name",typeCol="type",typeList=typeList)
  
  TypeNamesLat <- lrSplit(inputOutputList$names,nameCol = "name",typeCol="type",typeList=typeList)
 
  outByTypesLat <- getTypeToTypeTable(outputsLat,typesTable = outputsRef,oldTable = inputOutputList$outputs)
  inByTypesLat <- getTypeToTypeTable(inputsLat,typesTable = TypeNamesLat,oldTable = inputOutputList$inputs)
  
  
  return(neuronBag(outputs = outByTypesLat,
                   inputs = inByTypesLat,
                   names = TypeNamesLat,
                   outputs_raw = outputsLat,
                   inputs_raw=inputsLat,
                   outputsTableRef=outputsRef))
  
}

cxRetyping <- function(neurons,verbose=TRUE){
  #' Convenience function to deal with the tricky left/right asymetries
  if (verbose) message("Renaming PFL3")
  neurons <- redefineTypeByNameInList(neurons,typeList = c("PFL3"),pattern = "(^.*_L(?!.*irreg))|(^.*_R.*irreg)",perl=TRUE,newPostFixes = c("_L*","_R*"))
  if (verbose) message("Renaming PFL1/PFR_a")
  neurons <- redefineTypeByNameInList(neurons,typeList = c("PFR_a","PFL1"),pattern = "_L[2-7]|_R1",newPostFixes = c("_L*","_R*"))
  if (verbose) message("Renaming PFR_b")
  neurons <- redefineTypeByNameInList(neurons,typeList = c("PFR_b"),pattern = "(^.*_L(?!.*C9))|(^.*C1.*)",perl=TRUE,newPostFixes = c("_L*","_R*"))
  if (verbose) message("All other L/R retyping")
  neurons <- lateralizeInputOutputList(neurons)
  return(neurons)
}

bind_InoutLists <- function(...){
  full <- list(...)
  out <- neuronBag(outputs = distinct(bind_rows(lapply(full,function(i) i$outputs))),
                   inputs = distinct(bind_rows(lapply(full,function(i) i$inputs))),
                   names = distinct(bind_rows(lapply(full,function(i) i$names))),
                   outputs_raw = distinct(bind_rows(lapply(full,function(i) i$outputs_raw))),
                   inputs_raw = distinct(bind_rows(lapply(full,function(i) i$inputs_raw))),
                   outputsTableRef = distinct(bind_rows(lapply(full,function(i) i$outputsTableRef)))
              )
  return(out)
}

filter.neuronBag <- function(.nbag,filterPartners=FALSE,...){
  #` Meant to filter on $names
  #' @param nbag a neuronBag object to filter
  #' @param ... to be passed to a filtering function applied to the names field
  #' @param filterPartners : Whether to apply the filter to input/output neurons to
  .nbag$names <- filter(.nbag$names,...)
  
    .nbag$outputs <- filter(.nbag$outputs,type.from %in% .nbag$names$type)
    .nbag$outputs_raw <- filter(.nbag$outputs_raw,type.from %in% .nbag$names$type)
      
    .nbag$inputs <- filter(.nbag$inputs,type.to %in% .nbag$names$type)
    .nbag$inputs_raw <- filter(.nbag$inputs_raw,type.to %in% .nbag$names$type)
    
    if (filterPartners == TRUE){
      .nbag$outputs <- filter(.nbag$outputs,type.to %in% .nbag$names$type)
      .nbag$outputs_raw <- filter(.nbag$outputs_raw,type.to %in% .nbag$names$type)
      
      .nbag$inputs <- filter(.nbag$inputs,type.from %in% .nbag$names$type)
      .nbag$inputs_raw <- filter(.nbag$inputs_raw,type.from %in% .nbag$names$type)
    }
    
    .nbag$outputsTableRef <- filter(.nbag$outputsTableRef,type %in% .nbag$outputs$type.to)
    .nbag
}

getIntraBag <- function(nBag){
  #' Convenience to extract just connections between the central neurons of a bag
  filter(nBag,filterPartners = TRUE,type %in% unique(nBag$names$type))
}

summarizeConnectionTable <- function(connTable,groupFrom,groupTo,refOuts){
  
  groupType <- as.name(sub("\\.to","",groupTo))
  typesCount <- refOuts %>% group_by(!!(groupType)) %>%
    summarise(n=n())
  
  connTable <- connTable %>% 
    mutate(n = typesCount[["n"]][match(!!as.name(groupTo),typesCount[[groupType]])])
  
  connTable <- connTable %>% group_by(!!(as.name(groupFrom)),
                                      !!(as.name(groupTo)),
                                      to,
                                      roi) %>%
                             summarize(weightRelative = sum(weightRelative),
                                       weightRelativeTotal = sum(weightRelativeTotal),
                                       weight = sum(ROIweight),
                                       n=n[1],
                                       outputContribution = mean(outputContribution)) %>% ungroup() %>%
                             group_by(!!(as.name(groupFrom)),
                                      !!(as.name(groupTo)),
                                      roi) %>%
                             summarize(missingV = ifelse(is.null(n),0,n[1]-n()),
                                       varWeight = var(c(weightRelative,rep(0,missingV))),
                                       weightRelative = mean(c(weightRelative,rep(0,missingV))),
                                       weightRelativeTotal = mean(c(weightRelativeTotal,rep(0,missingV))),
                                       absoluteWeight = sum(weight),
                                       weight = mean(c(weight,rep(0,missingV))),
                                       outputContribution = outputContribution[1],
                                       n_targets = n(),
                                       n_type = n[1]) %>% select(-missingV) %>% ungroup()
  connTable
  
}
substractRois <- function(connections,largeRoi,smallRoi,newRoi,...){UseMethod("substractRois")}

substractRois.data.frame <- function(connections,largeRoi,smallRoi,newRoi,...){
  ## CHECK IT'S A RAW CONNECTION TABLE
  newRegionTable <- connections %>% 
    filter(roi %in% c(largeRoi,smallRoi)) %>% 
    group_by(to) %>%
    mutate(totalROIweight=totalROIweight[match(largeRoi,roi)] - ifelse(!(smallRoi %in% roi),0,totalROIweight[match(smallRoi,roi)])) %>%
    group_by(from) %>%
    mutate(totalPreROIweight=totalPreROIweight[match(largeRoi,roi)] - ifelse(!(smallRoi %in% roi),0,totalPreROIweight[match(smallRoi,roi)])) %>%
    group_by_if(names(.) %in% c(paste0(c("","name.","type.","databaseType.","previous.type."),"to"),paste0(c("","name.","type.","databaseType.","previous.type."),"from"),paste0("supertype.to",1:3),paste0("supertype.from",1:3))) %>%
    summarize(weight=weight[1],
              totalPreWeight=totalPreWeight[1],
              weightRelativeTotal=weightRelativeTotal[1],
              outputContributionTotal=outputContributionTotal[1],
              totalROIweight = totalROIweight[1],
              totalPreROIweight=totalPreROIweight[1],
              ROIweight=ROIweight[match(largeRoi,roi)] - ifelse(!(smallRoi %in% roi),0,ROIweight[match(smallRoi,roi)]),
              weightROIRelativeTotal=weightROIRelativeTotal[match(largeRoi,roi)] - ifelse(!(smallRoi %in% roi),0,weightROIRelativeTotal[match(smallRoi,roi)]),
              roi=newRoi
    ) %>%
    ungroup() %>%
    mutate(weightRelative=ROIweight/totalROIweight,
           outputContribution=ROIweight/totalPreROIweight)
  
  newRegionTable
}

substractRois.neuronBag <- function(connections,largeRoi,smallRoi,newRoi,...){
  new_inputsR <- substractRois(connections$inputs_raw,largeRoi,smallRoi,newRoi)
  new_outputsR <- substractRois(connections$outputs_raw,largeRoi,smallRoi,newRoi)
  new_inputs <- getTypeToTypeTable(new_inputsR,typesTable = connections$names,...)
  new_outputs <- getTypeToTypeTable(new_outputsR,typesTable = connections$outputsTableRef,...)
  neuronBag(outputs = new_outputs,
            inputs = new_inputs,
            names = connections$names,
            inputs_raw = new_inputsR,
            outputs_raw = new_outputsR,
            outputsTableRef = connections$outputsTableRef) 
  
}

localized <- function(nBag,rois,localRef,...){
  ## Take the perspective of the region as containing "all" the neurons of the considered types. Useful for some comparisons.
  new_inputsR <- nBag$inputs_raw %>% filter(roi %in% rois & to %in% localRef$bodyid)
  new_outputsR <- nBag$outputs_raw %>% filter(roi %in% rois & to %in% localRef$bodyid)
  new_names <- nBag$names %>% filter(bodyid %in% localRef$bodyid)
  new_refs <- nBag$outputsTableRef %>% filter(bodyid %in% localRef$bodyid)
  new_inputs <- getTypeToTypeTable(new_inputsR,typesTable = localRef,...)
  new_outputs <- getTypeToTypeTable(new_outputsR,typesTable = localRef,...)
  neuronBag(outputs = new_outputs,
            inputs = new_inputs,
            names = new_names,
            inputs_raw = new_inputsR,
            outputs_raw = new_outputsR,
            outputsTableRef = new_refs) 
  
}

combineRois <- function(connections,rois,newRoi,...){UseMethod("combineRois")}

combineRois.data.frame <- function(connections,rois,newRoi){
  ## CHECK IT'S A RAW CONNECTION TABLE
  newRegionTable <- connections %>% 
          filter(roi %in% rois) %>% 
          group_by(to) %>%
          mutate(totalROIweight=sum(totalROIweight[match(rois,roi)],na.rm=TRUE)) %>%
          group_by(from) %>%
          mutate(totalPreROIweight=sum(totalPreROIweight[match(rois,roi)],na.rm=TRUE)) %>%
          group_by_if(names(.) %in% c(paste0(c("","name.","type.","databaseType.","previous.type."),"to"),paste0(c("","name.","type.","databaseType.","previous.type."),"from"),paste0("supertype.to",1:3),paste0("supertype.from",1:3))) %>%
                                    summarize(roi=newRoi,
                                              weight=weight[1],
                                              totalPreWeight=totalPreWeight[1],
                                              weightRelativeTotal=weightRelativeTotal[1],
                                              outputContributionTotal=outputContributionTotal[1],
                                              totalROIweight = totalROIweight[1],
                                              totalPreROIweight=totalPreROIweight[1],
                                              ROIweight=sum(ROIweight),
                                              weightROIRelativeTotal=sum(weightROIRelativeTotal)) %>%
                                    ungroup() %>%
                                    mutate(weightRelative=ROIweight/totalROIweight,
                                           outputContribution=ROIweight/totalPreROIweight)
  
  newRegionTable
}

combineRois.neuronBag <- function(connections,rois,newRoi,...){
  new_inputsR <- combineRois(connections$inputs_raw,rois,newRoi)
  new_outputsR <- combineRois(connections$outputs_raw,rois,newRoi)
  new_inputs <- getTypeToTypeTable(new_inputsR,typesTable = connections$names,...)
  new_outputs <- getTypeToTypeTable(new_outputsR,typesTable = connections$outputsTableRef,...)
  neuronBag(outputs = new_outputs,
            inputs = new_inputs,
            names = connections$names,
            inputs_raw = new_inputsR,
            outputs_raw = new_outputsR,
            outputsTableRef = connections$outputsTableRef) 
            
}

getROISummary <- function(InOutList,filter=TRUE,rois = NULL){
  #' Build a pre roi summary of innervation for neurons in a neuronBag
  #' @param InOutList : a neuronBag object
  #' @param filter : if TRUE (the default), only return results in ROIs where significant type to 
  #' type connections are found. Otherwise consider all connections
  #' @param rois : a roiset to consider (if NULL consider all rois)
  #' 
  ROIOutputs <- InOutList$outputs_raw %>% group_by(roi,type.from,databaseType.from)   %>%
    summarize(OutputWeight = sum(weight)) %>% rename(type = type.from,databaseType=databaseType.from) %>% ungroup()
  
  ROIInputs <- InOutList$inputs_raw %>% group_by(roi,type.to,databaseType.to)   %>%
    summarize(InputWeight = sum(weight))  %>% rename(type = type.to,databaseType=databaseType.to) %>% ungroup()
  
  if (filter){
  ROIOutputs <- ROIOutputs %>% 
    filter(paste0(roi,type) %in% paste0(InOutList$outputs$roi,InOutList$outputs$type.from)) 
  ROIInputs <- ROIInputs %>%
    filter(paste0(roi,type) %in% paste0(InOutList$inputs$roi,InOutList$inputs$type.to))
  }
  
  if (!(is.null(rois))){
    ROIOutputs <- ROIOutputs %>% 
      filter(roi %in% rois$roi)
    ROIInputs <- ROIInputs %>%
      filter(roi %in% rois$roi)
  }
  
  roiSummary <- 
    full_join(ROIInputs,ROIOutputs,by=c("roi","type","databaseType")) %>% replace_na(list(InputWeight=0,OutputWeight=0)) %>%
    mutate(fullWeight = OutputWeight+InputWeight,
           deltaWeight = (OutputWeight - InputWeight)/fullWeight,
           supertype1 = supertype(databaseType,level=1),
           supertype2 = supertype(databaseType,level=2),
           supertype3 = supertype(databaseType,level=3))
  
  if (!is.null(rois)){roiSummary <- left_join(roiSummary,rois,by=roi)}
  return(roiSummary)
}

compressROISummary <- function(roiSummary,stat=median,level=1){
  #' Summarise a roi summary by supertype
  #' @param roiSummary : a table as returned by \code{getROISummary}
  #' @param stat : what statistic to use to summarise accross supertypes?
  #' @param level : what level of supertype to group by (default 1)
  #' 
  
  spt <- paste0("supertype",level)
  roiSummary %>% group_by(roi,!!(as.name(spt))) %>%
    summarise(type = (!!(as.name(spt)))[1],
              fullWeight = stat(fullWeight),
              deltaWeight = mean(deltaWeight),
              supertype2 = supertype2[1],
              supertype3 = supertype3[1],
              databaseType=databaseType[1]
              ) %>% ungroup()
}

haneschPlot <- function(roiTable,
                        roiSelect=selectRoiSet(getRoiTree()),
                        grouping=NULL,flip=FALSE,
                        alphaG=1,
                        alphaRois=0.15,
                        roiLabel=selectRoiSet(getRoiTree(),default_level = 0),
                        regionOutlines=T,
                        theme=theme_minimal(),
                        interactive=FALSE){
  roiTable <- roiTable %>% filter(roi %in% unique(roiSelect$roi))  %>% 
                      mutate(roi = factor(roi,levels=levels(roiSelect$roi)),
                             l4 = roiSelect$level4[match(roi,roiSelect$roi)],
                             side = roiSelect$side2[match(roi,roiSelect$roi)],
                             superroi = roiLabel$roi[match(l4,roiLabel$level4)]) %>%
                             arrange(roi) %>%
                             mutate(roiX = match(roi,unique(roi)))
  
  roiPos <- roiTable %>% group_by(superroi,side) %>%
                         summarize(xmin=min(roiX)-0.45,xmax=max(roiX)+0.45) %>% 
                         ungroup() #%>%
                         #filter(xmax-xmin>1)
  
  hanesch <- ggplot(data=roiTable,aes(x=roi,y=type))
  roiP <- roisPalette()
  
  hanesch <- hanesch +
    geom_line(aes(group=type),alpha=alphaG) 
  if (regionOutlines==TRUE){hanesch <- hanesch +
    geom_rect(data=roiPos,aes(xmin=xmin,xmax=xmax,ymin=-Inf,ymax=Inf,fill=superroi),alpha=alphaRois,inherit.aes = F) + 
              scale_fill_manual(name="Brain region",values=roiP,guide = guide_legend(reverse = TRUE)) +
              new_scale_fill()}
  if (interactive){
    hanesch <- hanesch + geom_point_interactive(data=roiTable,aes(size=fullWeight,fill=deltaWeight,x=roi,y=type,tooltip=paste0(type," in ",roi,"\nOutputs: ",OutputWeight,"\nInputs: ",InputWeight)),shape=21,alpha=alphaG)}
  else{
  hanesch <- hanesch + 
    geom_point(data=roiTable,aes(size=fullWeight,fill=deltaWeight,x=roi,y=type),shape=21,alpha=alphaG)}
  hanesch <- hanesch +
    scale_fill_gradient(name="Polarity",breaks=c(-1,-0.5,0,0.5,1),labels=c("Receives inputs","","Mixed","","Sends outputs"),low = "white", high = "black",
                        space = "Lab") +
    guides(fill = guide_legend(override.aes = list(size=5))) +
    scale_size_continuous(name = "# Synapses") + labs(y="Neuron type",x="Neuropile") + theme 
  
  if (!(is.null(grouping))){
    if (flip==TRUE){fct <- paste(". ~",grouping)}else{fct <- paste(grouping,"~ .")}
    hanesch <- hanesch + facet_grid(as.formula(fct),scale="free",space="free") 
  }
  
  if (flip==TRUE){hanesch <- hanesch + coord_flip()}
  hanesch + theme(axis.text.x = element_text(angle = 90))
  
}
  
  