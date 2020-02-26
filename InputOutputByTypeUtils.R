source("neuprintQueryUtils.R")
source("supertypeUtils.R")
library(pbapply)
library(parallel)

buildInputsOutputsByType <- function(typeQuery,fixed=FALSE,big=FALSE,nc=5,...){
  UseMethod("buildInputsOutputsByType")}

buildInputsOutputsByType.character <- function(typeQuery,fixed=FALSE,big=FALSE,nc=5,...){
  TypeNames <- distinct(bind_rows(lapply(typeQuery,neuprint_search,field="type",fixed=fixed))) %>%
                  mutate(databaseType = type)
  buildInputsOutputsByType(TypeNames,fixed=FALSE,big=big,nc=nc,...)
}
  
buildInputsOutputsByType.data.frame <- function(typeQuery,fixed=FALSE,selfRef=FALSE,big=FALSE,nc=5,...){
  if (big == TRUE){
    inoutList <- pblapply(unique(typeQuery$type),
                          function(t) buildInputsOutputsByType(typeQuery %>% filter(type == t),selfRef=selfRef,big=FALSE),cl = nc)
                          
    return(do.call(bind_InoutLists,inoutList))
  }
  
  outputsR <- getConnectionTable(typeQuery,synapseType = "POST",by.roi=TRUE,...)
  inputsR <- getConnectionTable(typeQuery,synapseType = "PRE",by.roi=TRUE,...)
  if (nrow(outputsR)==0){OUTByTypes <- NULL
                         outputsTableRef <- NULL
                         unknowns <- NULL
                         }else{
    OUTByTypes <- getTypeToTypeTable(outputsR)
    outputsR <- retype.na(outputsR)
    outputsTableRef <- getTypesTable(unique(outputsR$type.to))
    unknowns <- retype.na_meta(neuprint_get_meta(outputsR$to[!(outputsR$to %in% outputsTableRef$bodyid)]))
    }
  if (nrow(inputsR)==0){INByTypes <- NULL}else{
    if (selfRef){
      INByTypes <- getTypeToTypeTable(inputsR,typesTable = typeQuery)
    }else{
      INByTypes <- getTypeToTypeTable(inputsR)
    }
    inputsR <- retype.na(inputsR)}
  
  return(list(outputs = OUTByTypes,
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
  
  
  return(list(outputs = outByTypesLat,
              inputs = inByTypesLat,
              names = TypeNamesLat,
              outputs_raw = outputsLat,
              inputs_raw=inputsLat,
              outputsTableRef=outputsRef))
  
}

bind_InoutLists <- function(...){
  full <- list(...)
  out <- list(outputs = distinct(bind_rows(lapply(full,function(i) i$outputs))),
              inputs = distinct(bind_rows(lapply(full,function(i) i$inputs))),
              names = distinct(bind_rows(lapply(full,function(i) i$names))),
              outputs_raw = distinct(bind_rows(lapply(full,function(i) i$outputs_raw))),
              inputs_raw = distinct(bind_rows(lapply(full,function(i) i$inputs_raw))),
              outputsTableRef = distinct(bind_rows(lapply(full,function(i) i$outputsTableRef)))
              )
  return(out)
}


getROISummary <- function(InOutList,filter=TRUE){
  ROIOutputs <- InOutList$outputs_raw %>% group_by(roi,type.from,databaseTypeFrom)   %>%
    summarize(OutputWeight = sum(weight)) %>% rename(type = type.from,databaseType=databaseTypeFrom)
  
  ROIInputs <- InOutList$inputs_raw %>% group_by(roi,type.to,databaseTypeTo)   %>%
    summarize(InputWeight = sum(weight))  %>% rename(type = type.to,databaseType=databaseTypeTo)
  
  if (filter){
  ROIOutputs <- ROIOutputs %>% 
    filter(paste0(roi,type) %in% paste0(InOutList$outputs$roi,InOutList$outputs$type.from)) 
  ROIInputs <- ROIInputs %>%
    filter(paste0(roi,type) %in% paste0(InOutList$inputs$roi,InOutList$inputs$type.to))
  }
  
  roiSummary <- 
    full_join(ROIInputs,ROIOutputs,by=c("roi","type","databaseType")) %>% replace_na(list(InputWeight=0,OutputWeight=0)) %>%
    mutate(fullWeight = OutputWeight+InputWeight,
           deltaWeight = (OutputWeight - InputWeight)/fullWeight,
           supertype = supertype(type),
           megatype = megatype(supertype))
  
  return(roiSummary)
}

haneschPlot <- function(roiTable,roiSelect=unique(roiTable(roi)),by.supertype=F){
  roiTable <- roiTable %>% filter(roi %in% roiSelect)
  
  hanesch <- ggplot(roiTable,aes(x=roi,y=type)) + 
    geom_line(aes(group=type)) +
    geom_point(aes(size=fullWeight,fill=deltaWeight),colour="black",shape=21)+
    scale_fill_gradient(name="Polarity",breaks=c(-1,-0.5,0,0.5,1),labels=c("Receives inputs","","Mixed","","Sends outputs"),low = "white", high = "black",
                        space = "Lab", na.value = "grey50", guide = "legend",
                        aesthetics = "fill") +
    guides(fill = guide_legend(override.aes = list(size=5))) +
    scale_size_continuous(name = "# Synapses") +
    theme_minimal()
  if (by.supertype){
    hanesch <- hanesch + facet_grid(supertype~.,scale="free_y",space="free_y") + theme_gray()
  }
  hanesch + theme(axis.text.x = element_text(angle = 90)) 
  
}
  
  