source("neuprintQueryUtils.R")

buildInputsOutputsByType <- function(typeQuery,fixed=FALSE,...){
  TypeNames <- neuprint_search(typeQuery,field="type",fixed=fixed)
  outputs <- getConnectionTable(TypeNames,synapseType = "POST",by.roi=TRUE,...)
  inputs <- getConnectionTable(TypeNames,synapseType = "PRE",by.roi=TRUE,...)
  OUTByTypes <- getTypeToTypeTable(outputs)
  INByTypes <- getTypeToTypeTable(inputs)
  return(list(outputs = OUTByTypes,
              inputs = INByTypes,
              names = TypeNames,
              outputs_raw = outputs,
              inputs_raw = inputs
              ))
}

redefineTypeByNameInList <- function(IOList,
                                     typeList,
                                     pattern,
                                     newPostFixes,
                                     perl=FALSE){

  
  
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
  }
  
  outputNames <- unique(IOList$outputs$databaseTypeTo)
  outputsLatNames <- getTypesTable(outputNames)
  for (t in typeList){
    outputsLatNames <- redefineTypeByName(outputsLatNames,
                                          type=t,
                                          pattern=pattern,
                                          newPostFixes=newPostFixes,
                                          type_col="type",
                                          name_col="name",
                                          perl=perl)
  }
     
  IOList$outputs <- getTypeToTypeTable(IOList$outputs_raw,typesTable = outputsLatNames)
  IOList$inputs <- getTypeToTypeTable(IOList$inputs_raw,typesTable = IOList$names)

  
  return(IOList)
}

lateralizeInputOutputList <- function(inputOutputList,typeList=NULL){
 
  outputsLat <- lrSplit(inputOutputList$outputs_raw,nameCol = "name.from",typeCol = "type.from",typeList=typeList)
  outputsLat <- lrSplit(outputsLat,typeList=typeList)
  
  inputsLat <- lrSplit(inputOutputList$inputs_raw,nameCol = "name.from",typeCol = "type.from",typeList=typeList)
  inputsLat <- lrSplit(inputsLat,typeList=typeList)
 
  outputsLatNames <- getTypesTable(unique(outputsLat$databaseTypeTo))
  outputsLatNames <- lrSplit(outputsLatNames,nameCol="name",typeCol="type",typeList=typeList)
  TypeNamesLat <- lrSplit(inputOutputList$names,nameCol = "name",typeCol="type",typeList=typeList)
 
  outByTypesLat <- getTypeToTypeTable(outputsLat,typesTable = outputsLatNames)
  inByTypesLat <- getTypeToTypeTable(inputsLat,typesTable = TypeNamesLat)

  return(list(outputs = outByTypesLat,
              inputs = inByTypesLat,
              names = TypeNamesLat,
              outputs_raw = outputsLat,
              inputs_raw=inputsLat))
  
}

bind_InoutLists <- function(...){
  full <- list(...)
  out <- list(outputs = distinct(bind_rows(lapply(full,function(i) i$outputs))),
              inputs = distinct(bind_rows(lapply(full,function(i) i$inputs))),
              names = distinct(bind_rows(lapply(full,function(i) i$names))),
              outputs_raw = distinct(bind_rows(lapply(full,function(i) i$outputs_raw))),
              inputs_raw = distinct(bind_rows(lapply(full,function(i) i$inputs_raw))))
  return(out)
  
}