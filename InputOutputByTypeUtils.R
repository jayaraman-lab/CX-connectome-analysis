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

lateralizeInputOutputList <- function(inputOutputList){
 
  outputsLat <- lrSplit(inputOutputList$outputs_raw,nameCol = "name.from",typeCol = "type.from")
  outputsLat <- lrSplit(outputsLat)
  
  inputsLat <- lrSplit(inputOutputList$inputs_raw,nameCol = "name.from",typeCol = "type.from")
  inputsLat <- lrSplit(inputsLat)
 
  outputsLatNames <- getTypesTable(unique(outputsLat$databaseTypeTo))
  outputsLatNames <- lrSplit(outputsLatNames,nameCol="name",typeCol="type")
  TypeNamesLat <- lrSplit(inputOutputList$names,nameCol = "name",typeCol="type")
 
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