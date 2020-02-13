########################################################################################################################
# Convenience functions for query results processing
########################################################################################################################

### Neuprint search
getBodyIdsForList = function (neuronList,prefix="",postfix=".*",...){
  #' Get one dataframe of bodyIDs for all search strings in neuronList
  #' @param neuronList: A list of search strings to be passed.
  #' @param prefix: String to be added before each query (default "")
  #' @param postfix: String to be added after each query (default ".*")
  #' @param ...: Parameters to be passed to neuprint_search. Note that meta=FALSE won't work for now.
  #' @return A data frame of metadata for the results of all the queries
  #' @examples
  #' \dontrun{
  #' # Both will return the same
  #' getBodyIdsForList(c("PFL1","PFL2"))
  #' getBodyIdsForList(c("PFL1","PFL2"),postfix="",field="type")
  #' }
  
  neuronList <-  paste0(prefix,neuronList,postfix)
  bodiesList <- lapply(neuronList,neuprint_search,...)
  return(bind_rows(bodiesList))
}

### Connection table
getConnectionTable <- function(bodyIDs,synapseType, slctROI,by.roi, synThresh = 3,...){
  #' Get connection table of inputs and add meta
  #' @return Returns a connection table as data frame. Added columns are \code{weightRelativeTotal} which is 
  #' the relative weight considering all the synapses (irrespective of the ROI), and if ROI are used (either if
  #' \code{slctROI} has a value or \code{by.roi} is \code{TRUE}), \code{weightRelative} is the relative weight in the 
  #' ROI and \code{totalROIweight} is the absolute number of inputs this neuron receives in that region and 
  #' \code{weightROIRelativeTotal} is the weight in the ROI normalized by the total number of inputs (in all ROIs)
  #' @param bodyIDs: The bodyids of neurons who's connections should be queried or a metadata data.frame
  #' @param synapseType: Choose "PRE" or "POST" to get inputs or outputs of the neurons in bodyIDs, respectivly. "PRE" is 
  #' usually slower as it requires listing all the outputs of the input neurons to get the contribution of the outputs.
  #' @param slctROI: String specifying the ROI where connections should be queried. By default all the ROIs.
  #' @param by.roi: Passed to neuprint_connection_table. If returning all ROIs, should results be broken down by ROI?
  #' @param synThresh: Minimum number of synapses to consider a connection (default 3)
  #' @param ...: Other arguments to be passed to neuprint_connection_table
  
  UseMethod("getConnectionTable")}


getConnectionTable.default = function(bodyIDs, synapseType, slctROI=NULL,by.roi=FALSE, synThresh=3,...){
  refMeta <- neuprint_get_meta(bodyIDs)
  return(getConnectionTable(refMeta,synapseType,slctROI,by.roi,synThresh,...))
}

getConnectionTable.data.frame <- function(bodyIDs,synapseType, slctROI=NULL,by.roi=FALSE,synThresh=3,...){
  refMeta <- bodyIDs
  bodyIDs <- bodyIDs$bodyid
  myConnections <- neuprint_connection_table(bodyIDs, synapseType, slctROI,by.roi=by.roi,...)
  myConnections <- myConnections %>% drop_na(ROIweight) %>% filter(ROIweight>synThresh)
  partnerMeta <- neuprint_get_meta(myConnections$partner)
  refMetaOrig <- neuprint_get_meta(myConnections$bodyid)  ## To get the database type name
  
  myConnections <- filter(myConnections,partnerMeta$status == "Traced")
  partnerMeta <- filter(partnerMeta,status == "Traced")
  
  refMeta <- slice(refMeta,sapply(myConnections$bodyid,function(b) match(b,refMeta$bodyid)))
  refMetaOrig <- slice(refMetaOrig,sapply(myConnections$bodyid,function(b) match(b,refMetaOrig$bodyid)))
  
  myConnections <-myConnections %>%
    mutate(partnerName = partnerMeta[["name"]],
           name = refMeta[["name"]],
           partnerType = partnerMeta[["type"]],
           type = refMeta[["type"]])
  
  myConnections <- simplifyConnectionTable(myConnections)
  ## Normalization is always from the perspective of the output (fraction of inputs to the output neuron)
  if (synapseType == "PRE"){
    outMeta <- refMeta
    myConnections <- myConnections %>% mutate(databaseTypeTo = refMetaOrig$type,
                                              databaseTypeFrom = type.from)
    } else {
    outMeta <- partnerMeta
    myConnections <- myConnections %>% mutate(databaseTypeTo = type.to,
                                              databaseTypeFrom = refMetaOrig$type)
  }
  myConnections[["weightRelativeTotal"]] <- myConnections[["weight"]]/outMeta[["post"]]
  
  if (by.roi | !is.null(slctROI)){
    myConnections[["weightROIRelativeTotal"]] <- myConnections[["ROIweight"]]/outMeta[["post"]]
    outInfo <- neuprint_get_roiInfo(myConnections$to)
    if (synapseType == "PRE"){
      inputsTable <- neuprint_connection_table(unique(myConnections$from),"POST",slctROI,by.roi=by.roi,...)
      inputsTable <- inputsTable %>% drop_na(ROIweight) %>% mutate(from = bodyid)
      inp <- "partner"
    }else{
      inputsTable <- myConnections
      inp <- "bodyid"
    }
    
    totalPre <- inputsTable %>% group_by(from) %>%
                    summarise(totalPreROIweight = sum(ROIweight))
    
    postVar <- paste0(myConnections[["roi"]],".post")
  
    myConnections <- myConnections %>%
      mutate(totalROIweight = sapply(1:length(postVar),function(v) outInfo[[postVar[v]]][v]),
             totalPreROIweight = totalPre[["totalPreROIweight"]][match(myConnections$from,totalPre$from)]) %>%
      mutate(weightRelative=ROIweight/totalROIweight,
             outputContribution=ROIweight/totalPreROIweight) ## This is how much this connection accounts for the outputs of the input neuron (not the standard measure)
  }
  
  return( myConnections %>% drop_na(weightRelative) ) ## NA values can occur in rare cases where
                                                      ## synapse (pre/post) is split between ROIs
  
}

simplifyConnectionTable <- function(connectionTable){
  #' Change from a bodyid/partner/prepost to a from/to format
  #' @param connectionTable: A data frame in the bodyid/partner/prepost format (as returned by getConnectionTable)
  #' @return A data frame in the from/to/name.from/name.to... format
  #' 
  #' 
  if ("from" %in% names(connectionTable)){return(connectionTable)}else{
  connectionTable <- connectionTable %>% mutate(from = ifelse(prepost==1,!!as.name("bodyid"),!!as.name("partner")),
                                                to = ifelse(prepost==1,!!as.name("partner"),!!as.name("bodyid")),
                                                name.from = ifelse(prepost==1,!!as.name("name"),!!as.name("partnerName")),
                                                name.to = ifelse(prepost==1,!!as.name("partnerName"),!!as.name("name")),
                                                type.from = ifelse(prepost==1,!!as.name("type"),!!as.name("partnerType")),
                                                type.to = ifelse(prepost==1,!!as.name("partnerType"),!!as.name("type"))
                                                ) %>%
                                         select(-bodyid,-partner,-name,-partnerName,-partnerType,-type,-prepost)
  return(connectionTable)
  }
}

## Returns a table with all neurons belonging to types
getTypesTable <- function(types){
  #' Get a table of all instances of certain types
  #' @param types: A vector of type names
  #' @return A data frame of instances of those types
  #' 
  #' 
  return(bind_rows(lapply(types,function(t) neuprint_search(t,field="type",fixed=TRUE))))
}

redefineType <- function(table,type,condition,newTypes,type_col="type"){
  #' Split a type in 2 categories according to a condition
  #' @param table: A table with a column of types to be renamed
  #' @param type: The name of the type to split
  #' @param condition: A logical vector indicating which rows fulfill the split factor
  #' @param newTypes: A vector of strings of length 2: first value is the new type name if 
  #' condition is TRUE, 2d value the other case
  #' @param type_col: which column in the dataframe we want to consider
  #' @return A data frame with the type_col column renamed accordingly
  #' @examples
  #' \dontrun{
  #' ## Splitting PFL2s according to left/right in a PFL table
  #' PFLs <- getTypesTable(c("PFL1","PFL2","PFL3"))
  #' PFL_renamed <- redefineType(table=PFLs,
  #'                             type="PFL2",
  #'                             condition=grepl("_L",PFLs$name),
  #'                             newTypes=c("PFL2_L","PFL2_R"))
  #' }
  #' 
  table[[type_col]][table[[type_col]] == type] <-  newTypes[2]
  table[[type_col]][table[[type_col]] == newTypes[2] & condition] <-  newTypes[1] 
  return(table)
}

lrSplit <- function(connectionTable,
                    nameCol="name.to",
                    typeCol="type.to",
                    typeList=NULL){
  #' Retype neurons in a table according to L/R
  #' @param connectionTable : connectivity table to modify
  #' @param nameCol : the name of the column containing the names of the neurons (to be used
  #' to determine laterality)
  #' @param typeCol : the name of the column containing the types to modify
  #' @param typeList : which types to lateralize (by default all the neurons which names
  #' contains L or R)
  #' @example 
  #' \dontrun{
  #' PFLNames <- getTypesTable(c("PFL1","PFL2","PFL3"))
  #' 
  #' ## Rename only PFL2
  #' PFLNames2 <- lrSplit(PFLNames,nameCol="name",typeCol="type",typeList=c("PFL2"))
  #' 
  #' ##Rename all PFLs
  #' PFLNames3 <- lrSplit(PFLNames,nameCol="name",typeCol="type")
  #' 
  #' } 
  if (is.null(typeList)){
    typeList <- distinct(connectionTable,(!!as.name(nameCol)),(!!as.name(typeCol)))
    typeList <- filter(typeList,grepl("_R|_L",(!!as.name(nameCol)))) %>% na.omit()
    typeList <- typeList[[typeCol]]
  }
  for (t in typeList){
    connectionTable <- redefineType(table=connectionTable,
                                    type=t,
                                    condition=grepl("_L",connectionTable[[nameCol]]),
                                    newTypes=c(paste0(t,"_L"),paste0(t,"_R")),
                                    type_col = typeCol)
  }
  return(connectionTable)
}

retype.na <- function(connectionTable){
  #' Fill in the type and name field in case they are NAs, using the name field if it exists
  #' (removing the _L/_R) or the neuron id. By default expects a table in to/from format.
  #'
  connectionTable <- connectionTable %>% 
                     mutate(name.from = ifelse(is.na(name.from),from,name.from),
                            name.to = ifelse(is.na(name.to),to,name.to),
                            type.from = ifelse(is.na(type.from),gsub("_L$|_R$","",name.from),type.from),
                            type.to = ifelse(is.na(type.to),gsub("_L$|_R$","",name.to),type.to)
                            )
          
  return(connectionTable)
}

getTypeToTypeTable <- function(connectionTable,
                               majorOutputThreshold=0.8,
                               singleNeuronThreshold=0.01,
                               singleNeuronThresholdN=5,
                               pThresh = 0.05,
                               typesTable = NULL){
  #' Generate a table of type to type connections, keeping only the significant links
  #' @param connectionTable: A table of neuron to neuron connections (as generated by \code{getConnectionTable})
  #' @param majorOutputThreshold: Threshold of the fraction of the outputs of a presynaptic type 
  #' to be accounted for by an output type for us to keep this connection regardless of other 
  #' criterions. Should be close to 1, the rationale being if a neuron has only one type of partners
  #' in the region, one ought to consider it significant 
  #' @param singleNeuronThreshold: If a neuron is the only representent of its type, what fraction 
  #' of it's input should a presynaptic type represent for the connection to be kept?
  #' @param singleNeuronThresholdN: If a neuron is the only representent of its type, how many synapses from
  #' the presynaptic type should it minimally receive?
  #' @param pThresh: Significance level of the t-test for a connection to be kept
  #' @param  typesTable: A table of all the instances of the output types (as generated by a search 	  
  #' for example). If NULL (the default), it will be computed from the unique types of post synaptic	 
  #' partners. Necessary to use if there are some user defined types
  #' @return A data frame with the columns:
  #'              - type.to
  #'              - type.from
  #'              - weightRelative : the mean over the outputs of the relative input contribution
  #'              of the input type to the output type
  #'              - varWeight : the variance over the outputs of the relative input contribution
  #'              of the input type to the output type
  #'              -ci_low : lower bound confidence interval of the relative weight
  #'              - outputContribution : what proportion of the outputs of type.from does this connection
  #'              accounts for
  #'              - n_links : how many individual neuron to neuron connections does this connection contain?
  #'              -n_type : number of instances of type.to 
  #'              - databaseTypeTo/databaseTypeFrom : the neuprint defined types that contain the 
  #'              type.to/type.from used here
  #' @examples
  #' \dontrun{
  #' ## Getting a PFL outputs in the LAL table, and summarize it by type
  #' PFLs <- getTypesTable(c("PFL1","PFL2","PFL3"))
  #' PFLConnections <- getConnectionTable(PFLs,"POST","LAL(-GA)(R)")
  #' PFLConnectionByType <- getTypeToTypeTable(PFLConnections)
  #' 
  #'
  #' ## Splitting one output type according to left/right in a PFL table
  #' PFLConnections_renamed <- redefineType(table=PFLConnections,
  #'                                        type="AVL01op_pct",
  #'                                        condition=grepl("_L",PFLConnections$name.to),
  #'                                        newTypes=c("AVL01op_pct_L","AVL01op_pct_R"),
  #'                                        type_col = "type.to")
  #' ## One now needs to create an ad-hoc table and rename it similarly
  #' outputTypes <- getTypesTable(unique(PFLConnections)[["type.to"]])
  #' outputTypes_renamed <- redefineType(table=outputTypes,
  #'                             type="AVL01op_pct",
  #'                             condition=grepl("_L",OutputTypes$name),
  #'                             newTypes=c("AVL01op_pct_L","AVL01op_pct_R"),
  #'                             type_col = "type")
  #'                             
  #'## This particular transforms can also be achieved with lrSplit
  #'PFLConnections_renamed <- lrSplit(PFLConnections,typeList=c("AVL01op_pct" ))
  #'outputTypes_renamed <- lrSplit(outputTypes,nameCol="name",typeCol="type",typeList=c("AVL01op_pct" )) 
  #'                            
  #'## One can then use the function with the table
  #' PFLTypeToType <- getTypeToTypeTable(PFLConnections_renamed,typesTable = outputTypes_renamed)                       
  #'                             
  #' }
  #' 
  
  ## Counting instances for each post type 
  if (is.null(typesTable) & any(connectionTable$type.to != connectionTable$databaseTypeTo,na.rm=T))
    {
      stop("Some types are custom defined. You need to provide a `typesTable` argument.")
    }
  
  if (is.null(typesTable)){
    typesTable <- getTypesTable(unique(connectionTable$databaseTypeTo))
  } 
  
  typesCount <- typesTable %>% group_by(type) %>%
                               summarise(n=n())
  
  connectionTable <- connectionTable %>% 
                      mutate(n = typesCount[["n"]][match(type.to,typesCount[["type"]])])
 
  ## Renaming the unnamed neurons and treating them as single examples
  connectionTable <- connectionTable %>% mutate(n = replace_na(n,1))
  connectionTable <- retype.na(connectionTable)
  
  ## Gather the outputContributions
  connectionTable <-  connectionTable %>% group_by(from,type.to) %>%
                                          mutate(outputContribution = sum(outputContribution)) %>%
                                          group_by(type.from,type.to) %>%
                                          mutate(outputContribution = mean(outputContribution))
  
  ## This contains the neurons unique in their type that reach our hard threshold
  loners <- connectionTable %>% filter(n==1) %>%
                                group_by(type.from,type.to) %>%
                                summarize(weightRelative = sum(weightRelative),
                                          weight = sum(ROIweight),
                                          outputContribution = outputContribution[1],
                                          n_type = 1,
                                          n_targets = n(),
                                          databaseTypeTo = databaseTypeTo[1],
                                          databaseTypeFrom = databaseTypeFrom[1]) %>%
                                filter((weightRelative > singleNeuronThreshold & weight > singleNeuronThresholdN)| outputContribution > majorOutputThreshold)
  
  ## Main filter
  sTable <- connectionTable %>% filter(n>1) %>%
                                group_by(type.from,to,type.to) %>%
                                summarise(weightRelative = sum(weightRelative),
                                          weight = sum(ROIweight),
                                          n = n[1],
                                          outputContribution = outputContribution[1],
                                          databaseTypeTo = databaseTypeTo[1],
                                          databaseTypeFrom = databaseTypeFrom[1]) %>%
                                group_by(type.from,type.to) %>%
                                summarize(pVal = ifelse((all(weightRelative == weightRelative[1]) & n()==n[1]),   ## t.test doesn't run if values are constant. Keep those.
                                                                  0,
                                                                  t.test(c(weightRelative,unlist(replicate(n[1]-n(),0))),
                                                                                    alternative="greater",exact=FALSE)[["p.value"]]),
                                          varWeight = var(c(weightRelative,unlist(replicate(n[1]-n(),0)))),
                                          weightRelative = mean(c(weightRelative,unlist(replicate(n[1]-n(),0)))),
                                          weight = mean(c(weight,unlist(replicate(n[1]-n(),0)))),
                                          outputContribution = outputContribution[1],
                                          n_targets = n(),
                                          n_type = n[1],
                                          databaseTypeTo = databaseTypeTo[1],
                                          databaseTypeFrom = databaseTypeFrom[1]
                                ) %>% filter(pVal < pThresh | outputContribution > majorOutputThreshold) %>%
                                select(-pVal)
                                
          return(bind_rows(sTable,loners))                  

}

getConnectionTable_forSubset = function(preBodyIDs,postBodyIDs, slctROI=NULL,...){
  #' Filter connection table to contain only connections from preBodyIDs to postBodyIDs
  #' @return Returns a connection table as data frame.
  #' @param preBodyIDs: The bodyids of presynaptic neurons who's connections should be queried
  #' @param postBodyIDs: The bodyids of postsynaptic neurons who's connections should be queried
  #' @param slctROI: String specifying the ROI where connections should be queried.
  
  myConnections = getConnectionTable(preBodyIDs, "POST", slctROI,...)
  #myConnections = getConnectionTable(postBodyIDs, "PRE", slctROI,...)
  
  myConnections = myConnections %>% filter(to %in% postBodyIDs) 
  
  return( myConnections )
}
