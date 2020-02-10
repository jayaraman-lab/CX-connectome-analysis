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
getConnectionTable <- function(bodyIDs,synapseType, slctROI,by.roi,...){
  #' Get connection table of inputs and add meta
  #' @return Returns a connection table as data frame. Added columns are \code{weightRelativeTotal} which is 
  #' the relative weight considering all the synapses (irrespective of the ROI), and if ROI are used (either if
  #' \code{slctROI} has a value or \code{by.roi} is \code{TRUE}), \code{weightRelative} is the relative weight in the 
  #' ROI and \code{totalROIweight} is the absolute number of inputs this neuron receives in that region and 
  #' \code{weightROIRelativeTotal} is the weight in the ROI normalized by the total number of inputs (in all ROIs)
  #' @param bodyIDs: The bodyids of neurons who's connections should be queried or a metadata data.frame
  #' @param synapseType: Choose "PRE" or "POST" to get inputs or outputs of the neurons in bodyIDs, respectivly.
  #' @param slctROI: String specifying the ROI where connections should be queried. By default all the ROIs.
  #' @param by.roi: Passed to neuprint_connection_table. If returning all ROIs, should results be broken down by ROI?
  #' @param ...: Other arguments to be passed to neuprint_connection_table
  
  UseMethod("getConnectionTable")}


getConnectionTable.default = function(bodyIDs, synapseType, slctROI=NULL,by.roi=FALSE,...){
  refMeta <- neuprint_get_meta(bodyIDs)
  return(getConnectionTable(refMeta,synapseType,slctROI,by.roi,...))
}

getConnectionTable.data.frame <- function(bodyIDs,synapseType, slctROI=NULL,by.roi=FALSE,...){
  refMeta <- bodyIDs
  bodyIDs <- bodyIDs$bodyid
  myConnections <- neuprint_connection_table(bodyIDs, synapseType, slctROI,by.roi=by.roi,...)
  partnerMeta <- neuprint_get_meta(myConnections$partner)
  refMeta <- slice(refMeta,sapply(myConnections$bodyid,function(b) match(b,refMeta$bodyid)))
  
  myConnections <-myConnections %>%
    mutate(partnerName = partnerMeta[["name"]],
           name = refMeta[["name"]],
           partnerType = partnerMeta[["type"]],
           type = refMeta[["type"]])
  
  ## Normalization is always from the perspective of the output (fraction of inputs to the output neuron)
  if (synapseType == "PRE"){
    outMeta <- refMeta
    } else {
    outMeta <- partnerMeta
  }
  myConnections[["weightRelativeTotal"]] <- myConnections[["weight"]]/outMeta[["post"]]
  
  if (by.roi | !is.null(slctROI)){
    myConnections[["weightROIRelativeTotal"]] <- myConnections[["ROIweight"]]/outMeta[["post"]]
    
    if (synapseType == "PRE"){
      outInfo <- neuprint_get_roiInfo(myConnections$bodyid)
      inInfo <- neuprint_get_roiInfo(myConnections$partner)
    }else{
      outInfo <- neuprint_get_roiInfo(myConnections$partner)
      inInfo <- neuprint_get_roiInfo(myConnections$bodyid)
    }
    
    postVar <- paste0(myConnections[["roi"]],".post")
    preVar <- paste0(myConnections[["roi"]],".pre")
    myConnections <- myConnections %>%
      mutate(totalROIweight = sapply(1:length(postVar),function(v) outInfo[[postVar[v]]][v]),
             totalPreROIweight = sapply(1:length(preVar),function(v) inInfo[[preVar[v]]][v])) %>%
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
  return(bind_rows(lapply(types,function(t) neuprint_search(t,field="type"))))
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

## Filter and generate type to type connection
getTypeToTypeTable <- function(connectionTable,
                               majorOutputThreshold=0.8,
                               singleNeuronThreshold=0.01,
                               pThresh = 0.05,
                               typesTable = NULL){
  #' Generate a table of type to type connections, keeping only the significant links
  #' @param table: A table of neuron to neuron connections (as generated by \code{getConnectionTable})
  #' @param majorOutputThreshold: Threshold of the fraction of the outputs of a presynaptic type 
  #' to be accounted for by an output type for us to keep this connection regardless of other 
  #' criterions. Should be close to 1, the rationale being if a neuron has only one type of partners
  #' in the region, one ought to consider it significant 
  #' @param singleNeuronThreshold: If a neuron is the only representent of its type, what fraction 
  #' of it's input should a presynaptic type represent for the connection to be kept?
  #' @param pThresh: Significance level of the t-test for a connection to be kept
  #' @param typesTable: A table of all the instances of the output types (as generated by a search 
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
  #' @examples
  #' \dontrun{
  #' ## Getting a PFL outputs in the LAL table, and summarize it by type
  #' PFLs <- getTypesTable(c("PFL1","PFL2","PFL3"))
  #' PFLConnections <- getConnectionTable(PFLs,"POST","LAL(-GA)(R)")
  #' PFLConnectionByType <- getTypeToTypeTable(PFLConnections)
  #' 
  #'
  #' ## Splitting one output type according to left/right in a PFL table
  #' ## For ease:
  #' PFLConnections <- simplifyConnectionTable(PFLConnections)
  #' PFLConnections_renamed <- redefineType(table=PFLConnections,
  #'                                        type="AVL01op_pct",
  #'                                        condition=grepl("_L",PFLConnections$name.to),
  #'                                        newTypes=c("AVL01op_pct_L","AVL01op_pct_R"),
  #'                                        type_col = "type.to")
  #' ## One now needs to create an ad-hoc table and rename it similarly
  #' outputTypes <- getTypesTable(unique(PFLConnections)[["type.to"]])
  #' outputTypes_renamed <- redefineType(table=OutputTypes,
  #'                             type="AVL01op_pct",
  #'                             condition=grepl("_L",OutputTypes$name),
  #'                             newTypes=c("AVL01op_pct_L","AVL01op_pct_R"),
  #'                             type_col = "type")
  #'## One can then use the function with the table
  #' PFLTypeToType <- getTypeToTypeTable(PFLConnections_renamed,typesTable = outputTypes_renamed)                       
  #'                             
  #' }
  #' 
  connectionTable <- simplifyConnectionTable(connectionTable)
 
  ## Counting instances for each post type 
  if (is.null(typesTable)){
    typesTable <- getTypesTable(unique(connectionTable$type.to))
  }
  
  typesCount <- typesTable %>% group_by(type) %>%
                               summarise(n=n())
  
  connectionTable <- connectionTable %>% 
                      mutate(n = typesCount[["n"]][match(type.to,typesCount[["type"]])])
 
  ## Renaming the unnamed neurons and treating them as single examples
  connectionTable <- connectionTable %>% mutate(name.to = replace_na(name.to,"Unlabeled"),
                                                type.to = replace_na(type.to,"Unlabeled"),
                                                n = replace_na(n,1))
  
  unknowns <- which(connectionTable$type.to == "Unlabeled")
  idx <- sapply(connectionTable$to[unknowns], function(p) which(unique(connectionTable$to[unknowns])==p))
  connectionTable[["type.to"]][unknowns] <- paste0("Unlabeled_",idx)
  
  unknownNames <- which(connectionTable$name.to == "Unlabeled")
  connectionTable[["name.to"]][unknownNames] <- connectionTable[["type.to"]][unknownNames]
  
  
  ## Gather the outputContributions
  connectionTable <-  connectionTable %>% group_by(from,type.to) %>%
                                          mutate(outputContribution = sum(outputContribution)) %>%
                                          group_by(type.from,type.to) %>%
                                          mutate(outputContribution = mean(outputContribution))
  
  ## This contains the neurons unique in their type that reach our hard threshold
  loners <- connectionTable %>% filter(n==1) %>%
                                group_by(type.from,type.to) %>%
                                summarize(weightRelative = sum(weightRelative),
                                          outputContribution = outputContribution[1],
                                          n_type = 1,
                                          n_links = n()) %>%
                                filter(weightRelative > singleNeuronThreshold | outputContribution > majorOutputThreshold)
  
  ## Main filter
  sTable <- connectionTable %>% filter(n>1) %>%
                                group_by(type.from,to,type.to) %>%
                                summarise(weightRelative = sum(weightRelative),
                                          n = n[1],
                                          outputContribution = outputContribution[1]) %>%
                                group_by(type.from,type.to) %>%
                                summarize(pVal = t.test(c(weightRelative,unlist(replicate(n[1]-n(),0))),
                                                      alternative="greater")["p.value"],
                                          ci_low = t.test(c(weightRelative,unlist(replicate(n[1]-n(),0))),
                                                        alternative="greater")[["conf.int"]][1],
                                          varWeight = var(c(weightRelative,unlist(replicate(n[1]-n(),0)))),
                                          weightRelative = mean(c(weightRelative,unlist(replicate(n[1]-n(),0)))),
                                          outputContribution = outputContribution[1],
                                          n_links = n(),
                                          n_type = n[1]
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
  
  myConnections = myConnections %>% filter(partner %in% postBodyIDs) 
  
  return( myConnections )
}
