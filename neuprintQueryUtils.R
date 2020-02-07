########################################################################################################################
# Convenience functions for query results processing
########################################################################################################################

### Neuprint search
getBodyIdsForList = function (neuronList,prefix="",postfix=".*",...){
  #' Get one dataframe of bodyIDs for all search strings in neuronList
  #' @param neuronList: A list of search strings to be passed.
  #' @param addStart: Should a '.*' be added at the end of the query strings?
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
  
  return( myConnections )
  
}

## Change from a bodyid/partner/prepost to a from/to format
simplifyConnectionTable <- function(connectionTable){
  
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

## Filter and generate type to type connection
getTypeToTypeTable <- function(connectionTable){
  
  ## Grouping summarizing filtering. Problem: depends on connectionTable polarity. Write intermediate functon?
  
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
