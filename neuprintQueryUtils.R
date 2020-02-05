########################################################################################################################
# Convenience functions for query results processing
########################################################################################################################

### Neuprint search
getBodyIdsForList = function (neuronList,addStar=TRUE,...){
  #' Get one dataframe of bodyIDs for all search strings in neuronList
  #' @param neuronList: A list of search strings to be passed.
  #' @param addStart: Should a '.*' be added at the end of the query strings?
  #' @param ...: Parameters to be passed to neuprint_search. Note that meta=FALSE won't work for now.
  #' @return A data frame of metadata for the results of all the queries
  #' @examples
  #' \dontrun{
  #' # Both will return the same
  #' getBodyIdsForList(c("PFL1","PFL2"))
  #' getBodyIdsForList(c("PFL1","PFL2"),addStar=FALSE,field="type")
  #' }
  
  if (addStar){ neuronList <-  paste0(neuronList,'.*') }
  bodiesList <- lapply(neuronList,neuprint_search,...)
  return(bind_rows(bodiesList))
}


### Connection table

getConnectionTable = function(bodyIDs, synapseType, slctROI=NULL,by.roi=FALSE,...){
  #' Get connection table of inputs and add meta
  #' @return Returns a connection table as data frame. Added columns are \code{weightRelativeTotal} which is 
  #' the relative weight considering all the synapses (irrespective of the ROI), and if ROI are used (either if
  #' \code{slctROI} has a value or \code{by.roi} is \code{TRUE}), \code{weightRelative} is the relative weight in the 
  #' ROI and \code{totalROIweight} is the absolute number of inputs this neuron receives in that region
  #' @param bodyIDs: The bodyids of neurons who's connections should be queried.
  #' @param synapseType: Choose "PRE" or "POST" to get inputs or outputs of the neurons in bodyIDs, respectivly.
  #' @param slctROI: String specifying the ROI where connections should be queried. By default all the ROIs.
  #' @param by.roi: Passed to neuprint_connection_table. If returning all ROIs, should results be broken down by ROI?
  #' @param ...: Other arguments to be passed to neuprint_connection_table
  
  
  myConnections <- neuprint_connection_table(bodyIDs, synapseType, slctROI,by.roi=by.roi,...)
  refMeta <- neuprint_get_meta(myConnections$bodyid)
  partnerMeta <- neuprint_get_meta(myConnections$partner)
  
  
  myConnections <-myConnections %>%
                  mutate(partnerName = partnerMeta[["name"]],
                         name = refMeta[["name"]],
                         partnerType = partnerMeta[["type"]],
                         type = refMeta[["type"]])
  
  ## Normalization is always from the perspective of the output (fraction of inputs to the output neuron)
  if (synapseType == "PRE"){
    outMeta <- refMeta} 
  else {
      outMeta <- partnerMeta
      }
  myConnections[["weightRelativeTotal"]] <- myConnections[["weight"]]/outMeta[["post"]]
  
  if (by.roi | !is.null(slctROI)){
    if (synapseType == "PRE"){
      outInfo <- neuprint_get_roiInfo(myConnections$bodyid)
    }else{
        outInfo <- neuprint_get_roiInfo(myConnections$partner)
      }
                                         
  postVar <- paste0(myConnections[["roi"]],".post")
  myConnections <- myConnections %>%
                   mutate(totalROIweight = sapply(1:length(postVar),function(v) outInfo[[postVar[v]]][v])) %>%
                   mutate(weightRelative=ROIweight/totalROIweight)
  }
  
  #myConnections$nameid = paste(as.character(myConnections$name), as.character(myConnections$bodyid), sep = "_")
  #myConnections$partnerid = paste(as.character(myConnections$partnerName), as.character(myConnections$partner), sep = "_")
  
  return( myConnections )
}

getConnectionTable_forSubset = function(preBodyIDs,postBodyIDs, slctROI){
  #' Filter connection table to contain only connections from preBodyIDs to postBodyIDs
  #' @return Returns a connection table as data frame.
  #' @param preBodyIDs: The bodyids of presynaptic neurons who's connections should be queried
  #' @param postBodyIDs: The bodyids of postsynaptic neurons who's connections should be queried
  #' @param slctROI: String specifying the ROI where connections should be queried.
  
  myConnections = getConnectionTable(postBodyIDs, "PRE", slctROI)
  
  myConnections = myConnections %>% filter(partner %in% preBodyIDs) 
  
  return( myConnections )
}
