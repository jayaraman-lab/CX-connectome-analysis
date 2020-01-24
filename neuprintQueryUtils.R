########################################################################################################################
# Convenience functions for query results processing
########################################################################################################################

### Neuprint search
getBodyIdsForList = function (neuronList){
  #' Get one dataframe of bodyIDs for all search strings in neuronList

  bodyIDs = neuprint_search(paste(neuronList[1],'.*',sep=""))
  
  if (length(neuronList) > 1){
    for (i in seq(2,length(neuronList))){
      bodyIDs = bind_rows(bodyIDs, neuprint_search(paste(neuronList[i],'.*',sep="")))
    }
  }
  
  return(bodyIDs)
}


### Connection table

getConnectionTable = function(bodyIDs, synapseType, slctROI){
  #' Get connection table of inputs and add meta
  #' @return Returns a connection table as data frame.
  #' @param bodyIDs: The bodyids of neurons who's connections should be queried.
  #' @param synapseType: Choose "PRE" or "POST" to get inputs or outputs of the neurons in bodyIDs, respectivly.
  #' @param slctROI: String specifying the ROI where connections should be queried.
  
  myConnections = neuprint_connection_table(bodyIDs, synapseType, slctROI)
  
  myConnections = myConnections %>%
    mutate(partnerName = neuprint_get_meta(as.numeric(partner))[["name"]],
           name = neuprint_get_meta(as.numeric(bodyid))[["name"]],
           partnerType = neuprint_get_meta(as.numeric(partner))[["type"]],
           type = neuprint_get_meta(as.numeric(bodyid))[["type"]]) %>% 
    group_by(name) %>% 
    mutate(weightRelative = weight/sum(weight)) %>%
    ungroup()
  
  myConnections$nameid = paste(as.character(myConnections$name), as.character(myConnections$bodyid), sep = "_")
  myConnections$partnerid = paste(as.character(myConnections$partnerName), as.character(myConnections$partner), sep = "_")
  
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
