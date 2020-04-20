if (!require("dtplyr")) install.packages("dtplyr")
library(dtplyr)
source("R/supertypeUtils.R")

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
getConnectionTable <- function(bodyIDs,synapseType, slctROI,by.roi, synThresh = 3,chunk_connections=TRUE,chunk_meta=TRUE,verbose=FALSE,computeKnownRatio=FALSE,...){
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
  #' @param chunk_meta : to be passed to metadata and roiInfo queries (default TRUE)
  #' @param chunk_connection : to be passed to neuprint_connection_table (default TRUE)
  #' @param verbose : should the function report on its progress?
  #' @param computeKnownRatio : should the function also compute ratios (weightRelative and outputContribution) relative to known synaptic partners rather
  #' than relative to the total number of synapses
  #' @param ...: Other arguments to be passed to neuprint queries (neuprint_connection_table, getRoiInfo and neuprint_get_meta)
  
  UseMethod("getConnectionTable")}


getConnectionTable.default = function(bodyIDs, synapseType, slctROI=NULL,by.roi=FALSE, synThresh=3,chunk_connections=TRUE,chunk_meta=TRUE,verbose=FALSE,computeKnownRatio=FALSE,...){
  refMeta <- neuprint_get_meta(bodyIDs,chunk=chunk_meta,...)
  return(getConnectionTable(refMeta,synapseType,slctROI,by.roi,synThresh,chunk_connections=chunk_connections,chunk_meta=chunk_meta,computeKnownRatio=computeKnownRatio,...))
}

getConnectionTable.data.frame <- function(bodyIDs,synapseType, slctROI=NULL,by.roi=FALSE,synThresh=3,chunk_connections=TRUE,chunk_meta=TRUE,verbose=FALSE,computeKnownRatio=FALSE,...){
  refMeta <- bodyIDs
  bodyIDs <- neuprint_ids(bodyIDs$bodyid)

  if (verbose) message("Pull connections")
  myConnections_raw <- neuprint_connection_table(bodyIDs, synapseType, slctROI,by.roi=by.roi,chunk=chunk_connections,...)

  if (by.roi | !is.null(slctROI)){
   myConnections_raw <- myConnections_raw %>% drop_na(ROIweight)
   myConnections <- myConnections_raw %>% filter(ROIweight>synThresh)
  }else{
    myConnections <- myConnections_raw %>% filter(weight>synThresh)
  }
  
  if (nrow(myConnections)==0){return(NULL)}
  
  if (verbose) message("Pull metadata")
  partnerMeta <- neuprint_get_meta(myConnections$partner,chunk=chunk_meta,...)
  refMetaOrig <- neuprint_get_meta(myConnections$bodyid,chunk=chunk_meta,...)  ## To get the database type name
  
  myConnections <- filter(myConnections,partnerMeta$status =="Traced")
  partnerMeta <- filter(partnerMeta,status == "Traced")
  
  processConnectionTable(myConnections,myConnections_raw,refMeta,partnerMeta,refMetaOrig,synapseType,by.roi,slctROI,verbose,chunk_meta,chunk_connections,computeKnownRatio,...)  
}

processConnectionTable <- function(myConnections,myConnections_raw,refMeta,partnerMeta,refMetaOrig,synapseType,by.roi,slctROI,verbose,chunk_meta,chunk_connections,computeKnownRatio,...){  
  
  
  refMeta <- slice(refMeta,as.integer(sapply(myConnections$bodyid,function(b) match(b,refMeta$bodyid))))
  refMetaOrig <- slice(refMetaOrig,as.integer(sapply(myConnections$bodyid,function(b) match(b,refMetaOrig$bodyid))))

  myConnections <-mutate(myConnections,
                         partnerName = partnerMeta[["name"]],
                         name = refMeta[["name"]],
                         partnerType = partnerMeta[["type"]],
                         type = refMeta[["type"]])
  
  myConnections <- simplifyConnectionTable(myConnections)
  ## Normalization is always from the perspective of the output (fraction of inputs to the output neuron)
  if (synapseType == "PRE"){
    if (computeKnownRatio){
      knownTablePost <- myConnections_raw  %>% mutate(from = ifelse(prepost==1,bodyid,partner),
                                                to = ifelse(prepost==1,partner,bodyid)) %>% select(-bodyid,-partner)
      knownTablePre <- neuprint_connection_table(unique(myConnections$from),"POST",slctROI,by.roi=by.roi,chunk=chunk_connections,...) %>% drop_na(ROIweight) %>% mutate(from = ifelse(prepost==1,bodyid,partner),
                                                                                                                                                                   to = ifelse(prepost==1,partner,bodyid)) %>% select(-bodyid,-partner)
    }
    outMeta <- refMeta
    inMeta <- partnerMeta
    myConnections <- mutate(myConnections,databaseType.to = refMetaOrig$type,
                                          databaseType.from = type.from)
    } else {
      if (computeKnownRatio){
        knownTablePre <- myConnections_raw  %>% mutate(from = ifelse(prepost==1,bodyid,partner),
                                                   to = ifelse(prepost==1,partner,bodyid)) %>% select(-bodyid,-partner)
        knownTablePost <- neuprint_connection_table(unique(myConnections$to),"PRE",slctROI,by.roi=by.roi,chunk=chunk_connections,...) %>% drop_na(ROIweight) %>% mutate(from = ifelse(prepost==1,bodyid,partner),
                                                                                                                                                                 to = ifelse(prepost==1,partner,bodyid)) %>% select(-bodyid,-partner)
      }
    inMeta <- refMeta
    outMeta <- partnerMeta
    myConnections <- mutate(myConnections,databaseType.to = type.to,
                                          databaseType.from = refMetaOrig$type)
  }
  
  if (computeKnownRatio){
    knownTablePostTotal <- knownTablePost %>% group_by(to) %>% distinct(from,weight) %>% mutate(knownPostWeight = sum(weight)) %>% ungroup()
    knownTablePreTotal <- knownTablePre %>% group_by(from) %>% distinct(to,weight) %>% mutate(knownPreWeight = sum(weight)) %>% ungroup()
  }
  
  myConnections <-mutate(myConnections,weightRelativeTotal = weight/outMeta[["post"]],
                                       totalPreWeight = inMeta[["downstream"]][match(myConnections$from,inMeta$bodyid)],
                                       outputContributionTotal = weight/totalPreWeight,
                                       previous.type.to = databaseType.to,
                                       previous.type.from = databaseType.from
                                            )
  if  (computeKnownRatio){
    myConnections <-mutate(myConnections,
                           knownWeightRelativeTotal = weight/knownTablePostTotal$knownPostWeight[match(myConnections$to,knownTablePostTotal$to)],
                           knownTotalPreWeight = knownTablePreTotal$knownPreWeight[match(myConnections$from,knownTablePreTotal$from)],
                           knownOutputContributionTotal = weight/knownTotalPreWeight
    )
    
  }
  
  if (by.roi | !is.null(slctROI)){
    if (verbose) message("Pull roiInfo")
    myConnections[["weightROIRelativeTotal"]] <- myConnections[["ROIweight"]]/outMeta[["post"]]
    outInfo <- getRoiInfo(unique(myConnections$to),chunk=chunk_meta,...) %>% select(bodyid,roi,post)
    inInfo <- getRoiInfo(unique(myConnections$from),chunk=chunk_meta,...) %>% select(bodyid,roi,downstream)
  
    myConnections <- left_join(myConnections,inInfo,by=c("from" = "bodyid","roi"="roi")) %>% rename(totalPreROIweight = downstream)
 
    myConnections <- left_join(myConnections,outInfo,by=c("to" = "bodyid","roi"="roi")) %>% rename(totalROIweight = post) %>%
      mutate(weightRelative=ROIweight/totalROIweight,
             outputContribution=ROIweight/totalPreROIweight)
    
    if (computeKnownRatio){
      knownTablePostROI <- knownTablePost %>% group_by(to,roi)  %>% summarize(knownPostWeight = sum(ROIweight)) %>% ungroup()
      knownTablePreROI <- knownTablePre %>% group_by(from,roi)  %>% summarize(knownPreWeight = sum(ROIweight)) %>% ungroup()
 ## This is how much this connection accounts for the outputs of the input neuron (not the standard measure)
      myConnections <- myConnections %>% mutate(knownWeightRelative = ROIweight/knownTablePostROI$knownPostWeight[match(paste0(myConnections$to,myConnections$roi),paste0(knownTablePostROI$to,knownTablePostROI$roi))],
                                              knownTotalPreROIweight = knownTablePreROI$knownPreWeight[match(paste0(myConnections$from,myConnections$roi),paste0(knownTablePreROI$from,knownTablePreROI$roi))],
                                              knownOutputContribution = ROIweight/knownTotalPreROIweight) %>% drop_na(weightRelative)  ## NA values can occur in rare cases where
    }
    ## synapse (pre/post) is split between ROIs
  }else{
  myConnections <- mutate(myConnections,roi = "All brain", 
                                        outputContribution = outputContributionTotal,
                                        weightRelative = weightRelativeTotal,
                                        ROIweight = weight)
  if (computeKnownRatio){
    myConnections <- mutate(myConnections,
                            knownOutputContribution = knownOutputContributionTotal,
                            knownWeightRelative = knownWeightRelativeTotal)
    
  }
  }
  return(supertype(myConnections))
}

simplifyConnectionTable <- function(connectionTable){
  #' Change from a bodyid/partner/prepost to a from/to format
  #' @param connectionTable: A data frame in the bodyid/partner/prepost format (as returned by getConnectionTable)
  #' @return A data frame in the from/to/name.from/name.to... format
  #' 
  #' 
  if ("from" %in% names(connectionTable)){return(connectionTable)}else{
  connectionTable <- connectionTable %>% mutate(from = ifelse(prepost==1,bodyid,partner),
                                                to = ifelse(prepost==1,partner,bodyid),
                                                name.from = as.character(ifelse(prepost==1,name,partnerName)),
                                                name.to = as.character(ifelse(prepost==1,partnerName,name)),
                                                type.from = as.character(ifelse(prepost==1,type,partnerType)),
                                                type.to = as.character(ifelse(prepost==1,partnerType,type))
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
  
  typesTable <- bind_rows(lapply(types,function(t) neuprint_search(t,field="type",fixed=TRUE,exact=TRUE)))
  
  if (length(typesTable)>0){typesTable <- typesTable %>%
                                            mutate(databaseType = type)}
  return(typesTable)
}

redefineTypeByName <- function(table,type,pattern,newPostFixes,type_col="type",name_col="name",perl=FALSE){
  condition <- grepl(pattern,table[[name_col]],perl=perl)
  redefineType(table=table,
               type=type,
               condition=condition,
               newTypes=paste0(type,newPostFixes),
               type_col=type_col)
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
  table[[paste0("previous.",type_col)]][table[[type_col]] == type] <- type ## keeping track of the last named types
  table[[type_col]][table[[type_col]] == type] <-  newTypes[2]
  table[[type_col]][(table[[type_col]] == newTypes[2]) & condition] <-  newTypes[1] 
  return(table)
}

lrSplit <- function(connectionTable,
                    nameCol="name.to",
                    typeCol="type.to",
                    databaseCol =paste0("databaseType",ifelse(grepl("[.]",typeCol),
                                                              paste0(".",c(unlist(strsplit(typeCol,"\\.")),"")[2]),
                                                       "")),
                    typeList=NULL){
  #' Retype neurons in a table according to L/R
  #' @param connectionTable : connectivity table to modify
  #' @param nameCol : the name of the column containing the names of the neurons (to be used
  #' to determine laterality)
  #' @param typeCol : the name of the column containing the types to modify
  #' @param databaseCol : the column of the table holding database type names (we only want to change
  #' type if it's not been modified before)
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
    typeList <- filter(connectionTable,((!!as.name(typeCol)) == (!!as.name(databaseCol))))
  }else{
    typeList <- filter(connectionTable,((!!as.name(typeCol)) %in% typeList) & (!!as.name(typeCol)) == (!!as.name(databaseCol)))
  }
  typeList <- distinct(typeList,(!!as.name(nameCol)),(!!as.name(typeCol)))
  typeList <- filter(typeList,grepl("_R|_L",(!!as.name(nameCol)))) %>% na.omit()
  typeList <- typeList[[typeCol]]
  condition <- grepl("_L",connectionTable[[nameCol]])
  connectionTable[[paste0("previous.",typeCol)]] <- connectionTable[[typeCol]] ## keeping track of the last named types
  for (t in typeList){
    rightType <- paste0(t,"_R")
    leftType <- paste0(t,"_L")
    connectionTable[[typeCol]][connectionTable[[typeCol]] == t] <-  rightType
    connectionTable[[typeCol]][(connectionTable[[typeCol]] == rightType) & condition] <-  leftType 
  }
  return(connectionTable)
}

retype.na <- function(connectionTable){
  #' Fill in the type and name field in case they are NAs, using the name field if it exists
  #' (removing the _L/_R) or the neuron id. By default expects a table in to/from format.
  #'
  connectionTable <- connectionTable %>% 
                     mutate(name.from = ifelse(is.na(name.from),as.character(from),name.from),
                            name.to = ifelse(is.na(name.to),as.character(to),name.to),
                            type.from = ifelse(is.na(type.from),gsub("_L$|_R$","",name.from),type.from),
                            type.to = ifelse(is.na(type.to),gsub("_L$|_R$","",name.to),type.to)
                            )
          
  return(connectionTable)
}

retype.na_meta <- function(metaTable){
  #' Fill in the type and name field in case they are NAs, using the name field if it exists
  #' (removing the _L/_R) or the neuron id. By default expects a table in to/from format.
  #'
  metaTable <- metaTable %>% 
    mutate(
           name = ifelse(is.na(name),as.character(bodyid),name),
           type = ifelse(is.na(type),gsub("_L$|_R$","",name),type)
    )
  
  return(metaTable)
}

getTypeToTypeTable <- function(connectionTable,
                               majorOutputThreshold=0.8,
                               singleNeuronThreshold=0.01,
                               singleNeuronThresholdN=3,
                               pThresh = 0.05,
                               typesTable = NULL,
                               oldTable = NULL){
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
  #' @param oldTable : only to be used in renaming functions. The type to type table before renaming.
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
  #'              - databaseType.to/databaseType.from : the neuprint defined types that contain the 
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
  if (is.null(typesTable) & any(connectionTable$type.to != connectionTable$databaseType.to,na.rm=T))
    {
      stop("Some types are custom defined. You need to provide a `typesTable` argument.")
    }
  
  if (is.null(typesTable)){
    typesTable <- getTypesTable(unique(connectionTable$databaseType.to))
  } 
  
  typesCount <- typesTable %>% group_by(type) %>%
                               summarise(n=n())
  connectionTable <- connectionTable %>% group_by(type.from) %>%
                           mutate(n_from = length(unique(from))) %>% ungroup()
  connectionTable <- connectionTable %>% 
                      mutate(n = typesCount[["n"]][match(type.to,typesCount[["type"]])])
 
  ## Renaming the unnamed neurons and treating them as single examples
  connectionTable <- connectionTable %>% mutate(n = replace_na(n,1))
  connectionTable <- retype.na(connectionTable)
  
  ## Gather the outputContributions
  connectionTable <-  connectionTable %>% group_by(from,type.to,roi) %>%
                                          mutate_at(vars(any_of(c("outputContribution","outputContributionTotal","knownOutputContribution","knownOutputContributionTotal"))),sum) %>%
                                          group_by(type.from,type.to,roi) %>%
                                          mutate_at(vars(any_of(c("outputContribution","outputContributionTotal","knownOutputContribution","knownOutputContributionTotal"))),~sum(.[match(unique(from),from)])/n_from) %>%
                                          ungroup()
  
  if (!is.null(oldTable)){
    connectionTableOld <- oldTable %>% filter(paste0(type.to,type.from) %in% paste0(connectionTable$type.to,connectionTable$type.from))
    connectionTable <- connectionTable %>% filter((type.to != previous.type.to) | (type.from != previous.type.from))
  }
  
  ## This contains the neurons unique in their type that reach our hard threshold
  loners <- connectionTable %>% filter(n==1) %>%
                                group_by_if(names(.) %in% c("type.from","type.to","roi","previous.type.from","previous.type.to","outputContribution","outputContributionTotal","knownOutputContribution","knownOutputContributionTotal",
                                                            "databaseType.to","databaseType.from",paste0("supertype.to",1:3),paste0("supertype.from",1:3))) %>%
                                mutate_at(vars(any_of(c("weightRelative","weightRelativeTotal","knownWeightRelative","knownWeightRelativeTotal"))),sum) %>%
                                mutate(weight = sum(ROIweight),
                                       absoluteWeight = sum(ROIweight),
                                       n_type = 1,
                                       n_targets = n()) %>%
                                 summarize_at(vars(any_of(c("weightRelative","weightRelativeTotal","knownWeightRelative","knownWeightRelativeTotal","weight","absoluteWeight","n_type","n_targets"))),first) %>% ungroup() 
  
  group_In <- names(connectionTable)[names(connectionTable) %in% c("type.from","to","type.to","roi","previous.type.from","previous.type.to","n","outputContribution","outputContributionTotal","knownOutputContribution","knownOutputContributionTotal",
                                            "databaseType.to","databaseType.from",paste0("supertype.to",1:3),paste0("supertype.from",1:3))]
  
  group_Out <- names(connectionTable)[names(connectionTable) %in% c("type.from","type.to","roi","previous.type.from","previous.type.to","outputContribution","outputContributionTotal","knownOutputContribution","knownOutputContributionTotal",
                 "databaseType.to","databaseType.from",paste0("supertype.to",1:3),paste0("supertype.from",1:3))]
  
  if ("knownWeightRelative" %in% names(connectionTable)){
    sTable <- connectionTable %>% filter(n>1) %>%
      group_by_at(group_In) %>%
      summarize(weight=sum(ROIweight),
                weightRelative=sum(weightRelative),
                weightRelativeTotal=sum(weightRelativeTotal),
                knownWeightRelative=sum(knownWeightRelative),
                knownWeightRelativeTotal=sum(knownWeightRelativeTotal)) %>%
      group_by_at(group_Out) %>%
      summarize(n_targets = n(),
                n_type = n[1],
                absoluteWeight = sum(weight),
                weightM = sum(weight)/n_type,
                weightRM = sum(weightRelative)/n_type,#mean(c(weightRelative,unlist(missingV))),
                weightRelativeTotal = sum(weightRelativeTotal)/n_type,
                weightRM2 = sum(weightRelative^2)/n_type,
                knownWeightRelative = sum(knownWeightRelative)/n_type,
                knownWeightRelativeTotal = sum(knownWeightRelativeTotal)/n_type
      ) %>% ungroup()
  }else{
    sTable <- connectionTable %>% filter(n>1) %>%
      group_by_at(group_In) %>%
      summarize(weight=sum(ROIweight),
                weightRelative=sum(weightRelative),
                weightRelativeTotal=sum(weightRelativeTotal)) %>%
      group_by_at(group_Out) %>%
      summarize(n_targets = n(),
                n_type = n[1],
                absoluteWeight = sum(weight),
                weightM = sum(weight)/n_type,
                weightRM = sum(weightRelative)/n_type,
                weightRelativeTotal = sum(weightRelativeTotal)/n_type,
                weightRM2 = sum(weightRelative^2)/n_type
      ) %>% ungroup()
    
  }
  
  sTable <-sTable %>% mutate(varWeight = (n_type/(n_type-1))*(weightRM2 - weightRM^2),
                             sdWeight = sqrt(varWeight),
                             tt = weightRM/(sdWeight/sqrt(n_type)),
                             pVal = pt(tt,n_type-1,lower.tail = FALSE)) %>% rename(weight=weightM,weightRelative=weightRM)
  
  if (is.null(oldTable)){
    loners <-  loners %>% filter((weightRelative > singleNeuronThreshold & weight > singleNeuronThresholdN)| outputContribution > majorOutputThreshold)
    sTable <- sTable %>% filter(pVal < pThresh | (outputContribution > majorOutputThreshold & weight > singleNeuronThresholdN)) %>%
      select(-pVal)
  }else{
    loners <-  loners %>% filter((weightRelative > singleNeuronThreshold & weight > singleNeuronThresholdN)| outputContribution > majorOutputThreshold | 
                                   (paste0(previous.type.to,previous.type.from) %in% paste0(oldTable$type.to,oldTable$type.from)))
    sTable <- sTable%>% filter(pVal < pThresh | (outputContribution > majorOutputThreshold & weight > singleNeuronThresholdN) | 
                                 (paste0(previous.type.to,previous.type.from) %in% paste0(oldTable$type.to,oldTable$type.from))) %>%
                                select(-pVal)
    sTable <- bind_rows(sTable,connectionTableOld)
  }               
  
  
  return(bind_rows(sTable,loners))                  

}

getConnectionTable_forSubset <- function(preBodyIDs,postBodyIDs, slctROI=NULL,...){
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

getRoiInfo <- function(bodyids,...){
  #' Formats the results of a roiInfo query in a nice table
  #' 
  if (length(bodyids)==0){return(data.frame(bodyid=double(),roi=character(),pre=integer(),post=integer(),downstream=integer(),stringsAsFactors = FALSE))}
  roiInfo <- neuprint_get_roiInfo(bodyids,...)
  roiInfo <-  pivot_longer(roiInfo,cols=-bodyid,names_to=c("roi","field"),names_sep="\\.",values_to="count")
  roiInfo <- pivot_wider(roiInfo,names_from = "field",values_from="count")
  roiInfo
}