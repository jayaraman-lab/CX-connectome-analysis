library(neuprintrExtra)
library(neuprintr)
library(dplyr)
library(tidyr)
library(broom)
library(purrr)

#' Gets a table of all neurons in a sub-ROI, with information about the innervation relative to an englobing ROI. 
#' Used for in the validation pipeline to select neurons of interest.
#' @param subroi The ROI we want the neurons to innervate
#' @param Roi The ROI we want to compute innervation ratios relative to
#' @details For example, we use this function to compute how much of the innervation of a given neuron is confined to one glomerulus (the subroi) relative to 
#' the all PB (Roi)
#' @return a table of neurons, with added `downstreamRatio` and `upstreamRatio columns`
subroiSet <- function(subroi,Roi){
  
  subroiNeurons <- neuprint_bodies_in_ROI(subroi,all_segments = FALSE)
  subroiInfo <- neuprint_get_roiInfo(subroiNeurons$bodyid,all_segments = FALSE)
  subroiMeta <- neuprint_get_meta(subroiNeurons$bodyid,all_segments = FALSE) %>% filter(status=="Traced")
  subroiInfo <- select(subroiInfo,bodyid,starts_with(paste0(Roi,".")),starts_with(paste0(subroi,".")))
  
  subroiNeurons <- left_join(subroiMeta,subroiInfo,by = "bodyid") %>% filter(!is.na(type))
  subroiNeurons[is.na(subroiNeurons)] <- 0

  subroiNeurons <- rename_at(subroiNeurons,vars(starts_with(paste0(subroi,"."))),str_replace,fixed(subroi),"subroi") %>%
    mutate(upstreamRatio = subroi.upstream/!!as.name(paste0(Roi,".upstream")),
                   downstreamRatio = subroi.downstream/!!as.name(paste0(Roi,".downstream")),
                   databaseType=type)
  subroiNeurons[is.na(subroiNeurons)] <- 0
  subroiNeurons <- supertype(subroiNeurons) %>% mutate(subregion = subroi)
}

#' Fits a linear model between columns `predictor` and `predicted` in a table of 
#' innervation comparisons `compTable` as those returned in the validation notebooks, grouped by
#' columns `groups`
getFit <- function(compTable,
                   groups=c("side","databaseType","supertype2"),
                   predicted="C3",
                   predictor="CX"){
  compFits <- compTable %>%
    nest(data = !(any_of(groups))) %>% 
    mutate(
      fitRes = map(data, ~ lm(as.formula(paste0(predicted,"~",predictor)), data = .x)),
      tidied = map(fitRes, tidy),
      fitStats = map(fitRes,glance)
    ) 
  compFitsCoeff <- compFits %>% 
    unnest(tidied) %>% filter(term==predictor) %>% select(-data,-fitStats,-fitRes)
  
  compFitsStats <-  compFits %>% 
    unnest(fitStats) %>% select(-data,-fitRes,-tidied)
  left_join(compFitsCoeff,compFitsStats,by=groups)
}

#' Create a table of connectivity comparisons between two versions of the dataset
#' @param bodyidsIn The neurons to compare the inputs of
#' @param bodyidsOut The neurons to compare the outputs of
#' @param roi The region to make the comparison in
#' @param newDataset The new dataset
#' @param oldDataset
#' @return A data frame of all connections found in at least one of the two data frames, with columns:
#' - \code{side} specify if the connection is to be used for inputs or outputs comparisons
#' - \code{known_synapses} refers to the number of synapses made onto (or received from) traced partners (rather than fragments)
#' - \code{knownWeightRelative} and \code{knownOutputContribution} are similar to \code{weightRelative} and \code{outputContribution} except
#'  that they are relative to the synapses made to traced partners rather than relative to all the neuron's synapses
#' - \code{old} and \code{new} are \code{knownWeightRelative.old/new} if \code{side} is "inputs", \code{knownOutputContribution.old/new} is it's "outputs".
#'  Similarly, \code{bodyid} is the postsynaptic partner if \code{side} is "inputs" and the presynaptic partner otherwise 
#' - \code{known_synapses_ratio} is the fraction increase in the number of synapses to traced partners between the old and the
#'  new dataset for the neuron defined in \code{bodyid}
#' The function also computes completedness statistics.
getCompConnectionTable <- function(bodyidsIn,bodyidsOut,ROI,newDataset,oldDataset){
  inputsComp <- getCompConnectionTableInternal(bodyidsIn,ROI,"inputs",newDataset,oldDataset)
  outputsComp <- getCompConnectionTableInternal(bodyidsOut,ROI,"outputs",newDataset,oldDataset)
  inputsComp <- mutate(inputsComp,
                       old=knownWeightRelative.old,
                       new=knownWeightRelative.new,
                       oldDatasetcompletedness=input_completedness.old,
                       newDatasetcompletedness=input_completedness.new,
                       supertype2=supertype2.to) %>% mutate(side="inputs")
  
  outputsComp <- mutate(outputsComp,
                        old=knownOutputContribution.old,
                        new=knownOutputContribution.new,
                        oldDatasetcompletedness=output_completedness.old,
                        newDatasetcompletedness=output_completedness.new,
                        supertype2=supertype2.from) %>% mutate(side="outputs")
  
  bind_rows(inputsComp,outputsComp)
  
}

#Internal function
getCompConnectionTableInternal <- function(bodyids,ROI,polarity=c("inputs","outputs"),newDataset,oldDataset){
  polarity <- match.arg(polarity)
  
  newRawTable <- neuprint_connection_table(bodyids,ifelse(polarity=="inputs","PRE","POST"),roi=ROI,dataset=newDataset) %>% 
    drop_na(ROIweight)
  
  oldRawTable <- neuprint_connection_table(bodyids,ifelse(polarity=="inputs","PRE","POST"),roi=ROI,dataset=oldDataset) %>% 
    drop_na(ROIweight)

  newRefTable <- neuprint_get_meta(bodyids,dataset=newDataset)
  ## Get the types from the new dataset
  oldRefTable <- neuprint_get_meta(bodyids,dataset=oldDataset) %>% 
    mutate(type=newRefTable$type[match(bodyid,newRefTable$bodyid)]) 

  newPartnerTable <- neuprint_get_meta(newRawTable$partner,dataset=newDataset)

  ## Filter for bodyids that exist and are traced in both datasets
  newPartnersOldDataset <- neuprint_get_meta(newRawTable$partner,dataset=oldDataset) %>% filter(status == "Traced")
  oldPartnersNewDataset <- neuprint_get_meta(oldRawTable$partner,dataset=newDataset) %>% filter(status == "Traced")

  oldPartnerTable <- neuprint_get_meta(oldRawTable$partner,dataset=oldDataset) %>% 
    mutate(type=oldPartnersNewDataset$type[match(bodyid,oldPartnersNewDataset$bodyid)]) 

  newRawTable <- filter(newRawTable,partner %in% newPartnersOldDataset$bodyid)
  oldRawTable <- filter(oldRawTable,partner %in% oldPartnersNewDataset$bodyid)

  newPartnerTable <- filter(newPartnerTable,bodyid %in% newPartnersOldDataset$bodyid)
  oldPartnerTable <- filter(oldPartnerTable,bodyid %in% oldPartnersNewDataset$bodyid)

  
  newDatasetonnTable <- processConnectionTableOld(newRawTable,3,newRefTable,newPartnerTable,newRefTable,ifelse(polarity=="inputs","PRE","POST"),
                                         slctROI=ROI,by.roi=FALSE,verbose=FALSE,chunk_meta=TRUE,dataset=newDataset,computeKnownRatio=TRUE,chunk_connections = TRUE)
  
  oldDatasetonnTable <- processConnectionTableOld(oldRawTable,3,oldRefTable,oldPartnerTable,oldRefTable,ifelse(polarity=="inputs","PRE","POST"),
                                         slctROI=ROI,by.roi=FALSE,verbose=FALSE,chunk_meta=TRUE,dataset=oldDataset,computeKnownRatio = TRUE,chunk_connections = TRUE)
  
  groupVars <- c("from","to","roi","type.from","type.to",paste0("supertype",1:3,".from"),paste0("supertype",1:3,".to"))
  measureVars <- c(paste0(c("ROIweight",ifelse(polarity=="inputs","knownWeightRelative","knownOutputContribution")),".old"),
                   paste0(c("ROIweight",ifelse(polarity=="inputs","knownWeightRelative","knownOutputContribution")),".new"))
  extraMeasure <- c(paste0(c("input_completedness","output_completedness","knownTotalROIweight","knownTotalPreROIweight"),".old"),
    c(paste0(c("input_completedness","output_completedness","knownTotalROIweight","knownTotalPreROIweight"),".new")))
  
  measureReplace <- rep(0,length(measureVars))
  names(measureReplace) <- measureVars
  compTable <- full_join(newDatasetonnTable,oldDatasetonnTable,by = groupVars,suffix=c(".new",".old")) %>% select_at(c(groupVars,measureVars,extraMeasure)) %>%
    replace_na(as.list(measureReplace))
  compTable <- na_matchValue(compTable,"input_completedness.old","to") %>% 
               na_matchValue("input_completedness.new","to") %>%
               na_matchValue("knownTotalROIweight.old","to") %>%
               na_matchValue("knownTotalROIweight.new","to") %>%
               na_matchValue("output_completedness.old","from") %>% 
               na_matchValue("output_completedness.new","from") %>%
               na_matchValue("knownTotalPreROIweight.old","from") %>%
               na_matchValue("knownTotalPreROIweight.new","from")
  if (polarity == "inputs"){
    compTable <- mutate(compTable,
                        completedness.new=input_completedness.new,
                        completedness.old=input_completedness.old,
                        known_synapses.old=knownTotalROIweight.old,
                        known_synapses.new=knownTotalROIweight.new,
                        databaseType=type.to,
                        bodyid=to)
  }else{
    compTable <- mutate(compTable,
                        completedness.new=output_completedness.new,
                        completedness.old=output_completedness.old,
                        known_synapses.old=knownTotalPreROIweight.old,
                        known_synapses.new=knownTotalPreROIweight.new,
                        databaseType=type.from,
                        bodyid=from)
    
  }
  
  compTable <- mutate(compTable,known_synapses_ratio=(known_synapses.new - known_synapses.old)/known_synapses.old)
  compTable
}

na_matchValue <- function(table,variable,variableMatch){
  table[[variable]][is.na(table[[variable]])] <- table[[variable]][!is.na(table[[variable]])][match(table[[variableMatch]][is.na(table[[variable]])],
                                                                                                    table[[variableMatch]][!is.na(table[[variable]])])]
  table
}

# A modified internal neuprintrExtra function to compute the "known" family of statistics
processConnectionTableOld <- function(myConnections,synThresh,refMeta,partnerMeta,refMetaOrig,synapseType,by.roi,slctROI,verbose,chunk_meta,chunk_connections,computeKnownRatio,...){
  
  
  refMeta <- slice(refMeta,match(myConnections$bodyid,bodyid))
  refMetaOrig <- slice(refMetaOrig,match(myConnections$bodyid,bodyid))
  
  myConnections <-mutate(myConnections,
                         partnerName = partnerMeta[["name"]],
                         name = refMeta[["name"]],
                         partnerType = partnerMeta[["type"]],
                         type = refMeta[["type"]])
  
  myConnections <- simplifyConnectionTableLocal(myConnections)
  ## Normalization is always from the perspective of the output (fraction of inputs to the output neuron)
  if (synapseType == "PRE"){
    if (computeKnownRatio){
      knownTablePost <- myConnections
      if (length(myConnections$from)==0){knownTablePre <- empty_connTableLocal(by.roi | !(is.null(slctROI)))}else{
        knownTablePre <- neuprint_connection_table(unique(myConnections$from),"POST",slctROI,by.roi=by.roi,chunk=chunk_connections,...)}
      knownTablePre <- knownTablePre %>% mutate(from = bodyid,to = partner) %>% select(-bodyid,-partner,-prepost)
      
      if (by.roi | !is.null(slctROI)){
        knownTablePre <- tidyr::drop_na(knownTablePre,ROIweight)
        knownTablePre <- filter(knownTablePre,ROIweight>synThresh)
      }else{
        knownTablePre <- filter(knownTablePre,weight>synThresh)
      }
      
      knownPreMeta <- neuprint_get_meta(knownTablePre$to,chunk=chunk_meta,...)
      knownTablePre <- filter(knownTablePre,knownPreMeta$status=="Traced")
    }
    outMeta <- refMeta
    inMeta <- partnerMeta
    myConnections <- mutate(myConnections,databaseType.to = refMetaOrig$type,
                            databaseType.from = type.from)
  } else {
    if (computeKnownRatio){
      knownTablePre <- myConnections
      if (length(myConnections$to)==0){knownTablePost <- empty_connTableLocal(by.roi | !(is.null(slctROI)))}else{
        knownTablePost <- neuprint_connection_table(unique(myConnections$to),"PRE",slctROI,by.roi=by.roi,chunk=chunk_connections,...)}
      knownTablePost <- knownTablePost %>% mutate(from = partner,to = bodyid) %>% select(-bodyid,-partner,-prepost)
      if (by.roi | !is.null(slctROI)){
        knownTablePost <- tidyr::drop_na(knownTablePost,ROIweight)
        knownTablePost <- filter(knownTablePost,ROIweight>synThresh)
      }else{
        knownTablePost <- filter(knownTablePost,weight>synThresh)
      }
      knownPostMeta <- neuprint_get_meta(knownTablePost$from,chunk=chunk_meta,...)
      knownTablePost <- filter(knownTablePost,knownPostMeta$status=="Traced")
    }
    inMeta <- refMeta
    outMeta <- partnerMeta
    myConnections <- mutate(myConnections,databaseType.to = type.to,
                            databaseType.from = refMetaOrig$type)
  }
  
  if (computeKnownRatio){
    if (by.roi | !(is.null(slctROI))){
      knownTablePre <- tidyr::drop_na(knownTablePre,ROIweight)
    }
    knownTablePostTotal <- knownTablePost %>% group_by(to) %>% distinct(from,weight) %>% mutate(knownPostWeight = sum(weight)) %>% ungroup()
    knownTablePreTotal <- knownTablePre %>% group_by(from) %>% distinct(to,weight) %>% mutate(knownPreWeight = sum(weight)) %>% ungroup()
  }
  
  myConnections <-mutate(totalWeight = outMeta[["post"]],
                         myConnections,weightRelativeTotal = weight/outMeta[["post"]],
                         totalPreWeight = inMeta[["downstream"]][match(myConnections$from,inMeta$bodyid)],
                         outputContributionTotal = weight/totalPreWeight,
                         previous.type.to = databaseType.to,
                         previous.type.from = databaseType.from
  )
  if  (computeKnownRatio){
    myConnections <-mutate(myConnections,
                           knownTotalWeight = knownTablePostTotal$knownPostWeight[match(myConnections$to,knownTablePostTotal$to)],
                           knownWeightRelativeTotal = weight/knownTotalWeight,
                           knownTotalPreWeight = knownTablePreTotal$knownPreWeight[match(myConnections$from,knownTablePreTotal$from)],
                           knownOutputContributionTotal = weight/knownTotalPreWeight,
                           output_completednessTotal = knownTotalPreWeight/totalPreWeight,
                           input_completednessTotal =  knownTotalWeight/totalWeight
    )
    
  }
  
  if (by.roi | !is.null(slctROI)){
    if (verbose) message("Pull roiInfo")
    myConnections[["weightROIRelativeTotal"]] <- myConnections[["ROIweight"]]/outMeta[["post"]]
    outInfo <- getRoiInfo(unique(myConnections$to),chunk=chunk_meta,...) %>% select(bodyid,roi,post)
    inInfo <- getRoiInfo(unique(myConnections$from),chunk=chunk_meta,...)
    inInfo <-  select(inInfo,bodyid,roi,downstream)
    
    myConnections <- left_join(myConnections,inInfo,by=c("from" = "bodyid","roi"="roi")) %>% rename(totalPreROIweight = downstream)
    
    myConnections <- left_join(myConnections,outInfo,by=c("to" = "bodyid","roi"="roi")) %>% rename(totalROIweight = post) %>%
      mutate(weightRelative=ROIweight/totalROIweight,
             outputContribution=ROIweight/totalPreROIweight)
    
    if (computeKnownRatio){
      knownTablePostROI <- knownTablePost %>% group_by(to,roi)  %>% summarize(knownPostWeight = sum(ROIweight)) %>% ungroup()
      knownTablePreROI <- knownTablePre %>% group_by(from,roi)  %>% summarize(knownPreWeight = sum(ROIweight)) %>% ungroup()
      ## This is how much this connection accounts for the outputs of the input neuron (not the standard measure)
      myConnections <- myConnections %>% mutate(knownTotalROIweight = knownTablePostROI$knownPostWeight[match(paste0(myConnections$to,myConnections$roi),paste0(knownTablePostROI$to,knownTablePostROI$roi))],
                                                knownWeightRelative = ROIweight/knownTotalROIweight,
                                                knownTotalPreROIweight = knownTablePreROI$knownPreWeight[match(paste0(myConnections$from,myConnections$roi),paste0(knownTablePreROI$from,knownTablePreROI$roi))],
                                                knownOutputContribution = ROIweight/knownTotalPreROIweight,
                                                output_completedness = knownTotalPreROIweight/totalPreROIweight,
                                                input_completedness = knownTotalROIweight/totalROIweight
      ) %>% tidyr::drop_na(weightRelative)  ## NA values can occur in rare cases where
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

simplifyConnectionTableLocal <- function(connectionTable){
  
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

empty_connTableLocal <- function(has.roi = FALSE){
  res <- data.frame(bodyid=numeric(),partner=numeric(),prepost=integer(),weight=integer())
  if (has.roi){res <- mutate(res,roi=character(),ROIweight=integer())}
  res
}