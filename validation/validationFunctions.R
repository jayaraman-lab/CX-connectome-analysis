##### GET A TABLE WITH RATIOS OF INNERVATION OF A SUBROI OF A ROI
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

#### GET ROI INFO FOR INPUT/OUTPUT SETS ACCROSS 2 VERSIONS OF THE DATABASE
getCompRoiInfoInternal <- function(bodyids,ROI,polarity=c("Inputs","Outputs"),newC,oldC){
  polarity <- match.arg(polarity)
  newRefTable <- neuprint_get_meta(bodyids,conn=newC)

  newRoiTable <- getRoiInfo(bodyids,conn=newC) %>% filter(roi == ROI) %>% mutate(version="new")
  oldRoiTable <- getRoiInfo(bodyids,conn=oldC) %>% filter(roi == ROI) %>% mutate(version="old")

  roiTable <- bind_rows(newRoiTable,oldRoiTable) 
  if (polarity == "Inputs")
   roiTable <- pivot_wider(roiTable,names_from = version,values_from=upstream,id_cols = bodyid)
  else
    roiTable <- pivot_wider(roiTable,names_from = version,values_from=downstream,id_cols = bodyid)
 
  mutate(roiTable,databaseType=newRefTable$type[match(bodyid,newRefTable$bodyid)]) %>% supertype() %>% mutate(side=polarity)
}

getCompRoiInfo <- function(bodyidsIn,bodyidsOut,roi,newC,oldC){
  roiTableIn <- getCompRoiInfoInternal(bodyidsIn,roi,polarity="Inputs",newC=newC,oldC=oldC)
  roiTableOut <- getCompRoiInfoInternal(bodyidsOut,roi,polarity="Outputs",newC=newC,oldC=oldC)
  
  bind_rows(roiTableIn,roiTableOut) 
}

#### MAKE COMPARABLE NEURON TO NEURON CONNECTION TABLES BETWEEN OLD (missing type) AND NEW DATASET
getCompConnectionTableInternal <- function(bodyids,ROI,polarity=c("Inputs","Outputs"),newC,oldC){
  polarity <- match.arg(polarity)
  
  newRawTable <- neuprint_connection_table(bodyids,ifelse(polarity=="Inputs","PRE","POST"),roi=ROI,conn=newC) %>% drop_na(ROIweight)
  oldRawTable <- neuprint_connection_table(bodyids,ifelse(polarity=="Inputs","PRE","POST"),roi=ROI,conn=oldC) %>% drop_na(ROIweight)

  newRefTable <- neuprint_get_meta(bodyids,conn=newC)
  ## Get the types from the new dataset
  oldRefTable <- neuprint_get_meta(bodyids,conn=oldC) %>% mutate(type=newRefTable$type[match(bodyid,newRefTable$bodyid)]) 

  newPartnerTable <- neuprint_get_meta(newRawTable$partner,conn=newC)

  ## Filter for bodyids that exist and are traced in both datasets
  newPartnersOldDataset <- neuprint_get_meta(newRawTable$partner,conn=oldC) %>% filter(status == "Traced")
  oldPartnersNewDataset <- neuprint_get_meta(oldRawTable$partner,conn=newC) %>% filter(status == "Traced")

  oldPartnerTable <- neuprint_get_meta(oldRawTable$partner,conn=oldC) %>% mutate(type=oldPartnersNewDataset$type[match(bodyid,oldPartnersNewDataset$bodyid)]) 

  newRawTable <- filter(newRawTable,partner %in% newPartnersOldDataset$bodyid)
  oldRawTable <- filter(oldRawTable,partner %in% oldPartnersNewDataset$bodyid)

  newPartnerTable <- filter(newPartnerTable,bodyid %in% newPartnersOldDataset$bodyid)
  oldPartnerTable <- filter(oldPartnerTable,bodyid %in% oldPartnersNewDataset$bodyid)

  
  newConnTable <- processConnectionTableOld(newRawTable,3,newRefTable,newPartnerTable,newRefTable,ifelse(polarity=="Inputs","PRE","POST"),
                                         slctROI=ROI,by.roi=FALSE,verbose=FALSE,chunk_meta=TRUE,conn=newC,computeKnownRatio=TRUE,chunk_connections = TRUE)
  
  oldConnTable <- processConnectionTableOld(oldRawTable,3,oldRefTable,oldPartnerTable,oldRefTable,ifelse(polarity=="Inputs","PRE","POST"),
                                         slctROI=ROI,by.roi=FALSE,verbose=FALSE,chunk_meta=TRUE,conn=oldC,computeKnownRatio = TRUE,chunk_connections = TRUE)
  
  groupVars <- c("from","to","roi","type.from","type.to",paste0("supertype",1:3,".from"),paste0("supertype",1:3,".to"))
  measureVars <- c(paste0(c("ROIweight",ifelse(polarity=="Inputs","knownWeightRelative","knownOutputContribution")),".old"),
                   paste0(c("ROIweight",ifelse(polarity=="Inputs","knownWeightRelative","knownOutputContribution")),".new"))
  extraMeasure <- c(paste0(c("input_completedness","output_completedness","knownTotalROIweight","knownTotalPreROIweight"),".old"),
    c(paste0(c("input_completedness","output_completedness","knownTotalROIweight","knownTotalPreROIweight"),".new")))
  
  measureReplace <- rep(0,length(measureVars))
  names(measureReplace) <- measureVars
  compTable <- full_join(newConnTable,oldConnTable,by = groupVars,suffix=c(".new",".old")) %>% select_at(c(groupVars,measureVars,extraMeasure)) %>%
    replace_na(as.list(measureReplace))
  compTable <- na_matchValue(compTable,"input_completedness.old","to") %>% 
               na_matchValue("input_completedness.new","to") %>%
               na_matchValue("knownTotalROIweight.old","to") %>%
               na_matchValue("knownTotalROIweight.new","to") %>%
               na_matchValue("output_completedness.old","from") %>% 
               na_matchValue("output_completedness.new","from") %>%
               na_matchValue("knownTotalPreROIweight.old","from") %>%
               na_matchValue("knownTotalPreROIweight.new","from")
  if (polarity == "Inputs"){
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

getCompConnectionTable <- function(bodyidsIn,bodyidsOut,ROI,newC,oldC){
  inputsComp <- getCompConnectionTableInternal(bodyidsIn,ROI,"Inputs",newC,oldC)
  outputsComp <- getCompConnectionTableInternal(bodyidsOut,ROI,"Outputs",newC,oldC)
  inputsComp <- mutate(inputsComp,
                       old=knownWeightRelative.old,
                       new=knownWeightRelative.new,
                       oldCompletedness=input_completedness.old,
                       newCompletedness=input_completedness.new,
                       supertype2=supertype2.to) %>% mutate(side="Inputs")
  
  outputsComp <- mutate(outputsComp,
                        old=knownOutputContribution.old,
                        new=knownOutputContribution.new,
                        oldCompletedness=output_completedness.old,
                        newCompletedness=output_completedness.new,
                        supertype2=supertype2.from) %>% mutate(side="Outputs")
  
  bind_rows(inputsComp,outputsComp)
  
}

getFit <- function(compTable,groups=c("side","databaseType","supertype2"),predicted="C3",predictor="CX"){
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
      if (length(myConnections$from)==0){knownTablePre <- empty_connTable(by.roi | !(is.null(slctROI)))}else{
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
      if (length(myConnections$to)==0){knownTablePost <- empty_connTable(by.roi | !(is.null(slctROI)))}else{
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