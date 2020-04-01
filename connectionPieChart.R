# Utility functions for generating hierarchical pie charts (similar to those shown on neuprint explorer)

prepData4pieChart <- function(myType_bag, partnerDirection, slctRoi = NULL,  splitLR = TRUE, myType = NULL, minweight = 0){
  #' Prepare neuron bag data for plotting of pie chart..
  #' @param myType_bag: Neuron bag object for type of interest
  #' @param partnerDirection: Chose "input" vs. "output" partners
  #' @param slctRoi: Select ROIs to display. Ideally, these should be non-overlapping.
  #' @param splitLR: Split types into left and right population.
  #' @param myType: Optional filtering for subtype in neuron bag.
  #' @param minweight: Option to filter by fraction of contribution to inputs/outputs per region.
  
  if(splitLR){
    myType_bag = lateralizeInputOutputList(myType_bag)
  }
  if (partnerDirection == "input"){
    myTypeData = myType_bag$inputs_raw %>%
      rename(target = to, partner = from,
             targetType = type.to, partnerType = type.from)
  }else if(partnerDirection == "output"){
    myTypeData = myType_bag$outputs_raw %>%
      rename(target = from, partner = to,
             targetType = type.from, partnerType = type.to)
  }else{
    print("Choose input/output direction.")
    return()
  }
  
  if(!is.null(myType)){
    myTypeData = myTypeData %>% filter(targetType == myType)
  }
  
  myTypeData = myTypeData %>% group_by(target, partnerType, roi) %>% 
    summarise(sumweight = sum(weightRelativeTotal)) %>% filter(sumweight >= minweight)
  
  if(!is.null(slctRoi)){
    myTypeData = myTypeData %>% filter(roi %in% slctRoi)
  }
  
  myTypeData = myTypeData %>% group_by(roi, partnerType) %>% 
    summarise(avWeight = mean(sumweight)) %>% 
    arrange(roi, partnerType) %>% ungroup() 
  
  myTypeData$roi = as.factor(myTypeData$roi)
  myTypeData$partnerType = as.character(myTypeData$partnerType)
  
  return(myTypeData)
}


generatePieChartData <- function(myTypeIn, myType, rootlab){
  roiFrac = myTypeIn %>% group_by(roi) %>% summarise(roifrac = sum(avWeight)) 
  
  pieData = data.frame(labels = c(rootlab,
                                  unique(as.character(myTypeIn$roi)) ,
                                  myTypeIn$partnerType),
                       
                       parents = c("",
                                   rep(rootlab,length(unique(as.character(myTypeIn$roi)))) ,
                                   paste(rootlab, myTypeIn$roi, sep=" - ")),
                       
                       ids = c(rootlab,
                               paste(rootlab, unique(as.character(myTypeIn$roi)), sep=" - ") ,
                               paste(rootlab, myTypeIn$roi,myTypeIn$partnerType, sep=" - ")),
                       
                       values = c(sum(roiFrac$roifrac),
                                  roiFrac$roifrac, 
                                  myTypeIn$avWeight)
  )
  return(pieData)
}
