########################################################################################################################
# Tools to look at pathways through regions
########################################################################################################################

library("jsonlite")

### Get the relative percentage of pre and postsynaptic sites within a given ROI

preNpostRatios <- function(nronTypes, ROI) {
  
  # get all of the ROIs in the neuprint dataset
  rois = neuprint_ROIs()
  
  # initialize a data frame to hold the count ratios
  preNpost = data.frame(type = character(),bodyId = integer(),PBPre = integer(),otherPre = integer(),allPre = integer(),PBPost = integer(),otherPost = integer(),allPost = integer())
  
  # find the % input and output for each neuron
  for (tp in 1:length(nronTypes)){
    nmNow = gsub(")","\\\\\\\\\\\\\\\\)",nronTypes[[tp]])
    nmNow = gsub("\\(","\\\\\\\\\\\\\\\\(",nmNow)
    queryNow = "MATCH (m:`hemibrain_Meta`) WITH m.superLevelRois AS rois MATCH (neuron :`hemibrain_Neuron`) WHERE (neuron.type =~\\\"" %>%
      paste0(nmNow) %>%
      paste0("\\\" OR neuron.instance =~\\\"") %>% paste0(nmNow) %>%
      paste0("\\\" AND neuron.status = \\\"Traced\\\" AND neuron.cropped = false) RETURN neuron.bodyId AS bodyid, neuron.roiInfo AS roiInfo, neuron.pre AS npre, neuron.post AS npost")
    
    roiQuery = neuprint_fetch_custom(cypher = queryNow)
    
    # pull out the relevant data from the query
    bodyIds = roiQuery$data %>% lapply(function(x){x[[1]]})
    roiNums = roiQuery$data %>% lapply(function(x){x[[2]]})
    totPre = roiQuery$data %>% lapply(function(x){x[[3]]})
    totPost = roiQuery$data %>% lapply(function(x){x[[4]]})
    
    # populate the pre and post synaptic columns
    for (bId in 1:length(roiNums)){
      roiDatNow = jsonlite::fromJSON(as.character(roiNums[bId]))
      
      synVec = c(0,0,0,0)
      
      for (nm in 1:length(names(roiDatNow))){
        if (names(roiDatNow)[nm] %in% ROI){
          if ("pre" %in% names(roiDatNow[[nm]]))
            synVec[1] = roiDatNow[[nm]]$pre
          if ("post" %in% names(roiDatNow[[nm]]))
            synVec[3] = roiDatNow[[nm]]$post
        } else if (grepl(ROI,names(roiDatNow)[nm])){
          break
        } else {
          if ("pre" %in% names(roiDatNow[[nm]]))
            synVec[2] =+ roiDatNow[[nm]]$pre
          if ("post" %in% names(roiDatNow[[nm]]))
            synVec[4] =+ roiDatNow[[nm]]$post
        }
      }
      
      preNpost <- preNpost %>%
        rbind(data.frame(type = nronTypes[[tp]],
                         bodyId = bodyIds[[bId]],
                         regPre = synVec[1],
                         otherPre = synVec[2],
                         allPre = totPre[[bId]],
                         regPost = synVec[3],
                         otherPost = synVec[4],
                         allPost = totPost[[bId]]))
    }
    
  }
  
  preNpost <- preNpost %>% mutate(preRelative = regPre/allPre)
  preNpost <- preNpost %>% mutate(postRelative = regPost/allPost)
  preNpost <- preNpost %>% mutate(prePostRatio = preRelative/postRelative)
  preNpost <- preNpost %>% mutate(prePostDiff = preRelative-postRelative)
  
  return(preNpost)
  
}