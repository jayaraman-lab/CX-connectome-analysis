# This file contain functions for:
# 1) Get_AllNeurons_InRoi gets the bodies from a given ROI. Must specify whether structure is paired or not.
# 2)SynapseStats_And_Threshold computes synapse statistics and thresholds the data (still need to agree on method here)


Get_AllNeurons_InRoi <- function(slctROI, pairedRegion) {
  
  # Get neurons from neuprint
  if (pairedRegion) {
    roiR_Connect <- neuprint_find_neurons( paste(slctROI, "(R)", sep=""), all_segments = FALSE)
    roiL_Connect <- neuprint_find_neurons( paste(slctROI, "(L)", sep=""), all_segments = FALSE)
  } else {
    roi_Connect <- neuprint_find_neurons(slctROI, all_segments = FALSE) #not currently working for "NO"
  }
  
  
  # Remove rows where bodytype is "NA", converte name and type to factors, join data from left and right side
  if (pairedRegion) {
    
    roiR_Connect <- cleanup(roiR_Connect)
    roiL_Connect <- cleanup(roiL_Connect)
    
    #rename ROI columns 
    roiR_Connect = renameRoiColumn(roiR_Connect, paste(slctROI, "(R)", sep=""),"ROIR")
    roiL_Connect = renameRoiColumn(roiL_Connect, paste(slctROI, "(L)", sep=""),"ROIL")
    
    # join data from paired ROI
    roi_Connect <- full_join(roiR_Connect,roiL_Connect, by = c("bodyid", "bodyname", "bodytype", "neuronStatus", "size", "npre", "npost"))
    
    # get full synaptic count
    roi_Connect = roi_Connect %>% rowwise() %>%  mutate(ROI_pre = sum(ROIR_pre, ROIL_pre, na.rm = TRUE),  
                                                        ROI_post = sum(ROIR_post, ROIL_post, na.rm = TRUE)) %>% ungroup()
    
  } else {
    roi_Connect <- cleanup(roi_Connect)
    roi_Connect = renameRoiColumn(roi_Connect, slctROI, "ROI")
  }

  roi_Connect[is.na(roi_Connect)] <- 0
  
  return(roi_Connect)
}




SynapseStats_And_Threshold <- function(NamedBodies, SaveDir, PreMaxThresh, PostMaxThresh, Name) {

# Get total synapses (Pre+Post) and the relative percent of Pre and Post synapses in the ROI. ---------------
  
  # Compute the total number of synapses in the ROI
  NamedBodies=mutate(NamedBodies, PrePost_Syns = (NamedBodies$ROI_pre + NamedBodies$ROI_post))
  
  # Compute the percent of each body's PRE synapses in the ROI. 
  NamedBodies=mutate(NamedBodies, Pre_Relative = NamedBodies$ROI_pre/unlist(NamedBodies$npre))
  NamedBodies$Pre_Relative[is.nan(NamedBodies$Pre_Relative)] = 0
  
  # Compute the percent of each body's POST synapses in the ROI.
  NamedBodies=mutate(NamedBodies, Post_Relative = NamedBodies$ROI_post/unlist(NamedBodies$npost) )
  NamedBodies$Post_Relative[is.nan(NamedBodies$Post_Relative)] = 0

  # Average the relative PRE and relative POST synapses. It's best to average after computing the relative PRE and relative POST
  # so that the measure isn't biased by there being more post synapses than pre synapses on average. 
  NamedBodies=mutate(NamedBodies, PrePost_Relative = (Pre_Relative+Post_Relative)/2)
  

# Get Pre-to-Post connection table. ---------------------------------------

  
  # Pre-to-Post connections for named bodies (return all segements and filter) --> subet at end is important or it returns bodyids not contained in NamedBodies
  NamedConnectTable = neuprint_connection_table(NamedBodies$bodyid,"PRE",ROI, all_segments=FALSE)   %>%
    mutate( PreName =neuprint_get_meta(bodyid)[["name"]], PostName = neuprint_get_meta(partner)[["name"]] )   %>% 
    mutate( PreType =neuprint_get_meta(bodyid)[["type"]], PostType = neuprint_get_meta(partner)[["type"]])   %>%  
    group_by(PreName) %>% 
    ungroup() %>%
    rename(Pre_bodyid = bodyid, Post_bodyid = partner)
  NamedConnectTable = subset(NamedConnectTable, Pre_bodyid %in% as.numeric(NamedBodies$bodyid) & Post_bodyid %in% as.numeric(NamedBodies$bodyid) ) 
  
  
# Calculate the max PRE weight and max POST weight for each bodyID. --------


  # Loop over pre-synaptic neurons and find max output weight
  PreMaxWeights=numeric()
  for (prepre in (1:length(NamedBodies$bodyid)) )
  {
    Temp=subset(NamedConnectTable, NamedConnectTable$Pre_bodyid==unlist(NamedBodies$bodyid[[prepre]]) )
    if (length(Temp$Pre_bodyid)==0){PreMaxWeights[prepre]=0
    }else{PreMaxWeights[prepre]=max(Temp$weight)}
  }
  
  # Loop over post-synaptic neurons and find max output weight
  PostMaxWeights=numeric()
  for (postpost in (1:length(NamedBodies$bodyid)) )
  {
    Temp=subset(NamedConnectTable, NamedConnectTable$Post_bodyid==unlist(NamedBodies$bodyid[[postpost]]) )
    if (length(Temp$Post_bodyid)==0){PostMaxWeights[postpost]=0
    }else{PostMaxWeights[postpost]=max(Temp$weight)}
  }
  remove(Temp)
  
  # Add data to NamedBodies dataframe
  NamedBodies=mutate(NamedBodies, PreMax = PreMaxWeights)
  NamedBodies=mutate(NamedBodies, PostMax = PostMaxWeights)
  NamedBodies=mutate(NamedBodies, PrePost_Max = (PreMaxWeights+PostMaxWeights))


# Get the average synapse statistics by neuron type -----------------------

  
  # Get the average number of Pre and Post synapses for each neuron type for both relative and absolute synapse counts.
  NamedBodies_byType <- aggregate(x=NamedBodies[c("ROI_pre","ROI_post","PrePost_Syns","Pre_Relative",
                                                  "Post_Relative","PrePost_Relative","PreMax","PostMax","PrePost_Max")],by=list(unlist(NamedBodies$bodytype[])), FUN=mean)
  
  names(NamedBodies_byType) <- c("Type","PreSyns", "PostSyns","PrePost_Syns","Pre_Relative","Post_Relative","PrePost_Relative","PreMax","PostMax","PrePost_Max")
  
  
# Make scatter plot of max pre and post synaptic weight. ------------------
  
  # Interactive plot of maximum pre vs post synapse weight for named neuron types, where cursor brings up neuron type name
  f <- list(size = 18, color = "black")
  p3a <- plot_ly(type = 'scatter', x = NamedBodies_byType$PreMax, y = NamedBodies_byType$PostMax, mode='markers',
                 text = paste("Type: ", NamedBodies_byType$Type )) %>% 
    layout(xaxis = list(title = paste("Max PRE Weight in", ROI),titlefont = f), yaxis = list(title = paste("Max POST Weight in",ROI),titlefont = f), 
           title=" \n \n  Each point is a neuron type")
  htmlwidgets::saveWidget(as_widget(p3a), paste(SaveDir,Name,"Pre_Vs_Post_MaxWeightsPerType_widget.html",sep=""))
  remove(p3a)
  
  # Interactive plot of maximum pre+post weight VS relative pre+post weight
  f <- list(size = 18, color = "black")
  p3b <- plot_ly(type = 'scatter', x = NamedBodies_byType$PrePost_Max, y = NamedBodies_byType$PrePost_Relative, mode='markers',
                 text = paste("Type: ", NamedBodies_byType$Type )) %>% 
    layout(xaxis = list(title = paste("Max PRE+POST Weight in", ROI),titlefont = f), yaxis = list(title = paste("Relative Pre+POST synapses in",ROI),titlefont = f), 
           title=" \n \n  Each point is a neuron type")
  htmlwidgets::saveWidget(as_widget(p3b), paste(SaveDir,Name,"PrePostMax_Vs_PrePostRelativeSyns_PerType.html",sep=""))
  remove(p3b,f)

# Exclude neurons that don't meet the criteria ----------------------------
 
  
  # First subset the type data frame
  GoodNeuronTypes = NamedBodies_byType$PreMax > PreMaxThresh | NamedBodies_byType$PostMax > PostMaxThresh 
  ExcludedNeurons_byType=subset(NamedBodies_byType,!GoodNeuronTypes)
  NamedBodies_byType=subset(NamedBodies_byType,GoodNeuronTypes)
  
  # Now subset the neurons themselves based on type
  GoodNeurons = NamedBodies$bodytype %in% NamedBodies_byType$Type
  ExcludedNeurons=subset(NamedBodies,!GoodNeurons)
  NamedBodies=subset(NamedBodies,GoodNeurons)
  
  
# Return variables of interest --------------------------------------------

  Output= list("NamedBodies"=NamedBodies,"ExcludedNeurons"=ExcludedNeurons)
  return(Output)
  
  
  
}

