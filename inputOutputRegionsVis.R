# Collections of functions for visualization of input/output regions per neuron type


get_rois_df <- function (){
  rois = neuprint_ROIs()
  return (data.frame(keyName=names(rois), value=rois, row.names=NULL))
}


renameRoiColumn <- function(df, slctROI, NewName) {
  names(df)[names(df) == paste(slctROI,".pre",sep="")] <- paste(NewName, "_pre",sep="")
  names(df)[names(df) == paste(slctROI,".post",sep="")] <- paste(NewName, "_post",sep="")
  df = df %>%rename_all(funs(str_replace(., "\\(", "")))
  df = df %>%rename_all(funs(str_replace(., "\\)", "")))
  df[[paste(NewName, "_pre",sep="")]]=as.numeric(df[[paste(NewName, "_pre",sep="")]])
  df[[paste(NewName, "_post",sep="")]]=as.numeric(df[[paste(NewName, "_post",sep="")]])
  df$npre=as.numeric(unlist(df$npre))
  df$npost=as.numeric(unlist(df$npost))
  df$bodyname=as.character(df$bodyname)
  df$bodytype=as.character(df$bodytype)
  return (df)
}


#remove missing values
cleanup <- function(df) {
  df <- unnest(df, colnames(df))
  df <- na.omit(df)
  df$bodyname <- as.factor(df$bodyname)
  df$bodytype <- as.factor(df$bodytype)
  
  return (df)
}


# Select neurons that are "significant" based on synapse statistics. This function should be adapted to common criterial we agreed on.
filterNeuronTypes = function(df,minSynapses,minSumSynapses){
  
  # -> Compute statistics across neurons of the same bodytype
  roi_synStats = df %>% group_by(bodytype) %>% 
    mutate(ROI_post = as.numeric(ROI_post), ROI_pre = as.numeric(ROI_pre), npre = as.numeric(npre), npost = as.numeric(npost)) %>%
    summarise(mean_pre = mean(ROI_pre,na.rm = TRUE),
              mean_post = mean(ROI_post,na.rm = TRUE),
              sum_pre = sum(ROI_pre,na.rm = TRUE),
              sum_post = sum(ROI_post,na.rm = TRUE),
              fract_pre = mean(ROI_pre/npre,na.rm = TRUE),
              fract_post = mean(ROI_post/npost,na.rm = TRUE),
              fract_all = mean((ROI_pre+ROI_post)/(npost + npre)),na.rm = TRUE)
  
  df = full_join(full_join(filter(df, bodytype %in% filter(roi_synStats, mean_pre >= minSynapses)$bodytype),
                           filter(df, bodytype %in% filter(roi_synStats, mean_post >= minSynapses)$bodytype)),
                 full_join(filter(df, bodytype %in% filter(roi_synStats, sum_pre >= minSumSynapses)$bodytype),
                           filter(df, bodytype %in% filter(roi_synStats, sum_post >= minSumSynapses)$bodytype)))
  
  return(df)
}






### Return a dataframe listing neuron types and synapse counts per ROI
InputOutput_ROI_PerNeuron = function(NamedBodies, minSynapses){
  

# For each neuron, get the number of pre- and post-synapses in the --------
  
  roi_Connect = neuprint_get_roiInfo(NamedBodies$bodyid)
  
  roi_Connect =  mutate(roi_Connect,
                        name = neuprint_get_meta(bodyid)$name,
                        type = neuprint_get_meta(bodyid)$type,
                        neuron_hemi = "C",
                        pre = neuprint_get_meta(bodyid)$pre,
                        post = neuprint_get_meta(bodyid)$post,
                        cropped = neuprint_get_meta(bodyid)$cropped) 
  
  # rearrange so that information worth viewing is in first few columns (HARD CODING OF COLULMNS)
  roi_Connect = roi_Connect[ , c(1,(ncol(roi_Connect)-5):ncol(roi_Connect),2:(ncol(roi_Connect)-6)) ]
  
  # Assign neuron_hemi to Left (L), Right (R), or Center/NA (C)
  roi_Connect$neuron_hemi[grepl("_R",as.character(roi_Connect$name))]='R'
  roi_Connect$neuron_hemi[grepl("_L",as.character(roi_Connect$name))]='L'
  
  
# Get neuron type averages, keeping Left/Right/Center neurons sepa --------

  # Set all NAs to zeros before averaging
  roi_Connect[is.na((roi_Connect))] <- 0
  
  # Set cropped to 1 if TRUE and 0 if FALSE, so they can be averaged
  roi_Connect$cropped[roi_Connect$cropped==TRUE]=1
  
  # Average by type (HARD CODING OF COLUMNS)
  roi_Connect_bytype=aggregate(roi_Connect[ ,5:length(colnames(roi_Connect))], by=list(unlist(roi_Connect$type[]), unlist(roi_Connect$neuron_hemi)), FUN=mean)
  colnames(roi_Connect_bytype)[c(1,2)] <- c("type","neuron_hemi")
  
  

# For the single neurons and type-averages, parse the ROIs and reorganize such that ROI, #pre, and #post become variables  --------

  # Define a function that can do this for both single neuron (roi_Connect) and type-averages (roi_Connect_bytype) 
  # COL_START is column in the data frame where the ROI columns begin
  ParseROIConnectionDataframe = function(roi_Connect_bytype,COL_START){
    
    roi_Connect_bytype_df = data.frame(type = character(),
                                       roi = character(),
                                       fullroi= character(),
                                       neuron_hemi=character(),
                                       roi_hemi=character(),
                                       prepost = character(),
                                       count = numeric() )
    
    
    # Loop over columns, parse ROI info, and insert into new dataframe (HARD CODING OF COLUMNS)
    for (roi in as.character(colnames(roi_Connect_bytype))[COL_START:(length(colnames(roi_Connect_bytype)))] ) {
      
      roiname = as.character(unlist(strsplit(roi, "[.]"))[1])
      
      # Get hemisphere- this assumes that ROI names with"(R" ALWAYS MEANS right and "(L" MEANS left. 
      if (grepl("(R", roi,fixed=TRUE )){
        roi_hemi="R"
        roiname_nohemi=as.character(unlist(strsplit(roiname, "[(]R[)]"))[1])
      } else if ( grepl("(L", roi,fixed=TRUE ) ) {
        roi_hemi="L"
        roiname_nohemi=as.character(unlist(strsplit(roiname, "[(]L[)]"))[1])
      } else {roi_hemi='C'}
      
      
      Temp_df = data.frame(type = as.character(roi_Connect_bytype$type),
                           roi = roiname_nohemi,
                           fullroi = as.character(unlist(roi)),
                           neuron_hemi =  as.character(roi_Connect_bytype$neuron_hemi),
                           roi_hemi = as.character(roi_hemi),
                           prepost = as.character(unlist(strsplit(roi, "[.]"))[2]),
                           count = roi_Connect_bytype[[roi]]  ) 
      
      roi_Connect_bytype_df = rbind(roi_Connect_bytype_df, Temp_df)
      remove(list=c("Temp_df","roi_hemi"))
    }
    
    return(roi_Connect_bytype_df)
  }
  
  # Run function of both roi_Connect and roi_Connect_bytype
  roi_Connect_bytype_df= ParseROIConnectionDataframe(roi_Connect_bytype,6)
  roi_Connect_df= ParseROIConnectionDataframe(roi_Connect,8)
  

# Filter out ROIs with low synapse counts ---------------------------------

  HistDat=data.frame(Type=roi_Connect_bytype_df$type, NumberOfSyns= roi_Connect_bytype_df$count)
  p1<-ggplot(HistDat, aes(x=NumberOfSyns)) + geom_histogram(binwidth=10) + xlim(0, 500) + ylim(0, length(HistDat$NumberOfSyns[HistDat$NumberOfSyns>2 & HistDat$NumberOfSyns<=11]) )
  print(p1)

  roi_Connect_bytype_df = filter(roi_Connect_bytype_df, count > minSynapses)
  roi_Connect_df = filter(roi_Connect_df, count > minSynapses)
  

# Plot input and output of all ROIs, both Left/Right/C --------------------

  # Convert some factor columns to character columns
  roi_Connect_bytype_df$type=as.character(roi_Connect_bytype_df$type)
  roi_Connect_bytype_df$roi=as.character(roi_Connect_bytype_df$roi)
  roi_Connect_bytype_df$prepost=as.character(roi_Connect_bytype_df$prepost)
  roi_Connect_bytype_df$neuron_hemi=as.character(roi_Connect_bytype_df$neuron_hemi)
  roi_Connect_bytype_df$roi_hemi=as.character(roi_Connect_bytype_df$roi_hemi)
  
  roi_Connect_df$type=as.character(roi_Connect_df$type)
  roi_Connect_df$roi=as.character(roi_Connect_df$roi)
  roi_Connect_df$prepost=as.character(roi_Connect_df$prepost)
  roi_Connect_df$neuron_hemi=as.character(roi_Connect_df$neuron_hemi)
  roi_Connect_df$roi_hemi=as.character(roi_Connect_df$roi_hemi)
  
  # Rename post and pre and input and output
  roi_Connect_bytype_df$prepost[as.character(roi_Connect_bytype_df$prepost)=="post"]="Input ROI"
  roi_Connect_bytype_df$prepost[as.character(roi_Connect_bytype_df$prepost)=="pre"]="Output ROI"
  
  roi_Connect_df$prepost[as.character(roi_Connect_df$prepost)=="post"]="Input ROI"
  roi_Connect_df$prepost[as.character(roi_Connect_df$prepost)=="pre"]="Output ROI"
  
  # Plot neuron type ROI innervation pattern
  p2<-ggplot(roi_Connect_bytype_df, aes(x=roi, y=type, size=count)) + 
    geom_point(aes(color=roi)) + facet_grid(cols=vars(prepost)) + 
    theme_bw() + theme(axis.text.x = element_text(angle = 90))  + guides(fill=FALSE)
  print(p2)
  

# Return variables of interest --------------------------------------------

  Output= list("roi_Connect_bytype_df"=roi_Connect_bytype_df,"roi_Connect_df"=roi_Connect_df)  
  
  return(Output)
}




















