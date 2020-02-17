# Collections of functions for visualization of input/output regions per neuron type

Get_SynapseCount_In_Rois <- function (ROI, Paired){

  
  ### Get neurons in the ROI and compute synapse statistics
  NamedBodies = Get_AllNeurons_InRoi(ROI, Paired)
  NamedBodies = ComputeRelativeSynapseCounts(NamedBodies)
  # NamedBodies = ComputeMaxWeights(NamedBodies)
  
  
  ### Compute type averages of named bodies to get an idea of which neurons belong in which compartments
  ### and where threshold on relative synapse counts should be
  ### Note, this averaging should not be done in cases where every instance of that neuron type is not thought to innervate the ROI.
  ### For example, if the ROI were glomerulus R4 in the PB, averaging across all PENs would dramatically underestimate the average synapse count in that ROI. 
  ### Similarly, care should be taken when averaging neurons that innervate the left or right portion of an ROI (like the PENs for the NO(R) and NO(L) )
  NamedBodies_TypeAverage= aggregate(NamedBodies[, 6:length(colnames(NamedBodies))], list(NamedBodies$bodytype), mean) ### try mean and median here
  colnames(NamedBodies_TypeAverage)[1]="bodytype"
  
  
  ### Get the number of synapses in each ROI for individual neurons and individual neuron types
  Output=InputOutput_ROI_PerNeuron(NamedBodies)
  list2env(Output ,.GlobalEnv)
  remove(Output)
 
  
  ### Add ROI and neuron type names that will be displayed on plots
  roi_Connect_bytype_df$roi_names = paste(roi_Connect_bytype_df$roi, '(', roi_Connect_bytype_df$roi_hemi, ')', sep="")
  roi_Connect_bytype_df$type_name = paste(roi_Connect_bytype_df$type, '(',roi_Connect_bytype_df$neuron_hemi, ')', sep="")
  
  
  Output= list("roi_Connect_bytype_df"=roi_Connect_bytype_df,"roi_Connect_df"=roi_Connect_df,"NamedBodies"=NamedBodies,"NamedBodies_TypeAverage"=NamedBodies_TypeAverage )  
  return(Output)
  
}


### Return a dataframe listing neuron types and synapse counts per ROI
InputOutput_ROI_PerNeuron = function(NamedBodies){
  
  
  # For each neuron, get the number of pre- and post-synapses in every ROI --------
  
  roi_Connect = neuprint_get_roiInfo(NamedBodies$bodyid)
  
  roi_Connect =  mutate(roi_Connect,
                        name = neuprint_get_meta(bodyid)$name,
                        type = neuprint_get_meta(bodyid)$type,
                        neuron_hemi = "C",
                        pre = neuprint_get_meta(bodyid)$pre,
                        post = neuprint_get_meta(bodyid)$post,
                        cropped = neuprint_get_meta(bodyid)$cropped)# Leaving this in for now for future when cropped, or similar field, might mean something.
  
  # rearrange so that information worth viewing is in first few columns (HARD CODING OF COLULMNS)
  roi_Connect = roi_Connect[ , c(1,(ncol(roi_Connect)-5):ncol(roi_Connect),2:(ncol(roi_Connect)-6)) ]
  
  # Assign neuron_hemi to Left (L), Right (R), or Center/NA (C)
  roi_Connect$neuron_hemi[grepl("_R",as.character(roi_Connect$name))]='R'
  roi_Connect$neuron_hemi[grepl("_L",as.character(roi_Connect$name))]='L'
  
  
  # Get neuron type averages, keeping Left/Right/Center neurons separate --------
  
  # Set all NAs to zeros before averaging (Note, need to be careful here. 
  # We are setting to 0 the number of synapses for some neurons of a type that don't
  # innervate sub ROIs, then averaging across type. This will be delt with right below.)
  roi_Connect[is.na((roi_Connect))] <- 0 
  
  
  # Define a list of PB columnar neurons (or any neuron that innervates only a portion of the PB) and set the synapse counts to NA
  # for any PB glomerulus where there are no synapses. This is to ensure that, when we average across neuron types, we don't underestimate
  # the number of synapses in the PB glomeruli. 
  ROI_Names_Temp=colnames(roi_Connect)
  PB_Columns=which(grepl("[[:digit:]]",ROI_Names_Temp) & grepl("PB[(]",ROI_Names_Temp))
  PB_ColumnarNeuron_Rows=which(grepl("EP",roi_Connect$type) | grepl("PE",roi_Connect$type) |  grepl("PF",roi_Connect$type) |
                               grepl("P6",roi_Connect$type) |  grepl("IbSps",roi_Connect$type))

  
  # Loop over PB columnar neurons (rows) and PB glomeruli ROI (columns) and set to NAN any 0 entries
  for (rrr in PB_ColumnarNeuron_Rows){
    for (ccc in PB_Columns){
      if (roi_Connect[rrr, ccc] == 0 ){
        roi_Connect[rrr, ccc] = NA
      }
    }
  }

  # Make new ROIs for PB(L) and PB(R) pre and post counts
  PB_Columns_L_Pre  =which(grepl("[[:digit:]]",ROI_Names_Temp) & grepl("PB[(]L",ROI_Names_Temp) & grepl("pre",ROI_Names_Temp))
  PB_Columns_L_Post =which(grepl("[[:digit:]]",ROI_Names_Temp) & grepl("PB[(]L",ROI_Names_Temp) & grepl("post",ROI_Names_Temp))
  PB_Columns_R_Pre  =which(grepl("[[:digit:]]",ROI_Names_Temp) & grepl("PB[(]R",ROI_Names_Temp) & grepl("pre",ROI_Names_Temp))
  PB_Columns_R_Post =which(grepl("[[:digit:]]",ROI_Names_Temp) & grepl("PB[(]R",ROI_Names_Temp) & grepl("post",ROI_Names_Temp))  

  
  # Populate new ROIs by NA summing across PB glomeruli for L/R and pre/post separately
  roi_Connect = mutate(roi_Connect, "PB(L).pre" =  rowSums(roi_Connect[,PB_Columns_L_Pre], na.rm=TRUE),
                                    "PB(L).post" =  rowSums(roi_Connect[,PB_Columns_L_Post], na.rm=TRUE),
                                    "PB(R).pre" =  rowSums(roi_Connect[,PB_Columns_R_Pre], na.rm=TRUE),
                                    "PB(R).post" =  rowSums(roi_Connect[,PB_Columns_R_Post], na.rm=TRUE))
  
  
  # Set cropped to 1 if TRUE and 0 if FALSE, so they can be averaged
  roi_Connect$cropped[roi_Connect$cropped==TRUE]=1
  
  
  # NA Average by neuron type, keeping L/R neurons separate (HARD CODING OF COLUMNS)
  roi_Connect_bytype=aggregate(roi_Connect[ ,5:length(colnames(roi_Connect))], by=list(unlist(roi_Connect$type[]), unlist(roi_Connect$neuron_hemi)), FUN=mean, na.rm = TRUE) ### try mean and median here
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
      
      # Get hemisphere- this assumes that ROI names with"(R" ALWAYS MEANS right, "(L" MEANS left, and nothing means Center 
      if (grepl("(R", roi,fixed=TRUE )){
        roi_hemi="R"
        roiname_nohemi=as.character(unlist(strsplit(roiname, "[(]R[)]"))[1])
      } else if ( grepl("(L", roi,fixed=TRUE ) ) {
        roi_hemi="L"
        roiname_nohemi=as.character(unlist(strsplit(roiname, "[(]L[)]"))[1])
      } else {
        roi_hemi='C'
        roiname_nohemi=as.character(unlist(roiname))
      }
      
      
      Temp_df = data.frame(type = as.character(roi_Connect_bytype$type),
                           roi = roiname_nohemi,
                           fullroi = as.character(unlist(roi)),
                           neuron_hemi =  as.character(roi_Connect_bytype$neuron_hemi),
                           roi_hemi = as.character(roi_hemi),
                           prepost = as.character(unlist(strsplit(roi, "[.]"))[2]),
                           count = roi_Connect_bytype[[roi]]  ) 
      
      roi_Connect_bytype_df = rbind(roi_Connect_bytype_df, Temp_df)
      remove(list=c("Temp_df","roi_hemi","roiname"))
    }
    
    return(roi_Connect_bytype_df)
  }
  
  # Run function of both roi_Connect and roi_Connect_bytype
  roi_Connect_bytype_df= ParseROIConnectionDataframe(roi_Connect_bytype,6)
  roi_Connect_df= ParseROIConnectionDataframe(roi_Connect,8)
  
  # Coerce data type before output --------------------
  
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
  
  # Rename post and pre as input and output
  roi_Connect_bytype_df$prepost[as.character(roi_Connect_bytype_df$prepost)=="post"]="Input ROI"
  roi_Connect_bytype_df$prepost[as.character(roi_Connect_bytype_df$prepost)=="pre"]="Output ROI"
  
  roi_Connect_df$prepost[as.character(roi_Connect_df$prepost)=="post"]="Input ROI"
  roi_Connect_df$prepost[as.character(roi_Connect_df$prepost)=="pre"]="Output ROI"
  
  
  # Return variables of interest --------------------------------------------
  
  Output= list("roi_Connect_bytype_df"=roi_Connect_bytype_df,"roi_Connect_df"=roi_Connect_df)  
  
  return(Output)
}


ROI_RelativeSynapses <- function (ROISubset){
  
  # Loop over each connection and compute relative Pre and relative Post counts
  ROISubset$RelativeCount <- NA
  for (nnn in 1:length(ROISubset$RelativeCount)){
    CellType<-as.character(ROISubset$type[nnn])
    PrePost_Temp=ROISubset$prepost[nnn]
    if (PrePost_Temp ==  "Input ROI"){
      ROISubset$RelativeCount[nnn] = ROISubset$count[nnn] / sum(ROISubset$count[ ROISubset$type==CellType & ROISubset$prepost==  "Input ROI"])
    } else if (PrePost_Temp ==  "Output ROI") {
      ROISubset$RelativeCount[nnn] = ROISubset$count[nnn] / sum(ROISubset$count[ ROISubset$type==CellType & ROISubset$prepost==  "Output ROI"])
    }
  }
  return(ROISubset)
}



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




