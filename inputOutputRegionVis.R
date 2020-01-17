# Collections of functions for visualization of input/output regions per neuron type


get_rois_df <- function (){
  rois = neuprint_ROIs()
  return (data.frame(keyName=names(rois), value=rois, row.names=NULL))
}

renameRoiColumn <- function(df, slctROI) {
  df = df %>% rename_at(vars(starts_with(slctROI)), funs(str_replace(., slctROI, "ROI")))
  df = df %>%rename_at(vars(ends_with(".pre")), funs(str_replace(., ".pre", "_pre")))
  df = df %>%rename_at(vars(ends_with(".post")), funs(str_replace(., ".post", "_post")))
  df = df %>%rename_all(funs(str_replace(., "\\(", "")))
  df = df %>%rename_all(funs(str_replace(., "\\)", "")))
  
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
  roi_synStats = df %>% group_by(bodytype) %>% summarise(mean_pre = mean(ROI_pre),
                                                         mean_post = mean(ROI_post),
                                                         sum_pre = sum(ROI_pre),
                                                         sum_post = sum(ROI_post),
                                                         fract_pre = mean(ROI_pre/npre),
                                                         fract_post = mean(ROI_post/npost),
                                                         fract_all = mean((ROI_pre+ROI_post)/(npost + npre)))
  
  df = full_join(full_join(filter(df, bodytype %in% filter(roi_synStats, mean_pre >= minSynapses)$bodytype),
                           filter(df, bodytype %in% filter(roi_synStats, mean_post >= minSynapses)$bodytype)),
                 full_join(filter(df, bodytype %in% filter(roi_synStats, sum_pre >= minSumSynapses)$bodytype),
                           filter(df, bodytype %in% filter(roi_synStats, sum_post >= minSumSynapses)$bodytype)))
  
  return(df)
}