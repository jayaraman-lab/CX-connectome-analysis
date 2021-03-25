################### A collection of functions for plotting FB connectivity information ##################


Assign_FBcol_PBglom <- function(DF, name_field, side_field, PBglom_field, FBcol_field){
  #' Function for parsing neuron names to extract PB glomeruli and FB column info.

  # Assign side
  DF[[side_field]]=NA
  DF[[side_field]][grepl("_R",DF[[name_field]])]<-"R"
  DF[[side_field]][grepl("_L",DF[[name_field]])]<-"L"
  
  # Assign PB glomeruli
  DF[[PBglom_field]]=NA
  RightColumns=str_extract(DF[[name_field]], "_R(\\d+)")
  LeftColumns=str_extract(DF[[name_field]], "_L(\\d+)")
  DF[[PBglom_field]][which(!is.na(RightColumns))]=RightColumns[which(!is.na(RightColumns))]
  DF[[PBglom_field]][which(!is.na(LeftColumns))]=LeftColumns[which(!is.na(LeftColumns))]
  DF[[PBglom_field]]=sapply(DF[[PBglom_field]], substring, 2, 3)
  
  # Assign FB columns
  DF[[FBcol_field]]=NA
  DF[[FBcol_field]]=str_extract(DF[[name_field]], "_C(\\d+)")
  DF[[FBcol_field]]=sapply(DF[[FBcol_field]], substring, 2, 4)
  
  return(DF)
}




Plot_ColCol_Matix <- function(PFLn_Input_Network){
  
  P1 <- ggplot(PFLn_Input_Network) +
    scale_fill_gradient2(low="thistle", mid="blueviolet", high="black", 
                         midpoint =0.5*max(PFLn_Input_Network$weightRelative), 
                         limits=c(0,max(PFLn_Input_Network$weightRelative)*1)) +
    geom_tile(aes(to_name,from_name,fill=weightRelative)) +
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(), 
          strip.placement = "outside", strip.background = element_rect(fill=NA, colour="grey50")) +
    facet_grid(reorder(type.from, desc(type.from)) ~ type.to, space="free", scales="free",switch="both")
  
  return(P1)
}
  
  
  
  
  
  
  
  
  
  
  
  
  


