################### A collection of functions for plotting FB connectivity information ##################


Assign_FBcol_PBglom <- function(DF, name_field, side_field, PBglom_field, FBcol_field, body_field){
  

  # Assign side (pb or fb column half)
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
  
  # Assign FB columns of FX, PFX, and Delta0
  DF[[FBcol_field]]=NA
  DF[[FBcol_field]]=str_extract(DF[[name_field]], "_C(\\d+)")
  DF[[FBcol_field]]=sapply(DF[[FBcol_field]], substring, 2, 3)
  
  # Assign column to Delta6s
  Delta6_Temp=subset(DF, startsWith(as.character(DF[[name_field]]), "Delta6"))[c(body_field,name_field,FBcol_field)]
  if (length(Delta6_Temp[[name_field]])>0){
    Delta6_Temp[[FBcol_field]]=str_extract(as.character(Delta6_Temp[[name_field]]), "_(\\d+)")
    Delta6_Temp[[FBcol_field]]=paste("DC",sapply(as.character(Delta6_Temp[[FBcol_field]]), substring, 2, 3),sep="")
    Delta6_Temp=distinct(Delta6_Temp)
    for (nnn in 1:length(Delta6_Temp[[body_field]])){
      DF[[FBcol_field]][DF[[body_field]] == Delta6_Temp[[body_field]][nnn]]=as.character(Delta6_Temp[[FBcol_field]][nnn])
    }
  }
  return(DF)
  
}


Plot_Connectivity_Mappings <- function(PFX_Types, PFX_FB_Outputs, PFX_FB_Inputs, PFX_Dir_Connectivity){
  
  
  # Loop over each PFX neuron type and plot upstream and downstream connections, arranged by PBglom and FBcol
  for (nnn in 1:length(PFX_Types)){
    Temp_Type=PFX_Types[nnn]
    
    ###### Output plots #########################################################################################################################
    # Output table for this neuron
    Temp_OutputTable=subset(PFX_FB_Outputs, databaseType.from == Temp_Type)
    
    # Filter output table by neuron types that innervate all columns
    Temp_GoodTypes_Output=Temp_OutputTable %>% group_by(FBcol.to, type.to) %>% summarize(n=n()) %>% group_by(type.to) %>% summarize(NumOfCol=n())
    Temp_GoodTypes_Output=subset(Temp_GoodTypes_Output, NumOfCol>8 )
    Temp_OutputTable=subset(Temp_OutputTable, type.to %in% Temp_GoodTypes_Output$type.to)
    
    # Compute mean by input/output cell type and input PBglom and output FBcol
    Temp_OutputTable_Sum=Temp_OutputTable %>% group_by(type.from, type.to, PBglom.from, FBcol.to) %>% summarise(weight=mean(weight))
    
    # Make output plots  
    if (length(Temp_OutputTable$weight)>0){
      
      # Plot all matrices
      OutputsPlot = plotConnectivityMatrix(Temp_OutputTable_Sum, byGroup = "Glom_to_Col", connectionMeasure = "weight")
      OutputsPlot = OutputsPlot + facet_grid(reorder(type.from, desc(type.from)) ~ type.to, switch="both") + 
        ggtitle(paste("Outputs from ", Temp_Type, " in FB")) +  coord_fixed(ratio = 1) 
      ggsave(paste(PFX_Dir_Connectivity, Temp_Type,"_Outputs_inFB.png",sep=""),
             plot = OutputsPlot, device='png', scale = 1, width = 10, height = 10, units ="in", dpi = 500, limitsize = TRUE)
      
      # Plot matrices one-by-one for each downstream cell type
      Temp_Downstream_Types=unique(Temp_OutputTable_Sum$type.to)
      for (ttt in 1:length(Temp_Downstream_Types)){
        Temp_OutputTable_Sum_Subset=subset(Temp_OutputTable_Sum, type.to == Temp_Downstream_Types[ttt])
        
        OutputsPlot_Sub = plotConnectivityMatrix(Temp_OutputTable_Sum_Subset, byGroup = "Glom_to_Col", connectionMeasure = "weight")
        OutputsPlot_Sub = OutputsPlot_Sub + facet_grid(reorder(type.from, desc(type.from)) ~ type.to, switch="both") + 
          ggtitle(paste("Outputs from ", Temp_Type, " in FB")) +  coord_fixed(ratio = 1) 
        ggsave(paste(PFX_Dir_Connectivity, Temp_Type,"_OutoutTo_", Temp_Downstream_Types[ttt], "_inFB.png",sep=""),
               plot = OutputsPlot_Sub, device='png', scale = 1, width = 4, height = 5, units ="in", dpi = 500, limitsize = TRUE)
      }
      
    }
    
    ###### Input plots #########################################################################################################################
    # Inpute table for this neuron
    Temp_InputTable=subset(PFX_FB_Inputs, databaseType.to == Temp_Type)
    Temp_InputTable=subset(Temp_InputTable, type.from %in% Columnar_Types) 
    
    # Filter input table by neuron types that innervate all columns
    Temp_GoodTypes_Input=Temp_InputTable %>% group_by(FBcol.from, type.from) %>% summarize(n=n()) %>% group_by(type.from) %>% summarize(NumOfCol=n())
    Temp_GoodTypes_Input=subset(Temp_GoodTypes_Input, NumOfCol>8)
    Temp_InputTable=subset(Temp_InputTable, type.from %in% Temp_GoodTypes_Input$type.from)
    
    # Compute mean by input/output cell type and input PBglom and output FBcol
    Temp_InputTable_Sum=Temp_InputTable %>% group_by(type.from, type.to, PBglom.to, FBcol.from) %>% summarise(weight=mean(weight))
    
    # Make input plots  
    if (length(Temp_InputTable$weight)>0){
      InputsPlot = plotConnectivityMatrix(Temp_InputTable_Sum, byGroup = "Col_to_Glom", connectionMeasure = "weight")
      InputsPlot = InputsPlot + facet_grid(reorder(type.from, desc(type.from)) ~ type.to, switch="both") +
        ggtitle(paste("Inputs to ", Temp_Type, " in FB")) +  coord_fixed(ratio = 1) 
      ggsave(paste(PFX_Dir_Connectivity, Temp_Type,"_Inputs_inFB.png",sep=""),
             plot = InputsPlot, device='png', scale = 1, width = 10, height = 10, units ="in", dpi = 500, limitsize = TRUE)
      
      # Plot matrices one-by-one for each upstream cell type
      Temp_Upstream_Types=unique(Temp_InputTable_Sum$type.from)
      for (ttt in 1:length(Temp_Upstream_Types)){
        Temp_InputTable_Sum_Subset=subset(Temp_InputTable_Sum, type.from == Temp_Upstream_Types[ttt])
        
        InputsPlot_Sub = plotConnectivityMatrix(Temp_InputTable_Sum_Subset, byGroup = "Col_to_Glom", connectionMeasure = "weight")
        InputsPlot_Sub = InputsPlot_Sub + facet_grid(reorder(type.from, desc(type.from)) ~ type.to, switch="both") + 
          ggtitle(paste("Inputs to ", Temp_Type, " in FB")) +  coord_fixed(ratio = 1) 
        ggsave(paste(PFX_Dir_Connectivity, Temp_Type,"_InputFrom_", Temp_Upstream_Types[ttt], "_inFB.png",sep=""),
               plot = InputsPlot_Sub, device='png', scale = 1, width = 7.5, height = 5, units ="in", dpi = 500, limitsize = TRUE)
      }
      
    }
    
  }
  
  
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
  
  
  
  
  
  
  
  
  
  
  
  
  


