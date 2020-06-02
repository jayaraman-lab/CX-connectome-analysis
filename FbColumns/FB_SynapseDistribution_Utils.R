################### A collection of functions for loading, processing, and plotting FB synapse distributions ##################



###############################################################################################################################
################## Function for loading FB synapses by layer ##################################################################



Get_FBlayer_ColumnarSyns <- function(BodyIDs, ROI, Layer) {
  
  
  SynLocs =  neuprint_get_synapses(BodyIDs, roi = ROI)
  
  SynLocs =  mutate(SynLocs, type=neuprint_get_meta(bodyid)$type, partner_type=neuprint_get_meta(partner)$type,
                    x=as.numeric(x),y=as.numeric(y),z=as.numeric(z))
  
  SynLocs$prepost[SynLocs$prepost==0]="Output"
  SynLocs$prepost[SynLocs$prepost==1]="Input"
  
  SynLocs$Layer=Layer
  return(SynLocs)
}



###############################################################################################################################
################# Function for plotting synapse distribution for all FB layers ################################################



PlotSyns_byLayer <- function(Plot_Syns, Type, PlotName, Grouping, DIR, TOPLOT_TF, TOPLOTLEGEND_TF){
  
  P1 = Plot_Layer(Plot_Syns, "L1", Outline_L1_XZ, paste("L1 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF, "XZ")
  P2 = Plot_Layer(Plot_Syns, "L2", Outline_L2_XZ, paste("L2 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF, "XZ")
  P3 = Plot_Layer(Plot_Syns, "L3", Outline_L3_XZ, paste("L3 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF, "XZ")
  P4 = Plot_Layer(Plot_Syns, "L4", Outline_L4_XZ, paste("L4 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF, "XZ")
  P5 = Plot_Layer(Plot_Syns, "L5", Outline_L5_XZ, paste("L5 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF, "XZ")
  P6 = Plot_Layer(Plot_Syns, "L6", Outline_L6_XZ, paste("L6 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF, "XZ")
  P7 = Plot_Layer(Plot_Syns, "L7", Outline_L7_XZ, paste("L7 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF, "XZ")
  P8 = Plot_Layer(Plot_Syns, "L8", Outline_L8_XZ, paste("L8 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF, "XZ")
  P9 = Plot_Layer(Plot_Syns, "L9", Outline_L9_XZ, paste("L9 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF, "XZ")
  
  P12<-ggarrange(plotlist=list(P9,P8,P7,P6,P5,P4,P3,P2,P1), ncol =3, nrow = 3)
  ggsave(paste(DIR, PlotName,"_XZ.pdf",sep=""),
         plot = P12, device='pdf', scale = 1, width = 10, height = 12, units ="in", dpi = 500, limitsize = TRUE)
  ggsave(paste(DIR, PlotName,"_XZ.png",sep=""),
         plot = P12, device='png', scale = 1, width = 10, height = 12, units ="in", dpi = 500, limitsize = TRUE)
  
  P1b = Plot_Layer(Plot_Syns, "L1", Outline_L1_XY, paste("L1 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF, "XY")
  P2b = Plot_Layer(Plot_Syns, "L2", Outline_L2_XY, paste("L2 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF, "XY")
  P3b = Plot_Layer(Plot_Syns, "L3", Outline_L3_XY, paste("L3 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF, "XY")
  P4b = Plot_Layer(Plot_Syns, "L4", Outline_L4_XY, paste("L4 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF, "XY")
  P5b = Plot_Layer(Plot_Syns, "L5", Outline_L5_XY, paste("L5 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF, "XY")
  P6b = Plot_Layer(Plot_Syns, "L6", Outline_L6_XY, paste("L6 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF, "XY")
  P7b = Plot_Layer(Plot_Syns, "L7", Outline_L7_XY, paste("L7 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF, "XY")
  P8b = Plot_Layer(Plot_Syns, "L8", Outline_L8_XY, paste("L8 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF, "XY")
  P9b = Plot_Layer(Plot_Syns, "L9", Outline_L9_XY, paste("L9 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF, "XY")
  
  P12b<-ggarrange(plotlist=list(P9b,P8b,P7b,P6b,P5b,P4b,P3b,P2b,P1b), ncol =3, nrow = 3)
  ggsave(paste(DIR, PlotName,"_XY.pdf",sep=""),
         plot = P12b, device='pdf', scale = 1, width = 10, height = 12, units ="in", dpi = 500, limitsize = TRUE,  useDingbats=FALSE)
  ggsave(paste(DIR, PlotName,"_XY.pdf",sep=""),
         plot = P12b, device='pdf', scale = 1, width = 10, height = 12, units ="in", dpi = 500, limitsize = TRUE,  useDingbats=FALSE)
}



###############################################################################################################################
################# Function for plotting synapse distribution of each individual FB layer ######################################



Plot_Layer <- function(Plot_Syns, Layer_To_Plot, Outline, Title, Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF, View){
  
  # Get data
  PlotData=subset(Plot_Syns, Layer==Layer_To_Plot)
  
  # Get color map
  col_vector=GetColorMap(Grouping)
  
  # Get color variable
  if (Grouping == "FBcol"){
    ColorVar=PlotData$FBcol
  } else if (Grouping == "PBglom") {
    ColorVar=PlotData$PBglom
  }
  
  # Make plot
  P1 <- Plot_Synapse_Distribution(PlotData, View, Outline, ColorVar, col_vector, DIR, PlotName, Layer_To_Plot, TOPLOT_TF, TOPLOTLEGEND_TF)

  return(P1)
}



###############################################################################################################################
################# Generic plotting function for synapse distributions  ################################################



Plot_Synapse_Distribution <- function(PlotData, VIEW, Outline, ColorVar, col_vector, DIR, PlotName, PlotName2, TOPLOT_TF, TOPLOTLEGEND_TF, Size=2){
  
  # Get view
  if (VIEW=="XY"){
    PlotData = PlotData[!colnames(PlotData)=="Z"]
    colnames(PlotData)[colnames(PlotData)=="Y"]="Z"
    PlotData$Z=-PlotData$Z
  }
  
  # Plot just the data for now
  P1 <- ggplot() + geom_point(data=PlotData, aes(x=X, y=Z, color=ColorVar, shape=Side) ,  size=Size, alpha = 1, stroke = 0) + 
    geom_path(data=Outline, aes(x=c1, y=c2), size = 1) + xlim(c(-9000, 9000)) +  coord_fixed(ratio = 1) 

  # Save with scaled visible if TOPLOTLEGEND_TF us TRUE
  if (TOPLOTLEGEND_TF==TRUE){
    ggsave(paste(DIR, PlotName,"_",PlotName2 , "_Legend.pdf",sep=""), plot = P1, device='pdf', scale = 1, width = 8, height = 5, units ="in", dpi = 500, limitsize = TRUE) 
  }  
  
  # Format the rest of the plot, get rid of axes
  P1 <- P1 + theme_void() + guides(fill=FALSE) + theme(legend.position = "none") +
      scale_color_manual(values=col_vector, drop=FALSE)  +  theme(
      panel.background = element_rect(fill = "transparent", color = NA), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      panel.grid.major = element_blank(), # get rid of major grid
      panel.grid.minor = element_blank(), # get rid of minor grid
      legend.background = element_rect(fill = "transparent", color = NA), # get rid of legend bg
      legend.box.background = element_rect(fill = "transparent", color = NA) # get rid of legend panel bg
    )
  
  # Save plot
  if (TOPLOT_TF==TRUE){
    ggsave(paste(DIR, PlotName,"_",PlotName2, VIEW, ".pdf",sep=""), plot = P1, device='pdf', scale = 1, width = 8, height = 5, units ="in", dpi = 500, limitsize = TRUE,  bg = "transparent")
  } 
  
  return(P1)
}



###############################################################################################################################
################# Function for getting color map of FB columns or PB glomeruli ################################################



GetColorMap <- function(Grouping){
  if (Grouping == "FBcol_Old"){
    col_vector=brewer.pal(8, "Paired")
    col_vector[9]=col_vector[1]
    col_vector=col_vector[c(1,4,7,2,5,8,3,6,9)]
  } else if (Grouping == "FBcol"){
    col_vector=color(c("#CD4F39", "#7D26CD", "#63B8FF", "#FF1493", "#CD96CD", "#00868B", "#ADFF2F", "#0000FF", "#FFB90F"))
  } else if (Grouping == "FBcol_Function"){
    col_vector=color(c("#FFF688", "#FFA010", "#FF3610", "#FBB5DE", "#7D24E7", "#5C89C7", "#C7FFEF", "#81FB35", "#FFF688"))
  } else if (Grouping == "PBglom") {
    L_col_vector=brewer.pal(8, "Paired")
    L_col_vector[9]=L_col_vector[1]
    L_col_vector=L_col_vector[c(1,4,7,2,5,8,3,6,9)]
    R_col_vector=rev(L_col_vector)
    col_vector=c(L_col_vector, R_col_vector)
  }
  return(col_vector)
}



###############################################################################################################################
################# Function for plotting synapse distribution for each FB layer ################################################



PlotSyns_byColumn <- function(Plot_Syns, Title, PlotName, DIR, TOPLOT_TF, TOPLOTLEGEND_TF){
  
  # Get color vector
  col_vector=GetColorMap("FBcol")
  
  ### XY View
  P1 <- Plot_Synapse_Distribution(Plot_Syns, "XY", Outline_FB_XY, Plot_Syns$FBcol, col_vector, DIR, PlotName, "XY", TOPLOT_TF, TOPLOTLEGEND_TF)
  
  ### XZ View
  P1 <- Plot_Synapse_Distribution(Plot_Syns, "XZ", Outline_FB_XZ, Plot_Syns$FBcol, col_vector, DIR, PlotName, "XZ", TOPLOT_TF, TOPLOTLEGEND_TF)
  
}


###############################################################################################################################
####### Functions for plotting column positions as defined by Fx neuron synapses locations ####################################


PlotColumn_Types<- function(Plot_Syns, PlotName, DIR, PLOT, ID){
  Types=unique(Plot_Syns$type)
  for (nnn in 1:length(Types)){
    Temp_Type=subset(Plot_Syns, type == as.character(Types[nnn]))
    T_ID=ID[Plot_Syns$type == as.character(Types[nnn]) ]
    PlotColumn_Aves(Temp_Type, paste(PlotName, "_", Types[nnn],sep="") , DIR, PLOT, T_ID)
  }
}



PlotColumn_Aves <- function(Plot_Syns, PlotName, DIR, PLOT, T_ID){
  
  
  # get layers
  Layers_All=sort(unique(as.character(Plot_Syns$Layer)))
  
  # Make plots
  l1=Plot_Crosses(Plot_Syns, 1, DIR, PlotName, PLOT, T_ID)
  l2=Plot_Crosses(Plot_Syns, 2, DIR, PlotName, PLOT, T_ID)
  l3=Plot_Crosses(Plot_Syns, 3, DIR, PlotName, PLOT, T_ID)
  l4=Plot_Crosses(Plot_Syns, 4, DIR, PlotName, PLOT, T_ID)
  l5=Plot_Crosses(Plot_Syns, 5, DIR, PlotName, PLOT, T_ID)
  l6=Plot_Crosses(Plot_Syns, 6, DIR, PlotName, PLOT, T_ID)
  l7=Plot_Crosses(Plot_Syns, 7, DIR, PlotName, PLOT, T_ID)
  l8=Plot_Crosses(Plot_Syns, 8, DIR, PlotName, PLOT, T_ID)
  
  
  # Save plot
  PP=ggarrange(l8,l7,l6, l5, l4, l3, l2, l1,  nrow = 3, ncol = 3)
  ggsave(paste(DIR, PlotName, "_AllLayer" , ".png",sep=""), plot = PP, device='png', scale = 1, width = 8, height = 5, units ="in", dpi = 500, limitsize = TRUE,  bg = "transparent") 

}


Plot_Crosses <- function(Plot_Syns, LayerToPlot, DIR, PlotName, PLOT, ID){
  
  # Get color vector
  col_vector=GetColorMap("FBcol")
  
  # Get mean positions for this layer
  Temp_Dist=subset(Plot_Syns, Layer == paste("L",as.character(LayerToPlot),sep="")) 
  colnames(Temp_Dist)[which(colnames(Temp_Dist) %in% c("X_Mean","Z_Mean"))] = c("X","Z")
  Temp_ID=ID[Plot_Syns$Layer == paste("L",as.character(LayerToPlot),sep="")]
  
  # Get layer outline 
  eval(parse( text= paste("Outline=Outline_L", as.character(LayerToPlot),"_XZ",sep="") ))
  
  # Get cross lines
  Horizontal=data.frame(X=c(Temp_Dist[,"Cross_XL_X"], Temp_Dist[,"Cross_XR_X"]), Y=c(Temp_Dist[,"Cross_XL_Z"], Temp_Dist[,"Cross_XR_Z"]),
                        ID=c(Temp_ID, Temp_ID), COLOR=c(as.character(Temp_Dist[,"FBcol"]), as.character(Temp_Dist[,"FBcol"])))
  Vertical=data.frame(X=c(Temp_Dist[,"Cross_YU_X"], Temp_Dist[,"Cross_YD_X"]), Y=c(Temp_Dist[,"Cross_YU_Z"], Temp_Dist[,"Cross_YD_Z"]),
                        ID=c(Temp_ID, Temp_ID), COLOR=c(as.character(Temp_Dist[,"FBcol"]), as.character(Temp_Dist[,"FBcol"])))
  

  # Plot data
  p1 <- ggplot() + geom_path(data=Outline, aes(x=c1, y=c2), size = 1) + xlim(c(-9000, 9000)) +  coord_fixed(ratio = 1) +
    geom_line(data=Horizontal, aes(x=X, y=Y, group=ID, color=COLOR)) + 
    geom_line(data=Vertical, aes(x=X, y=Y, group=ID, color=COLOR)) +
    geom_point(data=Temp_Dist, aes(x=X,y=Z, color= FBcol),size=2 ) +
    theme_void() + guides(fill=FALSE) + theme(legend.position = "none") +
    scale_color_manual(values=col_vector, drop=FALSE)  +  theme(
      panel.background = element_rect(fill = "transparent", color = NA), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      panel.grid.major = element_blank(), # get rid of major grid
      panel.grid.minor = element_blank(), # get rid of minor grid
      legend.background = element_rect(fill = "transparent", color = NA), # get rid of legend bg
      legend.box.background = element_rect(fill = "transparent", color = NA), # get rid of legend panel bg
      panel.border = element_rect( fill = NULL,
                                   colour = NULL,
                                   size = NULL,
                                   linetype = NULL,
                                   color = NULL,
                                   inherit.blank = TRUE)
    ) 
                     
  

  
  if (PLOT==TRUE){
  ggsave(paste(DIR, PlotName, "_L", as.character(LayerToPlot) , ".png",sep=""), plot = p1, device='png', scale = 1, width = 8, height = 5, units ="in", dpi = 500, limitsize = TRUE,  bg = "transparent")
  }
  return(p1)
}



###############################################################################################################################
### Function for assigning FB columns to D0 neurons, here based on input column ###############################################

D0_OutputPlot_Connectivity <- function(D0_FB_Outputs){
  
  OutputsPlot_Sub <- ggplot(D0_FB_Outputs) + theme_classic() + theme(axis.text.x = element_text(angle = 90)) +
    scale_fill_gradient2(low="thistle", mid="blueviolet", high="black",
                         midpoint =0.5*max(D0_FB_Outputs$weightRelative),
                         limits=c(0,max(D0_FB_Outputs$weightRelative))) + 
    geom_tile(aes(FBcol.to, InputName,  fill=weightRelative)) + 
    facet_grid(reorder(type.from, desc(type.from)) ~ type.to, space="free", scales="free",switch="both") +
    theme(axis.text.y = element_blank()) + xlab("FB Columnar Type (POST)") + ylab("Delta0 Type (PRE)")
  
  return(OutputsPlot_Sub)
}


D0_InputPlot_Connectivity <- function(D0_FB_Inputs){
  
  InputsPlot_Sub <- ggplot(D0_FB_Inputs) + theme_classic() + theme(axis.text.x = element_text(angle = 90)) +
    scale_fill_gradient2(low="thistle", mid="blueviolet", high="black",
                         midpoint =0.5*max(D0_FB_Inputs$weightRelative),
                         limits=c(0,max(D0_FB_Inputs$weightRelative))) + 
    geom_tile(aes(OutputName, FBcol.from, fill=weightRelative)) + 
    facet_grid(reorder(type.from, desc(type.from)) ~ type.to, space="free", scales="free",switch="both") +
    theme(axis.text.x = element_blank()) + ylab("FB Columnar Type (PRE)") + xlab("Delta0 Type (POST)")
  
  return(InputsPlot_Sub)
}


Assign_FB_Columns_ToD0 <- function(D0_Distribution_Filt, D12s){
  
  
  # Get input and output tables for D0 neurons in FB
  Columnar_Types=unique(NamedBodies$bodytype[(startsWith(NamedBodies$bodytype, "PF")   | 
                                              startsWith(NamedBodies$bodytype, "F")  ) &
                                             !startsWith(NamedBodies$bodytype, "FB") ])
  D0_Types=as.character(unique(D0_Distribution_Filt$type))
  D0_Bag=buildInputsOutputsByType(D0_Types, slctROI="FB")
  D0_FB_Outputs_All=D0_Bag[["outputs_raw"]]
  D0_FB_Inputs_All=D0_Bag[["inputs_raw"]]
  
  
  # Subset data on columnar types and remove Delta12s 
  D0_FB_Outputs_Col=subset(D0_FB_Outputs_All, databaseType.to %in% Columnar_Types & !from %in% D12s$bodyid )
  D0_FB_Inputs_Col=subset(D0_FB_Inputs_All, databaseType.from %in% Columnar_Types & !to   %in% D12s$bodyid )
  
  
  # Get FB column names
  D0_FB_Outputs_Col$FBcol.to=NA
  D0_FB_Outputs_Col$FBcol.to=str_extract(D0_FB_Outputs_Col$name.to, "_C(\\d+)")
  D0_FB_Outputs_Col$FBcol.to=sapply(D0_FB_Outputs_Col$FBcol.to, substring, 2, 3)
  D0_FB_Outputs_Col$FBcol.to=factor(D0_FB_Outputs_Col$FBcol.to, levels=c("C1","C2","C3","C4","C5","C6","C7","C8","C9"))
  D0_FB_Outputs_Col$FBcol.from=NA
  D0_FB_Outputs_Col$FBcol.from=str_extract(D0_FB_Outputs_Col$name.from, "_C(\\d+)")
  D0_FB_Outputs_Col$FBcol.from=sapply(D0_FB_Outputs_Col$FBcol.from, substring, 2, 3)
  D0_FB_Outputs_Col$FBcol.from=factor(D0_FB_Outputs_Col$FBcol.from, levels=c("C1","C2","C3","C4","C5","C6","C7","C8","C9"))
  
  D0_FB_Inputs_Col$FBcol.from=NA
  D0_FB_Inputs_Col$FBcol.from=str_extract(D0_FB_Inputs_Col$name.from, "_C(\\d+)")
  D0_FB_Inputs_Col$FBcol.from=sapply(D0_FB_Inputs_Col$FBcol.from, substring, 2, 3)   
  D0_FB_Inputs_Col$FBcol.from=factor(D0_FB_Inputs_Col$FBcol.from, levels=c("C1","C2","C3","C4","C5","C6","C7","C8","C9"))
  D0_FB_Inputs_Col$FBcol.to=NA
  D0_FB_Inputs_Col$FBcol.to=str_extract(D0_FB_Inputs_Col$name.to, "_C(\\d+)")
  D0_FB_Inputs_Col$FBcol.to=sapply(D0_FB_Inputs_Col$FBcol.to, substring, 2, 3)   
  D0_FB_Inputs_Col$FBcol.to=factor(D0_FB_Inputs_Col$FBcol.to, levels=c("C1","C2","C3","C4","C5","C6","C7","C8","C9"))   
  
  
  # Make plots of all inputs/outputs
  Output_Plot_1 = plotConnectivityMatrix(D0_FB_Outputs_Col, byGroup = "id", connectionMeasure="weightRelative") 
  Output_Plot_1 = structureMatrixPlotByType(Output_Plot_1)
  ggsave(paste(D0_Dir_Connect, "Delta0_Outputs_All.png",sep=""),
         plot = Output_Plot_1, device='png', scale = 1, width = 10, height = 10, units ="in", dpi = 500, limitsize = TRUE)
  
  Input_Plot_1 = plotConnectivityMatrix(D0_FB_Inputs_Col, byGroup = "id", connectionMeasure="weightRelative") 
  Input_Plot_1 = structureMatrixPlotByType(Input_Plot_1)
  ggsave(paste(D0_Dir_Connect, "Delta0_Inputs_All.png",sep=""),
         plot = Input_Plot_1, device='png', scale = 1, width = 10, height = 10, units ="in", dpi = 500, limitsize = TRUE)
  
  
  # Average by type and column
  D0_FB_Outputs_Mean=D0_FB_Outputs_Col %>% group_by(type.to, FBcol.to, from, type.from, name.from, FBcol.from) %>% summarise(weightRelative=mean(weightRelative), weight=mean(weight))
  D0_FB_Outputs_Mean$InputName =paste(as.character(D0_FB_Outputs_Mean$type.from), "_", as.character(D0_FB_Outputs_Mean$FBcol.from), "_", as.character(D0_FB_Outputs_Mean$from),sep="")
  D0_FB_Outputs_Mean$OutputName=paste(as.character(D0_FB_Outputs_Mean$type.to)  , "_", as.character(D0_FB_Outputs_Mean$FBcol.to),sep="")
  D0_FB_Outputs_Mean$from=factor(D0_FB_Outputs_Mean$from, levels= unique(D0_FB_Outputs_Col$from) )
  D0_FB_Outputs_Mean$type.from =  substring(D0_FB_Outputs_Mean$type.from, 6, 100)
  D0_FB_Outputs_Mean$InputName = factor(D0_FB_Outputs_Mean$InputName, levels = sort(unique(D0_FB_Outputs_Mean$InputName)))
  
  
  D0_FB_Inputs_Mean=D0_FB_Inputs_Col %>% group_by(type.from, FBcol.from, to, type.to, name.to, FBcol.to) %>% summarise(weightRelative=mean(weightRelative), weight=mean(weight))
  D0_FB_Inputs_Mean$InputName =paste(as.character(D0_FB_Inputs_Mean$type.from), "_", as.character(D0_FB_Inputs_Mean$FBcol.from),sep="")
  D0_FB_Inputs_Mean$OutputName=paste(as.character(D0_FB_Inputs_Mean$type.to)  , "_", as.character(D0_FB_Inputs_Mean$FBcol.to), "_", as.character(D0_FB_Inputs_Mean$to),sep="")
  D0_FB_Inputs_Mean$to=factor(D0_FB_Inputs_Mean$to, levels= unique(D0_FB_Inputs_Col$to) )
  D0_FB_Inputs_Mean$type.to =  substring(D0_FB_Inputs_Mean$type.to, 6, 100)
  D0_FB_Inputs_Mean$OutputName = factor(D0_FB_Inputs_Mean$OutputName, levels = sort(unique(D0_FB_Inputs_Mean$OutputName)))
  
  
  # Plot Delta0 connectivity matrices to see if typing is decent
  OutputsPlot_Sub=D0_OutputPlot_Connectivity(D0_FB_Outputs_Mean)
  ggsave(paste(D0_Dir_Connect, "Delta0_Outputs_ColAve.png",sep=""),
         plot = OutputsPlot_Sub, device='png', scale = 1, width = 10, height = 10, units ="in", dpi = 500, limitsize = TRUE)
  
  InputsPlot_Sub=D0_InputPlot_Connectivity(D0_FB_Inputs_Mean)
  ggsave(paste(D0_Dir_Connect, "Delta0_Inputs_ColAve.png",sep=""),
         plot = InputsPlot_Sub, device='png', scale = 1, width = 10, height = 10, units ="in", dpi = 500, limitsize = TRUE)
  
  
 
  # remove cell types that don't innervate all columns
  OutputColumnsNum = D0_FB_Outputs_Mean %>% group_by(type.to, type.from) %>% summarize(NumOfCols=length(unique(FBcol.to)))
  OutputColumnsNum = subset(OutputColumnsNum, NumOfCols>=7)
  D0_FB_Outputs_MeanSub=subset(D0_FB_Outputs_Mean, type.to %in% OutputColumnsNum$type.to & type.from %in% OutputColumnsNum$type.from)
  
  InputColumnsNum  = D0_FB_Inputs_Mean %>% group_by(type.to, type.from) %>% summarize(NumOfCols=length(unique(FBcol.from)))
  InputColumnsNum  = subset(InputColumnsNum, NumOfCols>=7)
  D0_FB_Inputs_MeanSub=subset(D0_FB_Inputs_Mean, type.from %in% InputColumnsNum$type.from & type.to %in% InputColumnsNum$type.to)
  
  
  # Replot the data
  OutputsPlot_Sub2=D0_OutputPlot_Connectivity(D0_FB_Outputs_MeanSub)
  ggsave(paste(D0_Dir_Connect, "Delta0_Outputs_ColAve_Sub.png",sep=""),
         plot = OutputsPlot_Sub2, device='png', scale = 1, width = 10, height = 10, units ="in", dpi = 500, limitsize = TRUE)
  
  InputsPlot_Sub2=D0_InputPlot_Connectivity(D0_FB_Inputs_MeanSub)
  ggsave(paste(D0_Dir_Connect, "Delta0_Inputs_ColAve_Sub.png",sep=""),
         plot = InputsPlot_Sub2, device='png', scale = 1, width = 10, height = 10, units ="in", dpi = 500, limitsize = TRUE)
  
  
  # Now average filtered data by column, independent of partner type
  D0_FB_Outputs_MeanSubCol = D0_FB_Outputs_MeanSub  %>% group_by(FBcol.to, from, type.from, name.from, FBcol.from, InputName) %>%  summarise(weightRelative=mean(weightRelative))
  D0_FB_Inputs_MeanSubCol  = D0_FB_Inputs_MeanSub   %>% group_by(FBcol.from, to, type.to, name.to, FBcol.to, OutputName)      %>%  summarise(weightRelative=mean(weightRelative))
  
  
  
  # Plot the per-column connectivity matrices
  OutputsPlot_Sub3 <- ggplot(D0_FB_Outputs_MeanSubCol) + theme_classic() + theme(axis.text.x = element_text(angle = 90)) +
    scale_fill_gradient2(low="thistle", mid="blueviolet", high="black",
                         midpoint =0.5*max(D0_FB_Outputs_MeanSubCol$weightRelative),
                         limits=c(0,max(D0_FB_Outputs_MeanSubCol$weightRelative))) + 
    geom_tile(aes(FBcol.to, InputName,  fill=weightRelative)) + 
    facet_grid(rows = vars(type.from), space="free", scales="free",switch="both") +
    theme(axis.text.y = element_blank()) + xlab("FB Columnar Type (POST)") + ylab("Delta0 Type (PRE)")
  ggsave(paste(D0_Dir_Connect, "Delta0_Outputs_ColFinal.png",sep=""),
         plot = OutputsPlot_Sub3, device='png', scale = 1, width = 4, height = 10, units ="in", dpi = 500, limitsize = TRUE)
  
  
  InputsPlot_Sub3 <- ggplot(D0_FB_Inputs_MeanSubCol) + theme_classic() + theme(axis.text.x = element_text(angle = 90)) +
    scale_fill_gradient2(low="thistle", mid="blueviolet", high="black",
                         midpoint =0.5*max(D0_FB_Inputs_MeanSubCol$weightRelative),
                         limits=c(0,max(D0_FB_Inputs_MeanSubCol$weightRelative))) + 
    geom_tile(aes(OutputName, FBcol.from, fill=weightRelative)) + 
    facet_grid(cols = vars(type.to), space="free", scales="free",switch="both") +
    theme(axis.text.x = element_blank()) + ylab("FB Columnar Type (PRE)") + xlab("Delta0 Type (POST)")
  ggsave(paste(D0_Dir_Connect, "Delta0_Inputs_ColFinal.png",sep=""),
         plot = InputsPlot_Sub3, device='png', scale = 1, width = 10, height = 3, units ="in", dpi = 500, limitsize = TRUE)
  
  
  # Loop over Delta0s that talk to columnar neurons and assign them columns based on their connectivity
  print(paste(as.character(length(unique(c(D0_FB_Outputs_Mean$from, D0_FB_Inputs_Mean$to)))), " of ", as.character(length(unique(D0_Distribution_Filt$bodyid))), 
              " Delta0 neurons (that are not D12) can be assigned columns based on connectivity"))
  
  IDID=unique(D0_FB_Outputs_Mean$from)
  Delta0_OutputColumnAssignments=data.frame(from=character(), FBcol.from_new=character(), type.from=character())
  Case=3
  for (nnn in 1:length(IDID)){
    Temp_output=subset(D0_FB_Outputs_Mean, from == IDID[nnn])
    Temp_output=Temp_output %>% group_by(FBcol.to, from, type.from) %>% summarize(weight=mean(weight), n=n())
    Temp_output$NumOfSyns=Temp_output$weight * Temp_output$n
    if (Case==1){
      Column= round(mean(as.numeric(substring(Temp_output$FBcol.to,2,3))))
      Column = paste("C",as.character(Column),sep="")
    } else if (Case==2){
      Column = round(sum(as.numeric(substring(Temp_output$FBcol.to,2,3))*Temp_output$NumOfSyns)/sum(Temp_output$NumOfSyns))
      Column = paste("C",as.character(Column),sep="")
    } else if (Case==3){
      Temp_Ind=which(Temp_output$NumOfSyns == max(Temp_output$NumOfSyns))
      if (length(Temp_Ind)==1){
        Column=Temp_output$FBcol.to[Temp_Ind]
      } else {
        Column = round(sum(as.numeric(substring(Temp_output$FBcol.to,2,3))*Temp_output$NumOfSyns)/sum(Temp_output$NumOfSyns))
        Column = paste("C",as.character(Column),sep="")
      }
    }
    Temp_output_DF=data.frame(from=as.character(Temp_output$from[1]), FBcol.from_new=Column, type.from=as.character(Temp_output$type.from[1]))
    Delta0_OutputColumnAssignments<-rbind(Delta0_OutputColumnAssignments,Temp_output_DF)
  }
  
  Output_Assign_Plot=merge(D0_FB_Outputs_Mean, Delta0_OutputColumnAssignments, by="from")
  Output_Assign_Plot=distinct(Output_Assign_Plot[c("from","FBcol.from","FBcol.from_new")])
  Output_Assign_Plot=Output_Assign_Plot %>% group_by(FBcol.from, FBcol.from_new) %>% summarize(n=n())
  Output_Assign_Plot$FBcol.from=factor(Output_Assign_Plot$FBcol.from, sort(unique(Output_Assign_Plot$FBcol.from)))
  Output_Assign_Plot$FBcol.from_new=factor(Output_Assign_Plot$FBcol.from_new, sort(unique(Output_Assign_Plot$FBcol.from)))
  P1=ggplot() + geom_point(data=Output_Assign_Plot, aes(x=FBcol.from, y=FBcol.from_new, size=n)) + 
    xlab("Column from synapse locations") + ylab("Column from connectivity") + ggtitle("Using Delta0 ouputs")
  
  ggsave(paste(D0_Dir_Connect, "ColumnsAssigned_by_Output.png",sep=""),
         plot = P1, device='png', scale = 1, width = 7, height = 5, units ="in", dpi = 500, limitsize = TRUE)
  

  Output_Types=unique(Delta0_OutputColumnAssignments$type.from)
  for (nnn in 1:length(Output_Types)){
    Output_Assign_Plot=merge(D0_FB_Outputs_Mean, Delta0_OutputColumnAssignments, by="from")
    
    Output_Assign_Plot=subset(Output_Assign_Plot, type.from.y == Output_Types[nnn])
    Output_Assign_Plot=distinct(Output_Assign_Plot[c("from","FBcol.from","FBcol.from_new")])
    
    Output_Assign_Plot=Output_Assign_Plot %>% group_by(FBcol.from, FBcol.from_new) %>% summarize(n=n())
    Output_Assign_Plot$FBcol.from=factor(Output_Assign_Plot$FBcol.from, sort(unique(Output_Assign_Plot$FBcol.from)))
    Output_Assign_Plot$FBcol.from_new=factor(Output_Assign_Plot$FBcol.from_new, sort(unique(Output_Assign_Plot$FBcol.from_new)))
    
    P1=ggplot() + geom_point(data=Output_Assign_Plot, aes(x=FBcol.from, y=FBcol.from_new, size=n)) + 
      xlab("Column from synapse locations") + ylab("Column from connectivity") + ggtitle(paste("Using ", as.character( Output_Types[nnn]) , " outputs",sep=""))
    ggsave(paste(D0_Dir_Connect, "ColumnsAssigned_by_Output_",  as.character( Output_Types[nnn]) ,".png",sep=""),
           plot = P1, device='png', scale = 1, width = 7, height = 5, units ="in", dpi = 500, limitsize = TRUE)
  }
  

  IDID=unique(D0_FB_Inputs_Mean$to)
  Delta0_InputColumnAssignments=data.frame(to=character(), FBcol.to_new=character(), type.to=character())
  for (nnn in 1:length(IDID)){
    Temp_Input=subset(D0_FB_Inputs_Mean, to == IDID[nnn])
    Temp_Input=Temp_Input %>% group_by(FBcol.from, to, type.to) %>% summarize(weight=mean(weight), n=n())
    Temp_Input$NumOfSyns=Temp_Input$weight * Temp_Input$n
    if (Case==1){
      Column= round(mean(as.numeric(substring(Temp_Input$FBcol.from,2,3))))
      Column = paste("C",as.character(Column),sep="")
    } else if (Case==2){
      Column = round(sum(as.numeric(substring(Temp_Input$FBcol.from,2,3))*Temp_Input$NumOfSyns)/sum(Temp_Input$NumOfSyns))
      Column = paste("C",as.character(Column),sep="")
    } else if (Case==3){
      Temp_Ind=which(Temp_Input$NumOfSyns == max(Temp_Input$NumOfSyns))
      if (length(Temp_Ind)==1){
        Column=Temp_Input$FBcol.from[Temp_Ind]
      } else {
        Column = round(sum(as.numeric(substring(Temp_Input$FBcol.from,2,3))*Temp_Input$NumOfSyns)/sum(Temp_Input$NumOfSyns))
        Column = paste("C",as.character(Column),sep="")
      }
    }
    Temp_Input_DF=data.frame(to=as.character(Temp_Input$to[1]), FBcol.to_new=Column, type.to=as.character(Temp_Input$type.to[1]))
    Delta0_InputColumnAssignments<-rbind(Delta0_InputColumnAssignments,Temp_Input_DF)
  }
  
  
  Input_Assign_Plot=merge(D0_FB_Inputs_Mean, Delta0_InputColumnAssignments, by="to")
  Input_Assign_Plot=distinct(Input_Assign_Plot[c("to","FBcol.to","FBcol.to_new")])
  Input_Assign_Plot=Input_Assign_Plot %>% group_by(FBcol.to, FBcol.to_new) %>% summarize(n=n())
  Input_Assign_Plot$FBcol.to=factor(Input_Assign_Plot$FBcol.to, sort(unique(Input_Assign_Plot$FBcol.to)))
  Input_Assign_Plot$FBcol.to_new=factor(Input_Assign_Plot$FBcol.to_new, sort(unique(Input_Assign_Plot$FBcol.to)))
  P2 = ggplot() + geom_point(data=Input_Assign_Plot, aes(x=FBcol.to, y=FBcol.to_new, size=n)) + 
    xlab("Column from synapse locations") + ylab("Column from connectivity") + ggtitle("Using Delta0 inputs")
  
  ggsave(paste(D0_Dir_Connect, "ColumnsAssigned_by_Input.png",sep=""),
         plot = P2, device='png', scale = 1, width = 7, height = 5, units ="in", dpi = 500, limitsize = TRUE)


  Input_Types=unique(Delta0_InputColumnAssignments$type.to)
  for (nnn in 1:length(Input_Types)){
    Input_Assign_Plot=merge(D0_FB_Inputs_Mean, Delta0_InputColumnAssignments, by="to")
    
    Input_Assign_Plot=subset(Input_Assign_Plot, type.to.y == Input_Types[nnn])
    Input_Assign_Plot=distinct(Input_Assign_Plot[c("to","FBcol.to","FBcol.to_new")])
    
    Input_Assign_Plot=Input_Assign_Plot %>% group_by(FBcol.to, FBcol.to_new) %>% summarize(n=n())
    Input_Assign_Plot$FBcol.to=factor(Input_Assign_Plot$FBcol.to, sort(unique(Input_Assign_Plot$FBcol.to)))
    Input_Assign_Plot$FBcol.to_new=factor(Input_Assign_Plot$FBcol.to_new, sort(unique(Input_Assign_Plot$FBcol.to_new)))
    
    P1=ggplot() + geom_point(data=Input_Assign_Plot, aes(x=FBcol.to, y=FBcol.to_new, size=n)) + 
      xlab("Column from synapse locations") + ylab("Column from connectivity") + ggtitle(paste("Using ", as.character( Input_Types[nnn]) , " inputs",sep=""))
    ggsave(paste(D0_Dir_Connect, "ColumnsAssigned_by_Input_",  as.character( Input_Types[nnn]) ,".png",sep=""),
           plot = P1, device='png', scale = 1, width = 7, height = 5, units ="in", dpi = 500, limitsize = TRUE)
  }
  

  
}

###############################################################################################################################
### Function for assigning FB columns to PB-FB-XX and D0 neurons, here based on closest type average  #########################


Assign_FB_Columns <- function(N_Distribution, FX_Column_Positions_Filt, DIR){
  
  N_Distribution$Distance2FxCol=NA
  N_Distribution$FBcol=as.character(N_Distribution$FBcol)
  Column_Names=colnames(N_Distribution)
  
  # Get unique body IDs
  Uni_BodyID=unique(N_Distribution$bodyid)
  
  # Dataframe mapping bodyid to FBcol
  Bodyid_To_FBcol=data.frame(bodyid=numeric(), FBcol=as.character(), Distance2FxCol=numeric(), FBcol_Con=numeric())
  
  
  # Loop over each neuron
  for (nnn in 1:length(Uni_BodyID)){
    
    # Get all synapses from all layers for this neuron
    Neuron_Data=subset(N_Distribution, bodyid==Uni_BodyID[nnn])
    
    # Get the synapses for this neuron from the layer with the most synapses for this neuron type
    MaxLayer=Neuron_Data$Layer[which(Neuron_Data$NumberOfSyns == max(Neuron_Data$NumberOfSyns))]
    
    
    # If type is Delta0 and Layer 1 is present use that
    if ( startsWith(as.character(Neuron_Data$type[1]), "Delta0") & "L1" %in% Neuron_Data$Layer ){
      MaxLayer="L1"
    }
    
    # Get layer data
    Layer_Data=subset(Neuron_Data, Layer==MaxLayer)
    
    # Subset FX data based on this layer
    TempCol_Positions=subset(FX_Column_Positions_Filt, Layer==Layer_Data$Layer)
    
    # Compute distance to every FXs neurons synapse mean in this layer
    Temp_Distances_All= data.frame(Col= TempCol_Positions$FBcol, Distance= ((TempCol_Positions$X_Mean - Layer_Data$X_Mean)^2 +
                                                                              (TempCol_Positions$Z_Mean - Layer_Data$Z_Mean)^2 +
                                                                              (TempCol_Positions$Y_Mean - Layer_Data$Y_Mean)^2)^0.5)  # could do some normalization here
    
    # Convert column number to numeric
    Temp_Distances_All$Col=as.numeric(sapply(Temp_Distances_All$Col, substring, 2, 2))
    
    # Compute average distance by column
    Temp_Distances_Mean=Temp_Distances_All %>% group_by(Col) %>% summarise(Distance=mean(Distance)) #quantile(Distance,0.1)
    
    # Interpolate data
    #Temp_Distances_Mean_UP=as.data.frame(approx(Temp_Distances_Mean$Col, y = Temp_Distances_Mean$Distance, xout=seq(1,9,0.1), method="linear"))
    Temp_Distances_Mean_UP=as.data.frame(spline(x = Temp_Distances_Mean$Col, y = Temp_Distances_Mean$Distance,  n = 91, method = "fmm", xmin = 0.6, xmax = 9.4))
    colnames(Temp_Distances_Mean_UP)=c("Col","Distance")
    
    # Take min of function
    Min_Distance=min(Temp_Distances_Mean_UP$Distance)
    Min_IND= which(Temp_Distances_Mean_UP$Distance == Min_Distance)
    Min_C=Temp_Distances_Mean_UP$Col[Min_IND]
    Min_Col=round(Min_C)
    
    # Take min distance as closest column
    Temp_df=data.frame(bodyid=Uni_BodyID[nnn], FBcol = paste("C",as.character(Min_Col),sep=""), Distance2FxCol= Min_Distance, FBcol_Con=Min_C)
    Bodyid_To_FBcol=rbind(Bodyid_To_FBcol, Temp_df)
    
    
    if (Neuron_Data$type[1] %in% c("PFGs","PFL1","PFL2","PFL3") | startsWith(as.character(Neuron_Data$type[1]),"Delta0")){
      P1<-ggplot() + geom_point(data=Temp_Distances_All, aes(x=Col, y=Distance), position = position_jitter(width = 0.1), size=0.2) +
        geom_point(data=Temp_Distances_Mean, aes(x=Col, y=Distance),color="red") +
        ggtitle(paste(Neuron_Data$type[1], " ", Neuron_Data$PBglom[1], " ", Neuron_Data$bodyid[1] )) +
        geom_line(data=Temp_Distances_Mean_UP, aes(x=Col,y=Distance)) + scale_x_continuous(breaks=seq(1,9,1)) + 
        geom_point(data=data.frame(XXX = Min_C, YYY = Min_Distance), aes(x=XXX,y=YYY), color="blue", size=2) + 
        geom_point(data=data.frame(XXX = Min_Col, YYY = Min_Distance), aes(x=XXX,y=YYY), color="green", size=2)
      
      
      OUTLINE = eval(parse( text= paste("Outline=Outline_", as.character(MaxLayer),"_XZ",sep="") ))
      col_vector=GetColorMap("FBcol")
      
      P2<-ggplot() + geom_point(data=Layer_Data, aes(x=X_Mean, y=Z_Mean),  size=5, alpha = 1, stroke = 0, shape=16) +
        geom_point(data=TempCol_Positions, aes(x=X_Mean, y=Z_Mean, color=FBcol),  size=2, alpha = 1, stroke = 0, shape=16) + 
        geom_path(data=OUTLINE, aes(x=c1, y=c2), size = 1) + xlim(c(-9000, 9000)) +  coord_fixed(ratio = 1) + theme_void() + guides(fill=FALSE) + theme(legend.position = "none") +
        scale_color_manual(values=col_vector, drop=FALSE)  +  theme(
          panel.background = element_rect(fill = "transparent"), # bg of the panel
          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          panel.grid.major = element_blank(), # get rid of major grid
          panel.grid.minor = element_blank(), # get rid of minor grid
          legend.background = element_rect(fill = "transparent"), # get rid of legend bg
          legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
        )  +
        ggtitle( paste("Layer ", as.character(MaxLayer),  "   PBglom ", as.character(Neuron_Data$PBglom[1]) ,
                       "   FBcol ", as.character(Temp_df$FBcol[1]), "   ", as.character(Neuron_Data$type[1]), "   ",  as.character(Neuron_Data$bodyid[1]),sep="") )
      
      PP=ggarrange(P1, P2, ncol=2, nrow=1)
      
      ggsave(paste(DIR, Neuron_Data$type[1], "_", Neuron_Data$PBglom[1], "_" ,Neuron_Data$bodyid[1] ,".png",sep=""),
             plot = PP, device='png', scale = 1, width = 10, height = 5, units ="in", dpi = 300, limitsize = TRUE)
      
    }
    
  }
  
  N_Distribution=N_Distribution[ , !colnames(N_Distribution) %in% c("FBcol","Distance2FxCol")]
  N_Distribution_New=merge(N_Distribution, Bodyid_To_FBcol, by="bodyid")
  N_Distribution_New=N_Distribution_New[c(Column_Names,"FBcol_Con")]
  
  N_Distribution_New$FBcol=factor(N_Distribution_New$FBcol, levels = c("C1","C2","C3","C4","C5","C6","C7","C8","C9") )
  
  return(N_Distribution_New)
}

###############################################################################################################################
### Function for plotting mapping between PB glom and FB cols #################################################################



Plot_PBglom_FBcol_Mapping <- function(PFX_Distribution, DIR, NAME){
  
  # Get number of neurons going from each glom to each column
  PFX_GlomToColumn=PFX_Distribution %>% group_by(bodyid, type, PBglom, FBcol) %>% summarise(Num_layers = n(), FBcol_Con= mean(FBcol_Con))
  PFX_GlomToColumn=PFX_GlomToColumn %>% group_by(type, PBglom, FBcol) %>% summarise(Num_Neurons = n(), FBcol_Con= mean(FBcol_Con))
  
  
  PB_FB_Types=unique(PFX_GlomToColumn$type)
  for (nnn in 1:length(PB_FB_Types)){
    
    
    # Data for mapping
    PB_FB_Mapping=subset(PFX_GlomToColumn, type==PB_FB_Types[nnn] )
    
    
    # Set up node names and positions
    nodes = data.frame(name = c("L9","L8","L7","L6","L5","L4","L3","L2","L1",
                                "R1","R2","R3","R4","R5","R6","R7","R8","R9",
                                "C1","C2","C3","C4","C5","C6","C7","C8","C9") )
    nodes$x <- c( seq(from= 9 ,to = 1), seq(from= -1 ,to = -9),  seq(from= -4 ,to = 4))
    nodes$y <- c( rep(1, 18), rep(0, 9))
    
    
    # Assign node from/to pairs (numerical)
    PB_FB_Mapping$from = match(PB_FB_Mapping$PBglom, nodes$name)  
    PB_FB_Mapping$to   = match(PB_FB_Mapping$FBcol,  nodes$name)
    
    
    # build graph object
    graph = tbl_graph(nodes, PB_FB_Mapping)
    
    
    # Get column vector
    col_vector=GetColorMap("FBcol_Function")
    PB_FB_Mapping$Color = order(PB_FB_Mapping$FBcol)
    
    
    P1<-ggraph(graph,layout="manual",x=nodes$x,y=nodes$y) +
      geom_edge_diagonal(aes(width=Num_Neurons,color=FBcol),alpha=0.5,strength=0.5) +
      geom_node_point(size=5)  + scale_edge_color_manual(values=col_vector)  +
      geom_node_text(aes(label=name),angle=40,size=4, nudge_y = c(rep(0.06,18),rep(-0.06,9)) ) +
      theme_classic() + 
      theme(legend.text=element_text(size=6),legend.title=element_text(size=6),
            axis.line=element_blank(),axis.text.x=element_blank(),
            axis.text.y=element_blank(),axis.ticks=element_blank(),
            axis.title.x=element_blank(),axis.title.y=element_blank()) + 
      ggtitle(PB_FB_Mapping$type[1])
    
    ggsave(paste(DIR, "Graph_", PB_FB_Mapping$type[1], "_", NAME, ".png",sep=""), plot = P1, device='png', scale = 1, width = 8, height = 4.5, units ="in", dpi = 500, limitsize = TRUE) 
    ggsave(paste(DIR, "Graph_", PB_FB_Mapping$type[1], "_", NAME, ".pdf",sep=""), plot = P1, device='pdf', scale = 1, width = 8, height = 4.5, units ="in", dpi = 500, limitsize = TRUE) 
    
    
    PB_FB_Mapping$COL=as.numeric(sapply(PB_FB_Mapping$FBcol, substring, 2, 2))

    
    R_PB=subset(PB_FB_Mapping, startsWith(as.character(PBglom), "R"))
    L_PB=subset(PB_FB_Mapping, startsWith(as.character(PBglom), "L"))
    L_PB$PBglom=factor(L_PB$PBglom, levels = sort(unique(L_PB$PBglom), decreasing = TRUE))
    
    P2 <- ggplot() + geom_point(data=R_PB, aes(x=FBcol_Con, y=PBglom, size=1)) + 
      geom_point(data=R_PB, aes(x=COL, y=PBglom, size=1), color="blue") +
      scale_x_continuous(breaks=seq(from=1,to=9), limits=c(0,10)) 
    
    P3 <- ggplot() + geom_point(data=L_PB, aes(x=FBcol_Con, y=PBglom, size=1)) + 
      geom_point(data=L_PB, aes(x=COL, y=PBglom, size=1), color="blue") +
      scale_x_reverse(breaks=seq(from=1,to=9), limits=c(10,0)) 
      
    
    P12=ggarrange(P2,P3,ncol=1,nrow=2)
    P123=ggarrange(P1, P12, ncol=2, nrow=1)
    
    
    ggsave(paste(DIR, "Mapping_", PB_FB_Mapping$type[1], "_", NAME, ".png",sep=""), plot = P123, device='png', scale = 1, width = 12, height = 5, units ="in", dpi = 500, limitsize = TRUE) 
    ggsave(paste(DIR, "Mapping_", PB_FB_Mapping$type[1], "_", NAME, ".pdf",sep=""), plot = P123, device='pdf', scale = 1, width = 12, height = 5, units ="in", dpi = 500, limitsize = TRUE) 
    
    
  }
  return(PFX_GlomToColumn)
}


###############################################################################################################################
### Function for getting the width of each neuron's synapse cloud in a direction tangent to the FB layer ######################


Get_SynLayerDistribution <- function(All_Neurons, Thresh, PlotDir){
  
  
  Syn_Distribution <- data.frame(bodyid=as.numeric(), type=as.character(), Layer=as.character(), 
                                 Side=as.character(), PBglom=as.character(), FBcol=as.character(),
                                 X_Mean=as.numeric(), Y_Mean=as.numeric(), Z_Mean=as.numeric(), 
                                 X_STD_Orth=as.numeric(), Z_STD_Orth=as.numeric(),
                                 X_HW_Orth=as.numeric(), Z_HW_Orth=as.numeric(),
                                 Cross_XL_X=as.numeric(),  Cross_XL_Z=as.numeric(),
                                 Cross_XR_X=as.numeric(),  Cross_XR_Z=as.numeric(),
                                 Cross_YU_X=as.numeric(),  Cross_YU_Z=as.numeric(),
                                 Cross_YD_X=as.numeric(),  Cross_YD_Z=as.numeric(),
                                 Box_UL_X=as.numeric(),    Box_UL_Z=as.numeric(),
                                 Box_UR_X=as.numeric(),    Box_UR_Z=as.numeric(),
                                 Box_LL_X=as.numeric(),    Box_LL_Z=as.numeric(),
                                 Box_LR_X=as.numeric(),    Box_LR_Z=as.numeric(),
                                 LayerLength=as.numeric(), NumberOfSyns=as.numeric())
  
  
  Neuron_Types=unique(All_Neurons$type)
  Layers_All=sort(unique(All_Neurons$Layer))
  
  # Loop over neuron types
  for (ttt in 1:length(Neuron_Types)){
    SubData=subset(All_Neurons, type == Neuron_Types[ttt])
    
    
    # Loop of layers
    for (lll in 1:8){
      LayerData=subset(SubData, Layer==Layers_All[lll])
      eval(parse( text= paste("MID=Mid_L", as.character(lll),sep="") ))
      eval(parse( text= paste("OUTLINE=Outline_L", as.character(lll),"_XZ",sep="") ))
      
      
      # Loop over individual neurons
      NeuronIDs=unique(LayerData$bodyid)
      for (nnn in 1:length(NeuronIDs)){
        NeuronData=subset(LayerData, bodyid==NeuronIDs[nnn])
        
        # Rename Z as Y or it gets confusing when looking at plots
        Mean_YY=mean(NeuronData$Y)
        NeuronData=NeuronData[colnames(NeuronData) != "Y"]
        colnames(NeuronData)[colnames(NeuronData) == "Z"] <-"Y"
        
        # Get X range
        X_RangeMin = quantile(NeuronData$X, probs = 0.1) 
        X_RangeMax = quantile(NeuronData$X, probs = 0.9) 
        
        if (length(NeuronData$X)>Thresh ){
          
        ############# Project data into new coordinate frame whose main axis is locally parallel to FB layer ############# 
          
 
          
          # Get slope of mid range 
          Mid_subset=subset(MID, X>X_RangeMin & X<X_RangeMax)
          Layer_Length_Temp=sum((diff(MID$X)^2 + diff(MID$Y)^2)^0.5)
          
          # If no values in range, use two closest values
          if (length(Mid_subset$X)<2){
            Center=mean(c(X_RangeMin,X_RangeMax))
            INDS= order(abs(MID$X - Center))[c(1,2)]
            Mid_subset=MID[INDS,]
          }
          
          # Compute slope
          TempSlope=median( diff(Mid_subset$Y)/diff(Mid_subset$X) )
          
          # Define new axis (vectors should be length 1 and orthogonal to keep shape of point cloud)
          # Also, keep Ax_1 pointing towards right and Ax_2 pointing up, so that axes are consistent across clouds
          Ax_1=c(1,TempSlope)
          Ax_1=Ax_1/sum(Ax_1^2)^0.5

          Ax_2=c(1, -Ax_1[1]/Ax_1[2])
          Ax_2=Ax_2/sum(Ax_2^2)^0.5
          if (Ax_2[2]<0){
            Ax_2[1]=-Ax_2[1]
            Ax_2[2]=-Ax_2[2]} 

          
          # Make sure vectors are orthogonal
          if (Ax_1 %*% Ax_2 > 10^-15){print(paste("Vectors not orthogonal ttt=", as.character(ttt), "   lll = ", as.character(lll), "   nnn = ", as.character(nnn)))}
          
          # Project onto new x and y axes
          NeuronData$X_New= as.matrix(NeuronData[c("X","Y")]) %*% Ax_1  # This
          NeuronData$Y_New= as.matrix(NeuronData[c("X","Y")]) %*% Ax_2  # This
          

          # Define vectors to get from new space back to old space
          Ax1_NewToOld= t(Ax_1 * as.matrix(c(1,-1)))
          Ax2_NewToOld= t(Ax_2 * as.matrix(c(-1,1)))
          
          
        ############# Compute values we want to save for this neuron ############# 
          
          # Mean in original space 
          #TempMean=as.data.frame( t(colMedians( as.matrix(NeuronData[c("X","Y")])) )) # This
          
          # Mean and STD along new axes
          TempMean_New=as.data.frame( t(colMedians( as.matrix(NeuronData[c("X_New","Y_New")]))))
          TempSTD_New=as.data.frame(t(sapply(NeuronData[c("X_New","Y_New")], sd)))
          
          # Ranges along new axes
          TempRange_New=data.frame(X = abs(quantile(NeuronData$X_New, probs = 0.05) - quantile(NeuronData$X_New, probs = 0.95))/2,
                                   Y = abs(quantile(NeuronData$Y_New, probs = 0.05) - quantile(NeuronData$Y_New, probs = 0.95))/2)
          
          # Mean in original space 
          TempMean=data.frame(X = as.matrix(TempMean_New[c("X_New","Y_New")]) %*% t(Ax1_NewToOld), Y = as.matrix(TempMean_New[c("X_New","Y_New")]) %*% t(Ax2_NewToOld))
          
          
          # Compute Cross points (4 points making up a cross around synapse mean position)
          Cross_XL= data.frame(X = as.numeric(as.matrix((TempMean_New -  c(TempRange_New$X,0))) %*% t(Ax1_NewToOld)) , 
                               Y = as.numeric(as.matrix((TempMean_New -  c(TempRange_New$X,0))) %*% t(Ax2_NewToOld)) )
          Cross_XR= data.frame(X = as.numeric(as.matrix((TempMean_New +  c(TempRange_New$X,0))) %*% t(Ax1_NewToOld)) , 
                               Y = as.numeric(as.matrix((TempMean_New +  c(TempRange_New$X,0))) %*% t(Ax2_NewToOld)) )
          Cross_YU= data.frame(X = as.numeric(as.matrix((TempMean_New +  c(0,TempRange_New$Y))) %*% t(Ax1_NewToOld)) , 
                               Y = as.numeric(as.matrix((TempMean_New +  c(0,TempRange_New$Y))) %*% t(Ax2_NewToOld)) )
          Cross_YD= data.frame(X = as.numeric(as.matrix((TempMean_New -  c(0,TempRange_New$Y))) %*% t(Ax1_NewToOld)) , 
                               Y = as.numeric(as.matrix((TempMean_New -  c(0,TempRange_New$Y))) %*% t(Ax2_NewToOld)) )
          
          # Compute box points (4 points making up a box around synapse mean position)
          Box_UL= data.frame(X = as.numeric(as.matrix((TempMean_New +  TempRange_New*c(-1,1) )) %*% t(Ax1_NewToOld)) , 
                             Y = as.numeric(as.matrix((TempMean_New +  TempRange_New*c(-1,1) )) %*% t(Ax2_NewToOld)) )
          Box_UR= data.frame(X = as.numeric(as.matrix((TempMean_New +  TempRange_New*c(1,1) )) %*% t(Ax1_NewToOld)) , 
                             Y = as.numeric(as.matrix((TempMean_New +  TempRange_New*c(1,1) )) %*% t(Ax2_NewToOld)) )
          Box_LL= data.frame(X = as.numeric(as.matrix((TempMean_New +  TempRange_New*c(-1,-1) )) %*% t(Ax1_NewToOld)) , 
                             Y = as.numeric(as.matrix((TempMean_New +  TempRange_New*c(-1,-1) )) %*% t(Ax2_NewToOld)) )
          Box_LR= data.frame(X = as.numeric(as.matrix((TempMean_New +  TempRange_New*c(1,-1) )) %*% t(Ax1_NewToOld)) , 
                             Y = as.numeric(as.matrix((TempMean_New +  TempRange_New*c(1,-1) )) %*% t(Ax2_NewToOld)) )
          
          
          
        ############# Plot some random examples ############# 
          
          # Debugging plot
          PLOT=0
          if ( mod(ttt*lll*nnn,50) ==0 | abs(TempMean$X)<150 | (TempRange_New$X*2/Layer_Length_Temp)>2000){
            p1<-ggplot() +  geom_path(data=OUTLINE, aes(x=c1, y=c2), size = 1, color="red") +
              geom_point(data=NeuronData, aes(x=X,y=Y)) + 
              geom_path(data=MID, aes(x=X, y=Y), size = 1, color="orange") + 
              geom_path(data=Mid_subset, aes(x=X, y=Y), size = 1, color="green") + 
              geom_point(data=TempMean, aes(x=X,y=Y,size=50),color="red") + 
              geom_point(data=Cross_XL, aes(x=X,y=Y,size=50),color="cornflowerblue") + 
              geom_point(data=Cross_XR, aes(x=X,y=Y,size=50),color="blueviolet") + 
              geom_point(data=Cross_YU, aes(x=X,y=Y,size=50),color="darkred") + 
              geom_point(data=Cross_YD, aes(x=X,y=Y,size=50),color="gold") + 
              geom_point(data=Box_UL, aes(x=X,y=Y,size=50),color="darkslateblue")  + 
              geom_point(data=Box_UR, aes(x=X,y=Y,size=50),color="darkorange") + 
              geom_point(data=Box_LL, aes(x=X,y=Y,size=50),color="deeppink") +
              geom_point(data=Box_LR, aes(x=X,y=Y,size=50),color="darkslategray") + coord_fixed(ratio = 1) +
              ggtitle(paste(NeuronData$type[1],"_Layer_",as.character(lll),sep=""))
            
            ggsave(paste(PlotDir, "Cloud_", NeuronData$type[1],"_Layer_",as.character(lll),"_", as.character(ttt*lll*nnn),".pdf",sep=""),
                   plot = p1, device='pdf', scale = 1, width = 8, height = 5, units ="in", dpi = 300, limitsize = TRUE)
            
            remove(p1)
          
            }
          
          
        ############# Add data to dataframe ############# 
          
          TempDF <- data.frame(bodyid=NeuronData$bodyid[1], type=NeuronData$type[1], Layer=NeuronData$Layer[1], 
                               Side=NeuronData$Side[1], PBglom=NeuronData$PBglom[1], FBcol=NeuronData$FBcol[1],
                                         X_Mean=TempMean$X, Y_Mean=Mean_YY, Z_Mean=TempMean$Y, 
                                         X_STD_Orth=TempSTD_New$X_New, Z_STD_Orth=TempSTD_New$Y_New,
                                         X_HW_Orth=TempRange_New$X, Z_HW_Orth=TempRange_New$Y,
                                         Cross_XL_X=Cross_XL$X,    Cross_XL_Z=Cross_XL$Y,
                                         Cross_XR_X=Cross_XR$X,    Cross_XR_Z=Cross_XR$Y,
                                         Cross_YU_X=Cross_YU$X,    Cross_YU_Z=Cross_YU$Y,
                                         Cross_YD_X=Cross_YD$X,    Cross_YD_Z=Cross_YD$Y,
                                         Box_UL_X=Box_UL$X,        Box_UL_Z=Box_UL$Y,
                                         Box_UR_X=Box_UR$X,        Box_UR_Z=Box_UR$Y,
                                         Box_LL_X=Box_LL$X,        Box_LL_Z=Box_LL$Y,
                                         Box_LR_X=Box_LR$X,        Box_LR_Z=Box_LR$Y,
                                         LayerLength = sum((diff(MID$X)^2 + diff(MID$Y)^2)^0.5),
                                         NumberOfSyns = length(NeuronData$X) )

          Syn_Distribution=rbind(Syn_Distribution,TempDF)
          
        }
      }
    }
  }
  return(Syn_Distribution)
}

