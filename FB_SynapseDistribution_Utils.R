################### A collection of functions for loading, processing, and plotting FB synapse distributions ##################



###############################################################################################################################
################## Function for loading FB synapses by layer ##################################################################



Get_FBlayer_ColumnarSyns <- function(ROI, Layer, IncludedTypes) {
  NamedBodies_LX=Get_AllNeurons_InRoi(ROI, FALSE)
  NamedBodies_LX=subset(NamedBodies_LX, NamedBodies_LX$bodytype %in% IncludedTypes)
  SynLocs=GetSynapseLocs(NamedBodies_LX$bodyid, ROI,TRUE)
  SynLocs$Layer=Layer
  return(SynLocs)
}



###############################################################################################################################
################# Function for plotting synapse distribution for all FB layers ################################################



PlotSyns_byLayer <- function(Plot_Syns, Type, PlotName, Grouping, DIR, TOPLOT_TF, TOPLOTLEGEND_TF){
  
  P1 = Plot_Layer(Plot_Syns, "L1", Outline_L1_XZ, paste("L1 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF)
  P2 = Plot_Layer(Plot_Syns, "L2", Outline_L2_XZ, paste("L2 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF)
  P3 = Plot_Layer(Plot_Syns, "L3", Outline_L3_XZ, paste("L3 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF)
  P4 = Plot_Layer(Plot_Syns, "L4", Outline_L4_XZ, paste("L4 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF)
  P5 = Plot_Layer(Plot_Syns, "L5", Outline_L5_XZ, paste("L5 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF)
  P6 = Plot_Layer(Plot_Syns, "L6", Outline_L6_XZ, paste("L6 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF)
  P7 = Plot_Layer(Plot_Syns, "L7", Outline_L7_XZ, paste("L7 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF)
  P8 = Plot_Layer(Plot_Syns, "L8", Outline_L8_XZ, paste("L8 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF)
  P9 = Plot_Layer(Plot_Syns, "L9", Outline_L9_XZ, paste("L9 ", Type,sep=""), Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF)
  
  P12<-ggarrange(plotlist=list(P9,P8,P7,P6,P5,P4,P3,P2,P1), ncol =3, nrow = 3)
  ggsave(paste(DIR, PlotName,".png",sep=""),
         plot = P12, device='png', scale = 1, width = 12, height = 10, units ="in", dpi = 500, limitsize = TRUE)
  
}



###############################################################################################################################
################# Function for plotting synapse distribution of each individual FB layer ######################################



Plot_Layer <- function(Plot_Syns, Layer_To_Plot, Outline, Title, Grouping, DIR, PlotName, TOPLOT_TF, TOPLOTLEGEND_TF){
  
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
  P1 <- Plot_Synapse_Distribution(PlotData, "XZ", Outline, ColorVar, col_vector, DIR, PlotName, Layer_To_Plot, TOPLOT_TF, TOPLOTLEGEND_TF)

  return(P1)
}



###############################################################################################################################
################# Generic plotting function for synapse distributions  ################################################



Plot_Synapse_Distribution <- function(PlotData, VIEW, Outline, ColorVar, col_vector, DIR, PlotName, PlotName2, TOPLOT_TF, TOPLOTLEGEND_TF){
  
  # Get view
  if (VIEW=="XY"){
    PlotData = PlotData[!colnames(PlotData)=="Z"]
    colnames(PlotData)[colnames(PlotData)=="Y"]="Z"
    PlotData$Z=-PlotData$Z
  }
  
  # Plot just the data for now
  P1 <- ggplot() + geom_point(data=PlotData, aes(x=X, y=Z, color=ColorVar) ,  size=1, alpha = 1, stroke = 0, shape=16) + 
    geom_path(data=Outline, aes(x=c1, y=c2), size = 2) + xlim(c(-9000, 9000)) +  coord_fixed(ratio = 1) 

  # Save with scaled visible if TOPLOTLEGEND_TF us TRUE
  if (TOPLOTLEGEND_TF==TRUE){
    ggsave(paste(DIR, PlotName,"_",PlotName2 , "_Legend.png",sep=""), plot = P1, device='png', scale = 1, width = 8, height = 5, units ="in", dpi = 500, limitsize = TRUE) 
  }  
  
  # Format the rest of the plot, get rid of axes
  P1 <- P1 + theme_void() + guides(fill=FALSE) + theme(legend.position = "none") +
      scale_color_manual(values=col_vector)  +  theme(
      panel.background = element_rect(fill = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      panel.grid.major = element_blank(), # get rid of major grid
      panel.grid.minor = element_blank(), # get rid of minor grid
      legend.background = element_rect(fill = "transparent"), # get rid of legend bg
      legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
    )
  
  # Save plot
  if (TOPLOT_TF==TRUE){
    ggsave(paste(DIR, PlotName,"_",PlotName2 , ".png",sep=""), plot = P1, device='png', scale = 1, width = 8, height = 5, units ="in", dpi = 500, limitsize = TRUE,  bg = "transparent")
  } 
  
  return(P1)
}



###############################################################################################################################
################# Function for getting color map of FB columns or PB glomeruli ################################################



GetColorMap <- function(Grouping){
  if (Grouping == "FBcol"){
    col_vector=brewer.pal(8, "Paired")
    col_vector[9]=col_vector[1]
    col_vector=col_vector[c(1,4,7,2,5,8,3,6,9)]
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
  
  # Get column vector
  col_vector=GetColorMap("FBcol")
  
  ### XY View
  P1 <- Plot_Synapse_Distribution(Plot_Syns, "XY", Outline_FB_XY, Plot_Syns$FBcol, col_vector, DIR, PlotName, "XY", TOPLOT_TF, TOPLOTLEGEND_TF)
  
  ### XZ View
  P1 <- Plot_Synapse_Distribution(Plot_Syns, "XZ", Outline_FB_XZ, Plot_Syns$FBcol, col_vector, DIR, PlotName, "XZ", TOPLOT_TF, TOPLOTLEGEND_TF)
  
}







PlotSyns_byGlom <- function(Plot_Syns, Title, PlotName, DIR){
  
  # Get neurons that innervate the left or right pb
  Plot_Syns_Lpb=subset(Plot_Syns, startsWith(PBglom,"L") )
  Plot_Syns_Lpb$PBglom=factor( Plot_Syns_Lpb$PBglom, levels= sort(unique(Plot_Syns_Lpb$PBglom)) )
  
  Plot_Syns_Rpb=subset(Plot_Syns, startsWith(PBglom,"R") )
  Plot_Syns_Rpb$PBglom=factor( Plot_Syns_Rpb$PBglom, levels= sort(unique(Plot_Syns_Rpb$PBglom)) )
  
  
  L_col_vector=brewer.pal(8, "Paired")
  L_col_vector[9]=L_col_vector[1]
  L_col_vector=L_col_vector[c(1,4,7,2,5,8,3,6,9)]
  R_col_vector=rev(L_col_vector)
  RL_col_Vector=c(L_col_vector, R_col_vector)
  
  
  # Make color map for each that maps PB glomeruli to their appropriate columns
  P1<-ggplot() + geom_point(data=Plot_Syns_Lpb, aes(x=X, y=-Y, color=PBglom) ,  size=2, alpha = 0.05, stroke = 0, shape=16) + 
    coord_fixed(ratio = 1) +  theme_void() + guides(fill=FALSE) + theme(legend.position = "top") +
    geom_path(data=Outline_FB_XY, aes(x=c1, y=c2), size = 1) +  scale_color_manual(values=L_col_vector) + guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
    ggtitle(paste("LEFT ", Title)) + xlim(c( min(Outline_FB_XZ$c1)-200 , max(Outline_FB_XZ$c1)+200)) 
  
  
  P2<-ggplot() + geom_point(data=Plot_Syns_Lpb, aes(x=X, y=Z, color=PBglom) ,  size=2, alpha = 0.05, stroke = 0, shape=16) + 
    coord_fixed(ratio = 1) +  theme_void() + guides(fill=FALSE) + theme(legend.position = "top") +
    geom_path(data=Outline_FB_XZ, aes(x=c1, y=c2), size = 1) +  scale_color_manual(values=L_col_vector) + guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
    ggtitle(paste("LEFT ", Title)) + xlim(c( min(Outline_FB_XZ$c1)-200 , max(Outline_FB_XZ$c1)+200))
  
  
  
  P3<-ggplot() + geom_point(data=Plot_Syns_Rpb, aes(x=X, y=-Y, color=PBglom) ,  size=2, alpha = 0.05, stroke = 0, shape=16) + 
    coord_fixed(ratio = 1) +  theme_void() + guides(fill=FALSE) + theme(legend.position = "top") +
    geom_path(data=Outline_FB_XY, aes(x=c1, y=c2), size = 1) +  scale_color_manual(values=R_col_vector) + guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
    ggtitle(paste("RIGHT ", Title)) + xlim(c( min(Outline_FB_XZ$c1)-200 , max(Outline_FB_XZ$c1)+200))
  
  
  P4<-ggplot() + geom_point(data=Plot_Syns_Rpb, aes(x=X, y=Z, color=PBglom) ,  size=2, alpha = 0.05, stroke = 0, shape=16) + 
    coord_fixed(ratio = 1) +  theme_void() + guides(fill=FALSE) + theme(legend.position = "top") +
    geom_path(data=Outline_FB_XZ, aes(x=c1, y=c2), size = 1) +  scale_color_manual(values=R_col_vector) + guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
    ggtitle(paste("RIGHT ", Title)) + xlim(c( min(Outline_FB_XZ$c1)-200 , max(Outline_FB_XZ$c1)+200))
  
  
  
  P5<-ggplot() + geom_point(data=Plot_Syns, aes(x=X, y=-Y, color=PBglom) ,  size=2, alpha = 0.05, stroke = 0, shape=16) + 
    coord_fixed(ratio = 1) +  theme_void() + guides(fill=FALSE) + theme(legend.position = "top") +
    geom_path(data=Outline_FB_XY, aes(x=c1, y=c2), size = 1) +  scale_color_manual(values=RL_col_Vector) + guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
    ggtitle(paste("RIGHT ", Title)) + xlim(c( min(Outline_FB_XZ$c1)-200 , max(Outline_FB_XZ$c1)+200))
  
  
  P6<-ggplot() + geom_point(data=Plot_Syns, aes(x=X, y=Z, color=PBglom) ,  size=2, alpha = 0.05, stroke = 0, shape=16) + 
    coord_fixed(ratio = 1) +  theme_void() + guides(fill=FALSE) + theme(legend.position = "top") +
    geom_path(data=Outline_FB_XZ, aes(x=c1, y=c2), size = 1) +  scale_color_manual(values=RL_col_Vector) + guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
    ggtitle(paste("RIGHT ", Title)) + xlim(c( min(Outline_FB_XZ$c1)-200 , max(Outline_FB_XZ$c1)+200))
  
  
  
  P12<-ggarrange(plotlist=list(P1,P2,P3,P4,P5,P6), ncol =2, nrow = 3)
  ggsave(paste(DIR, PlotName,".png",sep=""),
         plot = P12, device='png', scale = 1, width = 8, height = 12, units ="in", dpi = 300, limitsize = TRUE)
  
}


###############################################################################################################################
### Function for assigning FB column to PB-FB-XX and D0 neurons based on closest FX neuron  ##################################


Assign_FB_Columns <- function(N_Distribution, FX_Distribution){
  
  N_Distribution$Distance2FxCol=NA
  N_Distribution$FBcol=as.character(N_Distribution$FBcol)
  
  Uni_BodyID=unique(N_Distribution$bodyid)
  # Loop over each neuron
  for (nnn in 1:length(N_Distribution$X_Mean)){
    
    # Subset FX data based on layer
    TempCol_Positions=subset(FX_Distribution, Layer==N_Distribution$Layer[nnn])
    
    # Compute distance to each column center
    Temp_Distances= ((TempCol_Positions$X_Mean - N_Distribution$X_Mean[nnn])^2 + (TempCol_Positions$Z_Mean - N_Distribution$Z_Mean[nnn])^2)^0.5 
    
    # Get index of min
    Temp_Ind=which(Temp_Distances == min(Temp_Distances))
    N_Distribution$Distance2FxCol[nnn]=min(Temp_Distances)
    
    # Assign column
    N_Distribution$FBcol[nnn]=as.character(TempCol_Positions$FBcol[Temp_Ind])
  }
  
  return(N_Distribution)
  
}


###############################################################################################################################
### Function for getting the width of each neuron's synapse cloud in a direction tangent to the FB layer ######################


Get_SynLayerDistribution <- function(All_Neurons, Thresh, PlotDir){
  
  
  Syn_Distribution <- data.frame(bodyid=as.numeric(), type=as.character(), Layer=as.character(), 
                                 Side=as.character(), PBglom=as.character(), FBcol=as.character(),
                                 X_Mean=as.numeric(), Z_Mean=as.numeric(), 
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
        NeuronData=NeuronData[colnames(NeuronData) != "Y"]
        colnames(NeuronData)[colnames(NeuronData) == "Z"] <-"Y"
        
        if (length(NeuronData$X)>Thresh){
          
        ############# Project data into new coordinate frame whose main axis is locally parallel to FB layer ############# 
          
          # Get X range
          X_RangeMin = quantile(NeuronData$X, probs = 0.1) 
          X_RangeMax = quantile(NeuronData$X, probs = 0.9) 
          
          # Get slope of mid range 
          Mid_subset=subset(MID, X>X_RangeMin & X<X_RangeMax)
          TempSlope=median( diff(Mid_subset$Y)/diff(Mid_subset$X) )
          
          # Define new axis (vectors should be length 1 and orthogonal to keep shape of point cloud)
          Ax_1=c(1,TempSlope)
          Ax_1=Ax_1/sum(Ax_1^2)^0.5
          Ax_2=c(1, -Ax_1[1]/Ax_1[2])
          Ax_2=Ax_2/sum(Ax_2^2)^0.5
          
          # Keep axes orientation the same
          if (mean(NeuronData$X)>0){
            Ax_1=Ax_1*-1;
          }
          
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
          TempMean=as.data.frame(t(colMeans(NeuronData[c("X","Y")]))) # This
          
          # Mean and STD along new axes
          TempMean_New=as.data.frame(t(colMeans(NeuronData[c("X_New","Y_New")])))
          TempSTD_New=as.data.frame(t(sapply(NeuronData[c("X_New","Y_New")], sd)))
          
          # Ranges along new axes
          TempRange_New=data.frame(X = abs(quantile(NeuronData$X_New, probs = 0.02) - quantile(NeuronData$X_New, probs = 0.98))/2,
                                   Y = abs(quantile(NeuronData$Y_New, probs = 0.02) - quantile(NeuronData$Y_New, probs = 0.98))/2)
          
          
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
          if ( mod(ttt*lll*nnn,50) ==0){
            p1<-ggplot() +  geom_path(data=OUTLINE, aes(x=c1, y=c2), size = 1, color="red") +
              geom_point(data=NeuronData, aes(x=X,y=Y)) + 
              geom_path(data=MID, aes(x=X, y=Y), size = 1, color="orange") + coord_fixed(ratio = 1) + 
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
             
            ggsave(paste(PlotDir, NeuronData$type[1],"_Layer_",as.character(lll),"_", as.character(ttt*lll*nnn),".png",sep=""),
                   plot = p1, device='png', scale = 1, width = 8, height = 5, units ="in", dpi = 300, limitsize = TRUE)
            
            #print(paste("Plotting ttt=", as.character(ttt), "   lll = ", as.character(lll), "   nnn = ", as.character(nnn)))
            
            p2<-ggplot() + geom_point(data=NeuronData, aes(x=X_New,y=Y_New))
            
            ggarrange(p1,p2, nrow = 1, ncol = 2)
          
            }
          
          
        ############# Add data to dataframe ############# 
          
          TempDF <- data.frame(bodyid=NeuronData$bodyid[1], type=NeuronData$type[1], Layer=NeuronData$Layer[1], 
                               Side=NeuronData$Side[1], PBglom=NeuronData$PBglom[1], FBcol=NeuronData$FBcol[1],
                                         X_Mean=TempMean$X, Z_Mean=TempMean$Y, 
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

