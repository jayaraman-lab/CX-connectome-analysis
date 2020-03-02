# Functions for plotting Fb synapse distributions




Get_SynLayerDistribution <- function(All_Neurons, Thresh){
  
  
  Syn_Distribution <- data.frame(bodyid=as.numeric(), type=as.character(), Layer=as.character(), 
                                Side=as.character(), PBglom=as.character(), FBcol=as.character(),
                                X_Mean=as.numeric(), X_STD=as.numeric(), 
                                Y_Mean=as.numeric(), Y_STD=as.numeric(), 
                                Z_Mean=as.numeric(), Z_STD=as.numeric(),
                                X_Layer_Range=as.numeric(),Y_Layer_Range=as.numeric(),Z_Layer_Range=as.numeric(),
                                X_Range=as.numeric(),Y_Range=as.numeric(),Z_Range=as.numeric()) 
  
  Neuron_Types=unique(All_Neurons$type)
  Layers_All=sort(unique(All_Neurons$Layer))
  
  # Loop over neuron types
  for (ttt in 1:length(Neuron_Types)){
    SubData=subset(All_Neurons, type == Neuron_Types[ttt])
    NeuronIDs=unique(SubData$bodyid)
    
    # Loop of layers
    for (lll in 1:9){
      LayerData=subset(SubData, Layer==Layers_All[lll])
      
      X_Layer_Range = abs(quantile(LayerData$X, probs = 0.01) - quantile(LayerData$X, probs = 0.99) )
      Y_Layer_Range = abs(quantile(LayerData$Y, probs = 0.01) - quantile(LayerData$Y, probs = 0.99) )
      Z_Layer_Range = abs(quantile(LayerData$Z, probs = 0.01) - quantile(LayerData$Z, probs = 0.99) )
      
      
      # Loop over individual neurons
      for (nnn in 1:length(NeuronIDs)){
        NeuronData=subset(LayerData, bodyid==NeuronIDs[nnn])
        
        
        if (length(NeuronData$X)>Thresh){
          
          # Mean and std
          TempMean=colMeans(NeuronData[c("X","Y","Z")])
          TempSTD=sapply(NeuronData[c("X","Y","Z")], sd)
          
          X_Range = abs(quantile(NeuronData$X, probs = 0.01) - quantile(NeuronData$X, probs = 0.99) )
          Y_Range = abs(quantile(NeuronData$Y, probs = 0.01) - quantile(NeuronData$Y, probs = 0.99) )
          Z_Range = abs(quantile(NeuronData$Z, probs = 0.01) - quantile(NeuronData$Z, probs = 0.99) )
          
          
          TempDF=data.frame(bodyid=NeuronData$bodyid[1], type=NeuronData$type[1], Layer=NeuronData$Layer[1], 
                            Side=NeuronData$Side[1], PBglom=NeuronData$PBglom[1], FBcol=NeuronData$FBcol[1],
                            X_Mean=TempMean["X"], X_STD=TempSTD["X"], 
                            Y_Mean=TempMean["Y"], Y_STD=TempSTD["Y"], 
                            Z_Mean=TempMean["Z"], Z_STD=TempSTD["Z"],
                            X_Layer_Range=X_Layer_Range,Y_Layer_Range=Y_Layer_Range,Z_Layer_Range=Z_Layer_Range,
                            X_Range=X_Range,Y_Range=Y_Range,Z_Range=Z_Range) 
          
          Syn_Distribution=rbind(Syn_Distribution,TempDF)
          
        }
      }
    }
  }
  return(Syn_Distribution)
}




Plot_Layer <- function(Plot_Syns, Layer_To_Plot, Outline, Title, Grouping){
  
  # Get data
  PlotData=subset(Plot_Syns, Layer==Layer_To_Plot)
  
  # Get color map
  if (Grouping == "FBcol"){
    
    col_vector=brewer.pal(8, "Paired")
    col_vector[9]=col_vector[1]
    col_vector=col_vector[c(1,4,7,2,5,8,3,6,9)]
    
    ColorVar=PlotData$FBcol
    
  } else if (Grouping == "PBglom") {
    
    L_col_vector=brewer.pal(8, "Paired")
    L_col_vector[9]=L_col_vector[1]
    L_col_vector=L_col_vector[c(1,4,7,2,5,8,3,6,9)]
    R_col_vector=rev(L_col_vector)
    col_vector=c(L_col_vector, R_col_vector)
    
    ColorVar=PlotData$PBglom
  }
  
  # Make plot
  P1<-ggplot() + geom_point(data=PlotData, aes(x=X, y=Z, color=ColorVar) ,  size=2, alpha = 0.05, stroke = 0, shape=16) + 
    coord_fixed(ratio = 1) +  theme_void() + guides(fill=FALSE) + theme(legend.position = "none") +
    geom_path(data=Outline, aes(x=c1, y=c2), size = 1) +  scale_color_manual(values=col_vector) + xlim(c( min(Outline_FB_XZ$c1)-200 , max(Outline_FB_XZ$c1)+200)) +
    ggtitle(Title)
  
  return(P1)
}






PlotSyns_byLayer <- function(Plot_Syns, Type, PlotName, Grouping){
  
  
  P1 = Plot_Layer(Plot_Syns, "L1", Outline_L1_XZ, paste("L1 ", Type,sep=""), Grouping )
  P2 = Plot_Layer(Plot_Syns, "L2", Outline_L2_XZ, paste("L2 ", Type,sep=""), Grouping )
  P3 = Plot_Layer(Plot_Syns, "L3", Outline_L3_XZ, paste("L3 ", Type,sep=""), Grouping )
  P4 = Plot_Layer(Plot_Syns, "L4", Outline_L4_XZ, paste("L4 ", Type,sep=""), Grouping )
  P5 = Plot_Layer(Plot_Syns, "L5", Outline_L5_XZ, paste("L5 ", Type,sep=""), Grouping )
  P6 = Plot_Layer(Plot_Syns, "L6", Outline_L6_XZ, paste("L6 ", Type,sep=""), Grouping )
  P7 = Plot_Layer(Plot_Syns, "L7", Outline_L7_XZ, paste("L7 ", Type,sep=""), Grouping )
  P8 = Plot_Layer(Plot_Syns, "L8", Outline_L8_XZ, paste("L8 ", Type,sep=""), Grouping )
  P9 = Plot_Layer(Plot_Syns, "L9", Outline_L9_XZ, paste("L9 ", Type,sep=""), Grouping )
  
  P12<-ggarrange(plotlist=list(P1,P2,P3,P4,P5,P6,P7,P8,P9), ncol =3, nrow = 3)
  ggsave(paste(PlotDir, PlotName,".png",sep=""),
         plot = P12, device='png', scale = 1, width = 8, height = 5, units ="in", dpi = 300, limitsize = TRUE)
  
}





SetColumnOrder <- function(TypeX_Syns){
  TypeX_MeanPos= aggregate(TypeX_Syns[c("X","Y","Z")], by = list(bodyid=TypeX_Syns$bodyid), FUN=mean)
  colnames(TypeX_MeanPos)[c(2,3,4)]<-c("mean_X","mean_Y","mean_Z")
  TypeX_MeanPos$Order_X=order(TypeX_MeanPos$mean_X)
  TypeX_Syns = merge(TypeX_Syns, TypeX_MeanPos)
  return(TypeX_Syns)
}








PlotSyns_byColumn <- function(Plot_Syns, Title, PlotName){
  
  
  col_vector=brewer.pal(8, "Paired")
  col_vector[9]=col_vector[1]
  col_vector=col_vector[c(1,4,7,2,5,8,3,6,9)]
  
  P1<-ggplot() + geom_point(data=Plot_Syns, aes(x=X, y=-Y, color=FBcol) ,  size=2, alpha = 0.05, stroke = 0, shape=16) + 
    coord_fixed(ratio = 1) +  theme_void() + guides(fill=FALSE) + theme(legend.position = "top") +
    geom_path(data=Outline_FB_XY, aes(x=c1, y=c2), size = 1) +  scale_color_manual(values=col_vector) + guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
    ggtitle(Title)
  
  
  P2<-ggplot() + geom_point(data=Plot_Syns, aes(x=X, y=Z, color=FBcol) ,  size=2, alpha = 0.05, stroke = 0, shape=16) + 
    coord_fixed(ratio = 1) +  theme_void() + guides(fill=FALSE) + theme(legend.position = "top") +
    geom_path(data=Outline_FB_XZ, aes(x=c1, y=c2), size = 1) +  scale_color_manual(values=col_vector) + guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
    ggtitle(Title)
  
  
  P12<-ggarrange(plotlist=list(P1,P2), ncol =2, nrow = 1)
  ggsave(paste(PlotDir, PlotName,".png",sep=""),
         plot = P12, device='png', scale = 1, width = 8, height = 5, units ="in", dpi = 300, limitsize = TRUE)
  
}







PlotSyns_byGlom <- function(Plot_Syns, Title, PlotName){
  
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
  ggsave(paste(PlotDir, PlotName,".png",sep=""),
         plot = P12, device='png', scale = 1, width = 8, height = 12, units ="in", dpi = 300, limitsize = TRUE)
  
}



























