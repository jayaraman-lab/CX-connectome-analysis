### This script contains functions for analysing and plotting PB-FB phase shifts
### and other phase-related analysis related to FB columns


Assign_Angles_ToNeurons <- function(Table_w_NeuronName){
  #' Function for assigning EPG and Delta7 angles to neurons
  Angle_Table=data.frame(PBglom.from=c("L9","L8","L7","L6","L5","L4","L3","L2","L1",
                                       "R1","R2","R3","R4","R5","R6","R7","R8","R9" ),
                         EPG_Angle.from=c(0, 45, 90, 135, 180, 225, 270, 315, 0,
                                          22.5, 67.5, 112.5, 157.5, 202.5, 247.5, 292.5, 337.5, 22.5),
                         D7_Angle.from=c(348.75, 33.75, 78.75, 123.75, 168.75, 213.75, 258.75, 303.75, 348.75,
                                         33.75, 78.75, 123.75, 168.75, 213.75, 258.75, 303.75, 348.75, 33.75))
  Table_w_NeuronName=inner_join(Table_w_NeuronName, Angle_Table, "PBglom.from")
  return(Table_w_NeuronName)
}



GetAllAngles <- function(TempConnectFun){
  #' Function for creating a distribution of angles according to the # of synapses at each angle
  if (nrow(TempConnectFun>0)){
    AllAngles=c()
    for (rrr in 1:nrow(TempConnectFun)){
      AllAngles=c(AllAngles, rep(TempConnectFun$EPG_Angle.from[rrr], TempConnectFun$ROIweight[rrr]) ) 
    }
  } else {
    AllAngles=NaN
  }
  return(AllAngles)
}


ComputeConsistency <- function(FB_Bodies_All, PB_FB_Outputs_N2N){
  #' Function for computing the propotation of postsynaptic neurons targets by a presynaptic PB-FB type
  Columnar_Types=sort(unique(FB_Bodies_All$type[startsWith(FB_Bodies_All$type, "PF" )     | startsWith(FB_Bodies_All$type, "FR" ) | 
                                                  startsWith(FB_Bodies_All$type, "FS" )     | startsWith(FB_Bodies_All$type, "FC" ) | 
                                                  startsWith(FB_Bodies_All$type, "vDelta" ) | startsWith(FB_Bodies_All$type, "hDelta")]))
  FBColumnarBodyCount=FB_Bodies_All %>% group_by(type) %>% summarize(NumofNeurons=n())
  FBColumnarBodyCount=subset(FBColumnarBodyCount, type %in% Columnar_Types)
  PrePostNum=PB_FB_Outputs_N2N %>% group_by(type.from,type.to) %>% summarise(PreNum=length(unique(from)), PostNum=length(unique(to)))
  colnames(FBColumnarBodyCount)=c("type.from","PreNumTotal")
  PrePostNum=merge(PrePostNum, FBColumnarBodyCount, by=c("type.from"))
  colnames(FBColumnarBodyCount)=c("type.to","PostNumTotal")
  PrePostNum=merge(PrePostNum, FBColumnarBodyCount, by=c("type.to"))
  PrePostNum$PostConsistency=PrePostNum$PostNum/PrePostNum$PostNumTotal
  return(PrePostNum)
}


PlotConnectMats <- function(PB_FB_Outputs_N2N, PrePostTypes, PlotDir){
  # Make directory just for these plots
  PlotDirMats=paste(PlotDir,"ConnectionMats/")
  if (!dir.exists(PlotDirMats)){dir.create(PlotDirMats)}
  
  # Loop over and plot N2N connection matrices to make sure no double-banded matrices remain
  PB_FB_Outputs_N2N$FBcol.from=factor(PB_FB_Outputs_N2N$FBcol.from, levels=c("C0","C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11","C12"))
  PB_FB_Outputs_N2N$FBcol.to=factor(PB_FB_Outputs_N2N$FBcol.to, levels=c("C0","C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11","C12"))
  
  for (ttt in 1:length(PrePostTypes$type.from)){
    Temp_Connect=subset(PB_FB_Outputs_N2N, type.from == PrePostTypes$type.from[ttt] & type.to == PrePostTypes$type.to[ttt])
    Pa=ggplot(Temp_Connect) + 

      scale_fill_gradient2(low="gray30", high="red", 
                           limits=c(0,max(Temp_Connect$ROIweight)), oob=squish)  +
      geom_tile(aes(name.to,name.from,fill=ROIweight)) + ylab(PrePostTypes$type.from[ttt]) + xlab(PrePostTypes$type.to[ttt]) +
      facet_grid(rows=vars(FBcol.from), cols=vars(FBcol.to), scales="free") +
      theme(strip.background = element_blank(), axis.text.x = element_text(angle = 90),
            panel.border = element_rect(colour = "grey50", fill = NA, size = 0.5), 
            panel.spacing = unit(0, "lines"), axis.ticks = element_blank(), aspect.ratio = 1) 
    ggsave(paste(PlotDirMats, "ConnectMat_",  as.character(ttt),"_", PrePostTypes$type.from[ttt], "_2_", PrePostTypes$type.to[ttt], ".png", sep=""),
           plot = Pa, device='png', scale = 1, width = 7, height = 5, units ="in", dpi = 200, limitsize = TRUE)
  }
  
}

