##### A collection of functions for analyzing dFB sleep circuit connectivity  ######################################


Connectivity_VectorsDistance <- function(connectTable_Vectors, Method, MetaNames){
  #' Function for computing distance between input/ouput connectivity vectors
  
  # Compute similarity matrix
  if (Method=='binary'){
    Sim=1-bin_dist(connectTable_Vectors,0)
  } else if (Method=='cosine') {
    Sim = 1-cos_dist(connectTable_Vectors)
  } 
  
  # Convert to data frame
  Sim <- as.matrix(Sim)
  Sim=melt(Sim)
  colnames(Sim) <- c("bodyid", "bodyid2", "Similarity")
  
  # Assign names
  Names=distinct(MetaNames[c("bodyid","name","type")])
  Sim=inner_join(Sim,Names, by="bodyid")
  colnames(Names)=c("bodyid2","name2","type2")
  Sim=inner_join(Sim,Names, by="bodyid2")
  return(Sim)
}


Plot_SleepWake_Graph <- function(GraphTable2Plot, Pdir){
  #' Function for making network graph of dFP sleep-wake type connectivity
  
  # Set node positions
  IDs_Type=sort(unique(c(unique(GraphTable2Plot$type.from),unique(GraphTable2Plot$type.to))))
  NODES=data.frame(name=IDs_Type)
  NODES$name=NODES$name[c(7,3,1,4,8,2,5,9,11,6,10)]
  NODES$angle=NA
  NODES$x=NA
  NODES$y=NA
  NODES$angle[c(1,2,3,4,5)]=seq(from = 0, to = pi , by = pi/(4))
  NODES$x=cos(NODES$angle)*200
  NODES$y=sin(NODES$angle)*100
  NODES$x[6]=-125
  NODES$y[6]=0
  NODES$x[7]=125
  NODES$y[7]=0
  NODES$x[8]=0
  NODES$y[8]=-220/1.2
  NODES$x[9]=100
  NODES$y[9]=-300/1.3
  NODES$x[10]=0
  NODES$y[10]=-100/1.2
  NODES$x[11]=-100
  NODES$y[11]=-300/1.3
  NODES$y=-NODES$y
  
  # Set node colors
  NODES$COLORZ="dodgerblue3"
  NODES$COLORZ[NODES$name == "FB6H" | NODES$name == "FB7B"]="magenta3"
  
  # Set edges
  EDGES=GraphTable2Plot[c("type.from","type.to","weightRelative")]
  colnames(EDGES)=c("from","to","weight")
  EDGES$from=as.character(EDGES$from)
  EDGES$to=as.character(EDGES$to)
  EDGES$COLORZ="Sleep to Sleep"
  EDGES$COLORZ[EDGES$from == "FB6H" | EDGES$from == "FB7B"]="Wake to sleep"
  EDGES$COLORZ[EDGES$to   == "FB6H" | EDGES$to   == "FB7B"]="Sleep to Wake"
  EDGES$COLORZ=factor(EDGES$COLORZ,levels=(rev(c("Sleep to Wake", "Wake to sleep", "Sleep to Sleep"))))
  
  # Make ggraph object
  graph = tbl_graph(NODES, EDGES)
  
  # Chose colors  
  pcCols <- paletteer_d("Polychrome::palette36")
  col_vector=rev(c("dodgerblue3","magenta3","gray50"))
  
  # plot network graph and save
  Px=ggraph(graph,layout="manual",x=NODES$x,y=NODES$y) + 
    geom_edge_fan( aes(width=weight, color=COLORZ), 
                   arrow = arrow(length = unit(2, 'mm'),ends = "last", type = "closed"), end_cap = circle(1, 'cm'), alpha=1,strength=1) +
    geom_node_point(color=NODES$COLORZ, size=18) +     theme_classic() + 
    geom_node_text(aes(label=name),angle=0,size=2, nudge_y = c(rep(0.06,18),rep(-0.06,9)), color='black' ) +
    theme(legend.text=element_text(size=6),legend.title=element_text(size=6),
          axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),axis.title.y=element_blank()) +
    scale_edge_width(range = c(1.25,3)) + scale_edge_color_manual(values=col_vector)
  ggsave(paste(Pdir, "SleepWake_T2T_Graph.png",sep=""), plot = Px, device='png', scale = 1, width = 8.2, height = 8, units ="in", dpi = 500, limitsize = TRUE)
  return(Px)
  
}


OutputTyping <- function(DF_ToType){
  #' Function for custom retyping of sleep-wake output types
  DF_ToType$SuperType_Custom=NA
  DF_ToType$SuperType_Custom[startsWith(DF_ToType$type.to,"5-HT") |
                               startsWith(DF_ToType$type.to,"SI") |
                               startsWith(DF_ToType$type.to,"SL") |
                               startsWith(DF_ToType$type.to,"SM") |
                               startsWith(DF_ToType$type.to,"DGI") |
                               startsWith(DF_ToType$type.to,"IB042") ]="SNP"
  DF_ToType$SuperType_Custom[startsWith(DF_ToType$type.to,"EL") |
                               startsWith(DF_ToType$type.to,"ER") |
                               startsWith(DF_ToType$type.to,"Ex") |
                               startsWith(DF_ToType$type.to,"TuB")|
                               startsWith(DF_ToType$type.to,"AOT")]="EB and BU"
  DF_ToType$SuperType_Custom[startsWith(DF_ToType$type.to,"PF") |
                               startsWith(DF_ToType$type.to,"hDelta") |
                               startsWith(DF_ToType$type.to,"vDelta") |
                               startsWith(DF_ToType$type.to,"FS") |
                               startsWith(DF_ToType$type.to,"FR") |
                               startsWith(DF_ToType$type.to,"FC")]="FB Columnar"
  DF_ToType$SuperType_Custom[startsWith(DF_ToType$type.to,"FB") | startsWith(DF_ToType$type.to,"OA")]="FB Tangential"
  DF_ToType$SuperType_Custom[startsWith(DF_ToType$type.to,"LAL")]="LAL"
  DF_ToType$SuperType_Custom[startsWith(DF_ToType$type.to,"CRE")]="CRE"
  DF_ToType$SuperType_Custom[startsWith(DF_ToType$type.to,"ATL") | startsWith(DF_ToType$type.to,"LH") |
                               startsWith(DF_ToType$type.to,"PAM") | startsWith(DF_ToType$type.to,"PPL204")]="Olfactory"
  DF_ToType$SuperType_Custom[is.na(DF_ToType$SuperType_Custom)]="Unknown Type"
  DF_ToType$SuperType_Custom=factor(DF_ToType$SuperType_Custom, 
                                    levels=rev(c("FB Tangential","FB Columnar","EB and BU","SNP","LAL","CRE","Olfactory","Unknown Type")))
  return(DF_ToType)
}


InputTyping <- function(DF_FromType){
  #' Function for custom retyping of sleep-wake input types
  DF_FromType$SuperType_Custom=NA
  DF_FromType$SuperType_Custom[startsWith(DF_FromType$type.from,"5-HT") |
                                 startsWith(DF_FromType$type.from,"SI") |
                                 startsWith(DF_FromType$type.from,"SL") |
                                 startsWith(DF_FromType$type.from,"SM") |
                                 startsWith(DF_FromType$type.from,"DGI") |
                                 startsWith(DF_FromType$type.from,"IB042") ]="SNP"
  DF_FromType$SuperType_Custom[startsWith(DF_FromType$type.from,"EL") |
                                 startsWith(DF_FromType$type.from,"ER") |
                                 startsWith(DF_FromType$type.from,"Ex") |
                                 startsWith(DF_FromType$type.from,"TuB")|
                                 startsWith(DF_FromType$type.from,"AOT")|
                                 startsWith(DF_FromType$type.from,"PEN")|
                                 startsWith(DF_FromType$type.from,"EPG")|
                                 startsWith(DF_FromType$type.from,"Delta7")]="EB + PB + BU"
  DF_FromType$SuperType_Custom[startsWith(DF_FromType$type.from,"PF") |
                                 startsWith(DF_FromType$type.from,"hDelta") |
                                 startsWith(DF_FromType$type.from,"vDelta") |
                                 startsWith(DF_FromType$type.from,"FS") |
                                 startsWith(DF_FromType$type.from,"FR") |
                                 startsWith(DF_FromType$type.from,"FC") |
                                 startsWith(DF_FromType$type.from,"PFN")]="FB Columnar"
  DF_FromType$SuperType_Custom[startsWith(DF_FromType$type.from,"FB") | startsWith(DF_FromType$type.from,"OA-VPM3")]="FB Tangential"
  DF_FromType$SuperType_Custom[startsWith(DF_FromType$type.from,"CRE")]="CRE"
  DF_FromType$SuperType_Custom[startsWith(DF_FromType$type.from,"LCN")]="LAL"
  DF_FromType$SuperType_Custom[startsWith(DF_FromType$type.from,"ATL") | startsWith(DF_FromType$type.from,"LH") |
                                 startsWith(DF_FromType$type.from,"PAM") | startsWith(DF_FromType$type.from,"PPL204")]="Olfactory"
  DF_FromType$SuperType_Custom[is.na(DF_FromType$SuperType_Custom)]="Unknown Type"
  DF_FromType$SuperType_Custom=factor(DF_FromType$SuperType_Custom, 
                                      levels=rev(c("FB Tangential","FB Columnar","EB + PB + BU","SNP","LAL","CRE", "Olfactory","Unknown Type")))
  return(DF_FromType)
}


DownstreamLayer <- function(SW_Downstream_Filt_Layer){
  #' Function for finding which FB layer each sleep-wake type projects most strongly to
  SW_Downstream_Layer_Sum=SW_Downstream_Filt_Layer %>% group_by(type.from, type.to, supertype3.to) %>%
    summarize(weightRelative_path=sum(weightRelative_path))
  SW_Downstream_Layer_Sum=subset(SW_Downstream_Layer_Sum, (startsWith(type.to,"FB") |  startsWith(type.to,"OA-") | 
                                                             startsWith(type.to,"ExR")) & !supertype3.to=="Unassigned")
  MaxOutputLayer= SW_Downstream_Layer_Sum %>% group_by(type.to) %>% slice(which.max(weightRelative_path))
  MaxOutputLayer$Layer=substr(MaxOutputLayer$type.from, 3, 3) 
  MaxOutputLayer=MaxOutputLayer[c("type.to","Layer")]
  return(MaxOutputLayer)
}


UpstreamLayer <- function(SW_Upstream_Filt_Layer){
  #' Function for finding which FB layer each sleep-wake types receives strongest input from
  SW_Upstream_Layer=subset(SW_Upstream_Filt_Layer, weightRelative_path>PathwayThresh)
  SW_Upstream_Layer_Sum=SW_Upstream_Layer %>% group_by(type.from, type.to, supertype3.from) %>%
    summarize(weightRelative_path=sum(weightRelative_path))
  SW_Upstream_Layer_Sum=subset(SW_Upstream_Layer_Sum, (startsWith(type.from,"FB") | startsWith(type.from,"OA") |
                                                         startsWith(type.to,"ExR")) & !supertype3.from=="Unassigned")
  MaxInputLayer= SW_Upstream_Layer_Sum %>% group_by(type.from) %>% slice(which.max(weightRelative_path))
  MaxInputLayer$Layer=substr(MaxInputLayer$type.to, 3, 3) 
  MaxInputLayer=MaxInputLayer[c("type.from","Layer")]
  return(MaxInputLayer)
}


Plot_DownStream_Pathways <- function(SW_Downstream, PathwayThresh, Step, PlotDir){
  #' Function for plotting pathways downstream of sleep-wake types

  # Filter out weak pathways
  SW_Downstream_Filt=subset(SW_Downstream, weightRelative_path > PathwayThresh) 
  
  # Get rid of pathways longer than "Step" length, but keep them as NANs so they show up in plots
  SW_Downstream_Filt$weightRelative_path[SW_Downstream_Filt$n_steps>Step]=NA 
  SW_Downstream_Filt_Sum=SW_Downstream_Filt %>% group_by(type.from, type.to, supertype3.to) %>% 
    summarize(weightRelative_path=sum(weightRelative_path, na.rm=TRUE)) 
  SW_Downstream_Filt_Sum$weightRelative_path[SW_Downstream_Filt_Sum$weightRelative_path==0]=NA 
  
  # Group downstream neurons into custom supertypes, and see which custom types are included in the unknown category to make sure none were missed.
  SW_Downstream_Filt_Sum=OutputTyping(SW_Downstream_Filt_Sum)
  ExcludiedDownTypes=unique(SW_Downstream_Filt_Sum$type.to[SW_Downstream_Filt_Sum$SuperType_Custom == "Unknown Type"])
  
  # Group downstream types by supertype and compute average pathway weight and the number of target types
  SW_Downstream_Filt_Sum_Bar=subset(SW_Downstream_Filt_Sum, !is.na(weightRelative_path)) %>% 
    group_by(type.from, SuperType_Custom) %>% summarize(weightRelative_path=mean(weightRelative_path), n=n())
  
  
  ##########################################################################################################
  #### Make bar plots of downstream supertypes' average pathway weight and the number of target types ######
  
  # Set palette
  BarPalette=rev(color(c("mediumpurple4","mediumpurple1","forestgreen","darkorange1","bisque3","gray50","thistle4","black")))
  
  # Plot average pathway weight (by supertype)
  Pa=ggplot(SW_Downstream_Filt_Sum_Bar, aes(type.from)) +   geom_bar( aes(weight=weightRelative_path, fill = SuperType_Custom)) +
    scale_fill_manual(values=BarPalette, drop=FALSE) +  theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    ylab("Mean pathway relative weight ") + ylim(0,0.2)
  ggsave(paste(PlotDir, "SleepWake_Downstream_PathwayWeight","_Step",as.character(Step),".png",sep=""),
         plot = Pa, device='png', scale = 1, width =8, height =4, units ="in", dpi = 500, limitsize = TRUE)
  print(Pa)
  
  # Plot number of downstream types targeted (by supertype)
  Pb=ggplot(SW_Downstream_Filt_Sum_Bar, aes(type.from)) +   geom_bar( aes(weight=n,fill = SuperType_Custom)) +
    scale_fill_manual(values=BarPalette, drop=FALSE) +  theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    ylab("# downstream types") + ylim(0,100)
  ggsave(paste(PlotDir, "SleepWake_Downstream_NumberOfTypes","_Step",as.character(Step),".png",sep=""),
         plot = Pb, device='png', scale = 1, width =8, height = 4, units ="in", dpi = 500, limitsize = TRUE)
  print(Pb)
  
  
  ########################################################################################################
  #### Plot pathway weight matrix, showing pathway strength from sleep-wake types to downstream types ####
  
  PlotData=subset(SW_Downstream_Filt_Sum, !(startsWith(type.to,"FB") | startsWith(type.to,"OA-")) & !supertype3.to=="Unassigned")
  PlotData$SuperType_Custom=factor(PlotData$SuperType_Custom, levels=rev(levels(PlotData$SuperType_Custom)))
  
  Pc=ggplot(PlotData) + geom_tile(aes(type.to,type.from,fill=weightRelative_path)) +
    scale_fill_paletteer_c("ggthemes::Blue", limits=c(0,0.1), breaks = c(0,0.05,0.1), oob=squish, na.value="gray50") +
    theme_classic() + theme(axis.text.x = element_text(angle = 90))  + 
    facet_grid(cols=vars(SuperType_Custom),space = "free",  scales = "free" ) + 
    theme(strip.background = element_blank(), panel.border = element_rect(colour = "grey50", fill = NA, size = 0.5), 
          panel.spacing = unit(0, "lines"), aspect.ratio = 1) 
  ggsave(paste(PlotDir, "SleepWake_DownstreamPaths_noFBt","_Step",as.character(Step),".png",sep=""),
         plot = Pc, device='png', scale = 1, width =12, height = 6, units ="in", dpi = 500, limitsize = TRUE)
  print(Pc)
  
  # Get the FB layer targeted most strongly be each sleep-wake type (for plotting)
  MaxOutputLayer=DownstreamLayer(SW_Downstream_Filt)
  
  # remake matrix plots but with faceting by region 
  PlotData=subset(SW_Downstream_Filt_Sum, (startsWith(type.to,"FB") | startsWith(type.to,"OA-")) & !supertype3.to=="Unassigned")
  PlotData=merge(PlotData, MaxOutputLayer, by="type.to")
  PlotData$Layer=factor(PlotData$Layer, levels=c("6","7"))
  
  Pd=ggplot(PlotData) + geom_tile(aes(type.to,type.from,fill=weightRelative_path)) +
    scale_fill_paletteer_c("ggthemes::Blue", limits=c(0,0.1), breaks = c(0,0.05,0.1), oob=squish, na.value="gray50") +
    theme_classic() + theme(axis.text.x = element_text(angle = 90))  + 
    facet_grid(cols=vars(Layer),space = "free",  scales = "free" ) + 
    theme(strip.background = element_blank(), panel.border = element_rect(colour = "grey50", fill = NA, size = 0.5), 
          panel.spacing = unit(0, "lines"), aspect.ratio = 1) 
  ggsave(paste(PlotDir, "SleepWake_DownstreamPaths_FBt","_Step",as.character(Step),".png",sep=""),
         plot = Pd, device='png', scale = 1, width =12, height = 6, units ="in", dpi = 500, limitsize = TRUE)
  print(Pd)
  
  return(ExcludiedDownTypes)
}


Plot_UpStream_Pathways <- function(SW_Upstream, PathwayThresh, Step, PlotDir){
  #' Function for plotting pathways upstream of sleep-wake types
  
  # Filter out weak pathways
  SW_Upstream_Filt=subset(SW_Upstream, weightRelative_path>PathwayThresh) 
  
  # Get rid of pathways longer than "Step" length, but keep them as NANs so they show up in plots.
  SW_Upstream_Filt$weightRelative_path[SW_Upstream_Filt$n_steps>Step]=NA 
  SW_Upstream_Filt_Sum=SW_Upstream_Filt %>% group_by(type.from, type.to, supertype3.from) %>%
    summarize(weightRelative_path=sum(weightRelative_path, na.rm=TRUE)) 
  SW_Upstream_Filt_Sum$weightRelative_path[SW_Upstream_Filt_Sum$weightRelative_path==0]=NA
  
  # Group Upstream neurons into custom supertypes 
  SW_Upstream_Filt_Sum=InputTyping(SW_Upstream_Filt_Sum)
  ExcludedTypes=unique(SW_Upstream_Filt_Sum$type.from[SW_Upstream_Filt_Sum$SuperType_Custom == "Unknown Type"])
  
  # Group Upstream types by supertype and compute average pathway weight and the number of target types
  SW_Upstream_Filt_Sum_Bar=subset(SW_Upstream_Filt_Sum, !is.na(weightRelative_path)) %>% 
    group_by(type.to, SuperType_Custom) %>% summarize(weightRelative_path=mean(weightRelative_path), n=n())
  
  
  ########################################################################################################
  #### Make Bar plots of upstream supertypes' average pathway weight and the number of target types ######
  
  # Set palette
  BarPalette=rev(color(c("mediumpurple4","mediumpurple1","forestgreen","darkorange1","bisque3","gray50","thistle4","black")))
  
  # Plot average pathway weight (by supertype)
  Pa=ggplot(SW_Upstream_Filt_Sum_Bar, aes(type.to)) +   geom_bar( aes(weight=weightRelative_path,fill = SuperType_Custom)) +
    scale_fill_manual(values=BarPalette, drop=FALSE) +  theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    ylab("Mean pathway relative weight ") + ylim(0,0.2)
  ggsave(paste(PlotDir, "SleepWake_Upstream_PathwayWeight","_Step",as.character(Step),".pdf",sep=""),
         plot = Pa, device='pdf', scale = 1, width =8, height =4, units ="in", dpi = 500, limitsize = TRUE)
  print(Pa)
  
  # Plot number of upstream types targeted (by supertype)
  Pb=ggplot(SW_Upstream_Filt_Sum_Bar, aes(type.to)) +   geom_bar( aes(weight=n,fill = SuperType_Custom)) +
    scale_fill_manual(values=BarPalette, drop=FALSE) +  theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    ylab("# Upstream types") + ylim(0,50)
  ggsave(paste(PlotDir, "SleepWake_Upstream_NumberOfTypes","_Step",as.character(Step),".pdf",sep=""),
         plot = Pb, device='pdf', scale = 1, width =8, height = 4, units ="in", dpi = 500, limitsize = TRUE)
  print(Pb)
  
  
  ########################################################################################################
  #### Plot pathway weight matrix, showing pathway strength from upstream types to sleep-wake types ######
  
  # Plot pathway weights to non FB tanengtial types 
  PlotData=subset(SW_Upstream_Filt_Sum, !startsWith(type.from,"FB") & !startsWith(type.from,"OA") & !supertype3.from=="Unassigned")
  PlotData$SuperType_Custom=factor(PlotData$SuperType_Custom, levels=rev(levels(PlotData$SuperType_Custom)))
  
  Pc=ggplot(PlotData) + geom_tile(aes(type.to,type.from,fill=weightRelative_path)) +
    scale_fill_paletteer_c("ggthemes::Blue", limits=c(0,0.1), breaks = c(0,0.05,0.1), oob=squish, na.value="gray50") +
    theme_classic() + theme(axis.text.x = element_text(angle = 90))  + 
    facet_grid(rows=vars(SuperType_Custom),space = "free",  scales = "free" ) + 
    theme(strip.background = element_blank(), panel.border = element_rect(colour = "grey50", fill = NA, size = 0.5),
          panel.spacing = unit(0, "lines"), aspect.ratio = 1) 
  ggsave(paste(PlotDir, "SleepWake_UpstreamPaths_noFBt","_Step",as.character(Step),".pdf",sep=""),
         plot = Pc, device='pdf', scale = 1, width =6, height = 12, units ="in", dpi = 500, limitsize = TRUE)
  print(Pc)
  
  # Get the FB layer that targets each sleep-wake type most strongly (for faceting)
  MaxInputLayer=UpstreamLayer(SW_Upstream_Filt)
  
  # Plot pathway weights to FB tanengtial types 
  PlotData=subset(SW_Upstream_Filt_Sum, (startsWith(type.from,"FB") | startsWith(type.from,"OA")) & !supertype3.from=="Unassigned")
  PlotData=merge(PlotData, MaxInputLayer, by="type.from")
  PlotData$Layer=factor(PlotData$Layer, levels=c("6","7"))
  
  Pd=ggplot(PlotData) + geom_tile(aes(type.to,type.from,fill=weightRelative_path)) +
    scale_fill_paletteer_c("ggthemes::Blue", limits=c(0,0.1), breaks = c(0,0.05,0.1), oob=squish, na.value="gray50") +
    theme_classic() + theme(axis.text.x = element_text(angle = 90))  + 
    facet_grid(rows=vars(Layer),space = "free",  scales = "free" ) + 
    theme(strip.background = element_blank(), panel.border = element_rect(colour = "grey50", fill = NA, size = 0.5), 
          panel.spacing = unit(0, "lines"), aspect.ratio = 1) 
  ggsave(paste(PlotDir, "SleepWake_UpstreamPaths_FBt","_Step",as.character(Step),".pdf",sep=""),
         plot = Pd, device='pdf', scale = 1, width =6, height = 12, units ="in", dpi = 500, limitsize = TRUE)
  print(Pd)
  
  return(ExcludedTypes)
}






