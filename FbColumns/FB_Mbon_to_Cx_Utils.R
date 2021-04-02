##### A collection of functions for analyzing MBON to FB connectivity  ######################################

Plot_Indirect_Matrices <- function(Matrix, OrderIn, OrderOut){
  # Function for plotting MBON-to-intermediate and intermediate-to-CX connection matrices
  Px=plotConnectivity(Matrix, slctROI = NULL, grouping = "type", connectionMeasure = "weightRelative",
                      xaxis = c("outputs"), facetInputs = NULL, facetOutputs = NULL, theme = theme_minimal(),
                      cmax = NULL, replacementLabels = NULL, orderIn = OrderIn, 
                      orderOut = OrderOut, legendName = NULL, showTable = "inputs")
  Px = Px + coord_fixed(ratio = 1) +
    theme(strip.background = element_blank(), 
          panel.border = element_rect(colour = "grey50", fill = NA), 
          panel.spacing = unit(0, "lines"),
          text = element_text(size=6)) 
  return(Px)
}


Compute_MBON_to_FBlayer <- function(Pathway, relativeWeighThresh){
  #' Function for compute the number of pathways from each MBON to each FB layer
  Pathways_Out=subset(Pathway, startsWith(databaseType.to,"FB") & weightRelative_N1>relativeWeighThresh & weightRelative_N2>relativeWeighThresh)
  Pathways_Out$FBlayer=paste("L", substr(Pathways_Out$databaseType.to,3,3), sep="")
  Pathways_Out=Pathways_Out %>% group_by(databaseType.from, FBlayer) %>%
    summarize(NumOfPaths=n(), NumOfFBt=length(unique(databaseType.to)), weightRelative_path=mean(weightRelative_path))
  Pathways_Out$FBlayer=factor(Pathways_Out$FBlayer, levels=c("L1","L2","L3","L4","L5","L6","L7","L8","L9"))
  return(Pathways_Out)
}


Plot_FbLayers_Targeted <- function(LayersTargeted){
  #' Function plotting the number of FB layers targetd by each MBON
  Px=ggplot(LayersTargeted) +
    geom_tile(aes(databaseType.from,FBlayer,fill=NumOfFBt)) + scale_y_discrete(drop=FALSE) + scale_x_discrete(drop=FALSE) +
    scale_fill_paletteer_c("ggthemes::Sunset-Sunrise Diverging", limits=c(0,14), breaks = c(0,7,14), oob=squish, na.value="gray50") +
    theme_classic() + theme(axis.text.x = element_text(angle = 90))  + 
    theme(strip.background = element_blank(), 
          panel.border = element_rect(colour = "grey50", fill = NA, size = 0.5), 
          panel.spacing = unit(0, "lines"))  + coord_fixed(ratio = 1)
  return(Px)
}


Get_Max_T2T <- function(MB2NonCX_Inp, Downstream_T2T_L2_Inp, relativeWeighThresh, WeightThresh){
  #' Function for getting  MBON-to-intermediate and intermediate-to-CX connection tables
  
  # Subset connection tables by weight thresholds and remove connections not belonging to MBON-->CX pathways
  L2_to_Cx=subset(Downstream_T2T_L2_Inp, databaseType.to %in% CXTypes  & weightRelative>relativeWeighThresh & weight>WeightThresh)
  L1_to_L2=subset(MB2NonCX_Inp, databaseType.to %in% L2_to_Cx$databaseType.from & weightRelative>relativeWeighThresh & weight>WeightThresh
                  & !startsWith(type.to,"MBON") & !databaseType.from=="MBON15-like")
  L2_to_Cx=subset(L2_to_Cx, databaseType.from %in% L1_to_L2$databaseType.to) 
  
  # Get type-to-type connections with max weights (as explained in main notebook)
  L1_to_L2_Max=L1_to_L2 %>% group_by(databaseType.from, databaseType.to) %>%  filter(weight == max(weight))
  L1_to_L2_Max$type.from=L1_to_L2_Max$databaseType.from
  L1_to_L2_Max$type.to=L1_to_L2_Max$databaseType.to
  
  L2_to_Cx_Max=L2_to_Cx %>% group_by(databaseType.from, databaseType.to) %>%  filter(weight == max(weight))
  L2_to_Cx_Max$type.from=L2_to_Cx_Max$databaseType.from
  L2_to_Cx_Max$type.to=L2_to_Cx_Max$databaseType.to
  
  OUT=list(L1_to_L2_Max, L2_to_Cx_Max)
  return(OUT)
}


Get_Direct_Graph <- function(Direct_Data, FigDir){
  #' Function containing low level code to build a network graph showing direct MBON to CX connections
  
  # Get node sizes (proportional to total outgoing/incoming weights)
  PostSize=Direct_Data %>% group_by(type.to) %>% summarize(Syns=sum(weight))
  colnames(PostSize)=c("name","Size")
  PreSize=Direct_Data %>% group_by(type.from) %>% summarize(Syns=sum(weight))
  colnames(PreSize)=c("name","Size")
  AllSize=rbind(PostSize, PreSize)
  AllSize$name[startsWith(AllSize$name,"MBON")]=substr(AllSize$name[startsWith(AllSize$name,"MBON")], 5, 6)
  
  # Pre Nodes
  MBON_Nodes=data.frame(name=(unique(Direct_Data$type.from) %>% as.character()))
  MBON_Nodes$name=substr(MBON_Nodes$name,5, 6)
  MBON_Nodes$x=0
  MBON_Nodes$y=(seq(from = 1, to = length(MBON_Nodes$name), by = 1))
  MBON_Nodes$y=(MBON_Nodes$y-mean(MBON_Nodes$y))*1.2
  
  # Post Nodes
  CX_Nodes=data.frame(name=rev((unique(Direct_Data$type.to) %>% as.character()))) 
  CX_Nodes$x=1
  CX_Nodes$y=rev(seq(from = 1, to = length(CX_Nodes$name), by = 1))
  CX_Nodes$y=CX_Nodes$y-mean(CX_Nodes$y)
  
  # Edges
  EDGES=Direct_Data[c("type.from","type.to","weightRelative")]
  colnames(EDGES)=c("from","to","weight")
  EDGES$from[startsWith(EDGES$from,"MBON")]=substr(EDGES$from[startsWith(EDGES$from,"MBON")], 5, 6)
  
  # Nodes
  NODES=rbind(MBON_Nodes, CX_Nodes)
  NODES=merge(NODES, AllSize, by="name")
  
  # Node Colors
  TransmitterColor=c("#F26629","#20B685","#3668AC") #ACH, GLU
  NODES$COLORZ="Black"
  NODES$COLORZ[NODES$name %in% c("12","13","15","21","22","23","24","26","27","29","33","35")]=TransmitterColor[1]
  NODES$COLORZ[NODES$name %in% c("01 ","03","04","05","06","30","34")]=TransmitterColor[2]
  NODES$COLORZ[NODES$name %in% c("09","11")]=TransmitterColor[3]
  
  # Node fills
  NODES$FILLZ=NODES$COLORZ
  NODES$FILLZ[NODES$name %in% c("21","22","23","24","26","27","29","33","35")]="white"
  
  #Edge colors
  EDGES$COLORZ="black"
  EDGES$COLORZ[EDGES$from %in% c("12","13","15","21","22","23","24","26","27","29","33","35")]=TransmitterColor[1]
  EDGES$COLORZ[EDGES$from %in% c("01 ","03","04","05","06","30","34")]=TransmitterColor[2]
  EDGES$COLORZ[EDGES$from %in% c("09","11")]=TransmitterColor[3]
  EDGES$COLORZ=factor(EDGES$COLORZ, levels=unique(EDGES$COLORZ))
  col_vector=levels(EDGES$COLORZ)
  
  # Make graph
  graph = tbl_graph(NODES, EDGES)
  
  # Plot graph
  P1=ggraph(graph,layout="manual",x=NODES$x,y=NODES$y) + 
    geom_edge_fan( aes(width=weight, color=COLORZ),
                   arrow = arrow(length = unit(2, 'mm'),ends = "last", type = "closed"),
                   end_cap = circle(0.66, 'cm'), alpha=1,strength=1, n=2, linemitre=1) +
    geom_node_point(aes(size = Size), color=NODES$COLORZ, fill=NODES$FILLZ) + scale_size(range = c(7,12)) +
    geom_node_text(aes(label=name),angle=0,size=3, nudge_x = NODES$x*0.08, color='black' ) +
    theme_classic() + 
    theme(legend.text=element_text(size=6),legend.title=element_text(size=6),
          axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),axis.title.y=element_blank()) +
    scale_edge_width(range = c(1,6)) + scale_edge_color_manual(values=col_vector)
  ggsave(paste(FigDir, "Direct_Graph.png",sep=""),
         plot = P1, device='png', scale = 1, width =8, height = 13, units ="in", dpi = 500, limitsize = TRUE)

}



Get_VisHygThermo_Graph <- function(Graph_Pathway, Int){
  #' Function containing low level code to build a network graph showing indirect MBON to CX connections
  
  # Convert pathway dataframe to connecitiy table (from-to pairs)
  Graph_Data_L1=Graph_Pathway[c("type.from","type_N1","weightRelative_N1","supertype1_N1")]
  colnames(Graph_Data_L1)=c("from","to","weight","supertype1_N1")
  Graph_Data_L1$Layer="L1"
  Graph_Data_L2=Graph_Pathway[c("type_N1","type.to","weightRelative_N2","supertype1_N1")]
  colnames(Graph_Data_L2)=c("from","to","weight","supertype1_N1")
  Graph_Data_L2$Layer="L2"
  Graph_Data=rbind(Graph_Data_L1,Graph_Data_L2)
  remove(Graph_Data_L1, Graph_Data_L2)
  
  # Get a list of intermediate and CX nodes
  IntermediateNodes=unique(Graph_Pathway$type_N1)
  CXNodes=unique(Graph_Pathway$type.to)
  
  # Make edges
  EDGES=Graph_Data
  
  # Make nodes
  NODES=data.frame(name=unique(c(EDGES$from,EDGES$to)))
  NODES$name=as.character(NODES$name)
  
  NODES$x=NA
  MBONlogica=startsWith(NODES$name,"MBON")
  NODES$x[MBONlogica]=-1
  NODES$y[MBONlogica]=seq(from=-1, to=1, by=2/(sum(MBONlogica)-1))*0.35
  
  Intermediatelogica=NODES$name %in% IntermediateNodes
  NODES$x[Intermediatelogica]=0
  NODES$y[Intermediatelogica]=seq(from=-1, to=1, by=2/(sum(Intermediatelogica)-1))*1
  
  CXlogica=NODES$name %in% CXNodes
  NODES$x[CXlogica]=1
  NODES$y[CXlogica]=seq(from=-1, to=1, by=2/(sum(CXlogica)-1))*1.3
  
  # Set node colors
  NODES$COLORZ="Black"
  NODES$FILLZ="Black"
  
  # Set edge color
  EDGES$supertype1_N1=factor(EDGES$supertype1_N1,levels=Int)
  col_vector=(color(c("navy","dodgerblue2","limegreen", "magenta3","hotpink1","darkviolet","slateblue1")))
  
  # Make graph
  graph = tbl_graph(NODES, EDGES)
  
  Gx=ggraph(graph,layout="manual",x=NODES$x,y=NODES$y) + 
    geom_edge_fan( aes(width=weight, color=EDGES$supertype1_N1),
                   arrow = arrow(length = unit(2, 'mm'),ends = "last", type = "closed"),
                   end_cap = circle(0.3, 'cm'), alpha=1,strength=1, n=2, linemitre=1) +
    geom_node_point(size=4, color=NODES$COLORZ, fill=NODES$FILLZ) +
    geom_node_text(aes(label=name),angle=0,size=3, nudge_x = NODES$x*0.18, color='red' ) +
    theme_classic() + 
    theme(legend.text=element_text(size=6),legend.title=element_text(size=6),
          axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),axis.title.y=element_blank()) +
    scale_edge_width(range = c(0.5,4)) + scale_edge_color_manual(values=col_vector, drop=FALSE)
  return(Gx)
}


