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


Plot_SleepWake_Graph <- function(GraphTable){
  
# Set node positions
IDs_Type=sort(unique(c(unique(GraphTable$type.from),unique(GraphTable$type.to))))
NODES=data.frame(name=IDs_Type)

NODES$name=NODES$name[c(7,3,1,4,8,2,5,9,11,6,10)]
NODES$angle=NA
NODES$x=NA
NODES$y=NA

Diff=pi/(4)
NODES$angle[c(1,2,3,4,5)]=seq(from = 0, to = pi , by = Diff)
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

# NODES
NODES$y=-NODES$y
NODES$COLORZ="dodgerblue3"
NODES$COLORZ[NODES$name == "FB6H" | NODES$name == "FB7B"]="magenta3"

# Edges
EDGES=GraphTable[c("type.from","type.to","weightRelative")]
colnames(EDGES)=c("from","to","weight")
EDGES$from=as.character(EDGES$from)
EDGES$to=as.character(EDGES$to)
EDGES$COLORZ="Sleep to Sleep"
EDGES$COLORZ[EDGES$from == "FB6H" | EDGES$from == "FB7B"]="Wake to sleep"
EDGES$COLORZ[EDGES$to   == "FB6H" | EDGES$to   == "FB7B"]="Sleep to Wake"
EDGES$COLORZ=factor(EDGES$COLORZ,levels=(rev(c("Sleep to Wake", "Wake to sleep", "Sleep to Sleep"))))

# Make ggraph object
graph = tbl_graph(NODES, EDGES)

# colors
pcCols <- paletteer_d("Polychrome::palette36")
col_vector=rev(c("dodgerblue3","magenta3","gray50"))

# Make plot and save
Px=ggraph(graph,layout="manual",x=NODES$x,y=NODES$y) + 
  geom_edge_fan( aes(width=weight, color=COLORZ), 
                 arrow = arrow(length = unit(2, 'mm'),ends = "last", type = "closed"), end_cap = circle(1, 'cm'), alpha=1,strength=1) +
  geom_node_point(color=NODES$COLORZ, size=18) +
  geom_node_text(aes(label=name),angle=0,size=2, nudge_y = c(rep(0.06,18),rep(-0.06,9)), color='black' ) +
  theme_classic() + 
  theme(legend.text=element_text(size=6),legend.title=element_text(size=6),
        axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank()) +
  scale_edge_width(range = c(1.25,3)) + scale_edge_color_manual(values=col_vector)
ggsave(paste(PlotDir, "SleepWake_T2T_Graph.png",sep=""), plot = Px, device='png', scale = 1, width = 8.2, height = 8, units ="in", dpi = 500, limitsize = TRUE)
return(Px)
}
