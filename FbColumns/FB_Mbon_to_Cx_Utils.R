##### A collection of functions for analyzing MBON to FB connectivity  ######################################


Get_Direct_Graph <- function(Graph_Data_Clust){
  #' Function containing low level code to build a network graph showing direct MBON to CX connections
  
  # Get node sizes
  PostSize=Graph_Data_Clust %>% group_by(type.to) %>% summarize(Syns=sum(weight))
  colnames(PostSize)=c("name","Size")
  PreSize=Graph_Data_Clust %>% group_by(type.from) %>% summarize(Syns=sum(weight))
  colnames(PreSize)=c("name","Size")
  AllSize=rbind(PostSize, PreSize)
  AllSize$name[startsWith(AllSize$name,"MBON")]=substr(AllSize$name[startsWith(AllSize$name,"MBON")], 5, 6)
  
  # Pre Nodes
  MBON_Nodes=data.frame(name=(unique(Graph_Data_Clust$type.from) %>% as.character()))
  MBON_Nodes$name=substr(MBON_Nodes$name,5, 6)
  MBON_Nodes$x=0
  MBON_Nodes$y=(seq(from = 1, to = length(MBON_Nodes$name), by = 1))
  MBON_Nodes$y=(MBON_Nodes$y-mean(MBON_Nodes$y))*1.2
  
  # Post Nodes
  CX_Nodes=data.frame(name=rev((unique(Graph_Data_Clust$type.to) %>% as.character()))) 
  CX_Nodes$x=1
  CX_Nodes$y=rev(seq(from = 1, to = length(CX_Nodes$name), by = 1))
  CX_Nodes$y=CX_Nodes$y-mean(CX_Nodes$y)
  
  # Edges
  EDGES=Graph_Data_Clust[c("type.from","type.to","weightRelative")]
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
  return(graph) 
}

