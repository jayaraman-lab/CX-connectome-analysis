#############################################################

nodesFromTypeTable <- function(type2typeTable){
  nodes <- data.frame(name = unique(c(type2typeTable$type.from,type2typeTable$type.to)),stringsAsFactors = F) 
}

edgesFromTypeTable <- function(type2typeTable,pathNodes = nodesFromTypeTable(type2typeTable)){
  type2typeTable %>% mutate(to = sapply(type.to, function(f) which(f == pathNodes$name)),
                            from = sapply(type.from, function(f) which(f == pathNodes$name)))
}


makePyramidGraph <- function(type2typeList,ROIs,polarity="inputs",plot=FALSE){
  type2typeTable <- type2typeList[[polarity]]
  sourceTable <- type2typeList[["names"]]
  type2typeTable <- type2typeTable %>% filter(roi %in% ROIs)
  nodes <- nodesFromTypeTable(type2typeTable)
  if (polarity == "inputs"){
    nodes <- nodes %>% mutate(layer=ifelse(name %in% sourceTable$type,2,1))
  }else{
      nodes <- nodes %>% mutate(layer=ifelse(name %in% sourceTable$type,1,2))
  }
  edges <- edgesFromTypeTable(type2typeTable)
  graph <- tbl_graph(nodes,edges)
  if (!plot){
      return(list(graph = graph,nodes = nodes,edges= edges))
  }else{
      return(pyramidGraph(graph,nodes,edges))
  }
}

pyramidGraph <- function(graphT,nodeT,edgeT){
    ggraph(graphT,layout="sugiyama",layers=nodeT$layer) + 
      geom_edge_fan(aes(width=weightRelative),colour="grey",alpha=0.5) + 
      geom_edge_loop(colour="grey",aes(direction=10,span=10),alpha=0.5) +
      geom_node_point(aes(color=name),size=5) + 
      geom_node_text(aes(label=name),angle=40,size=4) +
      guides(color="none") + 
      facet_edges(~roi)
}