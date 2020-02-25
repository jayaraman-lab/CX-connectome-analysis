#############################################################

nodesFromTypeTable <- function(type2typeTable){
  nodes <- data.frame(name = unique(c(type2typeTable$type.from,type2typeTable$type.to)),
                      stringsAsFactors = F) %>% mutate(supertype = supertype(name),
                                                       community = as.factor(match(supertype,unique(supertype))))
}

edgesFromTypeTable <- function(type2typeTable,pathNodes = nodesFromTypeTable(type2typeTable)){
  type2typeTable %>% mutate(to = sapply(type.to, function(f) which(f == pathNodes$name)),
                            from = sapply(type.from, function(f) which(f == pathNodes$name)))
}

graphFromInOut <- function(type2typeList,kind="intra"){
  nodesIn <- nodesFromTypeTable(type2typeList$inputs) 
  nodesOut <- nodesFromTypeTable(type2typeList$outputs)
  nodes <- distinct(bind_rows(nodesIn,nodesOut))
  
  inp <- type2typeList$inputs
  outp <- type2typeList$outputs
  if (kind=="intra"){
    nodes <- nodes %>% filter(name %in% type2typeList$names$type)
    inp <- type2typeList$inputs %>% filter(type.from %in% type2typeList$names$type)
    outp <- type2typeList$outputs %>% filter(type.to %in% type2typeList$names$type)
  }
  
  edges <- edgesFromTypeTable(distinct(bind_rows(inp,outp)),nodes)
  
  graph <- tbl_graph(nodes,edges)
}

makePyramidGraph <- function(type2typeList,ROIs=NULL,by.roi=FALSE,polarity="inputs",plot=FALSE){
  type2typeTable <- type2typeList[[polarity]]
  sourceTable <- type2typeList[["names"]]
  if (is.null(ROIs) & by.roi==FALSE){
    type2typeTable <- type2typeTable %>% group_by(type.from,type.to) %>%
                                     summarize_if(is.numeric,sum)
  }else{
    if (!is.null(ROIs)){
    type2typeTable <- type2typeTable %>% filter(roi %in% ROIs)
    }
  }
  nodes <- nodesFromTypeTable(type2typeTable,supertype=supertype)
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
      return(pyramidGraph(graph,nodes,edges,by.roi=by.roi))
  }
}

pyramidGraph <- function(graphT,nodeT,edgeT,by.roi=T){
    g <- ggraph(graphT,layout="sugiyama",layers=nodeT$layer) + 
      geom_edge_fan(aes(width=weightRelative),colour="grey",alpha=0.5) + 
      geom_edge_loop(colour="grey",aes(direction=10,span=10,width=weightRelative),alpha=0.5) +
      geom_node_point(aes(color=supertype),size=5) + 
      geom_node_text(aes(label=name),angle=40,size=4) +
      guides(color="none") 
    if (by.roi){
      g <- g+facet_edges(~roi)
    }
    g
}