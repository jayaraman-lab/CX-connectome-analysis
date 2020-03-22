#############################################################

nodesFromTypeTable <- function(type2typeTable){
  
  typesFrom <- type2typeTable %>% rename(type = type.from,databaseType = databaseType.from) %>% select(type,databaseType)
                                  
  
  typesTo <- type2typeTable %>% rename(type = type.to,databaseType = databaseType.to) %>% select(type,databaseType)
    

  typesTable <- distinct(bind_rows(typesFrom,typesTo))
  
  nodes <- typesTable %>% mutate(name = type) %>%
    supertype() %>% mutate(community = as.factor(match(supertype2,unique(supertype2))))
    
}

edgesFromTypeTable <- function(type2typeTable,pathNodes = nodesFromTypeTable(type2typeTable)){
  distinct(type2typeTable %>% mutate(to = sapply(type.to, function(f) which(f == pathNodes$name)),
                            from = sapply(type.from, function(f) which(f == pathNodes$name))))
}

makeGraph <- function(type2type,roiSelect=NULL,roiLabel=roiSelect,by.roi=FALSE,polarity="inputs"){UseMethod("makeGraph")}

makeGraph.data.frame <- function(type2type,roiSelect=NULL,roiLabel=roiSelect,by.roi=FALSE,polarity="inputs"){
  if (is.null(roiSelect) & by.roi==FALSE){
    type2type <- type2type %>% group_by(type.from,type.to) %>%
      summarize_if(is.numeric,sum) %>% ungroup()
  }else{
    if (!is.null(roiSelect)){
      type2type <- type2type %>% filter(roi %in% roiSelect$roi) %>% 
                    mutate(roi = factor(roi,levels=levels(roiSelect$roi)),
                    l4 = roiSelect$level4[match(roi,roiSelect$roi)],
                    side = roiSelect$side2[match(roi,roiSelect$roi)],
                    region = roiLabel$roi[match(l4,roiLabel$level4)])
    }
  }
  nodes <- nodesFromTypeTable(type2type)
  if (polarity == "inputs"){
    nodes <- nodes %>% mutate(layers=ifelse(name %in% unique(type2type$type.to),2,1))
  }else{
    nodes <- nodes %>% mutate(layers=ifelse(name %in% unique(type2type$type.to),1,2))
  }
  edges <- edgesFromTypeTable(type2type,nodes)
  graph <- tbl_graph(nodes,edges)
 
  return(list(graph = graph,nodes = nodes,edges= edges))
}

makeGraph.neuronBag <- function(type2type,roiSelect=NULL,roiLabel=roiSelect,by.roi=FALSE,polarity="inputs"){
  type2typeTable <- type2type[[polarity]]
  makeGraph(type2typeTable,roiSelect=roiSelect,roiLabel=roiLabel,by.roi=by.roi,polarity=polarity)
}

pyramidGraph <- function(graphT,nodeT,by.roi=T){
    g <- ggraph(graphT,layout="sugiyama")+#,layers=nodeT$layer) + 
      geom_edge_fan(aes(width=weightRelative),colour="grey",alpha=0.5) + 
      geom_edge_loop(colour="grey",aes(direction=10,span=10,width=weightRelative),alpha=0.5) +
      geom_node_point(aes(color=supertype2),size=5) + 
      geom_node_text(aes(label=name),angle=40,size=4)
    if (by.roi){
      g <- g+facet_edges(~roi)
    }
    g
}