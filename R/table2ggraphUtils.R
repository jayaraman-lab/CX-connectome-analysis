library(tidygraph)
# Utilities to transform a table of connections into a tidygraph object with metadata
#############################################################

nodesFromTypeTable <- function(type2typeTable){
  
  typesFrom <- type2typeTable %>% select(contains(".from")) %>% rename_with(~gsub(".from","",.x)) 
  typesTo <- type2typeTable %>% select(contains(".to")) %>% rename_with(~gsub(".to","",.x)) 
  
  typesTable <- distinct(rbind(typesFrom,typesTo))
  
  nodes <- typesTable %>% mutate(name = type) 
    
  nodes
}

edgesFromTypeTable <- function(type2typeTable,pathNodes = nodesFromTypeTable(type2typeTable)){
  distinct(type2typeTable %>% mutate(to = match(type.to,pathNodes$name),
                             from = match(type.from,pathNodes$name)))
}

makeGraph <- function(type2type,polarity="inputs"){UseMethod("makeGraph")}

makeGraph.data.frame <- function(type2type){
  
  nodes <- nodesFromTypeTable(type2type)
  edges <- edgesFromTypeTable(type2type,nodes)
  graph <- tbl_graph(nodes,edges)
  graph
}

makeGraph.neuronBag <- function(type2type,polarity="inputs"){
  type2typeTable <- type2type[[polarity]]
  makeGraph(type2typeTable,roiSelect=roiSelect,roiLabel=roiLabel,by.roi=by.roi,polarity=polarity)
}
