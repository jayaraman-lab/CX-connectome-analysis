########################################################################################
# Functions to analyze tangential/contextual inputs and their effects on CX networks   #
########################################################################################

# From DanTE: Get edges and nodes and plot graph from a given connection matrix
graphConTabPolyChrome <- function(conTab,xyLookup,textRepel,guideOnOff){
  
  # Get the table of nodes (types)
  nodes = data.frame(name = xyLookup$type)
  nodes$superType <- nodes$name %>% as.character %>% supertype()
  
  # Position the nodes according to the lookup table
  nodes$x <- sapply(nodes$name, function(x) xyLookup$x[match(x,xyLookup$type)])
  nodes$y <- sapply(nodes$name, function(x) xyLookup$y[match(x,xyLookup$type)])
  
  # Assign colors to the supertypes
  pal <- supertype2Palette()
  if (guideOnOff){
    sTScale <- scale_colour_manual(values = pal$pal, drop=TRUE,limits = pal$breaks)
  } else {
    sTScale <- scale_colour_manual(values = pal$pal, drop=TRUE,limits = pal$breaks,guide = FALSE)
  }
  sTScale_edge <- scale_edge_colour_manual(values = pal$pal, drop=TRUE,limits = pal$breaks,guide = FALSE)
  
  # Get the edges from the connection table
  edges <- conTab[which((conTab$type.from %in% nodes$name) & (conTab$type.to %in% nodes$name)),] %>%
    mutate(to = sapply(type.to, function(f) which(f == nodes$name)),
           from = sapply(type.from, function(f) which(f == nodes$name)))
  edges$superType <- edges$type.from %>% as.character %>% supertype()
  
  # Plot the network
  graph <- tbl_graph(nodes,edges)
  gg <-
    ggraph(graph,layout="manual",x=nodes$x,y=nodes$y) + 
    geom_edge_diagonal(aes(width=weightRelative,color=superType),alpha=0.5,
                       strength=0.2,
                       #arrow = arrow(length = unit(0.5, "cm")),
                       end_cap = circle(0.5, 'cm')) + 
    #geom_edge_loop(aes(direction=45,span=90,width=weightRelative,color=superType,strength=0.1),alpha=0.5) +
    geom_node_point(aes(color=superType),size=4) + 
    sTScale_edge +
    sTScale +
    #geom_node_text(aes(x = xyLookup$xTxt, y = xyLookup$yTxt, label=name),size=3,repel = textRepel) +
    geom_node_text(aes(label=name),size=3,repel = textRepel) +
    theme_classic() + theme(legend.text=element_text(size=6),legend.title=element_text(size=8),
                            axis.line=element_blank(),axis.text.x=element_blank(),
                            axis.text.y=element_blank(),axis.ticks=element_blank(),
                            axis.title.x=element_blank(),axis.title.y=element_blank(),) + 
    coord_fixed(ratio = 1,clip="off")
  return(gg)
}

