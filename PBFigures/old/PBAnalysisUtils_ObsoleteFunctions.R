# Build an xyLookup table for plotting a PB connectivity graph
xyLookupTable <- function(){
  nronGps <- list(c("EPG","Delta7","IbSpsP","SpsP"),
                  c("PEN_a(PEN1)","PEN_b(PEN2)","PEG","PFGs","PFR_a","PFR_b","PFL1","PFL2","PFL3",
                    "PFNa","PFNd","PFNm_a","PFNm_b","PFNp_a","PFNp_b","PFNp_c","PFNp_d","PFNp_e","PFNv"))
  
  xyLookup <- data.frame(type = c(), x = c(), y = c(), xTxt = c(), yTxt = c())
  for (g in 1:length(nronGps)){
    types <-  nronGps[[g]]
    
    xs <- numeric(length(types))
    ys <- numeric(length(types))
    if (g > 1){
      angs <- seq(-pi/3,pi/3,length.out=length(types))
      xsTxt <- xs + 1.1*sin(angs)
      xs <- xs + sin(angs)
      ysTxt <- ys - 1.1*cos(angs)
      ys <- ys - cos(angs)
    } else {
      xs <- xs + c(-0.25,0.25,-0.25,0.25)
      xsTxt <- xs
      ysTxt <- ys + 0.1
    }
    
    typeLookup <- data.frame(type = types, x = xs, y = ys, xTxt = xsTxt, yTxt = ysTxt)
    xyLookup <- rbind(xyLookup,typeLookup)
  }
  return(xyLookup)
  #ggplot(xyLookup,aes(x=x,y=y)) + geom_point() + scale_y_reverse() + geom_text(aes(label=type),hjust=1, vjust=1,size=3)
}

# Create a function for plotting a graph of PB connections
graphConTab <- function(conTab,xyLookup,textRepel,guideOnOff){
  
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
    #  sTScale <- scale_colour_manual(values = sTCols$color, drop=TRUE,limits = sTCols$supertype)
  } else {
    sTScale <- scale_colour_manual(values = pal$pal, drop=TRUE,limits = pal$breaks,guide = FALSE)
    #  sTScale <- scale_colour_manual(values = sTCols$color, drop=TRUE,limits = sTCols$supertype,guide = FALSE)
  }
  
  #sTScale_edge <- scale_edge_colour_manual(values = sTCols$color, drop=TRUE,limits = sTCols$supertype,guide = FALSE)
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
    geom_node_text(aes(x = xyLookup$xTxt, y = xyLookup$yTxt, label=name),size=3,repel = textRepel) +
    theme_classic() + theme(legend.text=element_text(size=6),legend.title=element_text(size=8),
                            axis.line=element_blank(),axis.text.x=element_blank(),
                            axis.text.y=element_blank(),axis.ticks=element_blank(),
                            axis.title.x=element_blank(),axis.title.y=element_blank(),) + 
    coord_fixed(ratio = 1,clip="off")
  
  return(gg)
}