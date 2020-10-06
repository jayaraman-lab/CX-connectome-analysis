allConnections <- function(connmat,eps=10^-3,maxIt=100){
  simple <- vector(mode="list")
  loops <- vector(mode="list")
  
  loops[[1]] <- diag(connmat)
  diag(connmat) <- 0
  simple[[1]] <- connmat
  
  i <- 1
  while(i<maxIt & norm(simple[[i]])>eps){
    simple[[i+1]] <- simple[[i]] %*% connmat
    loops[[i+1]] <- diag(simple[[i+1]])
    diag(simple[[i+1]]) <- 0
    i <- i+1  
  }
  return(list(simple=simple,loops=loops))
}

plotRecurrentGraph <- function(allCXRecurrence,focusType,nSt,focusOn=c("type.from","type.to"),centrality=TRUE,textAngle=0,localScale){
  focusOn <- match.arg(focusOn)
  
  pathType <- filter(allCXRecurrence,!!(sym(focusOn))==focusType & n_steps %in% nSt)
  
  pathSummaries <- group_by(pathType,type.from,type.to) %>% summarize(pathWeight=sum(weightRelative_path))
  
  typeFeedbackGraph <- pathDf2graphDf(pathType) %>% rename(type.from=from,type.to=to) %>% supertype()
  
  typeFeedbackGraph <- makeGraph(typeFeedbackGraph)
  typeFeedbackGraph$nodes <- mutate(typeFeedbackGraph$nodes,
                                    positionL=case_when(name %in% pathType$type.from ~ 1,
                                                        name %in% pathType$type_N1 ~ 2,
                                                        name %in% pathType$type_N2 ~ 3,
                                                        name %in% pathType$type_N3 ~ 4,
                                                        name %in% pathType$type.to ~ max(nSt)+1
                                    ))
  
  typeFeedbackGraph$graph <- typeFeedbackGraph$graph %>% mutate(position=case_when(name %in% pathType$type.from ~ max(nSt)+1,
                                                                                   name %in% pathType$type_N1 ~ max(nSt)+0,
                                                                                   name %in% pathType$type_N2 ~ max(nSt)-1,
                                                                                   name %in% pathType$type_N3 ~ max(nSt)-2,
                                                                                   name %in% pathType$type.to ~ 1
  ),
  
  importance=case_when(
    type %in% pathSummaries$type.to ~ pathSummaries$pathWeight[match(type,pathSummaries$type.to)],
    type %in% pathSummaries$type.to ~ 1,
    TRUE ~ 0.05))
  
  if(nrow(typeFeedbackGraph$nodes)>1){
    
    if(centrality){
      ggG <- ggraph(typeFeedbackGraph$graph,layout="centrality",cent=position) + geom_edge_fan0(aes(width=weightRelative,colour = supertype2.from)) + coord_fixed() + theme_paper_map()+ scale_edge_width(range=c(0.1,4),limits=c(0,1))+
        graphlayouts::draw_circle(col = "gray80", use = "cent")+ scale_color_manual(breaks=names(localScale),values=localScale,guide="none") + scale_edge_color_manual(breaks=names(localScale),values=localScale,guide="none")+
        geom_node_point(aes(color=supertype2,shape = type %in% allCXRecurrence$type.from,size=importance))+ geom_node_text(aes(label=name),size=2) + scale_size(range=c(0.1,6),limits=c(0,1))
    }else(
      ggG<- ggraph(typeFeedbackGraph$graph,layout="sugiyama",attributes="all",layers=typeFeedbackGraph$nodes$positionL) + geom_edge_fan0(aes(width=weightRelative,colour = supertype2.from)) + theme_paper_map()+ scale_edge_width(range=c(0.1,4),limits=c(0,1))+
        scale_color_manual(breaks=names(localScale),values=localScale,guide="none") + scale_edge_color_manual(breaks=names(localScale),values=localScale,guide="none")+
        geom_node_point(aes(color=supertype2,shape = type %in% allCXRecurrence$type.from,size=importance))+ geom_node_text(aes(label=name),size=2,angle=textAngle) + scale_size(range=c(0.1,6),limits=c(0,1)) +
        scale_y_continuous(expand = expansion(0.2)) + scale_x_continuous(expand = expansion(0.2)) + scale_shape_manual(values=c(15,16))
    )
    return(ggG)}
}  


plotOutAndCXGraph <- function(ggraphObject,colorScale=NULL){
  ebGG <- makeGraph(ggraphObject) 
  
  if (is.null(colorScale)){
    localScale <- paletteer_d("Polychrome::palette36")[c(3:4,6:(length(unique(ebGG$nodes$supertype2))+3))]
    names(localScale) <- unique(ebGG$nodes$supertype2)}else{localScale=colorScale}
  
  #  toyIdent <- as_tbl_graph(data.frame(from=c("CX output","CX output"),to=c("Target","Target"),roi=c("CX","Outside")) %>% mutate(name.from=from,name.to=to)) %>% mutate(x=c(0.4,0.8),y=c(0.3,0.7))
  #  toyRecip <- as_tbl_graph(data.frame(from=c("CX output","Target"),to=c("Target","CX output"),roi=c("Outside","CX")) %>% mutate(name.from=from,name.to=to)) %>% mutate(x=c(0.4,0.8),y=c(0.3,0.7))
  #  toyCommon <- as_tbl_graph(data.frame(from=c("CX output","CX output","Target1"),to=c("Target1","Target2","Target2"),roi=c("Outside","Outside","CX")) %>% mutate(name.from=from,name.to=to)) %>% mutate(x=c(0.6,0.4,0.8),y=c(0.3,0.7,0.7))
  
  #  toyRecipG <- ggraph(toyRecip,x=x,y=y) + geom_edge_fan(aes(colour=name.from,label=roi),label_size=2) +  geom_node_point(size=8,aes(colour=name)) +geom_node_text(aes(label=name),size=2) + 
  #    theme_paper_map(legend.position="none") + scale_color_paletteer_d("ggthemes::Color_Blind") + scale_edge_color_manual(values=paletteer_d("ggthemes::Color_Blind")[1:2]) + ylim(c(0.2,0.8)) + xlim(c(0,1)) + coord_fixed()
  
  ##  toyIdentG <- ggraph(toyIdent,x=x,y=y) + geom_edge_fan(aes(colour=name.from,label=roi),label_size=2) +  geom_node_point(size=8,aes(colour=name)) +geom_node_text(aes(label=name),size=2) +
  #    theme_paper_map(legend.position="none") + scale_color_paletteer_d("ggthemes::Color_Blind") + scale_edge_color_manual(values=paletteer_d("ggthemes::Color_Blind")[1:2]) + ylim(c(0.2,0.8)) + xlim(c(0,1)) + coord_fixed()
  
  #  toyCommonG <- ggraph(toyCommon,x=x,y=y) + geom_edge_fan(aes(colour=name.from,label=roi),label_size=2) +  geom_node_point(size=8,aes(colour=name)) +
  #    geom_node_text(aes(label=name),size=2) + theme_paper_map(legend.position="none") + scale_color_paletteer_d("ggthemes::Color_Blind") + 
  #    scale_edge_color_manual(values=paletteer_d("ggthemes::Color_Blind")[1:2]) + ylim(c(0,1)) + xlim(c(0,1)) 
  
  ggraph(ebGG$graph,layout="stress") +
    geom_edge_fan0(aes(colour = supertype2.from,width=weightRelative)) +
    geom_edge_loop(aes(colour = supertype2.from)) + scale_shape(name="") +
    geom_node_point(aes(color=supertype2,shape = category),size=5)+ geom_node_text(aes(label=name),size=2) + scale_color_manual(breaks=names(localScale),values=localScale) + scale_edge_color_manual(breaks=names(localScale),values=localScale)+
    facet_edges(motif~linkCat,ncol=2,drop=FALSE) + theme_paper_map(panel.background=element_rect(fill="gray95",linetype = "blank"))  + guides(edge_colour="none",colour = "none")+ scale_edge_width(range=c(0.1,4),limits=c(0,1))
  
  
}

collectInputsTable <- function(connTable,connTableRaw){
  inputNeurons <- getTypesTable(unique(connTable$databaseType.from))
  inputNeurons <- cxRetyping(inputNeurons)
  unknownInput <-cxRetyping(neuprint_get_meta(unique(connTableRaw$from[!(connTableRaw$from %in% inputNeurons$bodyid)]))  %>% mutate(
    name = ifelse(is.na(name),as.character(bodyid),name),
    type = ifelse(is.na(type),gsub("_L$|_R$","",name),type)
  ) %>% mutate(databaseType=NA))
  
  inputNeurons <- rbind(inputNeurons,unknownInput) %>% filter(type %in% connTable$type.from)
  
}