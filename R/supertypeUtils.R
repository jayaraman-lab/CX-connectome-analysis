supertype <- function(types,level=2){UseMethod("supertype")}

supertype.character <- function(types,level=2){
  #' Gets a supertype from a type name, with various levels of depth
  #' @param types : a type name or vector of type names
  #' @param level : depth of the supertype. Possible values are 1,2 or 3 (default 2), 1 being the finest 
  #' subdivision and 3 the coarsest
  #' @details For example, at level 1 delta0 neurons are just divided in DeltaA to K, at level 2 they are
  #' D0, and at level 3 they are FB Interneurons
  supertype <- types
  supertype[is.na(types)] <- "Other"
  
  supertype[grepl("FB.*",types)] <- str_extract(types,"FB[1-9]")[grepl("FB.*",types)]
  supertype[grepl("Delta0.*",types)] <- str_extract(types,"Delta0[A-N]")[grepl("Delta0.*",types)]
  supertype[grepl("Delta6.*",types)] <- str_extract(types,"Delta6[A-K]")[grepl("Delta6.*",types)]
  supertype[grepl("FC.*",types)] <- str_extract(types,"FC[1-9]")[grepl("FC.*",types)]
  supertype[grepl("FS.*",types)] <- str_extract(types,"FS[1-9]")[grepl("FS.*",types)]
  supertype[grepl("^FR.*",types)] <- "FR"
  supertype[grepl("^EL.*",types)] <- "EL"
  supertype[grepl("PFR.*",types)] <- "PFR"
  supertype[grepl("PFN.*",types)] <- str_extract(types,"PFN[a|d|m|p|v]")[grepl("PFN.*",types)]
  supertype[types %in% c("L-L(c)1","AoL-L(c)","L-LC(c)1","VeL-CLVe(c)1","L-Lic(c)","VeL-LVeC(c)",
                        "CL-LC(c)","VeL-LVe(c)2","CL-L(c)1","L-LC(c)2","VeLC-L(c)","CL-L(c)2",
                        "VeL-L(c)","CL-L(c)2","VeLC-LVe(c)3","VeL-LC(c)","L-L(c)2")] <- "L-L(c)"
  supertype[grepl("^C.*Sm$",types)] <- "C-Sm"
  supertype[types %in% c("CL-AtIb(c)1","CL-AtIb(c)2","CL-IbAtl(b)1","CL-SmIb","CL-AtIb(b)2","CL-AtIb(b)3")] <- "CL-Ib"
  supertype[types %in% c("LCPl-SmSi","WplL-Sm")] <- "L-Sm"  
  supertype[types %in% c("IbSpVeL-(LIb)(c)Sm(b)","VeLC-(LC)(c)Sm(b)")] <- "L-Sm(b)"
  supertype[types %in% c("L-VeSpIp","L-EpVeSp","L-VeSp","LVe-SpIp","LW-Sp")] <-  "L-Sp"
  supertype[startsWith(types,"L-Ve")] <- "L-Ve"
  supertype[startsWith(types,"L-WP")] <- "L-WP"
  supertype[startsWith(types,"LVe")] <- "LVe"
  supertype[startsWith(types,"IbSp")] <- "IbSp"
  supertype[startsWith(types,"Pl-C")] <- "Pl-C"
  supertype[types %in% c("SpVeL-(?)(c)","SpVeL-LSpIp(c)1","SpVeL-LSpIp(c)2")] <- "SpVeL-X(c)"
  supertype[types %in% c("VeLC-CLVe(c)","VeLC-LVe(c)")] <- "VeL-LVe"
  supertype[types %in% c("WL-(X)(c)","VeWL-VeX(c)","LW-X(c)")] <- "WL-X(c)"
  
  if (level == 1){return(supertype)}
  
  supertype[grepl("FB.*",types)] <- "FBt"
  supertype[grepl("Delta0.*",types)] <- "D0"
  supertype[grepl("Delta6.*",types)] <- "D6"
  supertype[grepl("PFL.*",types)] <- "PFL"
  supertype[grepl("PFN.*",types)] <- "PFN"
  supertype[grepl("PFR.*",types)] <- "PFR"
  supertype[grepl("PEN.*",types)] <- "PEN"
  supertype[grepl("LN.*",types)] <- "LN"
  supertype[grepl("ExR.*",types)] <- "ExR"
  supertype[grepl("FC.*",types)] <- "FC"
  supertype[grepl("^FR.*",types)] <- "FR"
  supertype[grepl("FS.*",types)] <- "FS"
  supertype[grepl("LCN.*",types)] <- "LN"
  supertype[grepl("^EL.*",types)] <- "EL"
  supertype[grepl("^R[1-6].*",types)] <- "Ring"
  supertype[grepl("SA.*",types)] <- "SA"
  supertype[grepl("SpsP.*",types)] <- "SPS-PB"
  supertype[grepl("OA_V.*",types)] <- "OA"
  supertype[grepl("P[1|6].*",types)] <- "P"
  
  if (level == 2){return(supertype)}

  supertype[grepl("Delta7|P[1|6].*",types)] <- "PB Interneurons"
  supertype[grepl("^PF.*",types)] <- "FB Columnar"
  supertype[grepl("EPG.*|PEG.*|PEN.*|^EL.*",types)] <- "EB Columnar"
  supertype[grepl("^FC.*|^FR.*|^FS.*",types)] <- "FB Output"
  supertype[grepl("FB[1-9].*",types)] <- "FB Tangential"
  supertype[grepl("Delta[0|6].*",types)] <- "FB Interneuron"
  
  supertype[types == supertype] <- "Other"
  supertype
}

supertype.neuronBag <- function(types){
  for (lev in 1:3){
    for (ty in c(".from",".to")){
      for (tab in c("inputs","outputs","inputs_raw","outputs_raw")){
        types[[tab]][[paste0("supertype",ty,lev)]] <- supertype(types[[tab]][[paste0("databaseType",ty)]],level=lev) 
      }
    }
    types$names[[paste0("supertype",lev)]] <-  supertype(types$names[[paste0("databaseType")]],level=lev)
    types$outputsTableRef[[paste0("supertype",lev)]] <-  supertype(types$outputsTableRef[[paste0("databaseType")]],level=lev)
  }
  types
}

supertype.data.frame <- function(types,level=1:3){
  renamable <- names(types)[names(types) %in% c("databaseType","databaseType.from","databaseType.to")]
  for (lev in level){
    for (ty in renamable){
        types[[paste0(sub("databaseType","supertype",ty),lev)]] <- supertype(types[[ty]],level=lev) 
    }
  }
  types
}

supertype.NULL <- function(types,level=NULL){
  return(NULL)
}