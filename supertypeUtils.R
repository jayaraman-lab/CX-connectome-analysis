supertype <- function(types,level=2){
  #' Supertypes can be defined at 3 levels of depth (default 2)
  supertype <- types
  if (level == 1){## Minimum change
    supertype[grepl("FB.*",types)] <- str_extract(types,"FB[1-9]")[grepl("FB.*",types)]
    supertype[grepl("Delta0.*",types)] <- str_extract(types,"Delta0[A-N]")[grepl("Delta0.*",types)]
    supertype[grepl("Delta6.*",types)] <- str_extract(types,"Delta6[A-K]")[grepl("Delta6.*",types)]
    supertype[grepl("FC.*",types)] <- str_extract(types,"FC[1-9]")[grepl("FC.*",types)]
    supertype[grepl("FS.*",types)] <- str_extract(types,"FS[1-9]")[grepl("FS.*",types)]
    supertype[grepl("^FR.*",types)] <- "FR"
    supertype[grepl("PFR.*",types)] <- "PFR"
    supertype[grepl("PFN.*",types)] <- str_extract(types,"PFN[a|d|m|p|v]")[grepl("PFN.*",types)]
  }
  if (level == 2){
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
    supertype[grepl("^R[1-6].*",types)] <- "Ring"
    supertype[grepl("SA.*",types)] <- "SA"
    supertype[grepl("SpsP.*",types)] <- "SPS-PB"
    supertype[grepl("AVM09v_.*|AVM09p_.*|AVM09x_.*",types)] <- "L-L(c)"
    supertype[grepl("AVL01op_.*|AVL01ln_.*",types)] <- "LW"
    supertype[grepl("PDL11f_.*|PVM03m_.*",types)] <- "L-Sm"
    supertype[grepl("1882031306|5813104472|PVM09r_.*",types)] <- "LW-(X)(c)"
    supertype[grepl("ADM05a_.*",types)] <- "CRE-LAL"
    supertype[grepl("OA_V.*",types)] <- "OA"
  }
  
  if (level ==3){
    supertype[grepl("^PF.*",types)] <- "FB Columnar"
    supertype[grepl("EPG.*|PEG.*|PEN.*",types)] <- "EB Columnar"
    supertype[grepl("^FC.*|^FR.*|^FS.*",types)] <- "FB Output"
    supertype[grepl("FB[1-9].*",types)] <- "FB Tangential"
    supertype[grepl("Delta[0|6].*",types)] <- "FB Interneuron"
  }
  supertype
}


#megatype <- function(supertypes){
#  megatype <- supertypes
#  megatype[supertypes %in% c("PFL","PFN","PFR","PFGs")] <- "FB Columnar"
#  megatype[supertypes %in% c("EPG","PEG","PEN")] <- "EB Columnar"
#  megatype[grepl("FB[1-9]",supertypes)] <- "FB tangential"
#  megatype[grepl("Delta0.*",supertypes)] <- "D0"
#  megatype[grepl("Delta6.*",supertypes)] <- "D6"
#  megatype
#}