supertype <- function(types){
  supertype <- types
  supertype[grepl("FB.*",types)] <- "FB tangential"
  supertype[grepl("PFL.*",types)] <- "PFL"
  supertype[grepl("PFN.*",types)] <- "PFN"
  supertype[grepl("PFR.*",types)] <- "PFR"
  supertype[grepl("PEN.*",types)] <- "PEN"
  supertype[grepl("LN.*",types)] <- "LN"
  supertype[grepl("Delta0.*",types)] <- "Delta0"
  supertype[grepl("Delta6.*",types)] <- "Delta6"
  supertype[grepl("ExR.*",types)] <- "ExR"
  supertype[grepl("FC.*",types)] <- "FC"
  supertype[grepl("^FR.*",types)] <- "FR"
  supertype[grepl("FS.*",types)] <- "FS"
  supertype[grepl("LCN.*",types)] <- "LN"
  supertype[grepl("^R[1-6].*",types)] <- "Ring neuron"
  supertype[grepl("SA.*",types)] <- "SA"
  supertype[grepl("SpsP.*",types)] <- "SPS-PB"
  supertype[grepl("AVM09v_.*|AVM09p_.*|AVM09x_.*",types)] <- "LAL-LAL contra"
  supertype[grepl("AVL01op_.*|AVL01ln_.*",types)] <- "LAL-WED"
  supertype[grepl("PDL11f_.*|PVM03m_.*",types)] <- "LAL-SMP"
  supertype[grepl("1882031306|5813104472|PVM09r_.*",types)] <- "LAL-WED-(X)contra"
  supertype[grepl("ADM05a_.*",types)] <- "CRE-LAL"
  supertype[grepl("OA_V.*",types)] <- "OA"
  supertype
}

megatype <- function(supertypes){
  megatype <- supertypes
  megatype[supertypes %in% c("PFL","PFN","PFR","PFGs")] <- "FB Columnar"
  megatype[supertypes %in% c("EPG","PEG","PEN")] <- "EB Columnar"
  megatype[supertypes %in% c("Delta6","Delta0")] <- "FB interneurons"
  megatype[supertypes %in% c("FC","FR","FS")] <- "FB outputs"
  megatype
}