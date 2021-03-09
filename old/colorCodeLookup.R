# color code
getSimpleTypeNames <- function(mydata){
  simpleTypes = gsub("_[R,L]", "", mydata)
  simpleTypes = gsub("PFR_a.*", "PFR", simpleTypes)
  simpleTypes = gsub("PFR_b.*", "PFR", simpleTypes)
  #simpleTypes = gsub("PFN.*", "PFN", simpleTypes)
  simpleTypes = gsub("PFNp_*", "PFNp", simpleTypes)
  simpleTypes = gsub("PFNm_.*", "PFNm", simpleTypes)
  simpleTypes = gsub("PFL.*", "PFL", simpleTypes)
  simpleTypes = gsub("FB[[:alnum:]]*", "FB", simpleTypes)
  simpleTypes = gsub("FB_[[:alnum:]]*", "FB", simpleTypes)
  simpleTypes = gsub("FC[[:alnum:]]*", "FC", simpleTypes)
  simpleTypes = gsub("EB-Q[[:alnum:]]*", "EB-Q", simpleTypes)
  #simpleTypes = gsub("ExR[[:alnum:]]*", "ExR", simpleTypes)
  simpleTypes = gsub("PB[[:alnum:]]*", "PB", simpleTypes)
  simpleTypes = gsub("SLP-AB[[:alnum:]]*", "SLP-AB", simpleTypes)
  simpleTypes = gsub("ER1[[:alnum:]]*_[[:alnum:]]*", "ER1", simpleTypes)
  simpleTypes = gsub("ER2[[:alnum:]]*_[[:alnum:]]*", "ER2", simpleTypes)
  simpleTypes = gsub("ER3a_[[:alnum:]]*", "ER3a", simpleTypes)
  simpleTypes = gsub("ER3d_[[:alnum:]]*", "ER3d", simpleTypes)
  simpleTypes = gsub("ER3m_[[:alnum:]]*", "ER3m", simpleTypes)
  simpleTypes = gsub("ER3p_[[:alnum:]]*", "ER3p", simpleTypes)
  simpleTypes = gsub("ER3w_[[:alnum:]]*", "ER3w", simpleTypes)
  simpleTypes = gsub("ER4d_[[:alnum:]]*", "ER4d", simpleTypes)
  simpleTypes = gsub("ER4m[[:alnum:]]*", "ER4m", simpleTypes)
  simpleTypes = gsub("ER5[[:alnum:]]*", "ER5", simpleTypes)
  simpleTypes = gsub("ER6[[:alnum:]]*", "ER6", simpleTypes)
  simpleTypes = gsub("TuBu09_[[:alnum:]]*", "TuBu09", simpleTypes)
  simpleTypes = gsub("AOTU[[:alnum:]]*", "AOTU", simpleTypes)
  simpleTypes = gsub("PEN_a\\(PEN1\\)", "PEN1", simpleTypes)
  simpleTypes = gsub("PEN_b\\(PEN2\\)", "PEN2", simpleTypes)
  simpleTypes = gsub("MC[[:alnum:]]*", "MC", simpleTypes)
  simpleTypes = gsub("TuTu.*", "TuTu", simpleTypes)
  simpleTypes = gsub("ADM06b_[[:alnum:]]*_[[:alnum:]]*", "ADM06b", simpleTypes)
  simpleTypes = gsub("MBON[[:alnum:]]*", "MBON", simpleTypes)
  simpleTypes = gsub("AVL[[:alnum:]]*_[[:alnum:]]*_pct", "AVL", simpleTypes)
  simpleTypes = gsub("AVL[[:alnum:]]*_pct", "AVL", simpleTypes)
  simpleTypes = gsub("AVM[[:alnum:]]*_pct", "AVM", simpleTypes)
  simpleTypes = gsub("ADL.*_pct", "ADL", simpleTypes)
  simpleTypes = gsub("PVL[[:alnum:]]*_pct", "PVL", simpleTypes)
  simpleTypes = gsub("PVM.*_pct", "PVM", simpleTypes)
  simpleTypes = gsub("LHPV6q1", "putWPN", simpleTypes)
  simpleTypes = gsub("PDL[12x,20s]{1}.*_pct", "PDLother", simpleTypes)
  simpleTypes = gsub("ADM11[[:alnum:]]*_pct", "ADM11", simpleTypes)
  simpleTypes = gsub("SAD.*", "SAD", simpleTypes)
  simpleTypes = gsub("WED.*", "WED", simpleTypes)
  simpleTypes = gsub("LAL138.*", "WLL", simpleTypes)
  simpleTypes = gsub("LAL.*", "LAL", simpleTypes)
  simpleTypes = gsub("olfactory multi .*", "other", simpleTypes)
  simpleTypes = gsub("LNO[[:alnum:]]{1}.*", "LNO", simpleTypes)
  simpleTypes = gsub("Delta7*", "Delta7", simpleTypes)
  simpleTypes = gsub("hDelta.*", "hDelta", simpleTypes)
  simpleTypes = gsub("vDelta.*", "vDelta", simpleTypes)
  simpleTypes = gsub("P6\\-8P9", "other", simpleTypes)
  simpleTypes = gsub("P1\\-9", "other", simpleTypes)
  simpleTypes = gsub("LC[[:alnum:]]{2}", "LC", simpleTypes)
  simpleTypes = gsub("LCm", "LC", simpleTypes)
  simpleTypes = gsub("OA-.*", "OA", simpleTypes)
  
  return(simpleTypes)
}

colorValueLookup = data.frame(
  type = c('ER1', 'ER2', 'ER3a', 'ER3d', 'ER3m', 'ER3p',  'ER3w', 'ER4d', 'ER4m', 'ER5', 'ER6',
           'ExR1','ExR2','ExR3','ExR4','ExR5','ExR6','ExR7','ExR8',
           'TuBu01', 'TuBu02', 'TuBu03', 'TuBu04', 'TuBu05', 'TuBu06', 'TuBu07', 'TuBu08', 'TuBu09', 'TuBu10',
           'AOTU', 'MC',  'TuTu', 'LC',
           'EPG', 'EPGt', 'PEN1', 'PEN2', 'PEG', 'EL',
           
           'PDM14j', 'WED','SAD', 'ADM06b','putWPN','ADM06s', 'PDM09','PDMother','FC','LAL',
           
           'AVL', 'AVM','ADL','ADM11', 'MBON', 'PVL', 'PDM', 'PDLother','PVM',
           'PFN', 'PFL', 'GLNO', 'Delta7', 'FB', 'LNO', 'LPsP', 'IbSpsP', 'OA', 'other', 'WLL',
           'hDelta', 'vDelta', 'PFR',
           'PFNa','PFNp','PFNm','PFNd','PFNv'),
  col = c( 367 ,   9,   34,    101,   32,    21,     58,    11,    12,    657,  517,
           468,   456,  467,   463,   464,   465,   466,    98,
           592,   591,   590,   589,   616,   617,   618,   619,   128,   130,
           259,   600,  26, 121,
           499,    499,    143,   144,    573,  640,
           
           430, 431, 630, 452,632, 103, 104, 651, 639, 105,
           
           420, 420, 420, 420, 76, 535, 535,535,535,
           657, 86, 656, 520, 468, 404,565, 447, 450, 651, 484,
           367, 9, 32,
           617,   618,   619,   128,   130)

)

supertype2Palette <- function(){
  CXTypes <- supertype(read_csv("C:/Users/turnerevansd/Dropbox (HHMI)/FIBSEM CX Paper Jan 2020/CX-cell-types042020",
                                col_types=cols(n.type = col_character())) %>% rename(databaseType=n.type))
  s2 <- unique(CXTypes$supertype2)
  pal <- paletteer::paletteer_d("Polychrome::palette36")[c(35,32,28,8,12,33,6,10,9,3,25,18,21,30,31,34,16,27,7,26,1,15,36)]
  names(pal) <- s2
  list(pal=pal,breaks=s2)
}

scale_color_CX_supertype <- function(...){
  pal <- supertype2Palette()
  scale_color_manual(values=pal$pal,breaks=pal$breaks,...)
}                                      

scale_fill_CX_supertype <- function(...){
  pal <- supertype2Palette()
  scale_fill_manual(values=pal$pal,breaks=pal$breaks,...)
}    
