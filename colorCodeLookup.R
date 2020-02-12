# color code
getSimpleTypeNames <- function(mydata){
  simpleTypes = gsub("Delta0[[:alnum:]]*", "Delta0", mydata)
  simpleTypes = gsub("Delta12[[:alnum:]]*", "Delta12", simpleTypes)
  simpleTypes = gsub("Delta6[[:alnum:]]*", "Delta6", simpleTypes)
  simpleTypes = gsub("Delta0-Delta12[(]*[[:alnum:]]*[)]*", "Delta0-12", simpleTypes)
  simpleTypes = gsub("FB[[:alnum:]]*", "FB", simpleTypes)
  simpleTypes = gsub("FB_[[:alnum:]]*", "FB", simpleTypes)
  simpleTypes = gsub("FB-Q[[:alnum:]]*", "FB-Q", simpleTypes)
  simpleTypes = gsub("EB-Q[[:alnum:]]*", "EB-Q", simpleTypes)
  #simpleTypes = gsub("ExR[[:alnum:]]*", "ExR", simpleTypes)
  simpleTypes = gsub("PB[[:alnum:]]*", "PB", simpleTypes)
  simpleTypes = gsub("SLP-AB[[:alnum:]]*", "SLP-AB", simpleTypes)
  simpleTypes = gsub("R1[[:alnum:]]*_[[:alnum:]]*", "R1", simpleTypes)
  simpleTypes = gsub("R2[[:alnum:]]*_[[:alnum:]]*", "R2", simpleTypes)
  simpleTypes = gsub("R3a_[[:alnum:]]*", "R3a", simpleTypes)
  simpleTypes = gsub("R3d_[[:alnum:]]*", "R3d", simpleTypes)
  simpleTypes = gsub("R3m_[[:alnum:]]*", "R3m", simpleTypes)
  simpleTypes = gsub("R3p_[[:alnum:]]*", "R3p", simpleTypes)
  simpleTypes = gsub("R3w_[[:alnum:]]*", "R3w", simpleTypes)
  simpleTypes = gsub("R4d_[[:alnum:]]*", "R4d", simpleTypes)
  simpleTypes = gsub("R4m[[:alnum:]]*", "R4m", simpleTypes)
  simpleTypes = gsub("R5[[:alnum:]]*", "R5", simpleTypes)
  simpleTypes = gsub("R6[[:alnum:]]*", "R6", simpleTypes)
  #simpleTypes = gsub("TuBu[[:alnum:]]*", "TuBu", simpleTypes)
  simpleTypes = gsub("TuBu09_[[:alnum:]]*", "TuBu09", simpleTypes)
  simpleTypes = gsub("PDM21a_[[:alnum:]]*", "PDM21a", simpleTypes)
  simpleTypes = gsub("PEN_a\\(PEN1\\)", "PEN1", simpleTypes)
  simpleTypes = gsub("PEN_b\\(PEN2\\)", "PEN2", simpleTypes)
  simpleTypes = gsub("MC[[:alnum:]]*", "MC", simpleTypes)
  simpleTypes = gsub("TuTu[[:alnum:]]*", "TuTu", simpleTypes)
  simpleTypes = gsub("ADM06b_[[:alnum:]]*_[[:alnum:]]*", "ADM06b", simpleTypes)
  simpleTypes = gsub("MBON[[:alnum:]]*", "MBON", simpleTypes)
  simpleTypes = gsub("AVL[[:alnum:]]*_[[:alnum:]]*_pct", "AVL", simpleTypes)
  simpleTypes = gsub("AVL[[:alnum:]]*_pct", "AVL", simpleTypes)
  simpleTypes = gsub("AVM[[:alnum:]]*_pct", "AVM", simpleTypes)
  simpleTypes = gsub("ADL.*_pct", "ADL", simpleTypes)
  simpleTypes = gsub("PVL[[:alnum:]]*_pct", "PVL", simpleTypes)
  simpleTypes = gsub("PVM.*_pct", "PVM", simpleTypes)
  simpleTypes = gsub("PDL27e_pct", "PDL27e", simpleTypes)
  simpleTypes = gsub("PDL[12x,20s]{1}.*_pct", "PDLother", simpleTypes)
  simpleTypes = gsub("ADM11[[:alnum:]]*_pct", "ADM11", simpleTypes)
  simpleTypes = gsub("ADM06p_pct", "ADM06p", simpleTypes)
  simpleTypes = gsub("ADM06d_pct", "ADM06d", simpleTypes)
  simpleTypes = gsub("ADM03.*_pct", "ADM03", simpleTypes)
  simpleTypes = gsub("ADM[02h, 04a, 05oh, 10h, 01w]{1}.*_pct", "ADMother", simpleTypes)
  simpleTypes = gsub("olfactory multi .*", "other", simpleTypes)
  simpleTypes = gsub("PDM09[[:alnum:]]*_pct", "PDM09", simpleTypes)
  simpleTypes = gsub("PDM14j.*_pct", "PDM14j", simpleTypes)
  simpleTypes = gsub("PDM14[m,r,d]{1}.*_pct", "PDM14other", simpleTypes)
  simpleTypes = gsub("PDM[26, 28, 08k]{1}[[:alnum:]]*_pct", "PDMother", simpleTypes)
  simpleTypes = gsub("LN[[:alnum:]]{1}.*", "LN", simpleTypes)
  simpleTypes = gsub("PFN.*", "PFN", simpleTypes)
  simpleTypes = gsub("Delta7*", "Delta7", simpleTypes)
  simpleTypes = gsub("P6\\-8P9", "other", simpleTypes)
  simpleTypes = gsub("P1\\-9", "other", simpleTypes)
  simpleTypes = gsub("PVL08j_a_pct", "other", simpleTypes)
  simpleTypes = gsub("LC[[:alnum:]]{2}", "LC", simpleTypes)
  simpleTypes = gsub("LCm", "LC", simpleTypes)
  simpleTypes = gsub("OA-.*", "OA", simpleTypes)
  return(simpleTypes)
}

colorValueLookup = data.frame(
  type = c('R1', 'R2', 'R3a', 'R3d', 'R3m', 'R3p',  'R3w', 'R4d', 'R4m', 'R5', 'R6',
           'ExR1','ExR2','ExR3','ExR4','ExR5','ExR6','ExR7','ExR8',
           'TuBu01', 'TuBu02', 'TuBu03', 'TuBu04', 'TuBu05', 'TuBu06', 'TuBu07', 'TuBu08', 'TuBu09', 'TuBu10',
           'PDM21a', 'MC',  'TuTu', 'LC',
           'EPG', 'EPGt', 'PEN1', 'PEN2', 'PEG', 'EQ5',
           'PDM14j', 'ADM06d','ADM06p', 'ADM06b','PDL27e','ADM06s', 'PDM09','PDMother','PDM14other','ADM03',
           'AVL', 'AVM','ADL','ADM11', 'MBON', 'PVL', 'PDM', 'PDLother','PVM',
           'PFN', 'GLN', 'Delta7', 'FB', 'LN', 'LPsP', 'IbSpsP', 'OA', 'other', 'ADMother' ),
  col = c( 367 ,   9,   34,    101,   32,    21,     58,    11,    12,    657,  517,
           468,   456,  467,   463,   464,   465,   466,    98,
           592,   591,   590,   589,   616,   617,   618,   619,   128,   130,
           259,   600,  26, 121,
           499,    499,    143,   144,    573,  640, 
           639, 430, 431, 630, 452,632, 103, 104,105, 651,
           420, 420, 420, 420, 76, 535, 535,535,535,
           657, 656, 520, 468, 404,565, 447, 450, 651, 484)
)