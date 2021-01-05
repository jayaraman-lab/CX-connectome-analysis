library(cowplot)

mm2pt <- 0.352777778
update_geom_defaults("text", list(size = 6*mm2pt*1.111111,font_family="arial",face="plain"))
theme_paper <- function(...){
  theme_cowplot(font_size=7*1.111111,font_family="arial",rel_small=6/7,rel_tiny = 6/7,rel_large = 12/7, line_size = 0.5*mm2pt) + theme(strip.background = element_blank(),plot.tag=element_text(size=12**1.111111,face="bold",family = "arial"),...)
  }

theme_paper_grid <- function(...){
  theme_minimal_grid(font_size=7*1.111111,font_family="arial",rel_small=6/7,rel_tiny = 6/7,rel_large = 12/7, line_size = 0.5*mm2pt) + theme(strip.background = element_blank(),plot.tag=element_text(size=12*1.111111,face="bold",family = "arial"),...)
}

theme_paper_hgrid <- function(...){
  theme_minimal_hgrid(font_size=7*1.111111,font_family="arial",rel_small=6/7,rel_tiny = 6/7,rel_large = 12/7, line_size = 0.5*mm2pt) + theme(strip.background = element_blank(),plot.tag=element_text(size=12*1.111111,face="bold",family = "arial"),...)
}

theme_paper_vgrid <- function(...){
  theme_minimal_vgrid(font_size=7*1.111111,font_family="arial",rel_small=6/7,rel_tiny = 6/7,rel_large = 12/7, line_size = 0.5*mm2pt) + theme(strip.background = element_blank(),plot.tag=element_text(size=12*1.111111,face="bold",family = "arial"),...)
}

theme_paper_map <- function(...){
  theme_map(font_size=7*1.111111,font_family="arial",rel_small=6/7,rel_tiny = 6/7,rel_large = 12/7, line_size = 0.5*mm2pt) + theme(strip.background = element_blank(),plot.tag=element_text(size=12*1.111111,face="bold",family = "arial"),...)
}