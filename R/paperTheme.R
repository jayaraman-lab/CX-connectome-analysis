library(cowplot)

update_geom_defaults("text", list(size = 6*0.35))
theme_paper <- function(...){
  theme_cowplot(font_size=7,font_family="Arial",rel_small=6/7,rel_tiny = 5/7,rel_large = 8/7) + theme(strip.background = element_blank(),...)
  }

