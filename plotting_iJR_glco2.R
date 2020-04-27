require('tidyverse')
library("RColorBrewer")
library("reshape2")
library("ggplot2")
library("hrbrthemes")

# Read the ECMs as row vectors (R convention), and add column names for each metabolite
ecms <- read.csv('data/iJR_hideallexceptglco2biomass.csv', header=TRUE)

glcox_nobm = ecms[(ecms$M_o2_e != 0) & (ecms$objective == 0),]
glcox_nobm = 10*glcox_nobm/sum(glcox_nobm)
glc_nobm = ecms[(ecms$M_glc__D_e!=0) &(ecms$M_o2_e == 0) & (ecms$objective == 0),]
glc_nobm = 10*glc_nobm/sum(glc_nobm)

norm_ecms = as.data.frame(t(apply(ecms, 1, function(x) {x / -x['objective']})))
norm_ecms$glc_yield = 1/norm_ecms$M_glc__D_e
norm_ecms$ecm_id = 1:nrow(norm_ecms)
norm_ecms[is.na(norm_ecms)] = 0

paper_colors = brewer.pal(3, 'RdYlBu')

ggplot(norm_ecms, aes(x=M_glc__D_e, y=M_o2_e, fill=glc_yield)) + 
  geom_point(size=6, stroke=1,pch=21, color='black') +
  coord_cartesian(xlim = c(0, 50), ylim = c(0, 40)) +
  scale_fill_gradient(low = paper_colors[1], high = paper_colors[3], space = "Lab",
                       guide=guide_colorbar(label.theme = element_text(family='Calibri',size=16), 
                                            title.theme = element_text(family='Calibri',size=22))) +
  theme(axis.title.x = element_text(family='Calibri', size=22,colour='#050C17FF'),
      axis.title.y = element_text(family='Calibri', size=22,colour='#050C17FF'),
      axis.text.y = element_text(family='Calibri',size=16,colour='#050C17FF'),
      axis.text.x = element_text(family='Calibri',size=16,colour='#050C17FF'),
      axis.line=element_line(colour='#050C17FF'),
      panel.background=element_blank(),
      panel.grid.major = element_line(colour='grey'))
