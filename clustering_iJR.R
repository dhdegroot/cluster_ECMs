require('tidyverse')
library("RColorBrewer")
library(extrafont)
font_import(pattern='calibri')
loadfonts(device = "win")

log_shift = 1.1
shift_pos = 6

# Read the ECMs as row vectors (R convention), and add column names for each metabolite
ecms <- read.csv('data/iJR_hideallexceptglco2biomass.csv', header=TRUE)

# Drop uninteresting ECMs
interesting_ecms <- ecms %>%
  filter(objective > 0)

# Drop unused metabolites, and empty ECMs (they are normally not present to begin with)
col_sums <- colSums(interesting_ecms)
row_sums <- rowSums(interesting_ecms)
col_indices <- col_sums != 0
row_indices <- row_sums != 0
filled_ecms <- interesting_ecms[row_indices,col_indices] %>%
  apply(1, function(x){x / x['objective']}) %>%
  t() %>%
  as.data.frame()
rownames(filled_ecms) <- 1:nrow(filled_ecms)

# Log-scale all coefficients except for objective coefficients. Then shift numbers such that they are all negative again
max_neg = max(filled_ecms[filled_ecms<0])
filled_ecms[filled_ecms<0] <- -log(-filled_ecms[filled_ecms<0]) + log_shift * log(-max_neg)

# Shift positive numbers too, to better use colourscale
filled_ecms[filled_ecms>0] <- filled_ecms[filled_ecms>0]*shift_pos

row_clust <- hclust(dist(filled_ecms,method='manhattan'), method = "complete") # clustering
row.order <- row_clust$order

# Convert hclust into a dendrogram and plot
hcd <- as.dendrogram(row_clust)
# Default plot
plot(hcd, type = "rectangle", ylab = "Height")

# Cluster metabolites
col.order <- hclust(dist(t(filled_ecms),method='manhattan'), method = "complete")$order
ordered_metabs <- attributes(filled_ecms)$names[col.order]

# Order ECMs according to clustering
clustered_ecms <- filled_ecms[row.order,] %>%
  as.data.frame() %>%
  mutate(ecm=1:n()) %>%
  gather('metabolite', 'stoich', -ecm)

# Order metabolites according to clustering
clustered_ecms$metabolite <- factor(clustered_ecms$metabolite,
                                    levels=ordered_metabs)

inv_get_labels <- function(orig){
  result = rep(NA,length(orig))
  for(i in 1:length(orig)){
    if (orig[i]<0){
      result[i] = -log(-orig[i]) + log_shift*log(-max_neg)
    }else if(orig[i]>0){
      result[i] = orig[i]*shift_pos
    }else{
      result[i] = 0
    }
  }
  result
}

paper_colors = brewer.pal(3, 'RdYlBu')

# Render clusters
clustered_ecms %>%
  ggplot(aes(x=ecm, y=metabolite, fill=stoich)) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 0, low = paper_colors[1], mid = "white",
                       high = paper_colors[3], space = "Lab", breaks=inv_get_labels(c(-10,-.001,0,1)),
                       labels=as.character(c(-10,-1e-3,0,1),format='e'),
                       na.value='white',
                       guide=guide_colorbar(label.theme = element_text(family='Calibri',size=16), 
                                            title.theme = element_text(family='Calibri',size=22))) +
  geom_raster() +
  theme(axis.title.x = element_text(family='Calibri', size=22),
        axis.title.y = element_text(family='Calibri', size=22),
        axis.text.y = element_text(angle = 0, hjust = 1, family='Calibri', size=16),
        axis.text.x = element_blank())
