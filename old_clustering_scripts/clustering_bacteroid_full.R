require('tidyverse')
library("RColorBrewer")
library(extrafont)
font_import(pattern='calibri.ttf')
loadfonts(device = "win")

log_shift <- 1.1 # Has to be greater than 1
log_scale <- TRUE
BIOMASS_ONE <- TRUE
CLUSTER_TERNARY <- FALSE

# Read the ECMs as row vectors (R convention), and add column names for each metabolite
ecms <- as.data.frame(read_csv('data/bacteroid_split_before_computation.csv', col_names=TRUE))

# Read in matching of metabolite ids and names
metab_info <- read_csv(file.path('data','metab_info_bacteroid.csv'),col_names=TRUE)
metab_names <- c()
for (col in colnames(ecms)) {
  metab_names <- c(metab_names,metab_info[metab_info$id==col,]$name)
}
colnames(ecms) <- metab_names

interesting_ecms <- ecms %>% filter(L_Malate_e0!= 0)

# Drop unused metabolites, and empty ECMs (they are normally not present to begin with)
col_sums <- colSums(abs(interesting_ecms))
row_sums <- rowSums(abs(interesting_ecms))
col_indices <- col_sums != 0
row_indices <- row_sums != 0
filled_ecms <- interesting_ecms[row_indices,col_indices] %>%
  apply(1, function(x){x / sum(abs(x))}) %>%
  t() %>%
  as.data.frame()

if(BIOMASS_ONE){
  filled_ecms[filled_ecms$Biomass!=0,] <- filled_ecms[filled_ecms$Biomass!=0,] %>%
    apply(1, function(x){x / x['Biomass']}) %>%
    t() %>%
    as.data.frame()
}

# Log-scale all coefficients
# First log-scale all positives, but add large number such that all numbers remain positive
# Get smallest positive number
if(log_scale){
  min_pos = min(filled_ecms[filled_ecms>0])
  filled_ecms[filled_ecms>0] <- log(filled_ecms[filled_ecms>0]) - log_shift * log(min_pos)
  # Then log-scale all negatives, but subtract large number such that all numbers remain negative
  max_neg = max(filled_ecms[filled_ecms<0])
  filled_ecms[filled_ecms<0] <- -log(-filled_ecms[filled_ecms<0]) + log_shift * log(-max_neg)
}

clust_weights <- setNames(as.list(rep(1,length(colnames(filled_ecms)))), colnames(filled_ecms))
clust_weights[c('NH3_e0','L_Malate_e0','O2_e0','Biomass')] <- c(60,80,100,100)

weighted_ecms <- filled_ecms * clust_weights
# Cluster ECMs
row.order <- hclust(dist(weighted_ecms,method='manhattan'), method = "average")$order # clustering

factor_diff_obj <- mean(filled_ecms[filled_ecms$Biomass!=0,]$'O2_e0') / mean(filled_ecms[filled_ecms$Biomass==0,]$'O2_e0')
filled_ecms[filled_ecms$Biomass==0,]<-filled_ecms[filled_ecms$Biomass==0,]*(factor_diff_obj)

# Cluster metabolites
sgn_ecms <- data.frame(filled_ecms)
sgn_ecms[sgn_ecms<0] = -1
sgn_ecms[sgn_ecms>0] = 1

col.order <- order(colSums(sgn_ecms[, names(sgn_ecms)!= 'Biomass']))
col.order <- c(col.order, which(names(filled_ecms)=='Biomass'))
if(CLUSTER_TERNARY){
  row.order <- hclust(dist(sgn_ecms,method='manhattan'), method = "average")$order # clustering
}
# col.order <- hclust(dist(t(filled_ecms),method='manhattan'), method = "average")$order
ordered_metabs <- attributes(filled_ecms)$names[col.order]

# Order ECMs according to clustering
clustered_ecms <- filled_ecms[row.order,] %>% 
  mutate(ecm=1:n()) %>%
  gather('metabolite', 'stoich', -ecm)

# Order metabolites according to clustering
clustered_ecms$metabolite <- factor(clustered_ecms$metabolite,
                                 levels=ordered_metabs)
if(log_scale){
  get_labels <- function(orig){
    result = rep(NA,length(orig))
    for(i in 1:length(orig)){
      if (orig[i]<0){
        result[i] = -exp(-(orig[i]-log_shift*log(-max_neg)))
      }else if(orig[i]>0){
        result[i] = exp(orig[i]+log_shift*log(min_pos))
      }else{
        result[i] = 0
      }
    }
    as.character(format(result,digits=3))
  }
  
  inv_get_labels <- function(orig){
    result = rep(NA,length(orig))
    for(i in 1:length(orig)){
      if (orig[i]<0){
        result[i] = -log(-orig[i]) + log_shift*log(-max_neg)
      }else if(orig[i]>0){
        result[i] = log(orig[i]) - log_shift*log(min_pos)
      }else{
        result[i] = 0
      }
    }
    result
  }
} else {
  inv_get_labels <- function(orig){
  result <- orig
  result
  }  
}

paper_colors = brewer.pal(3, 'RdYlBu')
clustered_ecms[clustered_ecms$stoich==0,]$stoich<-NA

# Render clusters
clustered_ecms %>%
  ggplot(aes(x=ecm, y=metabolite, fill=stoich)) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 0, low = paper_colors[1], mid = "grey90", na.value = 'white',
                       high = paper_colors[3], space = "Lab", breaks=inv_get_labels(c(-100,-.1,0,.1,100)),
                       labels=as.character(c(-100,-.1,0,.1,100)),
                       guide=guide_colorbar(label.theme = element_text(family='Calibri',size=16), 
                                          title.theme = element_text(family='Calibri',size=22))) +
  geom_raster() +
  theme(axis.title.x = element_text(family='Calibri', size=22),
        axis.title.y = element_text(family='Calibri', size=22),
        axis.text.y = element_text(angle = 0, hjust = 1, family='Calibri', size=16),
        axis.text.x = element_blank())

