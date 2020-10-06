require('tidyverse')
library("RColorBrewer")
library(extrafont)
font_import(pattern='calibri.ttf')
loadfonts(device = "win")

factor_diff_obj <- 1
log_shift_neg <- .2 # Has to be smaller than 1
log_shift_pos <- .05
BIOMASS_ONE <- TRUE
CLUSTER_TERNARY <- FALSE
log_scale <- TRUE

# Read the ECMs as row vectors (R convention), and add column names for each metabolite
ecms <- read.csv('data/conversions_ecolicore_hideoutputs_tagPDH.csv', header=TRUE)

# Read in matching of metabolite ids and names
metab_info <- read_csv(file.path('data','metab_info_ecolicore.csv'),col_names=TRUE)
metab_names <- c()
for (col in colnames(ecms)) {
  metab_names <- c(metab_names,metab_info[metab_info$id==col,]$name)
}
colnames(ecms) <- metab_names

interesting_ecms <- ecms

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
  filled_ecms[filled_ecms>0] <- log((filled_ecms[filled_ecms>0])/(log_shift_pos*min_pos))
  # Then log-scale all negatives, but subtract large number such that all numbers remain negative
  max_neg = max(filled_ecms[filled_ecms<0])
  filled_ecms[filled_ecms<0] <- -log((filled_ecms[filled_ecms<0])/(log_shift_neg*max_neg))
}

clust_weights = list("CO2"=50,"D-Glucose"=100,"H2O"=10,"Ammonium"=1,"O2"=4,"Phosphate"=1,"Biomass"=100,"PDH-flux"=100)

weighted_ecms <- filled_ecms * clust_weights
# Cluster ECMs
row.order <- hclust(dist(weighted_ecms,method='manhattan'), method = "average")$order # clustering

filled_ecms[filled_ecms$Biomass==0,]<-filled_ecms[filled_ecms$Biomass==0,]*(factor_diff_obj)

# Cluster metabolites
sgn_ecms <- data.frame(filled_ecms)
sgn_ecms[sgn_ecms<0] = -1
sgn_ecms[sgn_ecms>0] = 1

col.order <- order(colSums(sgn_ecms))
col.order <- col.order[col.order!=which(names(filled_ecms) == "Biomass")]
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

inv_get_labels <- function(orig){
  result = rep(NA,length(orig))
  for(i in 1:length(orig)){
    if (orig[i]<0){
      result[i] = -log((orig[i])/(log_shift_neg*max_neg))
    }else if(orig[i]>0){
      result[i] = log((orig[i])/(log_shift_pos*min_pos))
    }else{
      result[i] = 0
    }
  }
  result
}

paper_colors = brewer.pal(3, 'RdYlBu')
clustered_ecms[clustered_ecms$stoich==0,]$stoich<-NA

# Render clusters
clustered_ecms %>%
  ggplot(aes(x=ecm, y=metabolite, fill=stoich)) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 0, low = paper_colors[1], mid = "grey90",
                       high = paper_colors[3], space = "Lab", breaks=inv_get_labels(c(-16,-2,0,1)),
                       labels=as.character(c(-16,-2,0,1),format='e'),
                       na.value='white',
                       guide=FALSE) +
  geom_raster() +
  theme(axis.title.x = element_text(family='Calibri', size=22),
        axis.title.y = element_text(family='Calibri', size=22),
        axis.text.y = element_text(angle = 0, hjust = 1, family='Calibri', size=16),
        axis.text.x = element_blank())
