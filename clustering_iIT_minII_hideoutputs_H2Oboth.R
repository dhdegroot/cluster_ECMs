require('tidyverse')
library("RColorBrewer")
library(extrafont)
font_import(pattern='calibri.ttf')
loadfonts(device = "win")

factor_diff_obj <- 1
log_shift_neg <- .9 # Has to be smaller than 1
log_shift_pos <- 0.001
log_scale <- TRUE
BIOMASS_ONE <- TRUE
CLUSTER_TERNARY <- TRUE

# Read the ECMs as row vectors (R convention), and add column names for each metabolite
ecms <- read.csv('data/iIT_minII_H2O_both_hideoutputs_01102020.csv', header=TRUE)

# interesting_ecms <- ecms %>%
#   filter(objective > 0)
interesting_ecms<-ecms

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
  filled_ecms[filled_ecms$objective!=0,] <- filled_ecms[filled_ecms$objective!=0,] %>%
    apply(1, function(x){x / x['objective']}) %>%
    t() %>%
    as.data.frame()
}

# Read in matching of metabolite ids and names
metab_info <- read_csv(file.path('data','metab_info_iIT.csv'),col_names=TRUE)
metab_names <- c()
for (col in colnames(filled_ecms)) {
  metab_names <- c(metab_names,metab_info[metab_info$id==col,]$name)
}
colnames(filled_ecms) <- metab_names

# Keep only metabolites that show variability
BM_rows = filled_ecms$Biomass!=0
constant_metabs <- c()
for (col in colnames(filled_ecms)) {
  if (col!="Biomass") {
    if (max(filled_ecms[BM_rows,col]) - min(filled_ecms[BM_rows,col]) < 0.01*max(abs(filled_ecms[BM_rows,col]))){
      constant_metabs <- c(constant_metabs,col)
    }
  }
}
filled_ecms <- filled_ecms[ , !(names(filled_ecms) %in% constant_metabs)]

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


clust_weights <- setNames(as.list(rep(1,length(colnames(filled_ecms)))), colnames(filled_ecms))
clust_weights[c('Biomass','D-Glucose','Sulfate','D-Alanine','L-Alanine','O2','H2O','H+')] <- c(1000,800,300,300,300,300,80,1)

weighted_ecms <- filled_ecms * clust_weights

# Cluster ECMs
row.order <- hclust(dist(weighted_ecms,method='manhattan'), method = "average")$order # clustering

filled_ecms[filled_ecms$Biomass==0,]<-filled_ecms[filled_ecms$Biomass==0,]*(factor_diff_obj)

# Cluster metabolites
sgn_ecms <- data.frame(filled_ecms)
sgn_ecms[sgn_ecms<0] = -1
sgn_ecms[sgn_ecms>0] = 1

col.order <- order(colSums(sgn_ecms[, names(sgn_ecms)!= 'Biomass']))
col.order <- c(col.order, which(names(filled_ecms)=='Biomass'))
if(CLUSTER_TERNARY){
  clust_weights <- setNames(as.list(rep(1,length(colnames(filled_ecms)))), colnames(filled_ecms))
  clust_weights[c('Biomass','D-Glucose','Sulfate','D-Alanine','L-Alanine','O2','H2O','H+')] <- c(1000,800,300,300,300,300,80,1)
  
  weighted_sgn_ecms <- sgn_ecms * clust_weights
  row.order <- hclust(dist(weighted_sgn_ecms,method='manhattan'), method = "average")$order # clustering
}
# col.order <- hclust(dist(t(filled_ecms),method='manhattan'), method = "average")$order
ordered_metabs <- attributes(filled_ecms)$names[col.order]

man_ordered_metabs <- c("L-Arginine","O2","D-Glucose","L-Alanine","D-Alanine","Sulfate","H+","H2O","Biomass")

# Order ECMs according to clustering
clustered_ecms <- filled_ecms[row.order,] %>% 
  mutate(ecm=1:n()) %>%
  gather('metabolite', 'stoichiometry', -ecm)

# Order metabolites according to clustering
clustered_ecms$metabolite <- factor(clustered_ecms$metabolite,
                                 levels=man_ordered_metabs)
if(log_scale){
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
} else {
  inv_get_labels <- function(orig){
    result <- orig
    result
  }
}

paper_colors = brewer.pal(3, 'RdYlBu')
clustered_ecms[clustered_ecms$stoichiometry==0,]$stoichiometry<-NA

# Render clusters
clustered_ecms %>%
  ggplot(aes(x=ecm, y=metabolite, fill=stoichiometry)) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 0, low = paper_colors[1], mid = "grey95", na.value = 'white',
                       high = paper_colors[3], space = "Lab", breaks=inv_get_labels(c(-100,-.1,0,1)),
                       labels=as.character(c(-100,-.1,0,1)),
                       guide=guide_colorbar(label.theme = element_text(family='Calibri',size=16), 
                                          title.theme = element_text(family='Calibri',size=22))) +
  geom_raster() +
  theme(axis.title.x = element_text(family='Calibri', size=22),
        axis.title.y = element_text(family='Calibri', size=22),
        axis.text.y = element_text(angle = 0, hjust = 1, family='Calibri', size=16),
        axis.text.x = element_blank())

