require('tidyverse')
library("RColorBrewer")

log_shift <- 1
log_scale <- TRUE

# Read the ECMs as row vectors (R convention), and add column names for each metabolite
ecms <- read.csv('data/conversions_e_coli_core.csv', header=TRUE)
interesting_ecms <- ecms %>% filter(objective > 0)

# Drop unused metabolites, and empty ECMs (they are normally not present to begin with)
col_sums <- colSums(abs(interesting_ecms))
row_sums <- rowSums(abs(interesting_ecms))
col_indices <- col_sums != 0
row_indices <- row_sums != 0
filled_ecms <- interesting_ecms[row_indices,col_indices] %>%
  apply(1, function(x){x / x['objective']}) %>%
  t() %>%
  as.data.frame()

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

# Determine number of clusters
wss <- (nrow(filled_ecms)-1)*sum(apply(filled_ecms,2,var))
for (i in 11:50) wss[i-10] <- sum(kmeans(filled_ecms,
                                     centers=i)$withinss)
plot(1:40, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

# K-Means Cluster Analysis
fit <- kmeans(filled_ecms, 40) # 5 cluster solution
# get cluster means
aggregate(filled_ecms,by=list(fit$cluster),FUN=mean)
# append cluster assignment
filled_ecms <- data.frame(filled_ecms, fit$cluster)

row.order <- order(filled_ecms$fit.cluster)
  
# Cluster metabolites
sgn_ecms <- filled_ecms
sgn_ecms[sgn_ecms<0] = -1
sgn_ecms[sgn_ecms>0] = 1

col.order <- order(colSums(sgn_ecms))
# col.order <- hclust(dist(t(filled_ecms),method='manhattan'), method = "average")$order
ordered_metabs <- attributes(filled_ecms)$names[col.order]

# Order ECMs according to clustering
filled_ecms <- filled_ecms[row.order,] %>% 
  as.data.frame() %>%
  mutate(ecm=1:n()) %>%
  gather('metabolite', 'stoich', -ecm)

# filled_ecms <- filled_ecms %>%
#  as.data.frame() %>%
#  mutate(ecm=1:n()) %>%
#  gather('metabolite', 'stoich', -ecm)

# Order metabolites according to clustering
filled_ecms$metabolite <- factor(filled_ecms$metabolite,
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

# Render clusters
filled_ecms %>%
  ggplot(aes(x=ecm, y=metabolite, fill=stoich)) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 0, low = paper_colors[1], mid = paper_colors[2],
                       high = paper_colors[3], space = "Lab", breaks=inv_get_labels(c(-40,-1,0,1,40)),
                       labels=as.character(c(-40,-1,0,1,40)),
                       guide=guide_colorbar(label.theme = element_text(family='Calibri',size=12), 
                                          title.theme = element_text(family='Calibri',size=18))) +
  geom_raster() +
  theme(axis.title.x = element_text(family='Calibri', size=18),
        axis.title.y = element_text(family='Calibri', size=18),
        axis.text.y = element_text(angle = 0, hjust = 1, family='Calibri'),
        axis.text.x = element_blank())

