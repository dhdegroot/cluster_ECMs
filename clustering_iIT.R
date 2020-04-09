require('tidyverse')

# Read the ECMs as row vectors (R convention), and add column names for each metabolite
ecms <- read.csv('data/iIT341_allsubstrates_to_biomass.csv', header=TRUE)

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

# Log-scale all coefficients except for objective coefficients
# Remove the objective from dataframe, so that all numbers are negative
filled_ecms_substrates <- filled_ecms[!names(filled_ecms) %in% c("objective")]
log_filled_ecms_substrates <- log(-filled_ecms_substrates)

# Free memory
rm(ecms)

cluster <- hclust(dist(log_filled_ecms_substrates), method = "complete") # clustering

# Populate a list of absolute ECM indices, such that they are accompanied by their own
# cluster, and similar clusters are placed near theirs.

row.order <- c()
for(abs_cluster_index in cluster.order) {
  print(paste('Doing cluster', abs_cluster_index))
  kms_index <- floor((abs_cluster_index-1)/20) + 1
  rel_index <- ((abs_cluster_index-1) %% 20) + 1
  rows <- as.numeric(rownames(as.data.frame(kms[[kms_index]]$cluster[which(kms[[kms_index]]$cluster == rel_index)])))
  if (length(rows) > 1){
    #### This part is where clustering goes wrong. Replace with 'row.order <- append(row.order, rows)' to see difference with unordered ECMs within clusters.
    rows_ordered <- hclust(dist(filled_ecms[rows,]), method = "average")$order
    row.order <- append(row.order, rows_ordered)
  } else {
    row.order <- append(row.order, rows)
  }
}

# Cluster metabolites
col.order <- hclust(dist(t(filled_ecms)))$order

# Render plot
filled_ecms[rows, col.order] %>%
  # top_n(1000) %>%
  mutate(ecm=1:n()) %>%
  gather('metabolite', 'stoich', -ecm) %>%
  # group_by(ecm) %>%
  # mutate(norm_stoich=stoich/max(abs(stoich))) %>%
  ggplot(aes(x=ecm, y=metabolite, fill=stoich)) +
  geom_tile() + 
  scale_fill_gradient2(midpoint = 0, low = "magenta", mid = "black",
                       high = "cyan", space = "Lab" ) +
  geom_raster() +
  theme(axis.text.y = element_text(angle = 0, hjust = 1),
        axis.text.x = element_blank())

# Render clusters
norm_centers[cluster.order,] %>%
  as.data.frame() %>%
  mutate(center=1:n()) %>%
  gather('metabolite', 'stoich', -center) %>%
  ggplot(aes(x=center, y=metabolite, fill=stoich)) +
  geom_tile() + 
  scale_fill_gradient2(midpoint = 0, low = "magenta", mid = "black",
                       high = "cyan", space = "Lab" ) +
  geom_raster() +
  theme(axis.text.y = element_text(angle = 0, hjust = 1),
        axis.text.x = element_blank())
  