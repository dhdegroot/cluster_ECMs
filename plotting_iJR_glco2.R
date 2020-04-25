require('tidyverse')
library("RColorBrewer")
library("reshape2")
library("ggplot2")
library("hrbrthemes")

# Read the ECMs as row vectors (R convention), and add column names for each metabolite
ecms <- read.csv('data/iJR_hideallexceptglco2biomass.csv', header=TRUE)

norm_ecms = as.data.frame(t(apply(ecms, 1, function(x) {x / -x['objective']})))
norm_ecms$ecm_id = 1:nrow(norm_ecms)

norm_ecms_td <- melt(norm_ecms, id.vars='ecm_id', variable.name = "metabolite", value.name = "flux")

ggplot(norm_ecms, aes(x='M_glc__D_e', y='M_o2_e')) + geom_point()
ggplot(norm_ecms, aes(x=M_glc__D_e, y=M_o2_e)) + 
  geom_point(size=6) + theme_ipsum()


# Log-scale all coefficients except for objective coefficients
# Remove the objective from dataframe, so that all numbers are negative
filled_ecms_substrates <- filled_ecms[!names(filled_ecms) %in% c("objective")]
log_filled_ecms_substrates <- log(-filled_ecms_substrates)

# Replace Inf's by a low number (we want the difference in the clustering between zero and non-zero 
# to be large, but not infinite)
min_noninf = min(log_filled_ecms_substrates[log_filled_ecms_substrates!=-Inf])
log_filled_ecms_substrates[log_filled_ecms_substrates == -Inf] = 10* min_noninf

# Free memory
rm(ecms)

row.order <- hclust(dist(log_filled_ecms_substrates,method='manhattan'), method = "complete")$order # clustering

# Cluster metabolites
col.order <- hclust(dist(t(log_filled_ecms_substrates),method='manhattan'), method = "complete")$order
ordered_metabs <- attributes(log_filled_ecms_substrates)$names[col.order]
man_ordered_metabs <- c("M_ala__L_e","M_ala__D_e","M_arg__L_e","M_o2_e","M_pi_e","M_h_e","M_nh4_e","M_so4_e","M_pheme_e",
                        "M_fe2_e","M_his__L_e", "M_val__L_e","M_leu__L_e", "M_ile__L_e", "M_met__L_e","M_pime_e","M_thm_e")
log_filled_ecms_substrates[log_filled_ecms_substrates<min_noninf]<-NA

# Order ECMs according to clustering
log_filled_ecms_substrates <- log_filled_ecms_substrates[row.order,] %>%
  as.data.frame() %>%
  mutate(ecm=1:n()) %>%
  gather('metabolite', 'stoich', -ecm)

# Order metabolites according to clustering
log_filled_ecms_substrates$metabolite <- factor(log_filled_ecms_substrates$metabolite,
                                                levels=man_ordered_metabs)

# Render clusters
log_filled_ecms_substrates %>%
  ggplot(aes(x=ecm, y=metabolite, fill=stoich)) +
  geom_tile() + 
  scale_fill_gradient(guide=guide_colorbar(reverse=FALSE),
                      labels = function(orig){as.character(format(-exp(orig),digits=3))},
                      na.value='grey10') +
  geom_raster() +
  theme(axis.text.y = element_text(angle = 0, hjust = 1),
        axis.text.x = element_blank())
