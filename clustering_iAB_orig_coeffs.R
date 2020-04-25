require('tidyverse')
library("RColorBrewer")

# Read the ECMs as row vectors (R convention), and add column names for each metabolite
ecms <- read_csv('data/conversions_iAB.csv', col_names=FALSE)
names <- c('M_13dpg_c', 'M_23dpg_c', 'M_2kmb_c', 'M_2pg_c', 'M_35cgmp_c', 'M_3dhguln_c', 'M_3moxtyr_c', 'M_3pg_c', 'M_4pyrdx_c', 'M_5aop_c', 'M_5mdr1p_c', 'M_5mdru1p_c', 'M_5mta_c', 'M_5oxpro_c', 'M_6pgc_c', 'M_6pgl_c', 'M_ac_c', 'M_acald_c', 'M_acgam_c', 'M_acgam6p_c', 'M_acmana_c', 'M_acmanap_c', 'M_acnam_c', 'M_acnamp_c', 'M_ade_c', 'M_adn_c', 'M_adp_c', 'M_adprbp_c', 'M_adrnl_c', 'M_ahcys_c', 'M_akg_c', 'M_ala__L_c', 'M_alpa_hs_16_0_c', 'M_alpa_hs_18_1_c', 'M_alpa_hs_18_2_c', 'M_amet_c', 'M_ametam_c', 'M_amp_c', 'M_ap4a_c', 'M_arg__L_c', 'M_ascb__L_c', 'M_atp_c', 'M_band_c', 'M_bandmt_c', 'M_bilglcur_c', 'M_bilirub_c', 'M_biliverd_c', 'M_ca2_c', 'M_camp_c', 'M_cdp_c', 'M_cdpchol_c', 'M_cdpdag_hs_16_0_16_0_c', 'M_cdpdag_hs_16_0_18_1_c', 'M_cdpdag_hs_16_0_18_2_c', 'M_cdpdag_hs_18_1_18_1_c', 'M_cdpdag_hs_18_1_18_2_c', 'M_cdpdag_hs_18_2_16_0_c', 'M_cdpdag_hs_18_2_18_1_c', 'M_cdpea_c', 'M_chol_c', 'M_cholp_c', 'M_cl_c', 'M_cmp_c', 'M_co_c', 'M_co2_c', 'M_coa_c', 'M_cpppg3_c', 'M_crn_c', 'M_ctp_c', 'M_cys__L_c', 'M_dag_hs_16_0_16_0_c', 'M_dag_hs_16_0_18_1_c', 'M_dag_hs_16_0_18_2_c', 'M_dag_hs_18_1_18_1_c', 'M_dag_hs_18_1_18_2_c', 'M_dag_hs_18_2_16_0_c', 'M_dag_hs_18_2_18_1_c', 'M_dhap_c', 'M_dhdascb_c', 'M_dhmtp_c', 'M_dkmpp_c', 'M_dnad_c', 'M_dopa_c', 'M_e4p_c', 'M_etha_c', 'M_ethamp_c', 'M_f26bp_c', 'M_f6p_c', 'M_fad_c', 'M_fdp_c', 'M_fe2_c', 'M_fmn_c', 'M_for_c', 'M_fru_c', 'M_fum_c', 'M_g1p_c', 'M_g3p_c', 'M_g3pc_c', 'M_g6p_c', 'M_gal_c', 'M_gal1p_c', 'M_galt_c', 'M_gam_c', 'M_gam6p_c', 'M_gdp_c', 'M_glc__D_c', 'M_glcur_c', 'M_gln__L_c', 'M_glu__L_c', 'M_glucys_c', 'M_gly_c', 'M_glyc_c', 'M_glyc3p_c', 'M_gmp_c', 'M_gsn_c', 'M_gthox_c', 'M_gthrd_c', 'M_gtp_c', 'M_gua_c', 'M_guln__L_c', 'M_h_c', 'M_h2o_c', 'M_h2o2_c', 'M_hco3_c', 'M_hcys__L_c', 'M_hdca_c', 'M_hmbil_c', 'M_hxan_c', 'M_icit_c', 'M_imp_c', 'M_inost_c', 'M_ins_c', 'M_k_c', 'M_lac__D_c', 'M_lac__L_c', 'M_leuktrA4_c', 'M_leuktrB4_c', 'M_lgt__S_c', 'M_lnlc_c', 'M_lnlccoa_c', 'M_lnlccrn_c', 'M_lpchol_hs_16_0_c', 'M_lpchol_hs_18_1_c', 'M_lpchol_hs_18_2_c', 'M_mal__L_c', 'M_man_c', 'M_man1p_c', 'M_man6p_c', 'M_mepi_c', 'M_met__L_c', 'M_mi1345p_c', 'M_mi134p_c', 'M_mi145p_c', 'M_mi14p_c', 'M_mi1p__D_c', 'M_mthgxl_c', 'M_na1_c', 'M_nac_c', 'M_nad_c', 'M_nadh_c', 'M_nadp_c', 'M_nadph_c', 'M_ncam_c', 'M_nh4_c', 'M_nicrns_c', 'M_nicrnt_c', 'M_nmn_c', 'M_normete__L_c', 'M_nrpphr_c', 'M_o2_c', 'M_oaa_c')
names <- append(names, c('M_ocdcea_c', 'M_odecoa_c', 'M_odecrn_c', 'M_orn_c', 'M_orot_c', 'M_orot5p_c', 'M_pa_hs_16_0_16_0_c', 'M_pa_hs_16_0_18_1_c', 'M_pa_hs_16_0_18_2_c', 'M_pa_hs_18_1_18_1_c', 'M_pa_hs_18_1_18_2_c', 'M_pa_hs_18_2_16_0_c', 'M_pa_hs_18_2_18_1_c', 'M_pail45p_hs_16_0_16_0_c', 'M_pail45p_hs_16_0_18_1_c', 'M_pail45p_hs_16_0_18_2_c', 'M_pail45p_hs_18_1_18_1_c', 'M_pail45p_hs_18_1_18_2_c', 'M_pail45p_hs_18_2_16_0_c', 'M_pail45p_hs_18_2_18_1_c', 'M_pail4p_hs_16_0_16_0_c', 'M_pail4p_hs_16_0_18_1_c', 'M_pail4p_hs_16_0_18_2_c', 'M_pail4p_hs_18_1_18_1_c', 'M_pail4p_hs_18_1_18_2_c', 'M_pail4p_hs_18_2_16_0_c', 'M_pail4p_hs_18_2_18_1_c', 'M_pail_hs_16_0_16_0_c', 'M_pail_hs_16_0_18_1_c', 'M_pail_hs_16_0_18_2_c', 'M_pail_hs_18_1_18_1_c', 'M_pail_hs_18_1_18_2_c', 'M_pail_hs_18_2_16_0_c', 'M_pail_hs_18_2_18_1_c', 'M_pchol_hs_16_0_16_0_c', 'M_pchol_hs_16_0_18_1_c', 'M_pchol_hs_16_0_18_2_c', 'M_pchol_hs_18_1_18_1_c', 'M_pchol_hs_18_1_18_2_c', 'M_pchol_hs_18_2_16_0_c', 'M_pchol_hs_18_2_18_1_c', 'M_pdx5p_c', 'M_pe_hs_16_0_16_0_c', 'M_pe_hs_16_0_18_1_c', 'M_pe_hs_16_0_18_2_c', 'M_pe_hs_18_1_18_1_c', 'M_pe_hs_18_1_18_2_c', 'M_pe_hs_18_2_16_0_c', 'M_pe_hs_18_2_18_1_c', 'M_pep_c', 'M_phe__L_c', 'M_pheme_c', 'M_phpyr_c', 'M_pi_c', 'M_pmtcoa_c', 'M_pmtcrn_c', 'M_ppbng_c', 'M_ppi_c', 'M_ppp9_c', 'M_pppg9_c', 'M_prpp_c', 'M_ptrc_c', 'M_pyam5p_c', 'M_pydam_c', 'M_pydx_c', 'M_pydx5p_c', 'M_pydxn_c', 'M_pyr_c', 'M_r1p_c', 'M_r5p_c', 'M_ribflv_c', 'M_rnam_c', 'M_ru5p__D_c', 'M_s7p_c', 'M_sbt__D_c', 'M_spmd_c', 'M_sprm_c', 'M_thm_c', 'M_thmmp_c', 'M_thmpp_c', 'M_thmtp_c', 'M_udp_c', 'M_udpg_c', 'M_udpgal_c', 'M_udpglcur_c', 'M_ump_c', 'M_uppg3_c', 'M_urea_c', 'M_uri_c', 'M_utp_c', 'M_xmp_c', 'M_xu5p__D_c', 'M_xylt_c', 'M_xylu__D_c', 'M_xylu__L_c', 'M_35cgmp_e', 'M_3moxtyr_e', 'M_4pyrdx_e', 'M_5aop_e', 'M_ac_e', 'M_acald_e', 'M_acnam_e', 'M_ade_e', 'M_adn_e', 'M_adrnl_e', 'M_ala__L_e', 'M_arg__L_e', 'M_ascb__L_e', 'M_bilglcur_e', 'M_ca2_e', 'M_camp_e', 'M_chol_e', 'M_cl_e', 'M_co_e', 'M_co2_e', 'M_cys__L_e', 'M_dhdascb_e', 'M_dopa_e', 'M_etha_e', 'M_fe2_e', 'M_fru_e', 'M_fum_e', 'M_gal_e', 'M_gam_e', 'M_glc__D_e', 'M_gln__L_e', 'M_gluala_e', 'M_gly_e', 'M_glyc_e', 'M_gthox_e', 'M_h_e', 'M_h2o_e', 'M_h2o2_e', 'M_hco3_e', 'M_hcys__L_e', 'M_hdca_e', 'M_hxan_e', 'M_ins_e', 'M_k_e', 'M_lac__D_e', 'M_lac__L_e', 'M_leuktrA4_e', 'M_leuktrB4_e', 'M_lnlc_e', 'M_mal__L_e', 'M_man_e', 'M_mepi_e', 'M_met__L_e', 'M_na1_e', 'M_nac_e', 'M_ncam_e', 'M_nh4_e', 'M_normete__L_e', 'M_nrpphr_e', 'M_o2_e', 'M_ocdcea_e', 'M_orot_e', 'M_phe__L_e', 'M_pi_e', 'M_ptrc_e', 'M_pydam_e', 'M_pydx_e', 'M_pydxn_e', 'M_pyr_e', 'M_ribflv_e', 'M_spmd_e', 'M_sprm_e', 'M_thm_e', 'M_thmmp_e', 'M_urea_e', 'M_uri_e', 'objective'))
colnames(ecms) <- names

# Drop uninteresting ECMs
interesting_ecms <- ecms %>%
  filter(objective > 0 & M_gal_e == 0 & M_glc__D_e < 0 & M_lac__D_e == 0 & M_lac__L_e == 0 & M_man_e == 0 & M_fru_e == 0)

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
min_pos = min(filled_ecms[filled_ecms>0])
filled_ecms[filled_ecms>0] <- log(filled_ecms[filled_ecms>0]) - 1.1* log(min_pos)
# Then log-scale all negatives, but subtract large number such that all numbers remain negative
max_neg = max(filled_ecms[filled_ecms<0])
filled_ecms[filled_ecms<0] <- -log(-filled_ecms[filled_ecms<0]) + 1.1* log(-max_neg)

# Free memory
# rm(ecms)

row.order <- hclust(dist(filled_ecms,method='manhattan'), method = "complete")$order # clustering

# Cluster metabolites
col.order <- hclust(dist(t(filled_ecms),method='manhattan'), method = "complete")$order
ordered_metabs <- attributes(filled_ecms)$names[col.order]

# Order ECMs according to clustering
filled_ecms <- filled_ecms[row.order,] %>%
  as.data.frame() %>%
  mutate(ecm=1:n()) %>%
  gather('metabolite', 'stoich', -ecm)

# Order metabolites according to clustering
filled_ecms$metabolite <- factor(filled_ecms$metabolite,
                                 levels=ordered_metabs)

get_labels <- function(orig){
  result = rep(NA,length(orig))
  for(i in 1:length(orig)){
    if (orig[i]<0){
      result[i] = -exp(-(orig[i]-1.1*log(-max_neg)))
    }else if(orig[i]>0){
      result[i] = exp(orig[i]+1.1*log(min_pos))
    }else{
      result[i] = 0
    }
  }
  as.character(format(result,digits=3))
}

# Render clusters
filled_ecms %>%
  ggplot(aes(x=ecm, y=metabolite, fill=stoich)) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 0, low = "magenta", mid = "black",
                       high = "cyan", space = "Lab" ) +
  geom_raster() +
  theme(axis.title.x = element_text(family='Calibri', size=18),
        axis.title.y = element_text(family='Calibri', size=18),
        axis.text.y = element_text(angle = 0, hjust = 1, family='Calibri'),
        axis.text.x = element_blank())

