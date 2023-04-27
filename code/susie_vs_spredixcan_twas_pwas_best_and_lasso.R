# **AD-PWAS by PrediXcan (using LASSO model and best model (both)) for Xiaotong Sun;
# 	any threshold R2 to pre-select models: R2>0.01
#            Compare how different our results from Wingo-NC_2022-suppl table 11 (28 genes).
# Table 1. 28 (AD1) genes in the Supp-Table 1, what is the best z score, p-value, significance, what is the best FUSION model z-score, p-value, significance (any threshold R2 to pre-select models). Wingo’s z-score, p-value, PP4. (post to github, send to the protem-AD slack)
library(dplyr)
library(data.table)
library(tidyr)
library(readr)

# chrom id  pos type    region_tag1 region_tag2 cs_index    susie_pip   mu2 
susie877df <- fread("data/877.susieIrss.txt") %>%
  filter(type == "gene") %>%
  select(all_of(c("chrom", "id", "pos", "cs_index", "susie_pip", "mu2")))
names(susie877df) <- paste0("susie.", names(susie877df))

pwas_chr_pos_key <- fread("data/pwas_chr_pos_key.txt") %>%
  separate(ID, c("ENSG ID", "Gene Symbol"), "\\.", extra="merge", remove=T) %>%
  rename(pos.p0 = P0) %>%
  select(all_of(c("ENSG ID", "Gene Symbol", "CHR", "pos.p0")))


spxcan_best_df <- fread("data/AD_GWAS_GWAX_meta_888_pqtl_wingo_nc_2022_best_spredixcan.csv") %>%
  select(all_of(c("gene", "zscore", "pvalue", "pred_perf_r2")))
names(spxcan_best_df) <- paste0("spxcan.best.", names(spxcan_best_df))

spxcan_lasso_df <- fread("data/AD_GWAS_GWAX_meta_888_pqtl_wingo_nc_2022_lasso_spredixcan.csv") %>%
  select(all_of(c("gene", "zscore", "pvalue", "pred_perf_r2")))
names(spxcan_lasso_df) <- paste0("spxcan.lasso.", names(spxcan_lasso_df))

num_lasso_non_na <- spxcan_lasso_df %>%
  filter(!is.na(spxcan.lasso.pvalue)) %>%
  nrow()
num_best_non_na <- spxcan_best_df %>%
  filter(!is.na(spxcan.best.pvalue)) %>%
  nrow()

final_col_order <- c("ENSG ID", "Gene Symbol", "CHR", "pos.p0",
  "susie.cs_index", "susie.pip", "susie.mu2",

  "spxcan.best.zscore", "spxcan.lasso.zscore",
  "spxcan.best.pvalue", "spxcan.lasso.pvalue",
  "spxcan.best.bfr.thresh", "spxcan.lasso.bfr.thresh",
  "spxcan.best.pred_perf_r2", "spxcan.lasso.pred_perf_r2")



pwas_write_table_df <- susie877df %>%
  left_join(pwas_chr_pos_key, by = c("susie.id" = "ENSG ID")) %>%
  left_join(spxcan_best_df, by = c("susie.id" = "spxcan.best.gene")) %>%
  left_join(spxcan_lasso_df, by = c("susie.id" = "spxcan.lasso.gene")) %>%
  rename(`ENSG ID` = "susie.id",
         `susie.pip` = "susie.susie_pip") %>%
  mutate(`CHR` = `susie.chrom`,
         `spxcan.lasso.bfr.thresh` = .05 / num_lasso_non_na,
         `spxcan.best.bfr.thresh` = .05 / num_best_non_na)%>%
  select(all_of(c(final_col_order))) %>%
  filter(susie.pip > 0.8)

browser()


write_delim(pwas_write_table_df, "data/susie_vs_pwas_spxcan_best_and_lasso.tsv", delim = "\t")

# susie877df <- fread("data/877.susieIrss.txt") %>%
#   filter(type == "gene") %>%
#   select(all_of(c("id", "cs_index", "susie_pip", "mu2")))
# names(susie877df) <- paste0("susie.", names(susie877df))
# 
# 
# wingo_s4_df <- fread("data/wingo_nc_2022_s4.csv") %>%
#   filter(Trait %in% c("Alzheimer's disease (AD1)")) %>%
#   dplyr::select(all_of(c("Ensembl gene ID", "Gene symbol", "PWAS Z score", "PWAS p-value", "COLOC PP4")))
# names(wingo_s4_df) <- paste0("wingo.", names(wingo_s4_df))
# 
# spxcan_lasso_df <- fread("data/AD_GWAS_GWAX_meta_888_brain_PWAS_lasso_model.csv") %>%
#   select(all_of(c("gene", "zscore", "pvalue", "pred_perf_r2")))
# names(spxcan_lasso_df) <- paste0("spxcan.lasso.", names(spxcan_lasso_df))
# 
# # spxcan_best_df <- fread("data/AD_GWAS_GWAX_meta_888_brain_PWAS_best_model.csv")
# spxcan_best_df <- fread("data/AD_GWAS_GWAX_meta_888_brain_PWAS_pretend_model.csv") %>%
#   select(all_of(c("gene", "zscore", "pvalue", "pred_perf_r2")))
# names(spxcan_best_df) <- paste0("spxcan.best.", names(spxcan_best_df))
# 
# pwas_chr_pos_key <- fread("data/pwas_chr_pos_key.txt") %>%
#   separate(ID, c("ENSG ID", "Gene Symbol"), "\\.", extra="merge", remove=T) %>%
#   rename(pos.p0 = P0) %>%
#   select(all_of(c("ENSG ID", "Gene Symbol", "CHR", "pos.p0")))
# 
# final_col_order <- c("ENSG ID", "Gene Symbol", "CHR", "pos.p0",
#   "wingo.PWAS Z score", "wingo.PWAS p-value", "wingo.COLOC PP4",
#   "spxcan.lasso.zscore", "spxcan.lasso.pvalue", "spxcan.lasso.pred_perf_r2",
#   "spxcan.best.zscore", "spxcan.best.pvalue", "spxcan.best.pred_perf_r2",
#   "susie.cs_index", "susie.susie_pip", "susie.mu2")
# 
# twas_write_table_df <- wingo_s4_df %>%
#   left_join(chr_pos_key, by = c("wingo.Ensembl gene ID" = "ENSG ID")) %>%
#   left_join(spxcan_lasso_df, by = c("wingo.Ensembl gene ID" = "spxcan.lasso.gene")) %>%
#   left_join(spxcan_best_df, by = c("wingo.Ensembl gene ID" = "spxcan.best.gene")) %>%
#   left_join(susie877df, by = c("wingo.Ensembl gene ID" = "susie.id")) %>%
#   rename(`ENSG ID` = "wingo.Ensembl gene ID",
#          `Gene Symbol` = "wingo.Gene symbol")
#   select(all_of(c(final_col_order)))
# 
# 
# 
# 
# write_delim(twas_write_table_df, "data/wingo_s4_vs_twas_spxcan_best_and_lasso_and_susie.tsv", delim = "\t")