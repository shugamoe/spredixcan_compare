# **AD-PWAS by PrediXcan (using LASSO model and best model (both)) for Xiaotong Sun;
# 	any threshold R2 to pre-select models: R2>0.01
#            Compare how different our results from Wingo-NC_2022-suppl table 11 (28 genes).
# Table 1. 28 (AD1) genes in the Supp-Table 1, what is the best z score, p-value, significance, what is the best FUSION model z-score, p-value, significance (any threshold R2 to pre-select models). Wingo’s z-score, p-value, PP4. (post to github, send to the protem-AD slack)
library(dplyr)
library(data.table)
library(tidyr)
library(readr)

spxcan_best_df <- fread("data/AD_GWAS_GWAX_meta_888_pqtl_wingo_nc_2022_best_spredixcan.csv") %>%
  select(all_of(c("gene", "zscore", "pvalue", "pred_perf_r2")))

num_best_non_na <- spxcan_best_df %>%
  filter(!is.na(zscore)) %>%
  nrow()

best_thresh <- .05 / num_best_non_na

spxcan_lasso_df <- fread("data/AD_GWAS_GWAX_meta_888_pqtl_wingo_nc_2022_lasso_spredixcan.csv") %>%
  select(all_of(c("gene", "zscore", "pvalue", "pred_perf_r2")))

num_lasso_non_na <- spxcan_lasso_df %>%
  filter(!is.na(zscore)) %>%
  nrow()

lasso_thresh <- .05 / num_lasso_non_na

spxcan_lasso_df <- spxcan_lasso_df %>%
  mutate(bfr_sig = ifelse(pvalue < lasso_thresh, TRUE, FALSE),
         bfr_thresh = lasso_thresh,
         )

spxcan_best_df <- spxcan_best_df %>%
mutate(bfr_sig = ifelse(pvalue < best_thresh, TRUE, FALSE),
       bfr_thresh = best_thresh)

pwas_chr_pos_key <- fread("data/pwas_chr_pos_key.txt") %>%
  separate(ID, c("ENSG ID", "Gene Symbol"), "\\.", extra="merge", remove=T) %>%
  rename(pos.p0 = P0) %>%
  select(all_of(c("ENSG ID", "Gene Symbol", "CHR", "pos.p0")))

final_col_order <- c("ENSG ID", "Gene Symbol", "CHR", "pos.p0",
  "wingo.TWAS Z", "wingo.TWAS p-value", "wingo.COLOC PP4",

  "zscore", 
  "pvalue", 
  "pred_perf_r2", 
  "bfr_thresh",

  "susie.cs_index", "susie.susie_pip", "susie.mu2")

pwas_best_sig_write_table_df <- pwas_chr_pos_key %>%
  inner_join(spxcan_best_df %>% filter(bfr_sig == TRUE), by = c("ENSG ID" = "gene")) %>%
  select(any_of(c(final_col_order)))

pwas_lasso_sig_write_table_df <- pwas_chr_pos_key %>%
  inner_join(spxcan_lasso_df %>% filter(bfr_sig == TRUE), by = c("ENSG ID" = "gene")) %>%
  select(any_of(c(final_col_order)))

browser()


write_delim(pwas_best_sig_write_table_df, "data/sig_pwas_spxcan_best.tsv", delim = "\t")
write_delim(pwas_lasso_sig_write_table_df, "data/sig_pwas_spxcan_lasso.tsv", delim = "\t")

# susie877df <- fread("data/877.susieIrss.txt") %>%
#   filter(type == "gene") %>%
#   select(all_of(c("id", "cs_index", "susie_pip", "mu2")))
# names(susie877df) <- paste0("susie.", names(susie877df))
# 
# 
# wingo_s11_df <- fread("data/wingo_nc_2022_s11.csv") %>%
#   filter(Trait %in% c("Alzheimer's disease (AD1)")) %>%
#   dplyr::select(all_of(c("Ensembl gene ID", "Gene symbol", "TWAS Z", "TWAS p-value", "COLOC PP4")))
# names(wingo_s11_df) <- paste0("wingo.", names(wingo_s11_df))
# 
# spxcan_lasso_df <- fread("data/AD_GWAS_GWAX_meta_888_brain_TWAS_lasso_model.csv") %>%
#   select(all_of(c("gene", "zscore", "pvalue", "pred_perf_r2")))
# names(spxcan_lasso_df) <- paste0("spxcan.lasso.", names(spxcan_lasso_df))
# 
# # spxcan_best_df <- fread("data/AD_GWAS_GWAX_meta_888_brain_TWAS_best_model.csv")
# spxcan_best_df <- fread("data/AD_GWAS_GWAX_meta_888_brain_TWAS_pretend_model.csv") %>%
#   select(all_of(c("gene", "zscore", "pvalue", "pred_perf_r2")))
# names(spxcan_best_df) <- paste0("spxcan.best.", names(spxcan_best_df))
# 
# pwas_chr_pos_key <- fread("data/pwas_chr_pos_key.txt") %>%
#   separate(ID, c("ENSG ID", "Gene Symbol"), "\\.", extra="merge", remove=T) %>%
#   rename(pos.p0 = P0) %>%
#   select(all_of(c("ENSG ID", "Gene Symbol", "CHR", "pos.p0")))
# 
# final_col_order <- c("ENSG ID", "Gene Symbol", "CHR", "pos.p0",
#   "wingo.TWAS Z", "wingo.TWAS p-value", "wingo.COLOC PP4",
#   "spxcan.lasso.zscore", "spxcan.lasso.pvalue", "spxcan.lasso.pred_perf_r2",
#   "spxcan.best.zscore", "spxcan.best.pvalue", "spxcan.best.pred_perf_r2",
#   "susie.cs_index", "susie.susie_pip", "susie.mu2")
# 
# twas_write_table_df <- wingo_s11_df %>%
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
# write_delim(twas_write_table_df, "data/wingo_s11_vs_twas_spxcan_best_and_lasso_and_susie.tsv", delim = "\t")
