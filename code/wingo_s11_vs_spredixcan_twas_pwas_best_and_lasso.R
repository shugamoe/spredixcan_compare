# **AD-PWAS by PrediXcan (using LASSO model and best model (both)) for Xiaotong Sun;
# 	any threshold R2 to pre-select models: R2>0.01
#            Compare how different our results from Wingo-NC_2022-suppl table 11 (28 genes).
# Table 1. 28 (AD1) genes in the Supp-Table 1, what is the best z score, p-value, significance, what is the best FUSION model z-score, p-value, significance (any threshold R2 to pre-select models). Wingo’s z-score, p-value, PP4. (post to github, send to the protem-AD slack)
library(dplyr)
library(data.table)
library(tidyr)
library(readr)

susie877df <- fread("data/877.susieIrss.txt") %>%
  filter(type == "gene") %>%
  select(all_of(c("id", "cs_index", "susie_pip", "mu2")))
names(susie877df) <- paste0("susie.", names(susie877df))


wingo_s11_df <- fread("data/wingo_nc_2022_s11.csv") %>%
  filter(Trait %in% c("Alzheimer's disease (AD1)")) %>%
  dplyr::select(all_of(c("Ensembl gene ID", "Gene symbol", "TWAS Z", "TWAS p-value", "COLOC PP4", "SMR SNP chr")))
names(wingo_s11_df) <- paste0("wingo.", names(wingo_s11_df))

spxcan_best_df <- fread("data/AD_GWAS_GWAX_meta_888_pqtl_wingo_nc_2022_best_spredixcan.csv") %>%
  select(all_of(c("gene", "zscore", "pvalue", "pred_perf_r2")))
names(spxcan_best_df) <- paste0("spxcan.best.", names(spxcan_best_df))

spxcan_lasso_df <- fread("data/AD_GWAS_GWAX_meta_888_pqtl_wingo_nc_2022_lasso_spredixcan.csv") %>%
  select(all_of(c("gene", "zscore", "pvalue", "pred_perf_r2")))
names(spxcan_lasso_df) <- paste0("spxcan.lasso.", names(spxcan_lasso_df))

pwas_chr_pos_key <- fread("data/pwas_chr_pos_key.txt") %>%
  separate(ID, c("ENSG ID", "Gene Symbol"), "\\.", extra="merge", remove=T) %>%
  rename(pos.p0 = P0) %>%
  select(all_of(c("ENSG ID", "Gene Symbol", "CHR", "pos.p0")))

final_col_order <- c("ENSG ID", "Gene Symbol", "CHR", "pos.p0",
  "wingo.TWAS Z", "wingo.TWAS p-value", "wingo.COLOC PP4",

  "spxcan.best.zscore", "spxcan.lasso.zscore",
  "spxcan.best.pvalue", "spxcan.lasso.pvalue",
  "spxcan.best.pred_perf_r2", "spxcan.lasso.pred_perf_r2",

  "susie.cs_index", "susie.susie_pip", "susie.mu2")

pwas_write_table_df <- wingo_s11_df %>%
  left_join(pwas_chr_pos_key, by = c("wingo.Ensembl gene ID" = "ENSG ID")) %>%
  left_join(spxcan_best_df, by = c("wingo.Ensembl gene ID" = "spxcan.best.gene")) %>%
  left_join(spxcan_lasso_df, by = c("wingo.Ensembl gene ID" = "spxcan.lasso.gene")) %>%
  left_join(susie877df, by = c("wingo.Ensembl gene ID" = "susie.id")) %>%
  rename(`ENSG ID` = "wingo.Ensembl gene ID") %>%
  mutate(`Gene Symbol` = ifelse(is.na(`Gene Symbol`), `wingo.Gene symbol`, `Gene Symbol`),
         `CHR` = ifelse(is.na(`CHR`), `wingo.SMR SNP chr`, `CHR`)) %>%
  select(all_of(c(final_col_order)))


write_delim(pwas_write_table_df, "data/wingo_s11_vs_pwas_spxcan_best_and_lasso_with_susie.tsv", delim = "\t")

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
