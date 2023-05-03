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


# Schartz
schwartz_best_df <- fread("data/AD_GWAS_GWAX_meta_888_pqtl_wingo_nc_2022_best_spredixcan.csv") %>%
  select(all_of(c("gene", "zscore", "pvalue", "pred_perf_r2")))

num_best_non_na <- schwartz_best_df %>%
  filter(!is.na(zscore)) %>%
  nrow()

best_thresh <- .05 / num_best_non_na

schwartz_lasso_df <- fread("data/AD_GWAS_GWAX_meta_888_pqtl_wingo_nc_2022_lasso_spredixcan.csv") %>%
  select(all_of(c("gene", "zscore", "pvalue", "pred_perf_r2")))

num_lasso_non_na <- schwartz_lasso_df %>%
  filter(!is.na(zscore)) %>%
  nrow()

lasso_thresh <- .05 / num_lasso_non_na

schwartz_lasso_df <- schwartz_lasso_df %>%
  mutate(bfr_sig = ifelse(pvalue < lasso_thresh, TRUE, FALSE),
         bfr_thresh = lasso_thresh,
         )

schwartz_best_df <- schwartz_best_df %>%
mutate(bfr_sig = ifelse(pvalue < best_thresh, TRUE, FALSE),
       bfr_thresh = best_thresh)

names(schwartz_best_df) <- paste0("schwartz.best.", names(schwartz_best_df))
names(schwartz_lasso_df) <- paste0("schwartz.lasso.", names(schwartz_lasso_df))

# Wightmen
wight_best_df <- fread("data/spxcan_pwas_output/AD_Wightmen_NG_2021_888_pqtl_wingo_nc_2022_best_spredixcan.csv") %>%
  select(all_of(c("gene", "zscore", "pvalue", "pred_perf_r2")))

num_best_non_na <- wight_best_df %>%
  filter(!is.na(zscore)) %>%
  nrow()

best_thresh <- .05 / num_best_non_na

wight_lasso_df <- fread("data/spxcan_pwas_output/AD_Wightmen_NG_2021_888_pqtl_wingo_nc_2022_lasso_spredixcan.csv") %>%
  select(all_of(c("gene", "zscore", "pvalue", "pred_perf_r2")))

num_lasso_non_na <- wight_lasso_df %>%
  filter(!is.na(zscore)) %>%
  nrow()

lasso_thresh <- .05 / num_lasso_non_na

wight_lasso_df <- wight_lasso_df %>%
  mutate(bfr_sig = ifelse(pvalue < lasso_thresh, TRUE, FALSE),
         bfr_thresh = lasso_thresh,
         )

wight_best_df <- wight_best_df %>%
mutate(bfr_sig = ifelse(pvalue < best_thresh, TRUE, FALSE),
       bfr_thresh = best_thresh)

names(wight_best_df) <- paste0("wight.best.", names(wight_best_df))
names(wight_lasso_df) <- paste0("wight.lasso.", names(wight_lasso_df))

# Bellenguez
bellen_best_df <- fread("data/spxcan_pwas_output/AD_Bellenguez_GWAS_NG_2022_888_pqtl_wingo_nc_2022_best_spredixcan.csv") %>%
  select(all_of(c("gene", "zscore", "pvalue", "pred_perf_r2")))

num_best_non_na <- bellen_best_df %>%
  filter(!is.na(zscore)) %>%
  nrow()

best_thresh <- .05 / num_best_non_na

bellen_lasso_df <- fread("data/spxcan_pwas_output/AD_Bellenguez_GWAS_NG_2022_888_pqtl_wingo_nc_2022_lasso_spredixcan.csv") %>%
  select(all_of(c("gene", "zscore", "pvalue", "pred_perf_r2")))

num_lasso_non_na <- bellen_lasso_df %>%
  filter(!is.na(zscore)) %>%
  nrow()

lasso_thresh <- .05 / num_lasso_non_na

bellen_lasso_df <- bellen_lasso_df %>%
  mutate(bfr_sig = ifelse(pvalue < lasso_thresh, TRUE, FALSE),
         bfr_thresh = lasso_thresh,
         )

bellen_best_df <- bellen_best_df %>%
mutate(bfr_sig = ifelse(pvalue < best_thresh, TRUE, FALSE),
       bfr_thresh = best_thresh)

names(bellen_best_df) <- paste0("bellen.best.", names(bellen_best_df))
names(bellen_lasso_df) <- paste0("bellen.lasso.", names(bellen_lasso_df))


final_col_order <- c("ENSG ID", "Gene Symbol", "CHR", "pos.p0",
  "susie.cs_index", "susie.pip", "susie.mu2",

  "schwartz.best.zscore", "schwartz.lasso.zscore",
  "schwartz.best.pvalue", "schwartz.lasso.pvalue",
  "schwartz.best.bfr_thresh", "schwartz.lasso.bfr_thresh",
  "schwartz.best.pred_perf_r2", "schwartz.lasso.pred_perf_r2",

  "bellen.best.zscore", "bellen.lasso.zscore",
  "bellen.best.pvalue", "bellen.lasso.pvalue",
  "bellen.best.bfr_thresh", "bellen.lasso.bfr_thresh",
  "bellen.best.pred_perf_r2", "bellen.lasso.pred_perf_r2",

  "wight.best.zscore", "wight.lasso.zscore",
  "wight.best.pvalue", "wight.lasso.pvalue",
  "wight.best.bfr_thresh", "wight.lasso.bfr_thresh",
  "wight.best.pred_perf_r2", "wight.lasso.pred_perf_r2"
)



pwas_write_table_df <- susie877df %>%
  left_join(pwas_chr_pos_key, by = c("susie.id" = "ENSG ID")) %>%
  left_join(schwartz_best_df, by = c("susie.id" = "schwartz.best.gene")) %>%
  left_join(schwartz_lasso_df, by = c("susie.id" = "schwartz.lasso.gene")) %>%
  left_join(bellen_best_df, by = c("susie.id" = "bellen.best.gene")) %>%
  left_join(bellen_lasso_df, by = c("susie.id" = "bellen.lasso.gene")) %>%
  left_join(wight_best_df, by = c("susie.id" = "wight.best.gene")) %>%
  left_join(wight_lasso_df, by = c("susie.id" = "wight.lasso.gene")) %>%
  rename(`ENSG ID` = "susie.id",
         `susie.pip` = "susie.susie_pip") %>%
  mutate(`CHR` = `susie.chrom`) %>%
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
