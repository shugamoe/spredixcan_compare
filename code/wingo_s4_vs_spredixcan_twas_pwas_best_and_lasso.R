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


wingo_s4_a1_df <- fread("data/wingo_nc_2022_s4.csv") %>%
  filter(Trait %in% c("Alzheimer's disease (AD1)")) %>%
  dplyr::select(all_of(c("Ensembl gene ID", "Gene symbol", "PWAS Z score", "PWAS p-value", "COLOC PP4", "SMR SNP chr")))
names(wingo_s4_a1_df) <- paste0("wingo.", names(wingo_s4_a1_df))

# AD2?
wingo_s4_a2_df <- fread("data/wingo_nc_2022_s4.csv") %>%
  filter(Trait %in% c("Alzheimer's disease (AD2)")) %>%
  dplyr::select(all_of(c("Ensembl gene ID", "Gene symbol", "PWAS Z score", "PWAS p-value", "COLOC PP4", "SMR SNP chr")))
names(wingo_s4_a2_df) <- paste0("wingo.", names(wingo_s4_a2_df))

# Schartz
schwartz_gwasx_best_df <- fread("data/spxcan_pwas_output/AD_GWAS_GWAX_meta_888_pqtl_wingo_nc_2022_best_spredixcan.csv") %>%
  select(all_of(c("gene", "zscore", "pvalue", "pred_perf_r2")))

num_best_non_na <- schwartz_gwasx_best_df %>%
  filter(!is.na(zscore)) %>%
  nrow()

best_thresh <- .05 / num_best_non_na

schwartz_gwasx_lasso_df <- fread("data/spxcan_pwas_output/AD_GWAS_GWAX_meta_888_pqtl_wingo_nc_2022_lasso_spredixcan.csv") %>%
  select(all_of(c("gene", "zscore", "pvalue", "pred_perf_r2")))

num_lasso_non_na <- schwartz_gwasx_lasso_df %>%
  filter(!is.na(zscore)) %>%
  nrow()

lasso_thresh <- .05 / num_lasso_non_na

schwartz_gwasx_lasso_df <- schwartz_gwasx_lasso_df %>%
  mutate(bfr_sig = ifelse(pvalue < lasso_thresh, TRUE, FALSE),
         bfr_thresh = lasso_thresh,
         )

schwartz_gwasx_best_df <- schwartz_gwasx_best_df %>%
mutate(bfr_sig = ifelse(pvalue < best_thresh, TRUE, FALSE),
       bfr_thresh = best_thresh)

names(schwartz_gwasx_best_df) <- paste0("schwartz_gwasx.best.", names(schwartz_gwasx_best_df))
names(schwartz_gwasx_lasso_df) <- paste0("schwartz_gwasx.lasso.", names(schwartz_gwasx_lasso_df))

# Schartz (GWAS only)
schwartz_gwas_best_df <- fread("data/spxcan_pwas_output/AD_GWAS_only_888_pqtl_wingo_nc_2022_best_spredixcan.csv") %>%
  select(all_of(c("gene", "zscore", "pvalue", "pred_perf_r2")))

num_best_non_na <- schwartz_gwas_best_df %>%
  filter(!is.na(zscore)) %>%
  nrow()

best_thresh <- .05 / num_best_non_na

schwartz_gwas_lasso_df <- fread("data/spxcan_pwas_output/AD_GWAS_only_888_pqtl_wingo_nc_2022_lasso_spredixcan.csv") %>%
  select(all_of(c("gene", "zscore", "pvalue", "pred_perf_r2")))

num_lasso_non_na <- schwartz_gwas_lasso_df %>%
  filter(!is.na(zscore)) %>%
  nrow()

lasso_thresh <- .05 / num_lasso_non_na

schwartz_gwas_lasso_df <- schwartz_gwas_lasso_df %>%
  mutate(bfr_sig = ifelse(pvalue < lasso_thresh, TRUE, FALSE),
         bfr_thresh = lasso_thresh,
         )

schwartz_gwas_best_df <- schwartz_gwas_best_df %>%
mutate(bfr_sig = ifelse(pvalue < best_thresh, TRUE, FALSE),
       bfr_thresh = best_thresh)

names(schwartz_gwas_best_df) <- paste0("schwartz_gwas.best.", names(schwartz_gwas_best_df))
names(schwartz_gwas_lasso_df) <- paste0("schwartz_gwas.lasso.", names(schwartz_gwas_lasso_df))

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

# Jansen
jansen_best_df <- fread("data/spxcan_pwas_output/Jansen_meta_GWAS_2019_NG_888_pqtl_wingo_nc_2022_best_spredixcan.csv") %>%
  select(all_of(c("gene", "zscore", "pvalue", "pred_perf_r2")))

num_best_non_na <- jansen_best_df %>%
  filter(!is.na(zscore)) %>%
  nrow()

best_thresh <- .05 / num_best_non_na

jansen_lasso_df <- fread("data/spxcan_pwas_output/Jansen_meta_GWAS_2019_NG_888_pqtl_wingo_nc_2022_lasso_spredixcan.csv") %>%
  select(all_of(c("gene", "zscore", "pvalue", "pred_perf_r2")))

num_lasso_non_na <- jansen_lasso_df %>%
  filter(!is.na(zscore)) %>%
  nrow()

lasso_thresh <- .05 / num_lasso_non_na

jansen_lasso_df <- jansen_lasso_df %>%
  mutate(bfr_sig = ifelse(pvalue < lasso_thresh, TRUE, FALSE),
         bfr_thresh = lasso_thresh,
         )

jansen_best_df <- jansen_best_df %>%
mutate(bfr_sig = ifelse(pvalue < best_thresh, TRUE, FALSE),
       bfr_thresh = best_thresh)

names(jansen_best_df) <- paste0("jansen.best.", names(jansen_best_df))
names(jansen_lasso_df) <- paste0("jansen.lasso.", names(jansen_lasso_df))

pwas_chr_pos_key <- fread("data/pwas_chr_pos_key.txt") %>%
  separate(ID, c("ENSG ID", "Gene Symbol"), "\\.", extra="merge", remove=T) %>%
  rename(pos.p0 = P0) %>%
  select(all_of(c("ENSG ID", "Gene Symbol", "CHR", "pos.p0")))

final_col_order <- c("ENSG ID", "Gene Symbol", "CHR", "pos.p0",
  "wingo.PWAS Z score", "wingo.PWAS p-value", "wingo.COLOC PP4",

  "schwartz_gwas.best.zscore", "schwartz_gwas.lasso.zscore",
  "schwartz_gwas.best.pvalue", "schwartz_gwas.lasso.pvalue",
  "schwartz_gwas.best.pred_perf_r2", "schwartz_gwas.lasso.pred_perf_r2",

  "schwartz_gwasx.best.zscore", "schwartz_gwasx.lasso.zscore",
  "schwartz_gwasx.best.pvalue", "schwartz_gwasx.lasso.pvalue",
  "schwartz_gwasx.best.pred_perf_r2", "schwartz_gwasx.lasso.pred_perf_r2",

  "bellen.best.zscore", "bellen.lasso.zscore",
  "bellen.best.pvalue", "bellen.lasso.pvalue",
  "bellen.best.pred_perf_r2", "bellen.lasso.pred_perf_r2",

  "wight.best.zscore", "wight.lasso.zscore",
  "wight.best.pvalue", "wight.lasso.pvalue",
  "wight.best.pred_perf_r2", "wight.lasso.pred_perf_r2",

  "jansen.best.zscore", "jansen.lasso.zscore",
  "jansen.best.pvalue", "jansen.lasso.pvalue",
  "jansen.best.pred_perf_r2", "jansen.lasso.pred_perf_r2",

  "susie.cs_index", "susie.susie_pip", "susie.mu2")

a1_pwas_write_table_df <- wingo_s4_a1_df %>%
  left_join(pwas_chr_pos_key, by = c("wingo.Ensembl gene ID" = "ENSG ID")) %>%
  left_join(schwartz_gwasx_best_df, by = c("wingo.Ensembl gene ID" = "schwartz_gwasx.best.gene")) %>%
  left_join(schwartz_gwasx_lasso_df, by = c("wingo.Ensembl gene ID" = "schwartz_gwasx.lasso.gene")) %>%
  left_join(schwartz_gwas_best_df, by = c("wingo.Ensembl gene ID" = "schwartz_gwas.best.gene")) %>%
  left_join(schwartz_gwas_lasso_df, by = c("wingo.Ensembl gene ID" = "schwartz_gwas.lasso.gene")) %>%
  left_join(bellen_best_df, by = c("wingo.Ensembl gene ID" = "bellen.best.gene")) %>%
  left_join(bellen_lasso_df, by = c("wingo.Ensembl gene ID" = "bellen.lasso.gene")) %>%
  left_join(wight_best_df, by = c("wingo.Ensembl gene ID" = "wight.best.gene")) %>%
  left_join(wight_lasso_df, by = c("wingo.Ensembl gene ID" = "wight.lasso.gene")) %>%
  left_join(jansen_best_df, by = c("wingo.Ensembl gene ID" = "jansen.best.gene")) %>%
  left_join(jansen_lasso_df, by = c("wingo.Ensembl gene ID" = "jansen.lasso.gene")) %>%
  left_join(susie877df, by = c("wingo.Ensembl gene ID" = "susie.id")) %>%
  rename(`ENSG ID` = "wingo.Ensembl gene ID") %>%
  mutate(`Gene Symbol` = ifelse(is.na(`Gene Symbol`), `wingo.Gene symbol`, `Gene Symbol`),
         `CHR` = ifelse(is.na(`CHR`), `wingo.SMR SNP chr`, `CHR`)) %>%
  select(all_of(c(final_col_order)))

a2_pwas_write_table_df <- wingo_s4_a2_df %>%
  left_join(pwas_chr_pos_key, by = c("wingo.Ensembl gene ID" = "ENSG ID")) %>%
  left_join(schwartz_gwasx_best_df, by = c("wingo.Ensembl gene ID" = "schwartz_gwasx.best.gene")) %>%
  left_join(schwartz_gwasx_lasso_df, by = c("wingo.Ensembl gene ID" = "schwartz_gwasx.lasso.gene")) %>%
  left_join(schwartz_gwas_best_df, by = c("wingo.Ensembl gene ID" = "schwartz_gwas.best.gene")) %>%
  left_join(schwartz_gwas_lasso_df, by = c("wingo.Ensembl gene ID" = "schwartz_gwas.lasso.gene")) %>%
  left_join(bellen_best_df, by = c("wingo.Ensembl gene ID" = "bellen.best.gene")) %>%
  left_join(bellen_lasso_df, by = c("wingo.Ensembl gene ID" = "bellen.lasso.gene")) %>%
  left_join(wight_best_df, by = c("wingo.Ensembl gene ID" = "wight.best.gene")) %>%
  left_join(wight_lasso_df, by = c("wingo.Ensembl gene ID" = "wight.lasso.gene")) %>%
  left_join(jansen_best_df, by = c("wingo.Ensembl gene ID" = "jansen.best.gene")) %>%
  left_join(jansen_lasso_df, by = c("wingo.Ensembl gene ID" = "jansen.lasso.gene")) %>%
  left_join(susie877df, by = c("wingo.Ensembl gene ID" = "susie.id")) %>%
  rename(`ENSG ID` = "wingo.Ensembl gene ID") %>%
  mutate(`Gene Symbol` = ifelse(is.na(`Gene Symbol`), `wingo.Gene symbol`, `Gene Symbol`),
         `CHR` = ifelse(is.na(`CHR`), `wingo.SMR SNP chr`, `CHR`)) %>%
  select(all_of(c(final_col_order)))


write_delim(a1_pwas_write_table_df, "data/wingo_s4_ad1_vs_pwas_spxcan_best_and_lasso_with_susie.tsv", delim = "\t")
write_delim(a2_pwas_write_table_df, "data/wingo_s4_ad2_vs_pwas_spxcan_best_and_lasso_with_susie.tsv", delim = "\t")
