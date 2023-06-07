# **AD-PWAS by PrediXcan (using LASSO model and best model (both)) for Xiaotong Sun;
# 	any threshold R2 to pre-select models: R2>0.01
#            Compare how different our results from Wingo-NC_2022-suppl table 11 (28 genes).
# Table 1. 28 (AD1) genes in the Supp-Table 1, what is the best z score, p-value, significance, what is the best FUSION model z-score, p-value, significance (any threshold R2 to pre-select models). Wingo’s z-score, p-value, PP4. (post to github, send to the protem-AD slack)
library(dplyr)
library(data.table)
library(tidyr)
library(readr)

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

# Schartz
schwartz_gwasx_best_df <- fread("data/AD_GWAS_GWAX_meta_888_pqtl_wingo_nc_2022_best_spredixcan.csv") %>%
  select(all_of(c("gene", "zscore", "pvalue", "pred_perf_r2")))

num_best_non_na <- schwartz_gwasx_best_df %>%
  filter(!is.na(zscore)) %>%
  nrow()

best_thresh <- .05 / num_best_non_na

schwartz_gwasx_lasso_df <- fread("data/AD_GWAS_GWAX_meta_888_pqtl_wingo_nc_2022_lasso_spredixcan.csv") %>%
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

# Prep output
pwas_chr_pos_key <- fread("data/pwas_chr_pos_key.txt") %>%
  separate(ID, c("ENSG ID", "Gene Symbol"), "\\.", extra="merge", remove=T) %>%
  rename(pos.p0 = P0) %>%
  select(all_of(c("ENSG ID", "Gene Symbol", "CHR", "pos.p0")))

final_col_order <- c("ENSG ID", "Gene Symbol", "CHR", "pos.p0",
  "wingo.TWAS Z", "wingo.TWAS p-value", "wingo.COLOC PP4",

  # Schwartz
  "schwartz_gwasx.best.zscore", "schwartz_gwasx.lasso.zscore",
  "schwartz_gwasx.best.pvalue", "schwartz_gwasx.lasso.pvalue",
  "schwartz_gwasx.best.pred_perf_r2", "schwartz_gwasx.lasso.pred_perf_r2",
  "schwartz_gwasx.best.bfr_thresh",
  "schwartz_gwasx.lasso.bfr_thresh",
  "schwartz_gwasx.best.bfr_sig",
  "schwartz_gwasx.lasso.bfr_sig",

  # Bellenguez
  "bellen.best.zscore", "bellen.lasso.zscore",
  "bellen.best.pvalue", "bellen.lasso.pvalue",
  "bellen.best.pred_perf_r2", "bellen.lasso.pred_perf_r2",
  "bellen.best.bfr_thresh",
  "bellen.lasso.bfr_thresh",
  "bellen.best.bfr_sig",
  "bellen.lasso.bfr_sig",

  # Wightmen 
  "wight.best.zscore", "wight.lasso.zscore",
  "wight.best.pvalue", "wight.lasso.pvalue",
  "wight.best.pred_perf_r2", "wight.lasso.pred_perf_r2",
  "wight.best.bfr_thresh",
  "wight.lasso.bfr_thresh",
  "wight.best.bfr_sig",
  "wight.lasso.bfr_sig",

  "susie.cs_index", "susie.susie_pip", "susie.mu2")

all_final_col_order <- c("ENSG ID", "Gene Symbol", "CHR", "pos.p0",
  "wingo.TWAS Z", "wingo.TWAS p-value", "wingo.COLOC PP4",

  # Schwartz (GWAS only)
  "schwartz_gwas.best.zscore", "schwartz_gwas.lasso.zscore",
  "schwartz_gwas.best.pvalue", "schwartz_gwas.lasso.pvalue",
  "schwartz_gwas.best.pred_perf_r2", "schwartz_gwas.lasso.pred_perf_r2",
  "schwartz_gwas.best.bfr_thresh",
  "schwartz_gwas.lasso.bfr_thresh",
  "schwartz_gwas.best.bfr_sig",
  "schwartz_gwas.lasso.bfr_sig",

  # Schwartz
  "schwartz_gwasx.best.zscore", "schwartz_gwasx.lasso.zscore",
  "schwartz_gwasx.best.pvalue", "schwartz_gwasx.lasso.pvalue",
  "schwartz_gwasx.best.pred_perf_r2", "schwartz_gwasx.lasso.pred_perf_r2",
  "schwartz_gwasx.best.bfr_thresh",
  "schwartz_gwasx.lasso.bfr_thresh",
  "schwartz_gwasx.best.bfr_sig",
  "schwartz_gwasx.lasso.bfr_sig",

  # Bellenguez
  "bellen.best.zscore", "bellen.lasso.zscore",
  "bellen.best.pvalue", "bellen.lasso.pvalue",
  "bellen.best.pred_perf_r2", "bellen.lasso.pred_perf_r2",
  "bellen.best.bfr_thresh",
  "bellen.lasso.bfr_thresh",
  "bellen.best.bfr_sig",
  "bellen.lasso.bfr_sig",

  # Wightmen 
  "wight.best.zscore", "wight.lasso.zscore",
  "wight.best.pvalue", "wight.lasso.pvalue",
  "wight.best.pred_perf_r2", "wight.lasso.pred_perf_r2",
  "wight.best.bfr_thresh",
  "wight.lasso.bfr_thresh",
  "wight.best.bfr_sig",
  "wight.lasso.bfr_sig",

  # Jansen
  "jansen.best.zscore", "jansen.lasso.zscore",
  "jansen.best.pvalue", "jansen.lasso.pvalue",
  "jansen.best.pred_perf_r2", "jansen.lasso.pred_perf_r2",
  "jansen.best.bfr_thresh",
  "jansen.lasso.bfr_thresh",
  "jansen.best.bfr_sig",
  "jansen.lasso.bfr_sig",

  "susie.cs_index", "susie.susie_pip", "susie.mu2")


# Schwartz best (GWAS only)
pwas_best_sig_write_table_df <- pwas_chr_pos_key %>%
  inner_join(schwartz_gwas_best_df %>% filter(schwartz_gwas.best.bfr_sig == TRUE), by = c("ENSG ID" = "schwartz_gwas.best.gene")) %>%
  select(any_of(c(final_col_order)))

pwas_lasso_sig_write_table_df <- pwas_chr_pos_key %>%
  inner_join(schwartz_gwas_lasso_df %>% filter(schwartz_gwas.lasso.bfr_sig == TRUE), by = c("ENSG ID" = "schwartz_gwas.lasso.gene")) %>%
  select(any_of(c(final_col_order)))

write_delim(pwas_best_sig_write_table_df, "data/sig_schwartz_gwas_spxcan_best.tsv", delim = "\t")
write_delim(pwas_lasso_sig_write_table_df, "data/sig_schwartz_gwas_spxcan_lasso.tsv", delim = "\t")

# Schwartz best
pwas_best_sig_write_table_df <- pwas_chr_pos_key %>%
  inner_join(schwartz_gwasx_best_df %>% filter(schwartz_gwasx.best.bfr_sig == TRUE), by = c("ENSG ID" = "schwartz_gwasx.best.gene")) %>%
  select(any_of(c(final_col_order)))

pwas_lasso_sig_write_table_df <- pwas_chr_pos_key %>%
  inner_join(schwartz_gwasx_lasso_df %>% filter(schwartz_gwasx.lasso.bfr_sig == TRUE), by = c("ENSG ID" = "schwartz_gwasx.lasso.gene")) %>%
  select(any_of(c(final_col_order)))

write_delim(pwas_best_sig_write_table_df, "data/sig_schwartz_gwasx_spxcan_best.tsv", delim = "\t")
write_delim(pwas_lasso_sig_write_table_df, "data/sig_schwartz_gwasx_spxcan_lasso.tsv", delim = "\t")

# Bellen best
bellen_best_sig_write_table_df <- pwas_chr_pos_key %>%
  inner_join(bellen_best_df %>% filter(bellen.best.bfr_sig == TRUE), by = c("ENSG ID" = "bellen.best.gene")) %>%
  select(any_of(c(final_col_order)))

bellen_lasso_sig_write_table_df <- pwas_chr_pos_key %>%
  inner_join(bellen_lasso_df %>% filter(bellen.lasso.bfr_sig == TRUE), by = c("ENSG ID" = "bellen.lasso.gene")) %>%
  select(any_of(c(final_col_order)))

write_delim(bellen_best_sig_write_table_df, "data/sig_bellen_spxcan_best.tsv", delim = "\t")
write_delim(bellen_lasso_sig_write_table_df, "data/sig_bellen_spxcan_lasso.tsv", delim = "\t")

# Wight Best
wight_best_sig_write_table_df <- pwas_chr_pos_key %>%
  inner_join(wight_best_df %>% filter(wight.best.bfr_sig == TRUE), by = c("ENSG ID" = "wight.best.gene")) %>%
  select(any_of(c(final_col_order)))

wight_lasso_sig_write_table_df <- pwas_chr_pos_key %>%
  inner_join(wight_lasso_df %>% filter(wight.lasso.bfr_sig == TRUE), by = c("ENSG ID" = "wight.lasso.gene")) %>%
  select(any_of(c(final_col_order)))

write_delim(wight_best_sig_write_table_df, "data/sig_wight_spxcan_best.tsv", delim = "\t")
write_delim(wight_lasso_sig_write_table_df, "data/sig_wight_spxcan_lasso.tsv", delim = "\t")

# Jansen Best
jansen_best_sig_write_table_df <- pwas_chr_pos_key %>%
  inner_join(jansen_best_df %>% filter(jansen.best.bfr_sig == TRUE), by = c("ENSG ID" = "jansen.best.gene")) %>%
  select(any_of(c(final_col_order)))

jansen_lasso_sig_write_table_df <- pwas_chr_pos_key %>%
  inner_join(jansen_lasso_df %>% filter(jansen.lasso.bfr_sig == TRUE), by = c("ENSG ID" = "jansen.lasso.gene")) %>%
  select(any_of(c(final_col_order)))

write_delim(jansen_best_sig_write_table_df, "data/sig_jansen_spxcan_best.tsv", delim = "\t")
write_delim(jansen_lasso_sig_write_table_df, "data/sig_jansen_spxcan_lasso.tsv", delim = "\t")

all_pwas_lasso_sig_write_table_df <- pwas_chr_pos_key %>%
  full_join(schwartz_gwasx_lasso_df, by = c("ENSG ID" = "schwartz_gwasx.lasso.gene")) %>%
  full_join(schwartz_gwasx_best_df, by = c("ENSG ID" = "schwartz_gwasx.best.gene")) %>%
  full_join(schwartz_gwas_lasso_df, by = c("ENSG ID" = "schwartz_gwas.lasso.gene")) %>%
  full_join(schwartz_gwas_best_df, by = c("ENSG ID" = "schwartz_gwas.best.gene")) %>%
  full_join(bellen_best_df, by = c("ENSG ID" = "bellen.best.gene")) %>%
  full_join(bellen_lasso_df, by = c("ENSG ID" = "bellen.lasso.gene")) %>%
  full_join(wight_best_df, by = c("ENSG ID" = "wight.best.gene")) %>%
  full_join(wight_lasso_df, by = c("ENSG ID" = "wight.lasso.gene")) %>%
  full_join(jansen_best_df, by = c("ENSG ID" = "jansen.best.gene")) %>%
  full_join(jansen_lasso_df, by = c("ENSG ID" = "jansen.lasso.gene")) %>%
  select(any_of(c(all_final_col_order)))

write_delim(all_pwas_lasso_sig_write_table_df, "data/all_pwas_spxcan_lasso_vs_best_models.tsv", delim = "\t")
