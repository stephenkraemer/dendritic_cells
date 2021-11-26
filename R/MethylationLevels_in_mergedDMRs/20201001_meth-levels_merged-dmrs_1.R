# Methylation levels in merged DC DMRs ------------------------------------------------------
# Dr. Sina Stäble
# 01.10.2020

# load libraries
library(tidyverse)
#library(bsseq)
library(GenomicRanges)
library(gUtils)

# load merged DMRs ------------------------------------------------------

merged_dmrs.gr <- readRDS("/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/DMR calling/merged_dmrs/2020-10-01_merged-dmrs.rds")


# load methylation calls ------------------------------------------------------

input.path <- "/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/meth_calls/"


# load methylation calls of HSCs

hsc_1 <- read_delim(file = paste0(input.path, "mcalls_hsc_1_CG_chrom-merged_strands-merged.bed.gz"), 
                    "\t", escape_double = FALSE, trim_ws = TRUE)
hsc_2 <- read_delim(file = paste0(input.path, "mcalls_hsc_2_CG_chrom-merged_strands-merged.bed.gz"), 
                    "\t", escape_double = FALSE, trim_ws = TRUE)
hsc_3 <- read_delim(file = paste0(input.path, "mcalls_hsc_3_CG_chrom-merged_strands-merged.bed.gz"), 
                    "\t", escape_double = FALSE, trim_ws = TRUE)

hsc <- list(hsc_1, hsc_2, hsc_3) %>%
  map(~ select(., -motif, -score, -strand, -beta_value)) %>%
  purrr::reduce(., full_join, by = c("#chrom", "start", "end")) %>%
  mutate(., total = rowSums(select(., contains("n_total")))) %>%
  mutate(., meth = rowSums(select(., contains("n_meth")))) %>%
  select(., `#chrom`, start, end, total, meth)


# load methylation calls of CDPs

cdp_1 <- read_delim(file = paste0(input.path, "mcalls_cdp_1_CG_chrom-merged_strands-merged.bed.gz"), 
                    "\t", escape_double = FALSE, trim_ws = TRUE)
cdp_2 <- read_delim(file = paste0(input.path, "mcalls_cdp_2_CG_chrom-merged_strands-merged.bed.gz"), 
                    "\t", escape_double = FALSE, trim_ws = TRUE)
cdp_3 <- read_delim(file = paste0(input.path, "mcalls_cdp_3_CG_chrom-merged_strands-merged.bed.gz"), 
                    "\t", escape_double = FALSE, trim_ws = TRUE)
cdp_4 <- read_delim(file = paste0(input.path, "mcalls_cdp_4_CG_chrom-merged_strands-merged.bed.gz"), 
                    "\t", escape_double = FALSE, trim_ws = TRUE)

cdp <- list(cdp_1, cdp_2, cdp_3, cdp_4) %>%
  map(~ select(., -motif, -score, -strand, -beta_value)) %>%
  purrr::reduce(., full_join, by = c("#chrom", "start", "end")) %>%
  mutate(., total = rowSums(select(., contains("n_total")))) %>%
  mutate(., meth = rowSums(select(., contains("n_meth")))) %>%
  select(., `#chrom`, start, end, total, meth)


# load methylation calls of CMOP

cmop_1 <- read_delim(file = paste0(input.path, "mcalls_cmop_1_CG_chrom-merged_strands-merged.bed.gz"), 
                     "\t", escape_double = FALSE, trim_ws = TRUE)
cmop_2 <- read_delim(file = paste0(input.path, "mcalls_cmop_2_CG_chrom-merged_strands-merged.bed.gz"), 
                     "\t", escape_double = FALSE, trim_ws = TRUE)
    
cmop <- list(cmop_1, cmop_2) %>%
  map(~ select(., -motif, -score, -strand, -beta_value)) %>%
  purrr::reduce(., full_join, by = c("#chrom", "start", "end")) %>%
  mutate(., total = rowSums(select(., contains("n_total")))) %>%
  mutate(., meth = rowSums(select(., contains("n_meth")))) %>%
  select(., `#chrom`, start, end, total, meth)


# load methylation calls of DC-CD11b

dc_cd11b_1 <- read_delim(file = paste0(input.path, "mcalls_dc-cd11b_1_CG_chrom-merged_strands-merged.bed.gz"), 
                         "\t", escape_double = FALSE, trim_ws = TRUE)
dc_cd11b_2 <- read_delim(file = paste0(input.path, "mcalls_dc-cd11b_2_CG_chrom-merged_strands-merged.bed.gz"), 
                         "\t", escape_double = FALSE, trim_ws = TRUE)
dc_cd11b_3 <- read_delim(file = paste0(input.path, "mcalls_dc-cd11b_3_CG_chrom-merged_strands-merged.bed.gz"), 
                         "\t", escape_double = FALSE, trim_ws = TRUE)

dc_cd11b <- list(dc_cd11b_1, dc_cd11b_2, dc_cd11b_3) %>%
  map(~ select(., -motif, -score, -strand, -beta_value)) %>%
  purrr::reduce(., full_join, by = c("#chrom", "start", "end")) %>%
  mutate(., total = rowSums(select(., contains("n_total")))) %>%
  mutate(., meth = rowSums(select(., contains("n_meth")))) %>%
  select(., `#chrom`, start, end, total, meth)


# load methylation calls of DC-CD8a

dc_cd8a_1 <- read_delim(file = paste0(input.path, "mcalls_dc-cd8a_1_CG_chrom-merged_strands-merged.bed.gz"), 
                        "\t", escape_double = FALSE, trim_ws = TRUE)
dc_cd8a_2 <- read_delim(file = paste0(input.path, "mcalls_dc-cd8a_2_CG_chrom-merged_strands-merged.bed.gz"), 
                        "\t", escape_double = FALSE, trim_ws = TRUE)
dc_cd8a_3 <- read_delim(file = paste0(input.path, "mcalls_dc-cd8a_3_CG_chrom-merged_strands-merged.bed.gz"), 
                        "\t", escape_double = FALSE, trim_ws = TRUE)

dc_cd8a <- list(dc_cd8a_1, dc_cd8a_2, dc_cd8a_3) %>%
  map(~ select(., -motif, -score, -strand, -beta_value)) %>%
  purrr::reduce(., full_join, by = c("#chrom", "start", "end")) %>%
  mutate(., total = rowSums(select(., contains("n_total")))) %>%
  mutate(., meth = rowSums(select(., contains("n_meth")))) %>%
  select(., `#chrom`, start, end, total, meth)


# load methylation calls of Monos

mono_1 <- read_delim(file = paste0(input.path, "mcalls_monos-new_1_CG_chrom-merged_strands-merged.bed.gz"), 
                     "\t", escape_double = FALSE, trim_ws = TRUE)
mono_2 <- read_delim(file = paste0(input.path, "mcalls_monos-new_2_CG_chrom-merged_strands-merged.bed.gz"), 
                     "\t", escape_double = FALSE, trim_ws = TRUE)
mono_3 <- read_delim(file = paste0(input.path, "mcalls_monos-new_3_CG_chrom-merged_strands-merged.bed.gz"), 
                     "\t", escape_double = FALSE, trim_ws = TRUE)

mono <- list(mono_1, mono_2, mono_3) %>%
  map(~ select(., -motif, -score, -strand, -beta_value)) %>%
  purrr::reduce(., full_join, by = c("#chrom", "start", "end")) %>%
  mutate(., total = rowSums(select(., contains("n_total")))) %>%
  mutate(., meth = rowSums(select(., contains("n_meth")))) %>%
  select(., `#chrom`, start, end, total, meth)


# load methylation calls of MDPs

mdp_1 <- read_delim(file = paste0(input.path, "mcalls_mdp_1_CG_chrom-merged_strands-merged.bed.gz"), 
                    "\t", escape_double = FALSE, trim_ws = TRUE)
mdp_2 <- read_delim(file = paste0(input.path, "mcalls_mdp_2_CG_chrom-merged_strands-merged.bed.gz"), 
                    "\t", escape_double = FALSE, trim_ws = TRUE)

mdp <- list(mdp_1, mdp_2) %>%
  map(~ select(., -motif, -score, -strand, -beta_value)) %>%
  purrr::reduce(., full_join, by = c("#chrom", "start", "end")) %>%
  mutate(., total = rowSums(select(., contains("n_total")))) %>%
  mutate(., meth = rowSums(select(., contains("n_meth")))) %>%
  select(., `#chrom`, start, end, total, meth)


# load methylation calls of PDCs

pdc_1 <- read_delim(file = paste0(input.path, "mcalls_pdc_1_CG_chrom-merged_strands-merged.bed.gz"), 
                    "\t", escape_double = FALSE, trim_ws = TRUE)
pdc_2 <- read_delim(file = paste0(input.path, "mcalls_pdc_2_CG_chrom-merged_strands-merged.bed.gz"), 
                    "\t", escape_double = FALSE, trim_ws = TRUE)
pdc_3 <- read_delim(file = paste0(input.path, "mcalls_pdc_3_CG_chrom-merged_strands-merged.bed.gz"), 
                    "\t", escape_double = FALSE, trim_ws = TRUE)

pdc <- list(pdc_1, pdc_2, pdc_3) %>%
  map(~ select(., -motif, -score, -strand, -beta_value)) %>%
  purrr::reduce(., full_join, by = c("#chrom", "start", "end")) %>%
  mutate(., total = rowSums(select(., contains("n_total")))) %>%
  mutate(., meth = rowSums(select(., contains("n_meth")))) %>%
  select(., `#chrom`, start, end, total, meth)


# make list of all samples ------------------------------------------------------
    
  sample.list <- list(hsc, cdp, cmop, dc_cd11b, dc_cd8a, mono, mdp, pdc)
  
  names(sample.list) <- c("hsc", "cdp", "cmop", "dc_cd11b", "dc_cd8a", "mono", "mdp", "pdc")

# generate methylation mm, coverage mm and gr objects ------------------------------------------------------
 
  chr = c(1:19) 
  
  for (i in chr)
  {
    
    # filter for chr, add 1 bp for start and select columns ------------------------------------------------------
    
    x <- sample.list %>% 
          map(~filter(.x, `#chrom` == i)) %>%
          map(~mutate(.x, start = start + 1)) 
    assign(paste0("sample.list.chr", i), x)
    
    # make methylation matrix  --------------------------------------------------
    
    meth_mat <- x %>%
      map_dfc(~select(.x, 'meth')) %>%
      as.matrix(.)
    colnames(meth_mat) <- names(sample.list)
    assign(paste0("meth_mat.chr", i), meth_mat)
    
    # make coverage matrix  --------------------------------------------------
    
    cov_mat <- x %>%
      map_dfc(~select(.x, 'total')) %>%
      as.matrix(.)
    colnames(cov_mat) <- names(sample.list)
    assign(paste0("cov_mat.chr", i), cov_mat)
    
    # make GR object  --------------------------------------------------
    
    gr <- x[[1]] %>%
      select(`#chrom`, start, end) %>%
      dplyr::rename(seqnames = `#chrom`) %>%
      makeGRangesFromDataFrame(.) 
    assign(paste0("gr.chr", i), gr)
    
  }

# generate lists of methylation mm, coverage mm and gr objects for all chrom --------------------------------------------------    
  
  meth_mat.list <- mget(ls(pattern = "^meth_mat.chr"))
  cov_mat.list <- mget(ls(pattern = "^cov_mat.chr"))
  gr.list <- mget(ls(pattern = "^gr.chr"))


# methylation level per replicate per DMR using gUtils package   --------------------------------------------------
  
  out.dir <- "/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/DMR calling/meth-levels_merged-dmrs/"
  
  meth.all <- mget(ls(pattern = "^meth_mat.chr")) %>%
    map(~as.data.frame(.x)) %>%
    map_dfr(~dplyr::bind_rows(.x)) 
  #saveRDS(meth.all, file = file.path(out.dir, paste0(Sys.Date(), "_meth-matrix_DC-DMRs_allChr_pop.rds")))
  
  cov.all <- mget(ls(pattern = "^cov_mat.chr")) %>%
    map(~as.data.frame(.x)) %>%
    map_dfr(~dplyr::bind_rows(.x)) 
  #saveRDS(cov.all, file = file.path(out.dir, paste0(Sys.Date(), "_cov-matrix_DC-DMRs_allChr_pop.rds")))
  
  gr.all <- mget(ls(pattern = "^gr.chr")) %>%
    map(~as.data.frame(.x)) %>%
    map_dfr(~dplyr::bind_rows(.x)) %>%
    makeGRangesFromDataFrame(.)
  #saveRDS(gr.all, file = file.path(out.dir, paste0(Sys.Date(), "_gr-object_DC-DMRs_allChr_pop.rds")))
  
  # calculate methylation levels per replicate   --------------------------------------------------

  mcols(gr.all) <- as.data.frame(meth.all)
  dmrs_sum.meth <- merged_dmrs.gr %$% gr.all
  methsum <- as.data.frame(mcols(dmrs_sum.meth))
  
  mcols(gr.all) <- as.data.frame(cov.all)
  dmrs_sum.cov <- merged_dmrs.gr %$% gr.all
  covsum <- as.data.frame(mcols(dmrs_sum.cov))
  
  beta_value.df <- as.data.frame(as.matrix(methsum)/as.matrix(covsum))
  
  beta_value_allDMRs.df <- as.data.frame(merged_dmrs.gr) %>%
    select("seqnames", "start", "end") %>%
    dplyr::rename(., chr = "seqnames") %>%
    cbind(., beta_value.df) %>%
    arrange(chr, start)
  
  out.dir <- "/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/DMR calling/meth-levels_merged-dmrs/"
  saveRDS(beta_value_allDMRs.df, file = file.path(out.dir, paste0(Sys.Date(), "_beta-value_DC-DMRs_allChr_pop.rds")))
  write_csv(beta_value_allDMRs.df, path = file.path(out.dir, paste0(Sys.Date(), "_beta-value_DC-DMRs_allChr_pop.csv")))
  
  
  sessionInfo()
  
  # R version 3.6.2 (2019-12-12)
  # Platform: x86_64-pc-linux-gnu (64-bit)
  # Running under: CentOS Linux 7 (Core)
  # 
  # Matrix products: default
  # BLAS:   /usr/lib64/libblas.so.3.4.2
  # LAPACK: /usr/lib64/liblapack.so.3.4.2
  # 
  # Random number generation:
  #   RNG:     Mersenne-Twister 
  # Normal:  Inversion 
  # Sample:  Rounding 
  # 
  # locale:
  #   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
  # [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
  # [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
  # [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
  # 
  # attached base packages:
  #   [1] parallel  stats4    grid      stats     graphics  grDevices utils     datasets  methods   base     
  # 
  # other attached packages:
  #   [1] gUtils_0.2.0          data.table_1.13.0     GenomicRanges_1.38.0  GenomeInfoDb_1.22.1  
  # [5] IRanges_2.20.2        S4Vectors_0.24.4      BiocGenerics_0.32.0   RColorBrewer_1.1-2   
  # [9] circlize_0.4.10       ComplexHeatmap_2.2.0  dynamicTreeCut_1.63-1 forcats_0.5.0        
  # [13] stringr_1.4.0         dplyr_0.8.5           purrr_0.3.4           readr_1.3.1          
  # [17] tidyr_1.1.0           tibble_3.0.1          ggplot2_3.3.0         tidyverse_1.3.0      
  # 
  # loaded via a namespace (and not attached):
  #   [1] Rcpp_1.0.5             lubridate_1.7.9        lattice_0.20-41        png_0.1-7             
  # [5] assertthat_0.2.1       utf8_1.1.4             R6_2.4.1               cellranger_1.1.0      
  # [9] backports_1.1.10       reprex_0.3.0           httr_1.4.2             pillar_1.4.4          
  # [13] zlibbioc_1.32.0        GlobalOptions_0.1.2    rlang_0.4.7            readxl_1.3.1          
  # [17] rstudioapi_0.11        GetoptLong_1.0.2       RCurl_1.98-1.2         munsell_0.5.0         
  # [21] broom_0.5.6            compiler_3.6.2         modelr_0.1.8           pkgconfig_2.0.3       
  # [25] shape_1.4.5            tidyselect_1.1.0       GenomeInfoDbData_1.2.2 fansi_0.4.1           
  # [29] crayon_1.3.4           dbplyr_1.4.3           withr_2.3.0            bitops_1.0-6          
  # [33] nlme_3.1-149           jsonlite_1.7.1         gtable_0.3.0           lifecycle_0.2.0       
  # [37] DBI_1.1.0              magrittr_1.5           scales_1.1.1           cli_2.0.2             
  # [41] stringi_1.5.3          XVector_0.26.0         fs_1.5.0               xml2_1.3.2            
  # [45] ellipsis_0.3.1         generics_0.0.2         vctrs_0.3.0            rjson_0.2.20          
  # [49] tools_3.6.2            glue_1.4.2             hms_0.5.3              clue_0.3-57           
  # [53] colorspace_1.4-1       cluster_2.1.0          rvest_0.3.5            haven_2.2.0   