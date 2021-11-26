# Methylation levels in merged DC DMRs - DC CD8a ------------------------------------------------------
# Dr. Sina Stäble
# 01.10.2020

# load libraries
library(tidyverse)
library(GenomicRanges)
library(gUtils)

# load merged DMRs ------------------------------------------------------

merged_dmrs.gr <- readRDS("/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/DMR calling/merged_dmrs/2020-10-01_merged-dmrs.rds")


# load methylation calls ------------------------------------------------------

input.path <- "/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/meth_calls/"


# load methylation calls of DC-CD11b

dc_cd8a_1 <- read_delim(file = paste0(input.path, "mcalls_dc-cd8a_1_CG_chrom-merged_strands-merged.bed.gz"), 
                         "\t", escape_double = FALSE, trim_ws = TRUE)
dc_cd8a_2 <- read_delim(file = paste0(input.path, "mcalls_dc-cd8a_2_CG_chrom-merged_strands-merged.bed.gz"), 
                         "\t", escape_double = FALSE, trim_ws = TRUE)
dc_cd8a_3 <- read_delim(file = paste0(input.path, "mcalls_dc-cd8a_3_CG_chrom-merged_strands-merged.bed.gz"), 
                         "\t", escape_double = FALSE, trim_ws = TRUE)

# make list of all samples and replicates ------------------------------------------------------
    
  sample.list <- list(dc_cd8a_1, dc_cd8a_2, dc_cd8a_3)
  
  names(sample.list) <- c("dc_cd8a_1", "dc_cd8a_2", "dc_cd8a_3")

# generate methylation mm, coverage mm and gr objects ------------------------------------------------------
 
  chr = c(1:19) 
  
  for (i in chr)
  {
    
    # filter for chr, add 1 bp for start and select columns ------------------------------------------------------
    
    x <- sample.list %>% 
          map(~filter(.x, `#chrom` == i)) %>%
          map(~mutate(.x, start = start + 1)) %>%
          map(~select(.x, `#chrom`, start, end, n_total, n_meth))
    assign(paste0("sample.list.chr", i), x)
    
    # make methylation matrix  --------------------------------------------------
    
    meth_mat <- x %>%
      map_dfc(~select(.x, 'n_meth')) %>%
      as.matrix(.)
    colnames(meth_mat) <- names(sample.list)
    assign(paste0("meth_mat.chr", i), meth_mat)
    
    # make coverage matrix  --------------------------------------------------
    
    cov_mat <- x %>%
      map_dfc(~select(.x, 'n_total')) %>%
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
  #saveRDS(meth.all, file = file.path(out.dir, paste0(Sys.Date(), "_meth-matrix_DC-DMRs_allChr_CD-CD8a.rds")))
  
  cov.all <- mget(ls(pattern = "^cov_mat.chr")) %>%
    map(~as.data.frame(.x)) %>%
    map_dfr(~dplyr::bind_rows(.x)) 
  #saveRDS(cov.all, file = file.path(out.dir, paste0(Sys.Date(), "_cov-matrix_DC-DMRs_allChr_CD-CD8a.rds")))
  
  gr.all <- mget(ls(pattern = "^gr.chr")) %>%
    map(~as.data.frame(.x)) %>%
    map_dfr(~dplyr::bind_rows(.x)) %>%
    makeGRangesFromDataFrame(.)
  #saveRDS(gr.all, file = file.path(out.dir, paste0(Sys.Date(), "_gr-object_DC-DMRs_allChr_CD-CD8a.rds")))
  
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
  saveRDS(beta_value_allDMRs.df, file = file.path(out.dir, paste0(Sys.Date(), "_beta-value_DC-DMRs_allChr_CD-CD8a.rds")))
  write_delim(beta_value_allDMRs.df, path = file.path(out.dir, paste0(Sys.Date(), "_beta-value_DC-DMRs_allChr_CD-CD8a.csv")))
  
  
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
  #   [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
  # 
  # other attached packages:
  #   [1] gUtils_0.2.0         data.table_1.13.0    GenomicRanges_1.38.0 GenomeInfoDb_1.22.1 
  # [5] IRanges_2.20.2       S4Vectors_0.24.4     BiocGenerics_0.32.0  forcats_0.5.0       
  # [9] stringr_1.4.0        dplyr_0.8.5          purrr_0.3.4          readr_1.3.1         
  # [13] tidyr_1.1.0          tibble_3.0.1         ggplot2_3.3.0        tidyverse_1.3.0     
  # 
  # loaded via a namespace (and not attached):
  #   [1] tidyselect_1.1.0       haven_2.2.0            lattice_0.20-41        colorspace_1.4-1      
  # [5] vctrs_0.3.0            generics_0.0.2         rlang_0.4.7            pillar_1.4.4          
  # [9] glue_1.4.2             withr_2.3.0            DBI_1.1.0              dbplyr_1.4.3          
  # [13] modelr_0.1.8           readxl_1.3.1           GenomeInfoDbData_1.2.2 lifecycle_0.2.0       
  # [17] zlibbioc_1.32.0        munsell_0.5.0          gtable_0.3.0           cellranger_1.1.0      
  # [21] rvest_0.3.5            fansi_0.4.1            broom_0.5.6            Rcpp_1.0.5            
  # [25] scales_1.1.1           backports_1.1.10       XVector_0.26.0         jsonlite_1.7.1        
  # [29] fs_1.5.0               hms_0.5.3              stringi_1.5.3          grid_3.6.2            
  # [33] cli_2.0.2              tools_3.6.2            bitops_1.0-6           magrittr_1.5          
  # [37] RCurl_1.98-1.2         crayon_1.3.4           pkgconfig_2.0.3        ellipsis_0.3.1        
  # [41] xml2_1.3.2             reprex_0.3.0           lubridate_1.7.9        assertthat_0.2.1      
  # [45] httr_1.4.2             rstudioapi_0.11        R6_2.4.1               nlme_3.1-149          
  # [49] compiler_3.6.2 