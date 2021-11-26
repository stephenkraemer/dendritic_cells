# Methylation levels in merged DC DMRs - DC CD11b ------------------------------------------------------
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

dc_cd11b_1 <- read_delim(file = paste0(input.path, "mcalls_dc-cd11b_1_CG_chrom-merged_strands-merged.bed.gz"), 
                         "\t", escape_double = FALSE, trim_ws = TRUE)
dc_cd11b_2 <- read_delim(file = paste0(input.path, "mcalls_dc-cd11b_2_CG_chrom-merged_strands-merged.bed.gz"), 
                         "\t", escape_double = FALSE, trim_ws = TRUE)
dc_cd11b_3 <- read_delim(file = paste0(input.path, "mcalls_dc-cd11b_3_CG_chrom-merged_strands-merged.bed.gz"), 
                         "\t", escape_double = FALSE, trim_ws = TRUE)

# make list of all samples and replicates ------------------------------------------------------
    
  sample.list <- list(dc_cd11b_1, dc_cd11b_2, dc_cd11b_3)
  
  names(sample.list) <- c("dc_cd11b_1", "dc_cd11b_2", "dc_cd11b_3")

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
  #saveRDS(meth.all, file = file.path(out.dir, paste0(Sys.Date(), "_meth-matrix_DC-DMRs_allChr_CD-CD11b.rds")))
  
  cov.all <- mget(ls(pattern = "^cov_mat.chr")) %>%
    map(~as.data.frame(.x)) %>%
    map_dfr(~dplyr::bind_rows(.x)) 
  #saveRDS(cov.all, file = file.path(out.dir, paste0(Sys.Date(), "_cov-matrix_DC-DMRs_allChr_CD-CD11b.rds")))
  
  gr.all <- mget(ls(pattern = "^gr.chr")) %>%
    map(~as.data.frame(.x)) %>%
    map_dfr(~dplyr::bind_rows(.x)) %>%
    makeGRangesFromDataFrame(.)
  #saveRDS(gr.all, file = file.path(out.dir, paste0(Sys.Date(), "_gr-object_DC-DMRs_allChr_CD-CD11b.rds")))
  
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
  saveRDS(beta_value_allDMRs.df, file = file.path(out.dir, paste0(Sys.Date(), "_beta-value_DC-DMRs_allChr_CD-CD11b.rds")))
  write_delim(beta_value_allDMRs.df, path = file.path(out.dir, paste0(Sys.Date(), "_beta-value_DC-DMRs_allChr_CD-CD11b.csv")))
  
  
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
  #   [1] bsseq_1.22.0                SummarizedExperiment_1.16.1 DelayedArray_0.12.3        
  # [4] BiocParallel_1.20.1         matrixStats_0.57.0          Biobase_2.46.0             
  # [7] gUtils_0.2.0                data.table_1.13.0           GenomicRanges_1.38.0       
  # [10] GenomeInfoDb_1.22.1         IRanges_2.20.2              S4Vectors_0.24.4           
  # [13] BiocGenerics_0.32.0         forcats_0.5.0               stringr_1.4.0              
  # [16] dplyr_0.8.5                 purrr_0.3.4                 readr_1.3.1                
  # [19] tidyr_1.1.0                 tibble_3.0.1                ggplot2_3.3.0              
  # [22] tidyverse_1.3.0            
  # 
  # loaded via a namespace (and not attached):
  #   [1] httr_1.4.2               jsonlite_1.7.1           R.utils_2.10.1          
  # [4] DelayedMatrixStats_1.8.0 modelr_0.1.8             gtools_3.8.2            
  # [7] assertthat_0.2.1         Rsamtools_2.2.3          BSgenome_1.54.0         
  # [10] GenomeInfoDbData_1.2.2   cellranger_1.1.0         pillar_1.4.4            
  # [13] backports_1.1.10         lattice_0.20-41          glue_1.4.2              
  # [16] limma_3.42.2             XVector_0.26.0           rvest_0.3.5             
  # [19] colorspace_1.4-1         R.oo_1.24.0              Matrix_1.2-18           
  # [22] XML_3.99-0.3             pkgconfig_2.0.3          broom_0.5.6             
  # [25] haven_2.2.0              zlibbioc_1.32.0          scales_1.1.1            
  # [28] HDF5Array_1.14.4         generics_0.0.2           ellipsis_0.3.1          
  # [31] withr_2.3.0              cli_2.0.2                magrittr_1.5            
  # [34] crayon_1.3.4             readxl_1.3.1             R.methodsS3_1.8.1       
  # [37] fs_1.5.0                 fansi_0.4.1              nlme_3.1-149            
  # [40] xml2_1.3.2               tools_3.6.2              hms_0.5.3               
  # [43] lifecycle_0.2.0          Rhdf5lib_1.8.0           munsell_0.5.0           
  # [46] reprex_0.3.0             locfit_1.5-9.4           Biostrings_2.54.0       
  # [49] compiler_3.6.2           rlang_0.4.7              rhdf5_2.30.1            
  # [52] grid_3.6.2               RCurl_1.98-1.2           rstudioapi_0.11         
  # [55] bitops_1.0-6             gtable_0.3.0             DBI_1.1.0               
  # [58] R6_2.4.1                 GenomicAlignments_1.22.1 lubridate_1.7.9         
  # [61] rtracklayer_1.46.0       utf8_1.1.4               permute_0.9-5           
  # [64] stringi_1.5.3            Rcpp_1.0.5               vctrs_0.3.0             
  # [67] dbplyr_1.4.3             tidyselect_1.1.0   