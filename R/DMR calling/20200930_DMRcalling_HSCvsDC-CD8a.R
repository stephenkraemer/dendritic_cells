# DMR calling with bsseq - HSC vs DC CD8a ------------------------------------------------------
# Dr. Sina Stäble
# 30.09.2020

# load libraries
library(tidyverse)
library(bsseq)
library(DSS)
library(GenomicRanges)
library(gUtils)
library(data.table)

# load methylation calls ------------------------------------------------------

  # load methylation calls of HSCs

    input.path <- "/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/meth_calls/"
    
    hsc_1 <- read_delim(file = paste0(input.path, "mcalls_hsc_1_CG_chrom-merged_strands-merged.bed.gz"), 
                              "\t", escape_double = FALSE, trim_ws = TRUE)
    hsc_2 <- read_delim(file = paste0(input.path, "mcalls_hsc_2_CG_chrom-merged_strands-merged.bed.gz"), 
                        "\t", escape_double = FALSE, trim_ws = TRUE)
    hsc_3 <- read_delim(file = paste0(input.path, "mcalls_hsc_3_CG_chrom-merged_strands-merged.bed.gz"), 
                        "\t", escape_double = FALSE, trim_ws = TRUE)

  # load methylation calls of DC-CD8a
  
    dc_cd8a_1 <- read_delim(file = paste0(input.path, "mcalls_dc-cd8a_1_CG_chrom-merged_strands-merged.bed.gz"), 
                             "\t", escape_double = FALSE, trim_ws = TRUE)
    dc_cd8a_2 <- read_delim(file = paste0(input.path, "mcalls_dc-cd8a_2_CG_chrom-merged_strands-merged.bed.gz"), 
                             "\t", escape_double = FALSE, trim_ws = TRUE)
    dc_cd8a_3 <- read_delim(file = paste0(input.path, "mcalls_dc-cd8a_3_CG_chrom-merged_strands-merged.bed.gz"), 
                             "\t", escape_double = FALSE, trim_ws = TRUE)
    
# make list of all samples and replicates ------------------------------------------------------
    
  sample.list <- list(hsc_1, hsc_2, hsc_3,
                      dc_cd8a_1, dc_cd8a_2, dc_cd8a_3)
  
  names(sample.list) <- c("hsc_1", "hsc_2", "hsc_3", 
                          "dc_cd8a_1", "dc_cd8a_2", "dc_cd8a_3")

# generate Bsseq object ------------------------------------------------------
 
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
    
    # make BSseq object --------------------------------------------------
    
    BSobj <- BSseq(gr = gr, M = meth_mat, Cov = cov_mat, sampleNames = names(sample.list))
    assign(paste0("BSseq.chr", i), BSobj)
    rm(x, meth_mat, gr, cov_mat, BSobj)
  }

# generate lists of methylation mm, coverage mm and gr objects for all chrom --------------------------------------------------    
  
  #meth_mat.list <- mget(ls(pattern = "^meth_mat.chr"))
  #cov_mat.list <- mget(ls(pattern = "^cov_mat.chr"))
  #gr.list <- mget(ls(pattern = "^gr.chr"))
  
# perform DML Test   --------------------------------------------------

  BSobj.list <- mget(ls(pattern = "^BSseq.chr"))
  
  for(i in 1:length(BSobj.list))
  {
    dmlTest <- DMLtest(BSobj.list[[i]], 
                        group1 = c("hsc_1", "hsc_2", "hsc_3"), 
                        group2 = c("dc_cd8a_1", "dc_cd8a_2", "dc_cd8a_3"),
                        equal.disp = FALSE, smoothing = FALSE, smoothing.span = 500)
    x = str_extract(names(BSobj.list[i]), pattern = "chr[:digit:]+")
    out.dir <- "/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/DMR calling/HSC_vs_DC-CD8a/dmltest/"
    saveRDS(dmlTest, file = file.path(out.dir, paste0(Sys.Date(), "_dmltest_", "hsc_vs_dc-cd8a_", x, ".rds")))
    assign(paste0("dmlTest.", x), dmlTest)
    rm(dmlTest)
  }
    
  
# call DMRs   --------------------------------------------------
    
  dmlTest.list <- mget(ls(pattern = "^dmlTest.chr"))  

  for(i in 1:length(dmlTest.list))
  {  
    dmrs <- callDMR(dmlTest.list[[i]], 
                    delta = 0.1, 
                    p.threshold = 0.01, 
                    minlen = 50, 
                    minCG = 2, 
                    dis.merge = 50, 
                    pct.sig = 0.5)
    x = str_extract(names(dmlTest.list[i]), pattern = "chr[:digit:]+")
    out.dir <- "//icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/DMR calling/HSC_vs_DC-CD8a/dmrs/"
    saveRDS(dmrs, file = file.path(out.dir, paste0(Sys.Date(), "_dmrs_", "hsc_vs_dc-cd8a_", x, "_0.01.rds")))
    assign(paste0("dmrs.", x), dmrs)
    rm(dmrs)
  }
 
# make df with all DMRs (all Chr) --------------------------------------------------
    
  dmrs.list <- mget(ls(pattern = "^dmrs.chr"))   
  dmrs.all.df <- rbindlist(dmrs.list)
  
# save dmrs dataframe
  
  out.dir <- "/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/DMR calling/HSC_vs_DC-CD8a/"
  saveRDS(dmrs.all.df, file = file.path(out.dir, paste0(Sys.Date(), "_dmrs_", "hsc_vs_dc-cd8a_0.01.rds")))

  
  sessionInfo(package = NULL)
  
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
  #   [1] splines   stats4    parallel  stats     graphics  grDevices utils     datasets  methods  
  # [10] base     
  # 
  # other attached packages:
  #   [1] gUtils_0.2.0                data.table_1.13.0           DSS_2.34.0                 
  # [4] bsseq_1.22.0                SummarizedExperiment_1.16.1 DelayedArray_0.12.3        
  # [7] BiocParallel_1.20.1         matrixStats_0.57.0          Biobase_2.46.0             
  # [10] GenomicRanges_1.38.0        GenomeInfoDb_1.22.1         IRanges_2.20.2             
  # [13] S4Vectors_0.24.4            BiocGenerics_0.32.0         forcats_0.5.0              
  # [16] stringr_1.4.0               dplyr_0.8.5                 purrr_0.3.4                
  # [19] readr_1.3.1                 tidyr_1.1.0                 tibble_3.0.1               
  # [22] ggplot2_3.3.0               tidyverse_1.3.0            
  # 
  # loaded via a namespace (and not attached):
  #   [1] httr_1.4.2               jsonlite_1.7.1           R.utils_2.10.1          
  # [4] DelayedMatrixStats_1.8.0 modelr_0.1.8             gtools_3.8.2            
  # [7] assertthat_0.2.1         BSgenome_1.54.0          GenomeInfoDbData_1.2.2  
  # [10] cellranger_1.1.0         Rsamtools_2.2.3          pillar_1.4.4            
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
  # [61] rtracklayer_1.46.0       permute_0.9-5            stringi_1.5.3           
  # [64] Rcpp_1.0.5               vctrs_0.3.0              dbplyr_1.4.3            
  # [67] tidyselect_1.1.0 