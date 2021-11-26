# Methylation levels in merged DC DMRs - all DC populations ------------------------------------------------------
# Dr. Sina Stäble
# 01.10.2020

# load libraries
library(tidyverse)

# load meth matrix DMRs per population ------------------------------------------------------

DC_CD11b <- readRDS("/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/DMR calling/meth-levels_merged-dmrs/2020-10-01_beta-value_DC-DMRs_allChr_DC-CD11b.rds")
DC_CD8a <- readRDS("/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/DMR calling/meth-levels_merged-dmrs/2020-10-01_beta-value_DC-DMRs_allChr_DC-CD8a.rds")
CDP <- readRDS("/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/DMR calling/meth-levels_merged-dmrs/2020-10-01_beta-value_DC-DMRs_allChr_CDP.rds")
CMOP <- readRDS("/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/DMR calling/meth-levels_merged-dmrs/2020-10-01_beta-value_DC-DMRs_allChr_CMOP.rds")
HSC <- readRDS("/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/DMR calling/meth-levels_merged-dmrs/2020-10-01_beta-value_DC-DMRs_allChr_HSC.rds")
MDP <- readRDS("/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/DMR calling/meth-levels_merged-dmrs/2020-10-01_beta-value_DC-DMRs_allChr_MDP.rds")
MONO <- readRDS("/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/DMR calling/meth-levels_merged-dmrs/2020-10-01_beta-value_DC-DMRs_allChr_MONO.rds")
PDC <- readRDS("/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/DMR calling/meth-levels_merged-dmrs/2020-10-01_beta-value_DC-DMRs_allChr_PDC.rds")

DC_pop.list <- list(HSC, MDP, CDP, CMOP, DC_CD11b, DC_CD8a, PDC, MONO)

# merge all pop into one df  ------------------------------------------------------

DC_pop.df <- DC_pop.list %>%
  reduce(., full_join, by = c("chr", "start", "end"))
  
out.dir <- "/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/DMR calling/meth-levels_merged-dmrs/"
saveRDS(DC_pop.df, file = file.path(out.dir, paste0(Sys.Date(), "_beta-value_DC-DMRs_all-DC-pop.rds")))
write_csv(DC_pop.df, path = file.path(out.dir, paste0(Sys.Date(), "_beta-value_DC-DMRs_all-DC-pop.csv")))

sessionInfo()
# 
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
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] forcats_0.5.0   stringr_1.4.0   dplyr_0.8.5     purrr_0.3.4     readr_1.3.1     tidyr_1.1.0    
# [7] tibble_3.0.1    ggplot2_3.3.0   tidyverse_1.3.0
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.5       cellranger_1.1.0 pillar_1.4.4     compiler_3.6.2   dbplyr_1.4.3     tools_3.6.2     
# [7] jsonlite_1.7.1   lubridate_1.7.9  lifecycle_0.2.0  nlme_3.1-149     gtable_0.3.0     lattice_0.20-41 
# [13] pkgconfig_2.0.3  rlang_0.4.7      reprex_0.3.0     cli_2.0.2        DBI_1.1.0        rstudioapi_0.11 
# [19] haven_2.2.0      withr_2.3.0      xml2_1.3.2       httr_1.4.2       fs_1.5.0         generics_0.0.2  
# [25] vctrs_0.3.0      hms_0.5.3        grid_3.6.2       tidyselect_1.1.0 glue_1.4.2       R6_2.4.1        
# [31] fansi_0.4.1      readxl_1.3.1     modelr_0.1.8     magrittr_1.5     backports_1.1.10 scales_1.1.1    
# [37] ellipsis_0.3.1   rvest_0.3.5      assertthat_0.2.1 colorspace_1.4-1 stringi_1.5.3    munsell_0.5.0   
# [43] broom_0.5.6      crayon_1.3.4 