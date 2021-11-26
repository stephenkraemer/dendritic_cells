# Merge DMRs from multiple pairwise DMR calls ------------------------------------------------------
# Dr. Sina Stäble
# 30.09.2020

# load libraries
library(tidyverse)
library(GenomicRanges)
#library(gUtils)

# load DMRs

HSC_CDP <- readRDS("/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/DMR calling/HSC_vs_CDP/2020-10-01_dmrs_hsc_vs_cdp_0.01.rds") %>%
  makeGRangesFromDataFrame(.)

HSC_MDP <- readRDS("/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/DMR calling/HSC_vs_MDP/2020-10-01_dmrs_hsc_vs_mdp_0.01.rds") %>%
  makeGRangesFromDataFrame(.)

HSC_CMOP <- readRDS("/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/DMR calling/HSC_vs_CMOP/2020-10-01_dmrs_hsc_vs_cmop_0.01.rds") %>%
  makeGRangesFromDataFrame(.)

HSC_DCCD11b <- readRDS("/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/DMR calling/HSC_vs_DC-CD11b/2020-10-01_dmrs_hsc_vs_dc-cd11b_0.01.rds") %>%
  makeGRangesFromDataFrame(.)

HSC_DCCD8a <- readRDS("/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/DMR calling/HSC_vs_DC-CD8a/2020-10-01_dmrs_hsc_vs_dc-cd8a_0.01.rds") %>%
  makeGRangesFromDataFrame(.)

HSC_MONO <- readRDS("/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/DMR calling/HSC_vs_MONO/2020-10-01_dmrs_hsc_vs_mono_0.01.rds") %>%
  makeGRangesFromDataFrame(.)

HSC_PDC <- readRDS("/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/DMR calling/HSC_vs_PDC/2020-10-01_dmrs_hsc_vs_pdc_0.01.rds") %>%
  makeGRangesFromDataFrame(.)

# merge DMRs

all_dmrs <- c(HSC_CDP, HSC_MDP, HSC_CMOP, HSC_DCCD11b, HSC_DCCD8a, HSC_MONO, HSC_PDC)

merged_dmrs <- reduce(all_dmrs, 
                   #min.gapwidth=500
                   )
out.dir <- "/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/DMR calling/"
saveRDS(merged_dmrs, file = file.path(out.dir, paste0(Sys.Date(), "_merged-dmrs.rds")))

merged_dmrs.df <- as.data.frame(merged_dmrs)
write_delim(merged_dmrs.df, path = file.path(out.dir, paste0(Sys.Date(), "_merged-dmrs.csv")))

