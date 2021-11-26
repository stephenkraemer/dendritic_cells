#### 20201006_PCA (all hierarchy DMRs)

# Author: Dr. Sina Stäble

# load libraries

library(tidyverse)
library(FactoMineR)
library(factoextra)
library(RColorBrewer)


#load beta values in DC dmrs ----------------------------------------------------
dmrs <- readRDS("/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/DMR calling/meth-levels_merged-dmrs/2020-10-01_beta-value_DC-DMRs_all-DC-rep.rds")

# transform data table  ----------------------------------------------------

dmrs.t <- dmrs %>%
  select(., -chr, -start, -end) %>%
  t(.)

# calculate PCA ----------------------------------------------------
PCA <- PCA(dmrs.t, scale.unit = FALSE)

# plot PCA  ----------------------------------------------------

df <- PCA$ind$coord 
colnames(df) <- paste0("PC",1:5)


# PC1 on x-axis

setwd("/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/PCA")

for (a in c(1:4))
{
  pdf(file=paste0(Sys.Date(), "_PCA_all_dmrs", "PC1vsPC", (a+1), "_labled.pdf"), paper="a4r")
  population <- c(rep("HSC", times = 3), 
                    rep("MDP", times = 2), 
                    rep("CDP", times = 4),
                    rep("cMoP", times = 2),
                    rep("DC_CD11b", times = 3),
                    rep("DC_CD8a", times = 3), 
                    rep("PDC", times = 3),
                    rep("Mono", times = 3))
  
  
  population <- factor(population, levels = c("HSC", "MDP", "CDP", "cMoP", "DC_CD11b", "DC_CD8a", "PDC", "Mono"))
  
  #a = 2
  x <- ggplot(as.data.frame(df), aes(PC1, df[,a+1], colour = population)) +
    geom_point() +
    xlab(paste0("PC1 (", round(PCA$eig[1,2], digits = 1), " %)")) +
    ylab(paste0("PC", a+1 , " (", round(PCA$eig[a+1, 2], digits = 1), " %)")) +
    theme_bw() +
    scale_color_manual(values = c(HSC = "#1b9e77", 
                                  MDP="#d95f02", 
                                  CDP="#7570b3", 
                                  cMoP="#e7298a", 
                                  DC_CD11b="#66a61e",
                                  DC_CD8a="#e6ab02",
                                  PDC="#a6761d", 
                                  Mono="#666666")) +
    theme(legend.title=element_blank()) +
    ggtitle(paste0("PCA - all DC dmrs ", "PC", 1, " vs PC", a+1)) +
    geom_text(aes(label=base::rownames(df)), hjust = -0.1, vjust = -0.1, size = 3) + 
    theme(legend.position="none")
  
  print(x)
  dev.off()
}


# PC2 on x-axis

setwd("/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/PCA")

for (a in c(2:4))
{
  pdf(file=paste0(Sys.Date(), "_PCA_all_dmrs", "PC2vsPC", (a+1), "_labled.pdf"), paper="a4r")
  population <- c(rep("HSC", times = 3), 
                  rep("MDP", times = 2), 
                  rep("CDP", times = 4),
                  rep("cMoP", times = 2),
                  rep("DC_CD11b", times = 3),
                  rep("DC_CD8a", times = 3), 
                  rep("PDC", times = 3),
                  rep("Mono", times = 3))
  
  
  population <- factor(population, levels = c("HSC", "MDP", "CDP", "cMoP", "DC_CD11b", "DC_CD8a", "PDC", "Mono"))
  
  #a = 2
  x <- ggplot(as.data.frame(df), aes(PC2, df[,a+1], colour = population)) +
    geom_point() +
    xlab(paste0("PC2 (", round(PCA$eig[2,2], digits = 1), " %)")) +
    ylab(paste0("PC", a+1 , " (", round(PCA$eig[a+1, 2], digits = 1), " %)")) +
    theme_minimal() +
    scale_color_manual(values = c(HSC = "#1b9e77", 
                                  MDP="#d95f02", 
                                  CDP="#7570b3", 
                                  cMoP="#e7298a", 
                                  DC_CD11b="#66a61e",
                                  DC_CD8a="#e6ab02",
                                  PDC="#a6761d", 
                                  Mono="#666666")) +
    theme(legend.title=element_blank()) +
    ggtitle(paste0("PCA - all DC dmrs ", "PC", 2, " vs PC", a+1)) +
  geom_text(aes(label=base::rownames(df)), hjust = -0.1, vjust = -0.1, size = 3) + 
    theme(legend.position="none")
  
  print(x)
  dev.off()
}


# PC3 on x-axis

setwd("/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/PCA")

for (a in c(3:4))
{
  pdf(file=paste0(Sys.Date(), "_PCA_all_dmrs", "PC3vsPC", (a+1), "_labled.pdf"), paper="a4r")
  population <- c(rep("HSC", times = 3), 
                  rep("MDP", times = 2), 
                  rep("CDP", times = 4),
                  rep("cMoP", times = 2),
                  rep("DC_CD11b", times = 3),
                  rep("DC_CD8a", times = 3), 
                  rep("PDC", times = 3),
                  rep("Mono", times = 3))
  
  population <- factor(population, levels = c("HSC", "MDP", "CDP", "cMoP", "DC_CD11b", "DC_CD8a", "PDC", "Mono"))
  
  #a = 2
  x <- ggplot(as.data.frame(df), aes(PC2, df[,a+1], colour = population)) +
    geom_point() +
    xlab(paste0("PC3 (", round(PCA$eig[3,2], digits = 1), " %)")) +
    ylab(paste0("PC", a+1 , " (", round(PCA$eig[a+1, 2], digits = 1), " %)")) +
    theme_minimal() +
    scale_color_manual(values = c(HSC = "#1b9e77", 
                                  MDP="#d95f02", 
                                  CDP="#7570b3", 
                                  cMoP="#e7298a", 
                                  DC_CD11b="#66a61e",
                                  DC_CD8a="#e6ab02",
                                  PDC="#a6761d", 
                                  Mono="#666666")) +
    theme(legend.title=element_blank()) +
    ggtitle(paste0("PCA - all hierarchy dmrs ", "PC", 3, " vs PC", a+1)) +
    geom_text(aes(label=base::rownames(df)), hjust = -0.1, vjust = -0.1, size = 3)
  
  print(x)
  dev.off()
}

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
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] FactoMineR_2.3       RColorBrewer_1.1-2   circlize_0.4.10      ComplexHeatmap_2.2.0 forcats_0.5.0       
# [6] stringr_1.4.0        dplyr_0.8.5          purrr_0.3.4          readr_1.3.1          tidyr_1.1.0         
# [11] tibble_3.0.1         ggplot2_3.3.0        tidyverse_1.3.0     
# 
# loaded via a namespace (and not attached):
#   [1] ggrepel_0.8.2        Rcpp_1.0.5           lubridate_1.7.9      lattice_0.20-41      png_0.1-7           
# [6] assertthat_0.2.1     digest_0.6.25        R6_2.4.1             cellranger_1.1.0     backports_1.1.10    
# [11] reprex_0.3.0         httr_1.4.2           pillar_1.4.4         GlobalOptions_0.1.2  rlang_0.4.7         
# [16] readxl_1.3.1         rstudioapi_0.11      GetoptLong_1.0.2     labeling_0.3         munsell_0.5.0       
# [21] broom_0.5.6          compiler_3.6.2       modelr_0.1.8         pkgconfig_2.0.3      shape_1.4.5         
# [26] flashClust_1.01-2    tidyselect_1.1.0     fansi_0.4.1          crayon_1.3.4         dbplyr_1.4.3        
# [31] withr_2.3.0          MASS_7.3-53          leaps_3.1            nlme_3.1-149         jsonlite_1.7.1      
# [36] gtable_0.3.0         lifecycle_0.2.0      DBI_1.1.0            magrittr_1.5         scales_1.1.1        
# [41] cli_2.0.2            stringi_1.5.3        farver_2.0.3         fs_1.5.0             scatterplot3d_0.3-41
# [46] xml2_1.3.2           ellipsis_0.3.1       generics_0.0.2       vctrs_0.3.0          rjson_0.2.20        
# [51] tools_3.6.2          glue_1.4.2           hms_0.5.3            parallel_3.6.2       clue_0.3-57         
# [56] colorspace_1.4-1     cluster_2.1.0        rvest_0.3.5          haven_2.2.0