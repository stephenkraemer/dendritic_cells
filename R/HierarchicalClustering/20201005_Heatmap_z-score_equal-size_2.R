# methylation heatmap DC project
# 2020-10-05
# author: Dr. Sina Stäble

# load libraries
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)


# load cluster results with beta value matrix ------------------------------------------

dmrs_cluster_beta <- readRDS("/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/Hierarchical_clustering/cluster-id__hclust_z-score_heatmap_euclidean_deepSplit=1_minGap=0.1.rds")


# calculate z-score  ------------------------------------------

dmrs_cluster <- dmrs_cluster_beta %>%
  select(chr, start, end, cluster)

z_score.df <- dmrs_cluster_beta %>%
  select(-chr, -start, -end, -cluster) %>%
  as.matrix(.) %>%
  t(.) %>%
  scale(.) %>%
  t(.) %>%
  as.data.frame(.) %>%
  cbind(dmrs_cluster, .)


# sample random 350 regions per cluster ------------------------------------------

set.seed(100)
sampled_number = 350

z_score.df_sampled <- z_score.df %>%
  dplyr::group_by(cluster) %>%
  dplyr::sample_n(size = sampled_number) %>%
  dplyr::arrange(., cluster)

z_score.m_sampled <- z_score.df_sampled %>%
  ungroup(cluster) %>%
  select(-chr, -start, -end, -cluster) %>%
  as.matrix(.) 

# plot heatmap with ComplexHeatmap  ------------------------------------------

max(z_score.m_sampled, na.rm=T)
min(z_score.m_sampled, na.rm=T)

col_fun = colorRamp2(c(-2.5, 0, 2.5), c("#0571b0", "#f7f7f7", "#ca0020"))
col_fun(seq(-3, 3))

anno.cluster <- data.frame(cluster = c(rep("1", sampled_number), 
                             rep("2", sampled_number), 
                             rep("3", sampled_number), 
                             rep("4", sampled_number), 
                             rep("5", sampled_number), 
                             rep("6", sampled_number), 
                             rep("7", sampled_number), 
                             rep("8", sampled_number), 
                             rep("9", sampled_number)#, 
                             #rep("10", sampled_number), 
                             #rep("11", sampled_number), 
                             #rep("12", sampled_number), 
                             #rep("13", sampled_number), 
                             #rep("14", sampled_number)
                             ))

ha.cluster <- HeatmapAnnotation(df = anno.cluster, which = "row", width = unit(0.4, "cm"), show_legend = T)

split <- as.factor(c(rep("1", sampled_number), 
                     rep("2", sampled_number), 
                     rep("3", sampled_number), 
                     rep("4", sampled_number), 
                     rep("5", sampled_number), 
                     rep("6", sampled_number), 
                     rep("7", sampled_number), 
                     rep("8", sampled_number), 
                     rep("9", sampled_number)#, 
                     #rep("10", sampled_number), 
                     #rep("11", sampled_number), 
                     #rep("12", sampled_number), 
                     #rep("13", sampled_number), 
                     #rep("14", sampled_number)
                     )) 

h <- Heatmap(z_score.m_sampled, 
             col = col_fun, 
             cluster_rows = F, cluster_columns = T, 
             #clustering_distance_columns = "pearson",
             name = "z-score",
             split = split, gap = unit(2, "mm"),
             use_raster = T, raster_device = "png",
             column_title = "deepSplit=1_minGap=0.1")
             
             h + ha.cluster
             
             setwd("/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/Hierarchical_clustering")
             pdf(file=paste0(Sys.Date(), "_DC-Heatmap_", sampled_number, "sampled_euclidean_deepSplit=1_minGap=0.1.pdf"), paper="a4r", width = 5, height = 8)
             h + ha.cluster
             
             dev.off()
             
             
             
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
             #   [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
             # 
             # other attached packages:
             #   [1] RColorBrewer_1.1-2   circlize_0.4.10      ComplexHeatmap_2.2.0 forcats_0.5.0       
             # [5] stringr_1.4.0        dplyr_0.8.5          purrr_0.3.4          readr_1.3.1         
             # [9] tidyr_1.1.0          tibble_3.0.1         ggplot2_3.3.0        tidyverse_1.3.0     
             # 
             # loaded via a namespace (and not attached):
             #   [1] shape_1.4.5         GetoptLong_1.0.2    tidyselect_1.1.0    haven_2.2.0        
             # [5] lattice_0.20-41     colorspace_1.4-1    vctrs_0.3.0         generics_0.0.2     
             # [9] rlang_0.4.7         pillar_1.4.4        withr_2.3.0         glue_1.4.2         
             # [13] DBI_1.1.0           dbplyr_1.4.3        modelr_0.1.8        readxl_1.3.1       
             # [17] lifecycle_0.2.0     munsell_0.5.0       gtable_0.3.0        cellranger_1.1.0   
             # [21] rvest_0.3.5         GlobalOptions_0.1.2 parallel_3.6.2      fansi_0.4.1        
             # [25] broom_0.5.6         Rcpp_1.0.5          backports_1.1.10    scales_1.1.1       
             # [29] jsonlite_1.7.1      fs_1.5.0            rjson_0.2.20        hms_0.5.3          
             # [33] png_0.1-7           stringi_1.5.3       clue_0.3-57         cli_2.0.2          
             # [37] tools_3.6.2         magrittr_1.5        cluster_2.1.0       crayon_1.3.4       
             # [41] pkgconfig_2.0.3     ellipsis_0.3.1      xml2_1.3.2          reprex_0.3.0       
             # [45] lubridate_1.7.9     assertthat_0.2.1    httr_1.4.2          rstudioapi_0.11    
             # [49] R6_2.4.1            nlme_3.1-149        compiler_3.6.2             
             # 