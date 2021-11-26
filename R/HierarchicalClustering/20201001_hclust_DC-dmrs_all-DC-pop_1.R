# methylation in DC pop - hierarchical clustering
# 20201001
# author: Dr. Sina Stäble

# load libraries
library(tidyverse)
library(dynamicTreeCut)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(grDevices)


# load methylation df all DC pop / DC DMRs ----------------------------------------------------------------------------
    meth <- readRDS("/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/DMR calling/meth-levels_merged-dmrs/2020-10-01_beta-value_DC-DMRs_allChr_pop.rds")
    
# calculate z-score  ----------------------------------------------------------------------------

    meth.df <- meth %>%
      mutate(dmr = paste0("chr", chr, ":", start, "-", end)) %>%
      select(., -chr, -start, -end) %>%
      column_to_rownames(., var = "dmr") 
    
    meth.levels.t <- scale(t(meth.df))
    
    z_score <- meth.levels.t %>%
      t(.)

    z_score.m <- z_score %>%
      as.matrix(.)
    
    #set.seed(123)
    
    #z_score.m <- as.data.frame(z_score.m) %>% 
    #   sample_n(., 10000) %>%
    #   as.matrix(.)

# Distance matrix ---------------------------------------------------------------------------
    
    distances <- dist(z_score.m, method = "euclidean")
    out.dir <- "/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/Hierarchical_clustering//"
    saveRDS(distances, file = file.path(out.dir, paste0(Sys.Date(), "_distances.rds")))
    
    
# Hierarchical clustering with hclust
    dendro <- hclust(distances, method = "ward.D")
    plot(dendro)
    out.dir <- "/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/Hierarchical_clustering//"
    saveRDS(dendro, file = file.path(out.dir, paste0(Sys.Date(), "_dendro.rds")))
    
    
# Identify clusters with dynamicTreeCut
    distM=as.matrix(distances)
    
    clusters <- cutreeHybrid(dendro, distM,
                             cutHeight = NULL, minClusterSize = 0.005*49601, deepSplit = 1,
                             maxCoreScatter = NULL, minGap = 0.1,
                             maxAbsCoreScatter = NULL, minAbsGap = NULL,
                             minSplitHeight = NULL, minAbsSplitHeight = NULL,
                             externalBranchSplitFnc = NULL, minExternalSplit = NULL,
                             externalSplitOptions = list(),
                             externalSplitFncNeedsDistance = NULL,
                             assumeSimpleExternalSpecification = TRUE,
                             pamStage = FALSE, pamRespectsDendro = TRUE,
                             useMedoids = FALSE,
                             respectSmallClusters = TRUE,
                             verbose = 2, indent = 0)
    
# cbind clusters and methylation levels per population
    
    cluster <- clusters$labels
    cluster
    rownames(z_score.m) <- cluster
    
# save dataframe with DMRs and cluster id    
    
    dmrs_plus_cluster_id <- meth %>%
      cbind(., cluster) %>%
      arrange(., cluster)
    
    setwd("/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/Hierarchical_clustering")
    saveRDS(dmrs_plus_cluster_id, file = "cluster-id__hclust_z-score_heatmap_euclidean_deepSplit=1_minGap=0.1.rds")

# save DMRs plus Cluster Anno as bed. file -----------------------------------------------------    
    
DC_dmrs <- dmrs_plus_cluster_id %>%
      select(chr, start, end, cluster) %>%
      mutate(chr = paste0("chr", chr))
    
out.dir <- "/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/Hierarchical_clustering/"
write.table(DC_dmrs, file = paste0(out.dir, Sys.Date(), "_DC-dmrs_cluster-ids.bed"), quote = F, sep = "\t", row.names = F, col.names = F)
    
    
    