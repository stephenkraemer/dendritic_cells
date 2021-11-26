### Methylation levels ###

# 24.09.2020
# Dr. Sina Stäble

## load libraries

library(tidyverse)
library(plyr)
library(tidyselect)
library(data.table)
library(readr)
library(ggforce) # only works with R 3.5.1!!!

### sample 1000000 CpGs

set.seed(123)
rows <- sample(1:20383641, 1000000)

### hsc

## load and merge mcall files

pop <- "hsc"

setwd("/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/meth_calls/")

HSC <- list.files(pattern = paste0("mcalls_",pop)) %>%
  lapply(read_delim, "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  lapply(function(s) s[rows,]) %>% 
  reduce(full_join, by = c("#chrom", "start", "end", "motif", "score", "strand")) %>%
  select(-starts_with("beta_value")) %>%
  mutate(n_meth_all = rowSums(select(., contains("n_meth")))) %>%
  mutate(n_total_all = rowSums(select(., contains("n_total")))) %>%
  select("#chrom", "start", "end", "motif", "score", "strand", ends_with("_all")) %>%
  mutate(beta_value_all = (n_meth_all/n_total_all)) %>%
  setnames("beta_value_all", paste0(pop)) %>%
  select(-starts_with("n_")) 

setwd("/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/methylation_levels/")
saveRDS(HSC, file = "global_methylation_HSC.rds")


## load and merge mcall files

pop <- "cdp"

setwd("/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/meth_calls/")

CDP <- list.files(pattern = paste0("mcalls_",pop)) %>%
  lapply(read_delim, "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  lapply(function(s) s[rows,]) %>% 
  reduce(full_join, by = c("#chrom", "start", "end", "motif", "score", "strand")) %>%
  select(-starts_with("beta_value")) %>%
  mutate(n_meth_all = rowSums(select(., contains("n_meth")))) %>%
  mutate(n_total_all = rowSums(select(., contains("n_total")))) %>%
  select("#chrom", "start", "end", "motif", "score", "strand", ends_with("_all")) %>%
  mutate(beta_value_all = (n_meth_all/n_total_all)) %>%
  setnames("beta_value_all", paste0(pop)) %>%
  select(-starts_with("n_")) 

setwd("/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/methylation_levels/")
saveRDS(CDP, file = "global_methylation_CDP.rds")


## load and merge mcall files

pop <- "cmop"

setwd("/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/meth_calls/")

CMOP <- list.files(pattern = paste0("mcalls_",pop)) %>%
  lapply(read_delim, "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  lapply(function(s) s[rows,]) %>% 
  reduce(full_join, by = c("#chrom", "start", "end", "motif", "score", "strand")) %>%
  select(-starts_with("beta_value")) %>%
  mutate(n_meth_all = rowSums(select(., contains("n_meth")))) %>%
  mutate(n_total_all = rowSums(select(., contains("n_total")))) %>%
  select("#chrom", "start", "end", "motif", "score", "strand", ends_with("_all")) %>%
  mutate(beta_value_all = (n_meth_all/n_total_all)) %>%
  setnames("beta_value_all", paste0(pop)) %>%
  select(-starts_with("n_")) 

setwd("/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/methylation_levels/")
saveRDS(CMOP, file = "global_methylation_CMOP.rds")


## load and merge mcall files

pop <- "dc-cd8a"

setwd("/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/meth_calls/")

DCCD8A <- list.files(pattern = paste0("mcalls_",pop)) %>%
  lapply(read_delim, "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  lapply(function(s) s[rows,]) %>% 
  reduce(full_join, by = c("#chrom", "start", "end", "motif", "score", "strand")) %>%
  select(-starts_with("beta_value")) %>%
  mutate(n_meth_all = rowSums(select(., contains("n_meth")))) %>%
  mutate(n_total_all = rowSums(select(., contains("n_total")))) %>%
  select("#chrom", "start", "end", "motif", "score", "strand", ends_with("_all")) %>%
  mutate(beta_value_all = (n_meth_all/n_total_all)) %>%
  setnames("beta_value_all", paste0(pop)) %>%
  select(-starts_with("n_")) 

setwd("/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/methylation_levels/")
saveRDS(DCCD8A, file = "global_methylation_DCCD8A.rds")


## load and merge mcall files

pop <- "dc-cd11b"

setwd("/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/meth_calls/")

DCCD11B <- list.files(pattern = paste0("mcalls_",pop)) %>%
  lapply(read_delim, "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  lapply(function(s) s[rows,]) %>% 
  reduce(full_join, by = c("#chrom", "start", "end", "motif", "score", "strand")) %>%
  select(-starts_with("beta_value")) %>%
  mutate(n_meth_all = rowSums(select(., contains("n_meth")))) %>%
  mutate(n_total_all = rowSums(select(., contains("n_total")))) %>%
  select("#chrom", "start", "end", "motif", "score", "strand", ends_with("_all")) %>%
  mutate(beta_value_all = (n_meth_all/n_total_all)) %>%
  setnames("beta_value_all", paste0(pop)) %>%
  select(-starts_with("n_")) 

setwd("/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/methylation_levels/")
saveRDS(DCCD11B, file = "global_methylation_DCCD11B.rds")


## load and merge mcall files

pop <- "mdp"

setwd("/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/meth_calls/")

MDP <- list.files(pattern = paste0("mcalls_",pop)) %>%
  lapply(read_delim, "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  lapply(function(s) s[rows,]) %>% 
  reduce(full_join, by = c("#chrom", "start", "end", "motif", "score", "strand")) %>%
  select(-starts_with("beta_value")) %>%
  mutate(n_meth_all = rowSums(select(., contains("n_meth")))) %>%
  mutate(n_total_all = rowSums(select(., contains("n_total")))) %>%
  select("#chrom", "start", "end", "motif", "score", "strand", ends_with("_all")) %>%
  mutate(beta_value_all = (n_meth_all/n_total_all)) %>%
  setnames("beta_value_all", paste0(pop)) %>%
  select(-starts_with("n_")) 

setwd("/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/methylation_levels/")
saveRDS(MDP, file = "global_methylation_MDP.rds")


## load and merge mcall files

pop <- "monos"

setwd("/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/meth_calls/")

MONOS <- list.files(pattern = paste0("mcalls_",pop)) %>%
  lapply(read_delim, "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  lapply(function(s) s[rows,]) %>% 
  reduce(full_join, by = c("#chrom", "start", "end", "motif", "score", "strand")) %>%
  select(-starts_with("beta_value")) %>%
  mutate(n_meth_all = rowSums(select(., contains("n_meth")))) %>%
  mutate(n_total_all = rowSums(select(., contains("n_total")))) %>%
  select("#chrom", "start", "end", "motif", "score", "strand", ends_with("_all")) %>%
  mutate(beta_value_all = (n_meth_all/n_total_all)) %>%
  setnames("beta_value_all", paste0(pop)) %>%
  select(-starts_with("n_")) 

setwd("/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/methylation_levels/")
saveRDS(MONOS, file = "global_methylation_MONOS.rds")


## load and merge mcall files

pop <- "pdc"

setwd("/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/meth_calls/")

PDC <- list.files(pattern = paste0("mcalls_",pop)) %>%
  lapply(read_delim, "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  lapply(function(s) s[rows,]) %>% 
  reduce(full_join, by = c("#chrom", "start", "end", "motif", "score", "strand")) %>%
  select(-starts_with("beta_value")) %>%
  mutate(n_meth_all = rowSums(select(., contains("n_meth")))) %>%
  mutate(n_total_all = rowSums(select(., contains("n_total")))) %>%
  select("#chrom", "start", "end", "motif", "score", "strand", ends_with("_all")) %>%
  mutate(beta_value_all = (n_meth_all/n_total_all)) %>%
  setnames("beta_value_all", paste0(pop)) %>%
  select(-starts_with("n_")) 

setwd("/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/methylation_levels/")
saveRDS(PDC, file = "global_methylation_PDC.rds")


### merge coverage of all populations into one df

all_pop <- Reduce(full_join, list(HSC, MDP, CDP, CMOP, DCCD8A, DCCD11B, PDC, MONOS))

saveRDS(all_pop, file = "methylation-levels_all_pop.rds")

pop_gather <- all_pop %>%
  drop_na() %>%
  gather("hsc", "mdp", "cdp", "cmop", "dc-cd8a", "dc-cd11b", "pdc", "monos", 
        key = "pop", value = "meth_level") 

### make methylation violin plot 
setwd("/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/methylation_levels")
pdf(file = paste0(file = paste0(Sys.Date(), "_violin_methylation-level_all-pop_adjust-4.pdf")), width = 10, height = 6, paper = "a4r")

violin <- ggplot(data = pop_gather) +
  geom_violin(mapping = aes(x = factor(pop, levels=unique(pop)), y = meth_level), draw_quantiles = 0.5, adjust = 4) +
  #geom_violin(mapping = aes(x = factor(pop, levels=unique(pop)), y = meth_level), adjust = 3) +
  #geom_boxplot(mapping = aes(x = factor(pop, levels=unique(pop)), y = meth_level), width=0.1, outlier.shape = NA, coef = 0) +
  theme(axis.text.x = element_text(angle = -50, size = 10, hjust=-0.05), 
        axis.text.y = element_text(size = 10)) +
  ylab("CpG methylation") +
  xlab("") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5, fill=NA)) +
  scale_y_continuous(expand = c(0.01,0), breaks = seq(0, 1, 0.25)) 

violin

dev.off()


### make sina plot 

rows_sina <- sample(1:100000, 1000)

all_pop_sina <- all_pop %>%
  slice(rows_sina)

pop_gather_sina <- all_pop_sina %>%
  gather("hsc", "mdp", "cdp", "cmop", "dc-cd8a", "dc-cd11b", "pdc", "monos", 
         key = "pop", value = "meth_level") 

pdf(file = paste0(file = paste0(Sys.Date(), "_sina_methylation-level_all_pop_adjust=4.pdf")), width = 10, height = 6, paper = "a4r")

sina <- ggplot(data = pop_gather_sina) +
  geom_sina(mapping = aes(x = factor(pop, levels=unique(pop)), y = meth_level), size = 0.1, alpha = 0.1, adjust = 4, method = "density", maxwidth = NULL) +
  #geom_jitter(mapping = aes(x = factor(pop, levels=unique(pop)), y = meth_level), size = 0.1, alpha = 0.1) +
  theme(axis.text.x = element_text(angle = 90, size = 10, hjust=0.5), 
        axis.text.y = element_text(size = 10)) +
  ylab("CpG methylation") +
  xlab("") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5, fill=NA)) +
  scale_y_continuous(expand = c(0.01,0), breaks = seq(0, 1, 0.25)) 

sina

dev.off()


### make methylation violin plot and plot mean
setwd("/icgc/dkfzlsdf/analysis/C010/zimmerms/Collaboration_RosenbauerLab/methylation_levels")
pdf(file = paste0(file = paste0(Sys.Date(), "_violin_methylation-level_all-pop_adjust-3_mean.pdf")), width = 10, height = 6, paper = "a4r")

violin_mean <- ggplot(pop_gather, aes(x=factor(pop, levels=unique(pop)), y=meth_level)) + 
  geom_violin(trim=T, adjust = 3) +
  stat_summary(fun="mean",
               geom="point", color="black", size = 2) +
  theme(axis.text.x = element_text(angle = -50, size = 10, hjust=-0.05), 
        axis.text.y = element_text(size = 10)) +
  ylab("CpG methylation") +
  xlab("") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5, fill=NA)) +
  scale_y_continuous(expand = c(0.01,0), breaks = seq(0, 1, 0.25)) 

violin_mean

dev.off()


# calculate p-value with wilcoxon test

median(all_pop$hsc, na.rm = T) #0.9354839
median(all_pop$mdp, na.rm = T) #0.9285714
median(all_pop$cdp, na.rm = T) #0.925
median(all_pop$cmop, na.rm = T) #0.9090909
median(all_pop$`dc-cd8a`, na.rm = T) # 0.9090909
median(all_pop$`dc-cd11b`, na.rm = T) #0.9032258
median(all_pop$pdc, na.rm = T) # 0.9090909
median(all_pop$monos, na.rm = T) # 0.8888889


mean(all_pop$hsc, na.rm = T) #0.8205398
mean(all_pop$cmop, na.rm = T) #0.7881458
mean(all_pop$monos, na.rm = T) #0.7607818
mean(all_pop$mdp, na.rm = T) #0.8087739
mean(all_pop$cdp, na.rm = T) #0.8048729
mean(all_pop$`dc-cd11b`, na.rm = T) #0.7826421
mean(all_pop$`dc-cd8a`, na.rm = T) #0.7843942
mean(all_pop$pdc, na.rm = T) #0.7906624


wilcox.test(all_pop$hsc, all_pop$mdp, paired = T) # p-value < 2.2e-16
wilcox.test(all_pop$hsc, all_pop$cdp, paired = T) # p-value < 2.2e-16
wilcox.test(all_pop$hsc, all_pop$cmop, paired = T) # p-value < 2.2e-16
wilcox.test(all_pop$hsc, all_pop$`dc-cd11b`, paired = T) # p-value < 2.2e-16
wilcox.test(all_pop$hsc, all_pop$`dc-cd8a`, paired = T) # p-value < 2.2e-16
wilcox.test(all_pop$hsc, all_pop$pdc, paired = T) # p-value < 2.2e-16
wilcox.test(all_pop$hsc, all_pop$monos, paired = T) # p-value < 2.2e-16
