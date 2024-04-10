setwd("/data/public/fmosko/")
.libPaths("/data/public/fmosko/miniconda/envs/psupertime/lib/R/library")

library(Matrix)
library(SeuratObject)
library(Seurat)
library(psupertime)

sc_exp_obj = readRDS("/data/public/fmosko/RWorkspace/sc_exp_obj.rds")
time_points_0to24 = readRDS("/data/public/fmosko/RWorkspace/time_points_0to24.rds")

pseudo_obj <- psupertime(x = sc_exp_obj,
                         y = factor(time_points_0to24,
                                    levels = c("N2i", "N12", "N24")),
                         sel_genes = "all")

saveRDS(pseudo_obj, "/data/public/fmosko/RWorkspace/pseudo_obj.rds")