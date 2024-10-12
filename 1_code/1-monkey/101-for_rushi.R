no_source()
####This should be run in workstation

library(tidymass)
library(tidyverse)

rm(list = ls())

setwd(r4projects::get_project_wd())
source("1-code/tools.R")

setwd("3-data_analysis/monkey_rushi/")

result <-
  readr::read_rds("all_resLFCn.rds")

result1 <- result[[1]] %>%
  as.data.frame()

result2 <- result[[2]] %>%
  as.data.frame()

result3 <- result[[3]] %>%
  as.data.frame()

###GET gene id
library(clusterProfiler)
library(org.Hs.eg.db)
# ENSEMBL <-
#   bitr(
#     geneID = rownames(result1),
#     fromType = "SYMBOL",
#     toType = "ENSEMBL",
#     OrgDb = org.Hs.eg.db
#   )
#
# UNIPROT <-
#   bitr(
#     geneID = rownames(result1),
#     fromType = "SYMBOL",
#     toType = "UNIPROT",
#     OrgDb = org.Hs.eg.db
#   )
#
# ENTREZID <-
#   bitr(
#     geneID = rownames(result1),
#     fromType = "SYMBOL",
#     toType = "ENTREZID",
#     OrgDb = org.Hs.eg.db
#   )
#
# gene_id <-
#   ENSEMBL %>%
#   dplyr::distinct(SYMBOL, .keep_all = TRUE) %>%
#   dplyr::left_join(UNIPROT %>%
#                      dplyr::distinct(SYMBOL, .keep_all = TRUE),
#                    by = "SYMBOL") %>%
#   dplyr::left_join(ENTREZID %>%
#                      dplyr::distinct(SYMBOL, .keep_all = TRUE),
#                    by = "SYMBOL")
#
# save(gene_id, file = "gene_id")
load("gene_id")
dim(result1)

up1 <-
  result1 %>%
  dplyr::filter(padj < 0.05 & log2FoldChange > 0)
down1 <-
  result1 %>%
  dplyr::filter(padj < 0.05 & log2FoldChange < 0)
no1 <-
  result1 %>%
  dplyr::filter(padj > 0.05)

up2 <-
  result2 %>%
  dplyr::filter(padj < 0.05 & log2FoldChange > 0)
down2 <-
  result2 %>%
  dplyr::filter(padj < 0.05 & log2FoldChange < 0)
no2 <-
  result2 %>%
  dplyr::filter(padj > 0.05)

up3 <-
  result3 %>%
  dplyr::filter(padj < 0.05 & log2FoldChange > 0)
down3 <-
  result3 %>%
  dplyr::filter(padj < 0.05 & log2FoldChange < 0)
no3 <-
  result3 %>%
  dplyr::filter(padj > 0.05)

library(clusterProfiler)
library(org.Hs.eg.db)

# dim(up1)
# dim(up2)
# dim(up3)
# dim(down1)
# dim(down2)
# dim(down3)
# dim(no1)
# dim(no2)
# dim(no3)
#
# gene_list <-
#   gene_id$ENTREZID[match(rownames(up1), gene_id$SYMBOL)]
#
# gene_list <-
#   gene_list[!is.na(gene_list)]
#
# up_go_result1 <-
#   clusterProfiler::enrichGO(
#     gene = gene_list,
#     OrgDb = org.Hs.eg.db,
#     keyType = "ENTREZID",
#     ont = "ALL",
#     pAdjustMethod = "BH",
#     pvalueCutoff = 0.05,
#     qvalueCutoff = 0.05
#   )
#
# save(up_go_result1, file = "up_go_result1")
#
# gene_list <-
#   gene_id$ENTREZID[match(rownames(down1), gene_id$SYMBOL)]
#
# gene_list <-
#   gene_list[!is.na(gene_list)]
#
# down_go_result1 <-
#   clusterProfiler::enrichGO(
#     gene = gene_list,
#     OrgDb = org.Hs.eg.db,
#     keyType = "ENTREZID",
#     ont = "ALL",
#     pAdjustMethod = "BH",
#     pvalueCutoff = 0.05,
#     qvalueCutoff = 0.05
#   )
#
# save(down_go_result1, file = "down_go_result1")
#
#
# gene_list <-
#   gene_id$ENTREZID[match(rownames(no1), gene_id$SYMBOL)]
#
# gene_list <-
#   gene_list[!is.na(gene_list)]
#
# no_go_result1 <-
#   clusterProfiler::enrichGO(
#     gene = gene_list,
#     OrgDb = org.Hs.eg.db,
#     keyType = "ENTREZID",
#     ont = "ALL",
#     pAdjustMethod = "BH",
#     pvalueCutoff = 0.05,
#     qvalueCutoff = 0.05
#   )
#
# save(no_go_result1, file = "no_go_result1")
#
#
# gene_list <-
#   gene_id$ENTREZID[match(rownames(up2), gene_id$SYMBOL)]
#
# gene_list <-
#   gene_list[!is.na(gene_list)]
#
# up_go_result2 <-
#   clusterProfiler::enrichGO(
#     gene = gene_list,
#     OrgDb = org.Hs.eg.db,
#     keyType = "ENTREZID",
#     ont = "ALL",
#     pAdjustMethod = "BH",
#     pvalueCutoff = 0.05,
#     qvalueCutoff = 0.05
#   )
#
# save(up_go_result2, file = "up_go_result2")
#
# gene_list <-
#   gene_id$ENTREZID[match(rownames(down2), gene_id$SYMBOL)]
#
# gene_list <-
#   gene_list[!is.na(gene_list)]
#
# down_go_result2 <-
#   clusterProfiler::enrichGO(
#     gene = gene_list,
#     OrgDb = org.Hs.eg.db,
#     keyType = "ENTREZID",
#     ont = "ALL",
#     pAdjustMethod = "BH",
#     pvalueCutoff = 0.05,
#     qvalueCutoff = 0.05
#   )
#
# save(down_go_result2, file = "down_go_result2")
#
#
#
#
#
# gene_list <-
#   gene_id$ENTREZID[match(rownames(no2), gene_id$SYMBOL)]
#
# gene_list <-
#   gene_list[!is.na(gene_list)]
#
# no_go_result2 <-
#   clusterProfiler::enrichGO(
#     gene = gene_list,
#     OrgDb = org.Hs.eg.db,
#     keyType = "ENTREZID",
#     ont = "ALL",
#     pAdjustMethod = "BH",
#     pvalueCutoff = 0.05,
#     qvalueCutoff = 0.05
#   )
#
# save(no_go_result2, file = "no_go_result2")
#
#
# gene_list <-
#   gene_id$ENTREZID[match(rownames(up3), gene_id$SYMBOL)]
#
# gene_list <-
#   gene_list[!is.na(gene_list)]
#
# up_go_result3 <-
#   clusterProfiler::enrichGO(
#     gene = gene_list,
#     OrgDb = org.Hs.eg.db,
#     keyType = "ENTREZID",
#     ont = "ALL",
#     pAdjustMethod = "BH",
#     pvalueCutoff = 0.05,
#     qvalueCutoff = 0.05
#   )
#
# save(up_go_result3, file = "up_go_result3")
#
# gene_list <-
#   gene_id$ENTREZID[match(rownames(down3), gene_id$SYMBOL)]
#
# gene_list <-
#   gene_list[!is.na(gene_list)]
#
# down_go_result3 <-
#   clusterProfiler::enrichGO(
#     gene = gene_list,
#     OrgDb = org.Hs.eg.db,
#     keyType = "ENTREZID",
#     ont = "ALL",
#     pAdjustMethod = "BH",
#     pvalueCutoff = 0.05,
#     qvalueCutoff = 0.05
#   )
#
# save(down_go_result3, file = "down_go_result3")
#
#
#
# gene_list <-
#   gene_id$ENTREZID[match(rownames(no3), gene_id$SYMBOL)]
#
# gene_list <-
#   gene_list[!is.na(gene_list)]
#
# no_go_result3 <-
#   clusterProfiler::enrichGO(
#     gene = gene_list,
#     OrgDb = org.Hs.eg.db,
#     keyType = "ENTREZID",
#     ont = "ALL",
#     pAdjustMethod = "BH",
#     pvalueCutoff = 0.05,
#     qvalueCutoff = 0.05
#   )
#
# save(no_go_result3, file = "no_go_result3")

load("up_go_result1")
load("up_go_result2")
load("up_go_result3")

load("down_go_result1")
load("down_go_result2")
load("down_go_result3")

load("no_go_result1")
load("no_go_result2")
load("no_go_result3")

dim(up_go_result1@result)
dim(down_go_result1@result)
dim(no_go_result1@result)

dim(up_go_result2@result)
dim(down_go_result2@result)
dim(no_go_result2@result)

dim(up_go_result3@result)
dim(down_go_result3@result)
dim(no_go_result3@result)

# library(patchwork)
# clusterProfiler::dotplot(up_go_result1)
# clusterProfiler::dotplot(up_go_result2)
# clusterProfiler::dotplot(up_go_result3)
# 
# clusterProfiler::dotplot(down_go_result1)
# clusterProfiler::dotplot(down_go_result2)
# clusterProfiler::dotplot(down_go_result3)

# #####up_go_result1
# bp_sim <-
#   get_go_similarity_matrix(
#     go_id = up_go_result1@result$ID[up_go_result1@result$ONTOLOGY == "BP"],
#     ont = "BP",
#     measure = "Wang",
#     sim_cutoff = 0.7
#   )
# 
# cc_sim <-
#   get_go_similarity_matrix(
#     go_id = up_go_result1@result$ID[up_go_result1@result$ONTOLOGY == "CC"],
#     ont = "CC",
#     measure = "Wang",
#     sim_cutoff = 0.7
#   )
# 
# mf_sim <-
#   get_go_similarity_matrix(
#     go_id = up_go_result1@result$ID[up_go_result1@result$ONTOLOGY == "MF"],
#     ont = "MF",
#     measure = "Wang",
#     sim_cutoff = 0.7
#   )
# 
# go_information <-
#   up_go_result1@result
# 
# go_similarity_matrix <-
#   rbind(bp_sim,
#         cc_sim,
#         mf_sim)
# 
# merge_go_network(
#   go_information = go_information,
#   go_similarity_matrix = go_similarity_matrix,
#   path = "comparison1_up"
# )
# 
# load("comparison1_up/center_go_id")
# up_go_result1_new <-
#   up_go_result1
# 
# up_go_result1_new@result <-
#   up_go_result1_new@result %>%
#   dplyr::filter(ID %in% center_go_id)
# 
# plot <-
#   clusterProfiler::dotplot(up_go_result1_new)
# 
# ggsave(plot,
#        filename = "comparison1_up/dot_plot.pdf",
#        width = 7,
#        height = 7)
# 
# 
# #####down_go_result1
# bp_sim <-
#   get_go_similarity_matrix(
#     go_id = down_go_result1@result$ID[down_go_result1@result$ONTOLOGY == "BP"],
#     ont = "BP",
#     measure = "Wang",
#     sim_cutoff = 0.7
#   )
# 
# cc_sim <-
#   get_go_similarity_matrix(
#     go_id = down_go_result1@result$ID[down_go_result1@result$ONTOLOGY == "CC"],
#     ont = "CC",
#     measure = "Wang",
#     sim_cutoff = 0.7
#   )
# 
# mf_sim <-
#   get_go_similarity_matrix(
#     go_id = down_go_result1@result$ID[down_go_result1@result$ONTOLOGY == "MF"],
#     ont = "MF",
#     measure = "Wang",
#     sim_cutoff = 0.7
#   )
# 
# go_information <-
#   down_go_result1@result
# 
# go_similarity_matrix <-
#   rbind(bp_sim,
#         cc_sim,
#         mf_sim)
# 
# merge_go_network(
#   go_information = go_information,
#   go_similarity_matrix = go_similarity_matrix,
#   path = "comparison1_down"
# )
# 
# load("comparison1_down/center_go_id")
# down_go_result1_new <-
#   down_go_result1
# 
# down_go_result1_new@result <-
#   down_go_result1_new@result %>%
#   dplyr::filter(ID %in% center_go_id)
# 
# plot <-
#   clusterProfiler::dotplot(down_go_result1_new)
# 
# ggsave(plot,
#        filename = "comparison1_down/dot_plot.pdf",
#        width = 7,
#        height = 7)
# 
# 
# 
# #####up_go_result2
# bp_sim <-
#   get_go_similarity_matrix(
#     go_id = up_go_result2@result$ID[up_go_result2@result$ONTOLOGY == "BP"],
#     ont = "BP",
#     measure = "Wang",
#     sim_cutoff = 0.7
#   )
# 
# cc_sim <-
#   get_go_similarity_matrix(
#     go_id = up_go_result2@result$ID[up_go_result2@result$ONTOLOGY == "CC"],
#     ont = "CC",
#     measure = "Wang",
#     sim_cutoff = 0.7
#   )
# 
# mf_sim <-
#   get_go_similarity_matrix(
#     go_id = up_go_result2@result$ID[up_go_result2@result$ONTOLOGY == "MF"],
#     ont = "MF",
#     measure = "Wang",
#     sim_cutoff = 0.7
#   )
# 
# go_information <-
#   up_go_result2@result
# 
# go_similarity_matrix <-
#   rbind(bp_sim,
#         cc_sim,
#         mf_sim)
# 
# merge_go_network(
#   go_information = go_information,
#   go_similarity_matrix = go_similarity_matrix,
#   path = "comparison2_up"
# )
# 
# load("comparison2_up/center_go_id")
# up_go_result2_new <-
#   up_go_result2
# 
# up_go_result2_new@result <-
#   up_go_result2_new@result %>%
#   dplyr::filter(ID %in% center_go_id)
# 
# plot <-
#   clusterProfiler::dotplot(up_go_result2_new)
# 
# ggsave(plot,
#        filename = "comparison2_up/dot_plot.pdf",
#        width = 7,
#        height = 7)
# 
# 
# #####down_go_result2
# bp_sim <-
#   get_go_similarity_matrix(
#     go_id = down_go_result2@result$ID[down_go_result2@result$ONTOLOGY == "BP"],
#     ont = "BP",
#     measure = "Wang",
#     sim_cutoff = 0.7
#   )
# 
# cc_sim <-
#   get_go_similarity_matrix(
#     go_id = down_go_result2@result$ID[down_go_result2@result$ONTOLOGY == "CC"],
#     ont = "CC",
#     measure = "Wang",
#     sim_cutoff = 0.7
#   )
# 
# mf_sim <-
#   get_go_similarity_matrix(
#     go_id = down_go_result2@result$ID[down_go_result2@result$ONTOLOGY == "MF"],
#     ont = "MF",
#     measure = "Wang",
#     sim_cutoff = 0.7
#   )
# 
# go_information <-
#   down_go_result2@result
# 
# go_similarity_matrix <-
#   rbind(bp_sim,
#         cc_sim,
#         mf_sim)
# 
# merge_go_network(
#   go_information = go_information,
#   go_similarity_matrix = go_similarity_matrix,
#   path = "comparison2_down"
# )
# 
# load("comparison2_down/center_go_id")
# down_go_result2_new <-
#   down_go_result2
# 
# down_go_result2_new@result <-
#   down_go_result2_new@result %>%
#   dplyr::filter(ID %in% center_go_id)
# 
# plot <-
#   clusterProfiler::dotplot(down_go_result2_new)
# 
# ggsave(plot,
#        filename = "comparison2_down/dot_plot.pdf",
#        width = 7,
#        height = 7)
# 
# 
# 
# #####up_go_result3
# bp_sim <-
#   get_go_similarity_matrix(
#     go_id = up_go_result3@result$ID[up_go_result3@result$ONTOLOGY == "BP"],
#     ont = "BP",
#     measure = "Wang",
#     sim_cutoff = 0.7
#   )
# 
# cc_sim <-
#   get_go_similarity_matrix(
#     go_id = up_go_result3@result$ID[up_go_result3@result$ONTOLOGY == "CC"],
#     ont = "CC",
#     measure = "Wang",
#     sim_cutoff = 0.7
#   )
# 
# mf_sim <-
#   get_go_similarity_matrix(
#     go_id = up_go_result3@result$ID[up_go_result3@result$ONTOLOGY == "MF"],
#     ont = "MF",
#     measure = "Wang",
#     sim_cutoff = 0.7
#   )
# 
# go_information <-
#   up_go_result3@result
# 
# go_similarity_matrix <-
#   rbind(bp_sim,
#         cc_sim,
#         mf_sim)
# 
# merge_go_network(
#   go_information = go_information,
#   go_similarity_matrix = go_similarity_matrix,
#   path = "comparison3_up"
# )
# 
# load("comparison3_up/center_go_id")
# up_go_result3_new <-
#   up_go_result3
# 
# up_go_result3_new@result <-
#   up_go_result3_new@result %>%
#   dplyr::filter(ID %in% center_go_id)
# 
# plot <-
#   clusterProfiler::dotplot(up_go_result3_new)
# 
# ggsave(plot,
#        filename = "comparison3_up/dot_plot.pdf",
#        width = 7,
#        height = 7)
# 
# 
# #####down_go_result3
# bp_sim <-
#   get_go_similarity_matrix(
#     go_id = down_go_result3@result$ID[down_go_result3@result$ONTOLOGY == "BP"],
#     ont = "BP",
#     measure = "Wang",
#     sim_cutoff = 0.7
#   )
# 
# cc_sim <-
#   get_go_similarity_matrix(
#     go_id = down_go_result3@result$ID[down_go_result3@result$ONTOLOGY == "CC"],
#     ont = "CC",
#     measure = "Wang",
#     sim_cutoff = 0.7
#   )
# 
# mf_sim <-
#   get_go_similarity_matrix(
#     go_id = down_go_result3@result$ID[down_go_result3@result$ONTOLOGY == "MF"],
#     ont = "MF",
#     measure = "Wang",
#     sim_cutoff = 0.7
#   )
# 
# go_information <-
#   down_go_result3@result
# 
# go_similarity_matrix <-
#   rbind(bp_sim,
#         cc_sim,
#         mf_sim)
# 
# merge_go_network(
#   go_information = go_information,
#   go_similarity_matrix = go_similarity_matrix,
#   path = "comparison3_down"
# )
# 
# load("comparison3_down/center_go_id")
# down_go_result3_new <-
#   down_go_result3
# 
# down_go_result3_new@result <-
#   down_go_result3_new@result %>%
#   dplyr::filter(ID %in% center_go_id)
# 
# plot <-
#   clusterProfiler::dotplot(down_go_result3_new)
# 
# ggsave(plot,
#        filename = "comparison3_down/dot_plot.pdf",
#        width = 7,
#        height = 7)











#####no_go_result1
bp_sim <-
  get_go_similarity_matrix(
    go_id = no_go_result1@result$ID[no_go_result1@result$ONTOLOGY == "BP"],
    ont = "BP",
    measure = "Wang",
    sim_cutoff = 0.7
  )

cc_sim <-
  get_go_similarity_matrix(
    go_id = no_go_result1@result$ID[no_go_result1@result$ONTOLOGY == "CC"],
    ont = "CC",
    measure = "Wang",
    sim_cutoff = 0.7
  )

mf_sim <-
  get_go_similarity_matrix(
    go_id = no_go_result1@result$ID[no_go_result1@result$ONTOLOGY == "MF"],
    ont = "MF",
    measure = "Wang",
    sim_cutoff = 0.7
  )

go_information <-
  no_go_result1@result

go_similarity_matrix <-
  rbind(bp_sim,
        cc_sim,
        mf_sim)

merge_go_network(
  go_information = go_information,
  go_similarity_matrix = go_similarity_matrix,
  path = "comparison1_no"
)

load("comparison1_no/center_go_id")
no_go_result1_new <-
  no_go_result1

no_go_result1_new@result <-
  no_go_result1_new@result %>%
  dplyr::filter(ID %in% center_go_id)

plot <-
  clusterProfiler::dotplot(no_go_result1_new)
plot
ggsave(plot,
       filename = "comparison1_no/dot_plot.pdf",
       width = 7,
       height = 7)











#####no_go_result2
bp_sim <-
  get_go_similarity_matrix(
    go_id = no_go_result2@result$ID[no_go_result2@result$ONTOLOGY == "BP"],
    ont = "BP",
    measure = "Wang",
    sim_cutoff = 0.7
  )

cc_sim <-
  get_go_similarity_matrix(
    go_id = no_go_result2@result$ID[no_go_result2@result$ONTOLOGY == "CC"],
    ont = "CC",
    measure = "Wang",
    sim_cutoff = 0.7
  )

mf_sim <-
  get_go_similarity_matrix(
    go_id = no_go_result2@result$ID[no_go_result2@result$ONTOLOGY == "MF"],
    ont = "MF",
    measure = "Wang",
    sim_cutoff = 0.7
  )

go_information <-
  no_go_result2@result

go_similarity_matrix <-
  rbind(bp_sim,
        cc_sim,
        mf_sim)

merge_go_network(
  go_information = go_information,
  go_similarity_matrix = go_similarity_matrix,
  path = "comparison2_no"
)

load("comparison2_no/center_go_id")
no_go_result2_new <-
  no_go_result2

no_go_result2_new@result <-
  no_go_result2_new@result %>%
  dplyr::filter(ID %in% center_go_id)

plot <-
  clusterProfiler::dotplot(no_go_result2_new)
plot
ggsave(plot,
       filename = "comparison2_no/dot_plot.pdf",
       width = 7,
       height = 7)







#####no_go_result3
bp_sim <-
  get_go_similarity_matrix(
    go_id = no_go_result3@result$ID[no_go_result3@result$ONTOLOGY == "BP"],
    ont = "BP",
    measure = "Wang",
    sim_cutoff = 0.7
  )

cc_sim <-
  get_go_similarity_matrix(
    go_id = no_go_result3@result$ID[no_go_result3@result$ONTOLOGY == "CC"],
    ont = "CC",
    measure = "Wang",
    sim_cutoff = 0.7
  )

mf_sim <-
  get_go_similarity_matrix(
    go_id = no_go_result3@result$ID[no_go_result3@result$ONTOLOGY == "MF"],
    ont = "MF",
    measure = "Wang",
    sim_cutoff = 0.7
  )

go_information <-
  no_go_result3@result

go_similarity_matrix <-
  rbind(bp_sim,
        cc_sim,
        mf_sim)

merge_go_network(
  go_information = go_information,
  go_similarity_matrix = go_similarity_matrix,
  path = "comparison3_no"
)

load("comparison3_no/center_go_id")
no_go_result3_new <-
  no_go_result3

no_go_result3_new@result <-
  no_go_result3_new@result %>%
  dplyr::filter(ID %in% center_go_id)

plot <-
  clusterProfiler::dotplot(no_go_result3_new)
plot
ggsave(plot,
       filename = "comparison3_no/dot_plot.pdf",
       width = 7,
       height = 7)

