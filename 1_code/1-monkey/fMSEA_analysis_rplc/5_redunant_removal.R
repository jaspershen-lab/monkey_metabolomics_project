no_source()
####data preparation

library(tidymass)
library(tidyverse)

rm(list = ls())

setwd(r4projects::get_project_wd())
source("1-code/fMSEA/remove_redundancy4fmsea.R")
source("1-code/fMSEA/calculate_redundancy.R")
source("1-code/fMSEA/group_features.R")
source("1-code/tools.R")

setwd(r4projects::get_project_wd())

load("3-data_analysis/monkey_fMSEA_analysis_rplc/4_group_features/annotation_table")

dir.create("3-data_analysis/monkey_fMSEA_analysis_rplc/5_redundant_removal/")
setwd("3-data_analysis/monkey_fMSEA_analysis_rplc/5_redundant_removal/")

# ####The adduct distributation
# library(plyr)
# 
# compound_class <-
#   annotation_table %>%
#   plyr::dlply(.variables = .(compound_class))
# 
# length(compound_class)
# length(unique(annotation_table$HMDB.ID))
# 
# ####compound class size
# compound_class_size <-
#   compound_class %>%
#   purrr::map(nrow) %>%
#   unlist() %>%
#   data.frame(compound_class_size = .) %>%
#   dplyr::count(compound_class_size)
# 
# plot <-
#   compound_class_size %>%
#   ggplot(aes(compound_class_size, n)) +
#   geom_bar(stat = "identity") +
#   theme_bw() +
#   theme(
#     axis.title = element_text(size = 13),
#     axis.text = element_text(size = 12),
#     plot.background  = element_rect(fill = "transparent",
#                                     color = NA),
#     panel.background = element_rect(fill = "transparent",
#                                     color = NA),
#     legend.background = element_rect(fill = "transparent",
#                                      color = NA),
#     axis.text.x = element_text(
#       size = 12,
#       angle = 45,
#       hjust = 1,
#       vjust = 1
#     ),
#     panel.grid.minor = element_blank()
#   ) +
#   scale_x_continuous(
#     breaks = c(1:27),
#     name = c(1:27),
#     expand = expansion(mult = c(0.02, 0.02))
#   ) +
#   labs(x = "Score") +
#   ggforce::facet_zoom(
#     xlim = c(5.9, 20),
#     ylim = c(0, 1500),
#     horizontal = FALSE,
#     zoom.size = 1
#   )
# 
# plot
# 
# # ggsave(
# #   plot,
# #   filename = "compound_class_size_distributation.pdf",
# #   width = 10,
# #   height = 7,
# #   bg = "transparent"
# # )
# 
# 
# ###compound class polarity. positive, negative and mix
# compound_class_polarity <-
#   purrr::map(compound_class, function(x) {
#     polarity <- unique(x$polarity)
#     if (length(polarity) == 2) {
#       return("mix")
#     } else{
#       return(polarity)
#     }
#   }) %>%
#   unlist() %>%
#   data.frame(polarity = .) %>%
#   dplyr::count(polarity)
# 
# library(ggpubr)
# 
# compound_class_polarity$value <-
#   compound_class_polarity$n * 100 / sum(compound_class_polarity$n)
# 
# labs <- paste0(
#   compound_class_polarity$polarity,
#   " (",
#   round(compound_class_polarity$value, 2),
#   "%)"
# )
# 
# plot <-
#   ggpie(
#     data = compound_class_polarity,
#     label = labs,
#     lab.pos = "in",
#     x = "n",
#     lab.font = "white",
#     fill = "polarity",
#     color = "white",
#     palette = polarity_color
#   )
# plot
# 
# # ggsave(
# #   plot,
# #   filename = "compound_class_polarity.pdf",
# #   width = 7,
# #   height = 7,
# #   bg = "transparent"
# # )
# 
# ###adduct distributation
# ####only consider the Level 1 and Level 2
# temp_data <-
#   annotation_table %>%
#   dplyr::filter(Level != 3) %>%
#   dplyr::filter(isotopes == "[M]") %>%
#   select(Adduct, polarity) %>%
#   dplyr::count(Adduct, polarity)
# # dplyr::mutate(value = n * 100 / sum(n))
# 
# unique_adduct <-
#   unique(annotation_table$Adduct)
# 
# temp_data_pos <-
#   temp_data %>%
#   dplyr::filter(polarity == "positive") %>%
#   dplyr::mutate(value = n * 100 / sum(n)) %>%
#   dplyr::arrange(value) %>%
#   dplyr::mutate(Adduct = factor(Adduct, levels = Adduct))
# 
# temp_data_neg <-
#   temp_data %>%
#   dplyr::filter(polarity == "negative") %>%
#   dplyr::mutate(value = n * 100 / sum(n)) %>%
#   dplyr::arrange(value) %>%
#   dplyr::mutate(Adduct = factor(Adduct, levels = Adduct))
# 
# labs_pos <- paste0(temp_data_pos$Adduct,
#                    " (", round(temp_data_pos$value, 2), "%)")
# 
# plot_pos <-
#   ggpie(
#     data = temp_data_pos,
#     label = labs_pos,
#     x = "value",
#     fill = "Adduct",
#     color = "white"
#   ) +
#   ggsci::scale_fill_aaas()
# 
# plot_pos
# 
# ggsave(plot_pos,
#        filename = "adduct_positive_mode.pdf",
#        width = 7,
#        height = 7)
# 
# 
# labs_neg <- paste0(temp_data_neg$Adduct,
#                    " (", round(temp_data_neg$value, 2), "%)")
# 
# plot_neg <-
#   ggpie(
#     data = temp_data_neg,
#     label = labs_neg,
#     x = "value",
#     fill = "Adduct",
#     color = "white"
#   ) +
#   ggsci::scale_fill_nejm()
# 
# plot_neg
# 
# ggsave(plot_neg,
#        filename = "adduct_negative_mode.pdf",
#        width = 7,
#        height = 7)
# 
# 
# #####remove some compound classes
# remain_index <-
#   purrr::map(compound_class, function(x) {
#     if (any(
#       x$Adduct %in% c(
#         temp_data$Adduct,
#         "(M+H)+",
#         "(M+Na)+",
#         "(M+K)+",
#         "(M-H)-",
#         "(M+Cl)-"
#       )
#     ) | nrow(x) > 1) {
#       return(TRUE)
#     } else{
#       return(FALSE)
#     }
#   }) %>%
#   unlist() %>%
#   which()
# 
# # remove_index <-
# #   purrr::map(compound_class, function(x) {
# #     if (any(
# #       x$Adduct %in% c(
# #         temp_data$Adduct,
# #         "(M+H)+",
# #         "(M+Na)+",
# #         "(M+K)+",
# #         "(M-H)-",
# #         "(M+Cl)-"
# #       )
# #     ) | nrow(x) > 1) {
# #       return(FALSE)
# #     } else{
# #       return(TRUE)
# #     }
# #   }) %>%
# #   unlist() %>%
# #   which()
# 
# compound_class <-
#   compound_class[remain_index]
# 
# annotation_table <-
#   compound_class %>%
#   do.call(rbind, .) %>%
#   as.data.frame()
# 
# length(unique(annotation_table$HMDB.ID))
# length(compound_class)
# 
# ###score distributation
# temp_data <-
#   annotation_table %>%
#   dplyr::select(compound_class, score) %>%
#   dplyr::distinct(compound_class, score) %>%
#   dplyr::pull(score) %>%
#   table() %>%
#   as.data.frame()
# 
# colnames(temp_data) <- c("Score", "Freq")
# 
# temp_data <-
#   temp_data %>%
#   dplyr::mutate(Score = as.numeric(as.character(Score)))
# 
# plot <-
#   temp_data %>%
#   ggplot(aes(Score, Freq)) +
#   geom_bar(stat = "identity",
#            aes(x = Score,
#                y = Freq,
#                fill = Score),
#            show.legend = FALSE) +
#   theme_bw() +
#   # scale_y_continuous(expand = c(0, 10)) +
#   scale_x_continuous(
#     breaks = seq(0, 500, by = 20),
#     name = seq(0, 500, by = 20),
#     expand = expansion(mult = c(0.01, 0.01))
#   ) +
#   theme(
#     axis.title = element_text(size = 13),
#     axis.text = element_text(size = 12),
#     plot.background  = element_rect(fill = "transparent",
#                                     color = NA),
#     panel.background = element_rect(fill = "transparent",
#                                     color = NA),
#     legend.background = element_rect(fill = "transparent",
#                                      color = NA),
#     axis.text.x = element_text(
#       size = 12,
#       angle = 45,
#       hjust = 1,
#       vjust = 1
#     ),
#     panel.grid.minor = element_blank()
#   ) +
#   labs(x = "Score") +
#   ggforce::facet_zoom(
#     xlim = c(98, 500),
#     ylim = c(0, 4000),
#     horizontal = FALSE,
#     zoom.size = 1
#   )
# 
# plot
# 
# ggsave(
#   plot,
#   filename = "score_distributation.pdf",
#   width = 10,
#   height = 7,
#   bg = "transparent"
# )
# 
# 
# ######--------------------------------------------------------------------------
# ######gold standard from Level 1
# library(plyr)
# 
# golden_standard <-
#   annotation_table %>%
#   dplyr::filter(!is.na(Level)) %>%
#   dplyr::filter(Level == 1 & SS > 0.7) %>%
#   dplyr::select(
#     variable_id,
#     Compound.name2 = Compound.name,
#     CAS.ID2 = CAS.ID,
#     HMDB.ID2 = HMDB.ID,
#     KEGG.ID2 = KEGG.ID,
#     Adduct2 = Adduct,
#     mz.error2 = mz.error,
#     RT.error2 = RT.error,
#     SS2 = SS,
#     Total.score2 = Total.score,
#     Database2 = Database
#   )
# 
# temp_data <-
#   golden_standard %>%
#   dplyr::left_join(annotation_table, by = "variable_id") %>%
#   dplyr::select(
#     c(
#       compound_class,
#       variable_id,
#       Compound.name,
#       Compound.name2,
#       CAS.ID,
#       CAS.ID2,
#       HMDB.ID,
#       HMDB.ID2,
#       KEGG.ID,
#       KEGG.ID2,
#       Adduct,
#       Adduct2,
#       Database,
#       Database2,
#       Formula,
#       isotopes,
#       score
#     )
#   ) %>%
#   plyr::dlply(.variables = .(variable_id)) %>%
#   purrr::map(function(x) {
#     x %>%
#       dplyr::arrange(desc(score), compound_class)
#   })
# 
# length(temp_data)
# 
# temp_data[[1]]
# 
# temp_data %>%
#   lapply(function(x) {
#     max(x$score)
#   })  %>%
#   unlist() %>%
#   `>`(100) %>% which
# 
# temp_data[[19]]

###remove redundant annotation according to compound
# library(plyr)
###redundancy is
calculate_redundancy(annotation_table = annotation_table)

annotation_table <-
  remove_redundancy4fmsea(annotation_table = annotation_table, 
                          score_cutoff = 20)

calculate_redundancy(annotation_table = annotation_table)

library(data.table)

annotation_table

save(annotation_table, file = "annotation_table")
