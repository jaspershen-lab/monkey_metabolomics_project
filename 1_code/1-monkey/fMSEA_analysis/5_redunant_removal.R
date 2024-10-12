no_source()
####data preparation

library(tidymass)
library(tidyverse)

rm(list = ls())

setwd(r4projects::get_project_wd())
source("1-code/fMSEA/remove_redundancy.R")
source("1-code/fMSEA/calculate_redundancy.R")

#####load database
load(
  "3-data_analysis/monkey_fMSEA_analysis/2_metabolite_annotation/database/hmdb_ms1.rda"
)

setwd(r4projects::get_project_wd())

load("3-data_analysis/monkey_fMSEA_analysis/4_group_features/annotation_table")
# load(
#   "3-data_analysis/monkey_fMSEA_analysis/4_group_features/expression_data_hilic_neg"
# )
# load(
#   "3-data_analysis/monkey_fMSEA_analysis/4_group_features/expression_data_hilic_pos"
# )
# load(
#   "3-data_analysis/monkey_fMSEA_analysis/4_group_features/expression_data_rplc_neg"
# )
# load(
#   "3-data_analysis/monkey_fMSEA_analysis/4_group_features/expression_data_rplc_pos"
# )
# load("3-data_analysis/monkey_fMSEA_analysis/4_group_features/sample_info_hilic_neg")
# load("3-data_analysis/monkey_fMSEA_analysis/4_group_features/sample_info_hilic_pos")
# load("3-data_analysis/monkey_fMSEA_analysis/4_group_features/sample_info_rplc_neg")
# load("3-data_analysis/monkey_fMSEA_analysis/4_group_features/sample_info_rplc_pos")

dir.create("3-data_analysis/monkey_fMSEA_analysis/5_redundant_removal/")
setwd("3-data_analysis/monkey_fMSEA_analysis/5_redundant_removal/")

###score distributation
temp_data <-
  annotation_table %>%
  dplyr::select(compound_class, score) %>%
  dplyr::distinct(compound_class, score) %>%
  dplyr::pull(score) %>%
  table() %>%
  as.data.frame()

colnames(temp_data) <- c("Score", "Freq")

temp_data <-
  temp_data %>%
  dplyr::mutate(Score = as.numeric(as.character(Score)))

plot <-
  temp_data %>%
  ggplot(aes(Score, Freq)) +
  geom_bar(stat = "identity",
           aes(x = Score,
               y = Freq,
               fill = Score),
           show.legend = FALSE) +
  theme_bw() +
  # scale_y_continuous(expand = c(0, 10)) +
  scale_x_continuous(breaks = seq(0, 200, by = 20),
                     name = seq(0, 200, by = 20)) +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    plot.background  = element_rect(fill = "transparent",
                                    color = NA),
    panel.background = element_rect(fill = "transparent",
                                    color = NA),
    legend.background = element_rect(fill = "transparent",
                                     color = NA),
    axis.text.x = element_text(
      size = 12,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    panel.grid.minor = element_blank()
  ) +
  labs(x = "Score") +
  ggforce::facet_zoom(
    xlim = c(80, 200),
    ylim = c(0, 8500),
    horizontal = FALSE,
    zoom.size = 1
  )

plot

# ggsave(
#   plot,
#   filename = "score_distributation.pdf",
#   width = 10,
#   height = 7,
#   bg = "transparent"
# )

######--------------------------------------------------------------------------
######gold standard from Level 1
library(plyr)

golden_standard <-
  annotation_table %>%
  dplyr::filter(!is.na(Level)) %>%
  dplyr::filter(Level == 1 & SS > 0.7) %>%
  dplyr::select(
    variable_id,
    Compound.name2 = Compound.name,
    CAS.ID2 = CAS.ID,
    HMDB.ID2 = HMDB.ID,
    KEGG.ID2 = KEGG.ID,
    Adduct2 = Adduct,
    mz.error2 = mz.error,
    RT.error2 = RT.error,
    SS2 = SS,
    Total.score2 = Total.score,
    Database2 = Database
  )

temp_data <-
  golden_standard %>%
  dplyr::left_join(annotation_table, by = "variable_id") %>%
  dplyr::select(
    c(
      compound_class,
      variable_id,
      Compound.name,
      Compound.name2,
      CAS.ID,
      CAS.ID2,
      HMDB.ID,
      HMDB.ID2,
      KEGG.ID,
      KEGG.ID2,
      Adduct,
      Adduct2,
      Database,
      Database2,
      Formula,
      isotopes,
      score
    )
  ) %>%
  plyr::dlply(.variables = .(variable_id)) %>%
  purrr::map(function(x) {
    x %>%
      dplyr::arrange(desc(score), compound_class)
  })

length(temp_data)

temp_data[[2]]

temp_data %>% lapply(function(x)
  max(x$score))  %>% unlist() %>% `>`(180) %>% which()

temp_data[[196]]

#score distribution
score_rule <-
  data.frame(
    rule =
      c(
        "Adduct M+H",
        "Isotope (Adduct M+H)",
        "Other positive adduct",
        "Isotope (Other positive adduct)",
        "Adduct M-H",
        "Isotope (Adduct M-H)",
        "Other negative adduct",
        "Isotope (Other negative adduct)"
      ),
    score = c(50, 20, 20, 10, 50, 20, 20, 10),
    stringsAsFactors = FALSE
  )

3 * 3 * 3 * 3

test <- list(c(1, 1), c(1, 0), c(0, 0))

# comb <- vector(mode = "list", length = 3^4)
comb <- NULL

for (i in 1:3) {
  for (j in 1:3) {
    for (k in 1:3) {
      for (z in 1:3) {
        comb <- c(comb, list(c(test[[i]], test[[j]], test[[k]], test[[z]])))
      }
    }
  }
}

comb <-
  comb %>%
  do.call(cbind, .)

remove_idx <-
  apply(comb, 2, function(x) {
    all(x == 0)
  }) %>%
  which()

comb <- comb[, -remove_idx]

score <- score_rule$score * comb
score <- score %>% colSums()

idx <- order(score, decreasing = TRUE)

score <- score[idx]

comb <- comb[, idx]

score_rule <-
  data.frame(score_rule, comb, stringsAsFactors = FALSE)

plot1 <-
  data.frame(index = factor(paste('X', 1:length(score), sep = ""),
                            levels = paste('X', 1:length(score), sep = "")),
             score,
             stringsAsFactors = TRUE) %>%
  ggplot(aes(index, score)) +
  geom_point(
    stat = "identity",
    aes(color = as.character(score)),
    shape = 16,
    show.legend = FALSE
  ) +
  geom_segment(aes(
    x = index,
    y = 0,
    xend = index,
    yend = score,
    color = as.character(score)
  ),
  show.legend = FALSE) +
  theme_classic() +
  scale_y_continuous(expand = expansion(mult = c(0, .05))) +
  labs(x = "", y = "Score") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    plot.margin = unit(c(0, 0, 0, 0), "pt")
  )

plot1

plot2 <-
  score_rule %>%
  dplyr::select(-score) %>%
  tidyr::pivot_longer(cols = -rule,
                      names_to = "index",
                      values_to = "value") %>%
  dplyr::mutate(rule = factor(
    rule,
    levels =
      c(
        "Adduct M+H",
        "Isotope (Adduct M+H)",
        "Other positive adduct",
        "Isotope (Other positive adduct)",
        "Adduct M-H",
        "Isotope (Adduct M-H)",
        "Other negative adduct",
        "Isotope (Other negative adduct)"
      ) %>% rev()
  )) %>%
  ggplot(aes(index, rule)) +
  geom_tile(aes(fill = as.character(value)),
            color = "#FF6F00FF",
            show.legend = FALSE) +
  theme_bw() +
  labs(x = "", y = "") +
  scale_y_discrete(expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(values = c("0" = "white", "1" = "#008EA0FF")) +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 10),
    plot.margin = unit(c(0, 0, 0, 0), "pt")
  )

plot2

library(patchwork)

plot <-
  plot1 + plot2 + patchwork::plot_layout(ncol = 1, heights = c(1, 3))

plot

# ggsave(
#   plot,
#   filename = "score_rule.pdf",
#   width = 12,
#   height = 7,
#   bg = "transparent"
# )

###remove redundant annotation according to compound
library(plyr)
###redundancy is
calculate_redundancy(annotation_table = annotation_table)

annotation_table <-
  remove_redundancy(annotation_table = annotation_table)

calculate_redundancy(annotation_table = annotation_table)

save(annotation_table, file = "annotation_table")

# library(pmd)
# data("spmeinvivo")
# str(spmeinvivo)
# 
# pmd <- getpaired(spmeinvivo, rtcutoff = 5)
# plotrtg(pmd)
# plotpaired(pmd)
# 
# diff <- unique(pmd$paired$diff2)[1]
# index <- pmd$paired$diff2 == diff
# plotpaired(pmd,index)
# 
# std <- getstd(pmd)
# plotstd(std)
# 
# plotpca(std$data[std$stdmassindex,],
#         lv = as.numeric(as.factor(std$group$sample_group)),
#         main = paste(sum(std$stdmassindex),"independent peaks"))
# 
# x <-
#   data.frame(
#     mz = std$mz,
#     rt = std$rt,
#     rtcluster = std$rtcluster,
#     index = std$stdmassindex
#   ) %>%
#   dplyr::arrange(rt) %>%
#   plyr::dlply(.variables = .(rtcluster))
# 
