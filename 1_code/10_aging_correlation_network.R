no_source()
####This should be run in workstation

library(tidymass)
library(tidyverse)

rm(list = ls())

setwd(r4projects::get_project_wd())
source("1_code/tools.R")

load("3_data_analysis/monkey_metabolomics_data_preparation/metabolite/object")

dir.create("3_data_analysis/monkey_aging_correlation_network")
setwd("3_data_analysis/monkey_aging_correlation_network")

object <-
  object %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(!is.na(age)) %>%
  dplyr::mutate(age_class = case_when(age <= 4 ~ "young",
                                      age > 4 ~ "old"))

####
object %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::count(age_class)

young_object <-
  object %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(age_class == "young")

old_object <-
  object %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(age_class == "old")

###calculate correlation
young_cor <-
  young_object %>%
  `+`(1) %>%
  log(2) %>%
  massstat::cor_mass_dataset(
    margin = "variable",
    method = "spearman",
    p_adjust_method = "BH",
    data_type = "longer"
  )

old_cor <-
  old_object %>%
  `+`(1) %>%
  log(2) %>%
  massstat::cor_mass_dataset(
    margin = "variable",
    method = "spearman",
    p_adjust_method = "BH",
    data_type = "longer"
  )

young_cor <-
  young_cor %>%
  dplyr::filter(p_adjust < 0.05)

old_cor <-
  old_cor %>%
  dplyr::filter(p_adjust < 0.05)

unique(c(young_cor$from,
         young_cor$to))

unique(c(old_cor$from,
         old_cor$to))

load(
  here::here(
    "3_data_analysis/monkey_metabolomics_metabolite_annotation/hmdb_ms1.rda"
  )
)

variable_info <-
  extract_variable_info(object) %>%
  dplyr::left_join(hmdb_ms1@spectra.info[, c("HMDB.ID", "Kingdom", "Super_class", "Class", "Sub_class")],
                   by = "HMDB.ID")

young_cor <-
  young_cor %>%
  dplyr::left_join(variable_info[, c("variable_id", "Super_class", "Class")],
                   by = c("from" = "variable_id")) %>%
  dplyr::rename(from_super_class = Super_class,
                from_class = Class) %>%
  dplyr::left_join(variable_info[, c("variable_id", "Super_class", "Class")],
                   by = c("to" = "variable_id")) %>%
  dplyr::rename(to_super_class = Super_class,
                to_class = Class)

old_cor <-
  old_cor %>%
  dplyr::left_join(variable_info[, c("variable_id", "Super_class", "Class")],
                   by = c("from" = "variable_id")) %>%
  dplyr::rename(from_super_class = Super_class,
                from_class = Class) %>%
  dplyr::left_join(variable_info[, c("variable_id", "Super_class", "Class")],
                   by = c("to" = "variable_id")) %>%
  dplyr::rename(to_super_class = Super_class,
                to_class = Class)

grid_color <-
  metabolite_super_class_color

mat_young <-
  young_cor %>%
  dplyr::filter(!is.na(from_super_class) &
                  !is.na(to_super_class)) %>%
  dplyr::count(from_super_class, to_super_class) %>%
  dplyr::filter(n >= 5) %>%
  tidyr::pivot_wider(names_from = "to_super_class", values_from = "n") %>%
  tibble::column_to_rownames(var = "from_super_class") %>%
  as.data.frame()

mat_young[which(is.na(mat_young), arr.ind = TRUE)] <- 0

grid_color1 <-
  grid_color[names(grid_color) %in% c(colnames(mat_young), rownames(mat_young))]

circos.clear()
circos.par(gap.after = 0)
chordDiagram(as.matrix(mat_young),
             grid.col = grid_color1,
             order = names(grid_color1))
circos.clear()

mat_old <-
  old_cor %>%
  dplyr::filter(!is.na(from_super_class) &
                  !is.na(to_super_class)) %>%
  dplyr::count(from_super_class, to_super_class) %>%
  dplyr::filter(n >= 5) %>%
  tidyr::pivot_wider(names_from = "to_super_class", values_from = "n") %>%
  tibble::column_to_rownames(var = "from_super_class") %>%
  as.data.frame()

mat_old[which(is.na(mat_old), arr.ind = TRUE)] <- 0

grid_color2 <-
  grid_color[names(grid_color) %in% c(colnames(mat_old), rownames(mat_old))]

circos.clear()
circos.par(gap.after = 0)
chordDiagram(as.matrix(mat_old),
             grid.col = grid_color1,
             order = names(grid_color1))
circos.clear()

#####node distributation
edge_data_young <-
  young_cor

node_data_young <-
  data.frame(node = unique(c(
    edge_data_young$from, edge_data_young$to
  ))) %>%
  dplyr::left_join(variable_info, by = c("node" = "variable_id"))

library(igraph)
library(ggraph)
library(tidygraph)
young_graph <-
  tbl_graph(nodes = node_data_young, edges = edge_data_young)



edge_data_old <-
  old_cor

node_data_old <-
  data.frame(node = unique(c(edge_data_old$from, edge_data_old$to))) %>%
  dplyr::left_join(variable_info, by = c("node" = "variable_id"))

library(igraph)
library(ggraph)
library(tidygraph)
old_graph <-
  tbl_graph(nodes = node_data_old, edges = edge_data_old)

library(ggmosaic)

temp_data <-
  rbind(
    data.frame(Super_class = node_data_young[, c("Super_class")], age = 'young'),
    data.frame(Super_class = node_data_old[, c("Super_class")], age = 'old')
  ) %>%
  dplyr::filter(!is.na(Super_class))

plot <-
  temp_data %>%
  ggplot() +
  geom_mosaic(aes(x = product(age, Super_class), fill = Super_class),
              offset = 0.01,
              show.legend = FALSE) +
  theme_mosaic() +
  scale_fill_manual(values = metabolite_super_class_color[names(metabolite_super_class_color) %in%
                                                            unique(temp_data$Super_class)]) +
  labs(x = "", y = "") +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    panel.border = element_rect(color = "black",
                                fill = "transparent"),
    legend.position = "right"
  )

plot

ggsave(
  plot = plot,
  filename = "node_distributation.pdf",
  width = 12,
  height = 7
)



edge_data_young

mat_young <-
  young_cor %>%
  dplyr::filter(!is.na(from_super_class) &
                  !is.na(to_super_class)) %>%
  dplyr::count(from_super_class, to_super_class) %>%
  dplyr::filter(n >= 5) %>%
  dplyr::mutate(
    class = case_when(
      from_super_class == to_super_class ~ "same",
      from_super_class != to_super_class ~ "different"
    )
  )

mat_old <-
  old_cor %>%
  dplyr::filter(!is.na(from_super_class) &
                  !is.na(to_super_class)) %>%
  dplyr::count(from_super_class, to_super_class) %>%
  dplyr::filter(n >= 5) %>%
  dplyr::mutate(
    class = case_when(
      from_super_class == to_super_class ~ "same",
      from_super_class != to_super_class ~ "different"
    )
  )

range(mat_young$n)
range(mat_old$n)

# mat <-
# rbind(data.frame(mat_young, class = "young"),
#       data.frame(mat_old, class = "old"))

plot1 <-
  mat_young %>%
  ggplot(aes(from_super_class, to_super_class)) +
  geom_point(aes(size = n,
                 color = class)) +
  theme_base +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text = element_text(size = 10)
  ) +
  scale_size_continuous(limits = c(5, 3516), range = c(1, 8)) +
  labs(x = "", y = "") +
  scale_color_manual(values = c(
    "same" = ggsci::pal_nejm()(n = 9)[4],
    "different" = ggsci::pal_nejm()(n = 9)[2]
  ))

plot2 <-
  mat_old %>%
  ggplot(aes(from_super_class, to_super_class)) +
  geom_point(aes(size = n,
                 color = class)) +
  theme_base +
  theme(axis.text.x = element_text(
    angle = 45,
    vjust = 1,
    hjust = 1
  ),
  axis.text = element_text(size = 10)) +
  scale_size_continuous(limits = c(5, 3516), range = c(1, 8)) +
  labs(x = "", y = "") +
  scale_color_manual(values = c(
    "same" = ggsci::pal_nejm()(n = 9)[4],
    "different" = ggsci::pal_nejm()(n = 9)[2]
  ))


library(patchwork)

plot <-
  plot1 +
  plot2 +
  patchwork::plot_layout(ncol = 1)
plot
ggsave(plot,
       filename = "edge_young_old.pdf",
       width = 9.35,
       height = 8)

# mat %>%
#   ggplot(aes(from_super_class, to_super_class)) +
#   geom_point(aes(size = n,
#                  color = class),
#              alpha = 0.5) +
#   theme_base +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
#   scale_size_continuous(limits = c(5, 3516), range = c(1, 5)) +
#   labs(x = "", y = "") +
#   scale_color_manual(values = young_color)

plot <-
  young_graph %>%
  activate(what = "nodes") %>%
  filter(!is.na(Super_class)) %>%
  dplyr::mutate(Degree = centrality_degree(mode = 'all')) %>%
  ggraph(layout = "kk") +
  geom_edge_link(aes(color = correlation),
                 alpha = 0.5,
                 show.legend = TRUE) +
  scale_edge_color_gradient2(
    low = unname(marker_color["Down"]),
    mid = "white",
    high = unname(marker_color["Up"])
  ) +
  geom_node_point(aes(color = Super_class,
                      size = Degree)) +
  scale_color_manual(values = metabolite_super_class_color[names(metabolite_super_class_color) %in% node_data_young$Super_class]) +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot
ggsave(plot,
       filename = "yound_network.pdf",
       width = 10,
       height = 8)

plot <-
  old_graph %>%
  activate(what = "nodes") %>%
  filter(!is.na(Super_class)) %>%
  dplyr::mutate(Degree = centrality_degree(mode = 'all')) %>%
  ggraph(layout = "kk") +
  geom_edge_link(aes(color = correlation),
                 alpha = 0.5,
                 show.legend = TRUE) +
  scale_edge_color_gradient2(
    low = unname(marker_color["Down"]),
    mid = "white",
    high = unname(marker_color["Up"])
  ) +
  geom_node_point(aes(color = Super_class,
                      size = Degree)) +
  scale_color_manual(values = metabolite_super_class_color[names(metabolite_super_class_color) %in% node_data_young$Super_class]) +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot
ggsave(plot,
       filename = "old_network.pdf",
       width = 10,
       height = 8)
