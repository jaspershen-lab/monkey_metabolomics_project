no_source()
####This should be run in workstation

library(tidymass)
library(tidyverse)

rm(list = ls())

setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")

load("3_data_analysis/monkey_sex/object")

dir.create("3_data_analysis/monkey_sex_correlation_network")
setwd("3_data_analysis/monkey_sex_correlation_network")


####
object %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::count(sex)

male_object <-
  object %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(sex == "M")

famale_object <-
  object %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(sex == "F")

###calculate correlation
male_cor <-
  male_object %>%
  `+`(1) %>%
  log(2) %>%
  massstat::cor_mass_dataset(
    margin = "variable",
    method = "spearman",
    p_adjust_method = "BH",
    data_type = "longer"
  )

famale_cor <-
  famale_object %>%
  `+`(1) %>%
  log(2) %>%
  massstat::cor_mass_dataset(
    margin = "variable",
    method = "spearman",
    p_adjust_method = "BH",
    data_type = "longer"
  )

male_cor <-
  male_cor %>%
  dplyr::filter(p_adjust < 0.05)

famale_cor <-
  famale_cor %>%
  dplyr::filter(p_adjust < 0.05)

unique(c(male_cor$from,
         male_cor$to))

unique(c(famale_cor$from,
         famale_cor$to))

load(
  here::here(
    "3_data_analysis/monkey_metabolomics_metabolite_annotation/hmdb_ms1.rda"
  )
)

variable_info <-
  extract_variable_info(object) %>%
  dplyr::left_join(hmdb_ms1@spectra.info[, c("HMDB.ID", "Kingdom", "Super_class", "Class", "Sub_class")],
                   by = "HMDB.ID")

male_cor <-
  male_cor %>%
  dplyr::left_join(variable_info[, c("variable_id", "Super_class", "Class")],
                   by = c("from" = "variable_id")) %>%
  dplyr::rename(from_super_class = Super_class,
                from_class = Class) %>%
  dplyr::left_join(variable_info[, c("variable_id", "Super_class", "Class")],
                   by = c("to" = "variable_id")) %>%
  dplyr::rename(to_super_class = Super_class,
                to_class = Class)

famale_cor <-
  famale_cor %>%
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

mat_male <-
  male_cor %>%
  dplyr::filter(!is.na(from_super_class) &
                  !is.na(to_super_class)) %>%
  dplyr::count(from_super_class, to_super_class) %>%
  dplyr::filter(n >= 5) %>%
  tidyr::pivot_wider(names_from = "to_super_class", values_from = "n") %>%
  tibble::column_to_rownames(var = "from_super_class") %>%
  as.data.frame()

mat_male[which(is.na(mat_male), arr.ind = TRUE)] <- 0

grid_color1 <-
  grid_color[names(grid_color) %in% c(colnames(mat_male), rownames(mat_male))]

circos.clear()
circos.par(gap.after = 0)
chordDiagram(as.matrix(mat_male),
             grid.col = grid_color1,
             order = names(grid_color1))
circos.clear()

mat_female <-
  famale_cor %>%
  dplyr::filter(!is.na(from_super_class) &
                  !is.na(to_super_class)) %>%
  dplyr::count(from_super_class, to_super_class) %>%
  dplyr::filter(n >= 5) %>%
  tidyr::pivot_wider(names_from = "to_super_class", values_from = "n") %>%
  tibble::column_to_rownames(var = "from_super_class") %>%
  as.data.frame()

mat_female[which(is.na(mat_female), arr.ind = TRUE)] <- 0

grid_color2 <-
  grid_color[names(grid_color) %in% c(colnames(mat_female), rownames(mat_female))]

circos.clear()
circos.par(gap.after = 0)
chordDiagram(as.matrix(mat_female),
             grid.col = grid_color1,
             order = names(grid_color2))
circos.clear()

#####node distributation
edge_data_male <-
  male_cor

node_data_male <-
  data.frame(node = unique(c(
    edge_data_male$from, edge_data_male$to
  ))) %>%
  dplyr::left_join(variable_info, by = c("node" = "variable_id"))

library(igraph)
library(ggraph)
library(tidygraph)

male_graph <-
  tbl_graph(nodes = node_data_male, edges = edge_data_male)

edge_data_female <-
  famale_cor

node_data_female <-
  data.frame(node = unique(c(edge_data_female$from, edge_data_female$to))) %>%
  dplyr::left_join(variable_info, by = c("node" = "variable_id"))

library(igraph)
library(ggraph)
library(tidygraph)

female_graph <-
  tbl_graph(nodes = node_data_female, edges = edge_data_female)

library(ggmosaic)

temp_data <-
  rbind(
    data.frame(Super_class = node_data_male[, c("Super_class")], sex = 'M'),
    data.frame(Super_class = node_data_female[, c("Super_class")], sex = 'F')
  ) %>%
  dplyr::filter(!is.na(Super_class))

plot <-
  temp_data %>%
  ggplot() +
  geom_mosaic(aes(x = product(sex, Super_class), fill = Super_class),
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



edge_data_male

mat_male <-
  male_cor %>%
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

mat_female <-
  famale_cor %>%
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

range(mat_male$n)
range(mat_female$n)

# mat <-
# rbind(data.frame(mat_male, class = "young"),
#       data.frame(mat_female, class = "old"))

plot1 <-
  mat_male %>%
  ggplot(aes(from_super_class, to_super_class)) +
  geom_point(aes(size = n,
                 color = class)) +
  theme_base +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text = element_text(size = 10)
  ) +
  scale_size_continuous(limits = c(5, 525), range = c(1, 8)) +
  labs(x = "", y = "") +
  scale_color_manual(values = c(
    "same" = ggsci::pal_nejm()(n = 9)[4],
    "different" = ggsci::pal_nejm()(n = 9)[2]
  ))

plot2 <-
  mat_female %>%
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
  scale_size_continuous(limits = c(5, 525), range = c(1, 8)) +
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
       filename = "edge_male_female.pdf",
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
  male_graph %>%
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
  scale_color_manual(values = metabolite_super_class_color[names(metabolite_super_class_color) %in% node_data_male$Super_class]) +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot
ggsave(plot,
       filename = "male_network.pdf",
       width = 10,
       height = 8)

plot <-
  female_graph %>%
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
  scale_color_manual(values = metabolite_super_class_color[names(metabolite_super_class_color) %in% node_data_male$Super_class]) +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot
ggsave(plot,
       filename = "female_network.pdf",
       width = 10,
       height = 8)
