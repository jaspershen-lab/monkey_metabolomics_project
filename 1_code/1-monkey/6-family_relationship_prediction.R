no_source()
####This should be run in workstation

library(tidymass)
library(tidyverse)

rm(list = ls())

setwd(r4projects::get_project_wd())
source("1-code/tools.R")

load("3-data_analysis/monkey_metabolomics_data_preparation/metabolite/object")

dir.create("3-data_analysis/family_relationship_prediction")
setwd("3-data_analysis/family_relationship_prediction")

sample_info <-
  extract_sample_info(object) %>%
  dplyr::filter(group == "Subject")

####family tree using ggraph
####father
library(igraph)
library(ggraph)
library(tidygraph)

temp_data <-
  sample_info %>%
  dplyr::select(sample_id, sex, mother, Father, age) %>%
  dplyr::rename(father = Father)

node_data <-
  rbind(temp_data) %>%
  dplyr::rename(node = sample_id) %>%
  dplyr::distinct(node, .keep_all = TRUE) %>%
  dplyr::filter(node != "unknown")

edge_data1 <-
  temp_data %>%
  dplyr::select(from = father,
                to = sample_id) %>%
  dplyr::mutate(class = "father_child") %>%
  dplyr::filter(from != "unknown")

edge_data2 <-
  temp_data %>%
  dplyr::select(from = mother,
                to = sample_id) %>%
  dplyr::mutate(class = "mother_child") %>%
  dplyr::filter(from != "unknown")

edge_data <-
  rbind(edge_data1,
        edge_data2)

# node_data <-
#   node_data %>%
#   dplyr::filter(node %in% unique(c(edge_data$from, edge_data$to)))

edge_data <-
  edge_data %>%
  dplyr::filter(from %in% node_data$node & to %in% node_data$node)

# node_data <-
#   node_data %>%
#   dplyr::mutate(in_dataset =
#                   case_when(is.na(age) ~ "No",
#                             TRUE ~ "Yes"))
node_data$age[is.na(node_data$age)] <- 5

temp_graph <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = TRUE)

plot <-
  ggraph(temp_graph, layout = "fr") +
  geom_edge_link(arrow = arrow(length = unit(2, 'mm'),
                               type = "open"),
                 aes(color = class)) +
  geom_node_point(aes(fill = sex,
                      size = age),
                  shape = 21) +
  scale_fill_manual(values = sex_color) +
  scale_edge_color_manual(values = c(
    "father_child" = unname(sex_color["M"]),
    "mother_child" = unname(sex_color["F"])
  )) +
  scale_size_continuous(range = c(1, 5)) +
  ggraph::theme_graph() +
  shadowtext::geom_shadowtext(
    aes(x = x,
        y = y,
        label = node),
    bg.colour = 'white',
    color = "black",
    size = 2,
    check_overlap = TRUE
  )
extrafont::loadfonts()
plot
ggsave(plot,
       filename = "family_tree.pdf",
       width = 8,
       height = 7)

library(igraph)

###get the family distance for each two nodes
# family_distance <-
#   1:(nrow(node_data) - 1) %>%
#   purrr::map(function(idx1) {
#     cat(idx1, " ")
#     (idx1 + 1):nrow(node_data) %>%
#       purrr::map(function(idx2) {
#         family_score <-
#           get_family_score(graph = temp_graph,
#                            from = node_data$node[idx1],
#                            to = node_data$node[idx2])
#         data.frame(
#           from = node_data$node[idx1],
#           to = node_data$node[idx2],
#           family_distance = family_score
#         )
#       }) %>%
#       dplyr::bind_rows() %>%
#       as.data.frame()
#   })
#
# save(family_distance, file = "family_distance")
load("family_distance")
family_distance <-
  family_distance %>%
  dplyr::bind_rows() %>%
  as.data.frame()

family_distance %>%
  dplyr::filter(family_distance != Inf) %>%
  ggplot(aes(to, from)) +
  geom_tile(aes(fill = family_distance))

temp_data <-
  family_distance %>%
  dplyr::arrange(family_distance) %>%
  dplyr::filter(
    from %in% c("129F", "69U", "26B", "108F", "DBKW") &
      to %in% c("129F", "69U", "26B", "108F", "DBKW")
  )

sub_graph <-
  igraph::subgraph(graph = graph, 
                   v = match(c("129F", "69U", "26B", "108F", "DBKW"),
                             node_data$node))

plot <-
  ggraph(sub_graph, layout = "fr") +
  geom_edge_link(arrow = arrow(length = unit(10, 'mm'),
                               type = "open"),
                 aes(color = class)) +
  geom_node_point(aes(fill = sex),
                  size = 15,
                  shape = 21) +
  scale_fill_manual(values = sex_color) +
  scale_edge_color_manual(values = c(
    "father_child" = unname(sex_color["M"]),
    "mother_child" = unname(sex_color["F"])
  )) +
  scale_size_continuous(range = c(10, 15)) +
  ggraph::theme_graph() +
  shadowtext::geom_shadowtext(
    aes(x = x,
        y = y,
        label = node),
    bg.colour = 'white',
    color = "black",
    size = 10,
    check_overlap = TRUE
  )
extrafont::loadfonts()
plot
# ggsave(plot,
#        filename = "example_family_tree.pdf",
#        width = 8,
#        height = 7)


temp_data2 <-
  temp_data

temp_data2$from <- temp_data$to
temp_data2$to <- temp_data$from

temp_data <-
  rbind(temp_data,
        temp_data2)

plot <- 
temp_data %>%
  dplyr::mutate(family_distance = as.character(family_distance)) %>%
  ggplot(aes(to, from)) +
  geom_tile(aes(fill = family_distance), color = "black") +
  scale_fill_manual(values = family_distance_color) +
  theme_base +
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  scale_y_discrete(expand = expansion(mult = c(0, 0))) +
  theme(panel.grid = element_blank()) +
  geom_text(aes(label =family_distance), color = "white") +
  labs(x = "", y = "")

ggsave(plot, filename = "example_family_distance.pdf", width = 8.8, height = 7)

subject_id <- 
  object %>% 
  activate_mass_dataset(what = "sample_info") %>% 
  filter(class == "Subject") %>% 
  pull(sample_id)

temp_object <- 
  object  %>% 
  mutate_rsd(according_to_samples = subject_id) %>% 
  activate_mass_dataset(what = "variable_info") %>% 
  filter(rsd > 0)

# sample_cor_data <-
#   seq_len(nrow(family_distance)) %>%
#   purrr::map(function(idx) {
#     cat(idx, " ")
#     value1 <-
#       temp_object %>%
#       `+`(1) %>%
#       log(2) %>%
#       activate_mass_dataset(what = "expression_data") %>%
#       dplyr::pull(family_distance$from[idx])
#     
#     value2 <-
#       temp_object %>%
#       `+`(1) %>%
#       log(2) %>%
#       activate_mass_dataset(what = "expression_data") %>%
#       dplyr::pull(family_distance$to[idx])
#     result <-
#       cor.test(value1, value2)
#     data.frame(cor = as.numeric(result$estimate),
#                p_value = as.numeric(result$p.value))
#   }) %>%
#   dplyr::bind_rows()
# 
# save(sample_cor_data, file = "sample_cor_data")
load("sample_cor_data")

temp_data <-
  cbind(family_distance,
        sample_cor_data) %>%
  dplyr::mutate(family_distance = case_when(family_distance == Inf ~ 10,
                                            TRUE ~ family_distance)) %>% 
  # dplyr::mutate(family_distance = as.character(family_distance)) %>%
  dplyr::mutate(family_distance = factor(family_distance, levels = c(1:4,10)))

library(ggsignif)

plot <- 
temp_data %>%
  ggplot(aes(family_distance, cor, group=family_distance)) +
  theme_base +
  geom_jitter(aes(color = family_distance),
               shape = 16,
               alpha = 0.8,
              show.legend = FALSE) +
  geom_boxplot(outlier.shape = NA,
               fill = "transparent",
               color = "black") +
  scale_x_discrete(labels = c("1", "2", "3", "4", "Inf")) +
  scale_fill_manual(values = c(family_distance_color, "10" = "grey")) +
  scale_color_manual(values = c(family_distance_color, "10" = "grey")) +
  labs(x = "Family distance", y = "Spearman correlation") +
  ggsignif::geom_signif(
    test = "wilcox.test",
    comparisons = list(
    c("1", "2"),
    c("1", "3"),
    c("1", "4"),
    c("1", "10"),
    c("2", "3"),
    c("2", "4"),
    c("2", "10")
  ),
  y_position = c(1, 0.99, 0.98, 0.97, 0.96, 0.95, 0.94))
plot

ggsave(plot, filename = "family_distance_vs_correlation.pdf", width = 6, height = 8)
