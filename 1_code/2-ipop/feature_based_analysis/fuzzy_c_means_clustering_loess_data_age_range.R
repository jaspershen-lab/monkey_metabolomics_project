no_source()

rm(list = ls())

setwd(r4projects::get_project_wd())

library(tidyverse)
library(tidymass)

###load("data)
load("3-data_analysis/ipop_feature_analysis/loess_fit/object_corss_section_loess")

dir.create("3-data_analysis/ipop_feature_analysis/fuzzy_clustering_loess_data",
           recursive = TRUE)

setwd("3-data_analysis/ipop_feature_analysis/fuzzy_clustering_loess_data")

object_corss_section_loess

object_corss_section_loess %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::count(mode)

###remove samples without adjusted_age
object_corss_section_loess <-
  object_corss_section_loess %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(!is.na(adjusted_age))

####clustering
library(Mfuzz)

object_corss_section_loess[1, , drop = TRUE] %>%
  unlist() %>%
  density() %>%
  plot()

object_corss_section_loess@sample_info$adjusted_age %>% range

object_corss_section_loess@sample_info %>%
  ggplot(aes(adjusted_age)) +
  geom_histogram(binwidth = 5, color = "black")

sum(object_corss_section_loess@sample_info$adjusted_age <= 35)
sum(
  object_corss_section_loess@sample_info$adjusted_age > 35 &
    object_corss_section_loess@sample_info$adjusted_age <= 45
)
sum(
  object_corss_section_loess@sample_info$adjusted_age > 45 &
    object_corss_section_loess@sample_info$adjusted_age <= 50
)
sum(
  object_corss_section_loess@sample_info$adjusted_age > 50 &
    object_corss_section_loess@sample_info$adjusted_age <= 55
)
sum(
  object_corss_section_loess@sample_info$adjusted_age > 55 &
    object_corss_section_loess@sample_info$adjusted_age <= 60
)
sum(
  object_corss_section_loess@sample_info$adjusted_age > 60 &
    object_corss_section_loess@sample_info$adjusted_age <= 65
)
sum(object_corss_section_loess@sample_info$adjusted_age > 65)

age_index <-
  data.frame(
    from = c(25, 35, 45, 50, 55, 60, 65, 70),
    to = c(35, 45, 50, 55, 60, 65, 70, 76)
  )

temp_data <-
  log(object_corss_section_loess@expression_data + 1, 2)

temp_data_mean <-
  seq_len(nrow(age_index)) %>%
  purrr::map(function(i) {
    idx <- as.numeric(age_index[i,])
    temp_data[, which(
      object_corss_section_loess@sample_info$adjusted_age > idx[1] &
        object_corss_section_loess@sample_info$adjusted_age <= idx[2]
    )] %>%
      apply(1, mean)
  }) %>%
  dplyr::bind_cols() %>%
  as.data.frame()

colnames(temp_data_mean) <-
  paste(age_index$from, age_index$to, sep = "_")

rownames(temp_data_mean) <-
  object_corss_section_loess@variable_info$variable_id

idx = 3

object_corss_section_loess[idx, , drop = TRUE] %>%
  unlist() %>%
  density() %>%
  plot()

temp_data[idx, ] %>%
  as.numeric() %>%
  density() %>%
  plot()

time <- colnames(temp_data_mean)

temp_data <- rbind(time, temp_data_mean)

row.names(temp_data)[1] <- "time"
rownames(temp_data)

# write.table(
#   temp_data,
#   file = "temp_data.txt",
#   sep = '\t',
#   quote = FALSE,
#   col.names = NA
# )

#read it back in as an expression set
data <- table2eset(filename = "temp_data.txt")
data.s <- standardise(data)
m1 <- mestimate(data.s)
m1

# plot <-
#   Dmin(
#     data.s,
#     m = m1,
#     crange = seq(2, 20, 1),
#     repeats = 3,
#     visu = TRUE
#   )
#
# plot <-
#   plot %>%
#   data.frame(distance = plot,
#              k = seq(2, 20, 1)) %>%
#   ggplot(aes(k, distance)) +
#   geom_point(shape = 21, size = 4, fill = "black") +
#   # geom_smooth() +
#   geom_segment(aes(
#     x = k,
#     y = 0,
#     xend = k,
#     yend = distance
#   )) +
#   theme_bw() +
#   theme(
#     # legend.position = c(0, 1),
#     # legend.justification = c(0, 1),
#     panel.grid = element_blank(),
#     axis.title = element_text(size = 13),
#     axis.text = element_text(size = 12),
#     panel.background = element_rect(fill = "transparent", color = NA),
#     plot.background = element_rect(fill = "transparent", color = NA),
#     legend.background = element_rect(fill = "transparent", color = NA)
#   ) +
#   labs(x = "Cluster number",
#        y = "Min. centroid distance") +
#   scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
#
# plot
#
# ggsave(plot,
#        filename = "distance_k_number.pdf",
#        width = 7,
#        height = 7)

cluster <- 12

c <- mfuzz(data.s, c = cluster, m = m1)

# ####any two clusters with correlation > 0.8 should be considered as one
library(corrplot)
layout(1)
center <- c$centers
rownames(center) <- paste("Cluster", rownames(center), sep = ' ')

corrplot::corrplot(
  corr = cor(t(center)),
  type = "full",
  diag = TRUE,
  order = "hclust",
  hclust.method = "ward.D",
  # addrect = 5,
  col = colorRampPalette(colors = rev(
    RColorBrewer::brewer.pal(n = 11, name = "Spectral")
  ))(n = 100),
  number.cex = .7,
  addCoef.col = "black"
)

mfuzz.plot(
  eset = data.s,
  # min.mem = 0.6,
  cl = c,
  mfrow = c(3, 4),
  time.labels = time,
  new.window = FALSE
)

library(ComplexHeatmap)

###
cluster_color <-
  ggsci::pal_jama()(n = 7)[1:12]

names(cluster_color) <- as.character(1:12)

plot <-
  center %>%
  as.data.frame() %>%
  tibble::rowid_to_column(var = "cluster") %>%
  tidyr::pivot_longer(cols = -cluster,
                      names_to = "time",
                      values_to = "value") %>%
  dplyr::mutate(time = factor(time, levels = unique(time))) %>%
  dplyr::mutate(cluster = as.character(cluster)) %>%
  ggplot(aes(time, value)) +
  geom_hline(yintercept = 0) +
  geom_line(aes(group = cluster, color = cluster), show.legend = FALSE) +
  # geom_point(aes(group = cluster, fill = cluster), show.legend = FALSE, shape = 21, size = 3) +
  # geom_smooth(aes(color = cluster, group = cluster), se = FALSE) +
  scale_color_manual(values = cluster_color) +
  scale_fill_manual(values = cluster_color) +
  # facet_grid(rows = vars(class)) +
  theme_bw() +
  labs(x = "", y = "Z-score") +
  theme(
    legend.position = c(0, 1),
    legend.justification = c(0, 1),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 12
    ),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )

plot

ggsave(plot,
       filename = "cluster_center_plot.pdf",
       width = 7,
       height = 7)

centers <- c$centers

cluster_info <-
  data.frame(
    variable_id = names(c$cluster),
    c$membership,
    cluster = c$cluster,
    stringsAsFactors = FALSE
  ) %>%
  arrange(cluster)

###see the largest membership for all the features
c$membership %>%
  apply(1, function(x) {
    max(x)
  }) %>%
  density() %>%
  plot()

c$membership %>%
  apply(1, function(x) {
    max(x)
  }) %>%
  range()

max_membership <-
  c$membership %>%
  apply(1, function(x) {
    max(x)
  }) %>%
  data.frame(max_membership = .) %>%
  dplyr::mutate(index = "x")

library(gghalves)

plot <-
  max_membership %>%
  ggplot() +
  geom_half_boxplot(aes(x = index, y = max_membership), side = "l") +
  geom_half_violin(aes(x = index, y = max_membership), side = "r") +
  theme_bw() +
  labs(x = "", y = "Max membership") +
  theme(panel.grid.minor = element_blank())

plot

quantile(x = max_membership$max_membership)

membership_cutoff <-
  min(max_membership$max_membership)

# ggsave(plot,
#        filename = "max_membership.pdf",
#        width = 7,
#        height = 7)
# ggsave(plot,
#        filename = "max_membership.png",
#        width = 7,
#        height = 7)

plot <-
  c$membership %>%
  apply(1, function(x) {
    sum(x >= 0.135)
  }) %>%
  data.frame(number = .) %>%
  ggplot(aes(number)) +
  geom_bar(color = "black", fill = "#008EA0FF") +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  labs(x = "Cluster numer", y = "Metabolic feature count")
plot
# ggsave(plot,
#        filename = "cluster_number.pdf",
#        width = 7,
#        height = 7)
# ggsave(plot,
#        filename = "cluster_number.png",
#        width = 7,
#        height = 7)

####plot for each cluster
idx <- 1

temp_data <-
  data.s %>% as.data.frame() %>%
  t() %>% as.data.frame()

for (idx in 1:12) {
  cat(idx, " ")
  
  cluster_data <-
    cluster_info %>%
    # dplyr::filter(cluster == idx) %>%
    dplyr::select(1, 1 + idx, cluster)
  
  colnames(cluster_data)[2] <- c("membership")
  
  cluster_data <-
    cluster_data %>%
    dplyr::filter(membership > membership_cutoff)
  
  # cluster_data <-
  #   cluster_data %>%
  #   dplyr::filter(membership > 0.5)
  
  path <- paste("cluster", idx, sep = "_")
  dir.create(path)
  
  openxlsx::write.xlsx(
    cluster_data,
    file = file.path(path, paste("cluster", idx, ".xlsx", sep = "")),
    asTable = TRUE,
    overwrite = TRUE
  )
  
  temp_center <-
    centers[idx, , drop = TRUE] %>%
    data.frame(time = names(.),
               value = .,
               stringsAsFactors = FALSE) %>%
    dplyr::mutate(time = factor(time, levels = time)) %>%
    dplyr::mutate(time_point = time)
  
  temp_center$time <-
    temp_center$time %>%
    as.character() %>%
    stringr::str_split(pattern = "_") %>%
    lapply(function(x) {
      mean(as.numeric(x))
    }) %>%
    unlist()
  
  temp <-
    temp_data[cluster_data$variable_id,] %>%
    data.frame(
      membership = cluster_data$membership,
      .,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ) %>%
    tibble::rownames_to_column(var = "variable_id") %>%
    tidyr::pivot_longer(
      cols = -c(variable_id, membership),
      names_to = "time",
      values_to = "value"
    ) %>%
    dplyr::mutate(time = factor(time, levels = unique(time))) %>%
    dplyr::mutate(time_point = time)
  
  temp$time <-
    temp$time %>%
    as.character() %>%
    stringr::str_split(pattern = "_") %>%
    lapply(function(x) {
      mean(as.numeric(x))
    }) %>%
    unlist()
  
  plot <-
    temp %>%
    dplyr::arrange(membership, variable_id) %>%
    dplyr::mutate(variable_id = factor(variable_id, levels = unique(variable_id))) %>%
    ggplot(aes(time, value, group = variable_id)) +
    geom_line(aes(color = membership), alpha = 0.7) +
    theme_bw() +
    theme(
      legend.position = c(0, 1),
      legend.justification = c(0, 1),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size = 12
      ),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent", color = NA)
    ) +
    labs(
      x = "",
      y = "Z-score",
      title = paste(
        "Cluster ",
        idx,
        " (",
        nrow(cluster_data),
        " metabolic features)",
        sep = ""
      )
    ) +
    geom_line(
      mapping = aes(time, value, group = 1),
      data = temp_center,
      size = 2
    ) +
    geom_hline(yintercept = 0) +
    viridis::scale_color_viridis() +
    # scale_color_gradient(low = "white", high = "red") +
    scale_x_continuous(breaks = c(temp_center$time),
                       labels = temp_center$time_point)
  
  plot
  
  ggsave(
    plot,
    filename = file.path(path, paste("cluster", idx, ".pdf", sep = "")),
    width = 8,
    height = 7
  )
  
}

dim(cluster_data)

table(cluster_info$cluster)

cluster_info <-
  unique(cluster_info$cluster) %>%
  purrr::map(function(x) {
    temp <-
      cluster_info %>%
      # dplyr::filter(cluster == x) %>%
      dplyr::select(variable_id, paste0("X", x), cluster)
    colnames(temp)[2] <- "membership"
    temp <-
      temp %>%
      dplyr::filter(membership >= membership_cutoff)
    temp <-
      temp %>%
      dplyr::mutate(cluster_raw = cluster) %>%
      dplyr::mutate(cluster = x)
    temp
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()

cluster_info %>%
  dplyr::count(cluster)

cluster_info %>%
  dplyr::filter(membership > 0.5) %>%
  dplyr::count(cluster)

final_cluster_info <-
  cluster_info

save(final_cluster_info, file = "final_cluster_info")

openxlsx::write.xlsx(
  final_cluster_info,
  file = "final_cluster_info.xlsx",
  asTable = TRUE,
  overwrite = TRUE
)


######upsetplot to show the overlap between different clusters
library(ComplexHeatmap)
library(UpSetR)

final_cluster_info %>%
  dplyr::count(cluster)

temp_data <-
  unique(final_cluster_info$cluster) %>%
  purrr::map(function(x) {
    final_cluster_info %>%
      dplyr::filter(cluster == x) %>%
      pull(variable_id) %>%
      unique()
  })

names(temp_data) <-
  unique(final_cluster_info$cluster)

temp_data <-
  make_comb_mat(temp_data, mode = "intersect")

set_name(temp_data)
comb_name(temp_data)
comb_size(temp_data)
comb_degree(temp_data)

temp_data2 <-
  temp_data[comb_degree(temp_data) == 2]

plot <-
  UpSet(
    temp_data2,
    set_order = as.character(1:12),
    right_annotation = upset_right_annotation(temp_data2, add_numbers = TRUE),
    top_annotation = upset_top_annotation(temp_data2, add_numbers = TRUE)
  )

plot <-
  ggplotify::as.ggplot(plot)
plot
# ggsave(plot, filename = "cluster_overlap.pdf", width = 9, height = 5)


temp_data <-
  unique(final_cluster_info$cluster) %>%
  purrr::map(function(x) {
    unique(final_cluster_info$cluster) %>%
      purrr::map(function(y) {
        id1 <-
          final_cluster_info %>%
          dplyr::filter(cluster == x) %>%
          pull(variable_id)
        
        id2 <-
          final_cluster_info %>%
          dplyr::filter(cluster == y) %>%
          pull(variable_id)
        
        data.frame(
          from = x,
          to = y,
          jaccard_index =
            length(intersect(id1, id2)) / length(unique(c(id1, id2)))
        )
      }) %>%
      dplyr::bind_rows()
  }) %>%
  dplyr::bind_rows()

temp_data <-
  temp_data %>%
  tidyr::pivot_wider(names_from = to, values_from = jaccard_index) %>%
  tibble::column_to_rownames(var = "from") %>%
  as.matrix()

diag(temp_data) <- NA

library(ggheatmap)

plot <- ggheatmap(
  temp_data,
  cluster_rows = T,
  cluster_cols = T,
  border = "black"
) %>%
  ggheatmap_theme(1, theme = list(theme(axis.text.x = element_text())))
plot

# ggsave(plot, filename = "cluster_overlap_heatmap.pdf", width = 8, height = 7)

dim(object_corss_section_loess)

object_corss_section_loess@variable_info %>%
  dplyr::count(mode)

final_cluster_info %>%
  dplyr::left_join(object_corss_section_loess@variable_info, by = "variable_id") %>%
  dplyr::count(cluster, mode)

final_cluster_info %>%
  dplyr::left_join(object_corss_section_loess@variable_info, by = "variable_id") %>%
  dplyr::count(cluster, Level)
