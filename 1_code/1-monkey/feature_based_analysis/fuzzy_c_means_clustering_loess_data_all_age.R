no_source()

rm(list = ls())

setwd(r4projects::get_project_wd())

library(tidyverse)
library(tidymass)

###load("data)
load("3-data_analysis/monkey_feature_analysis/loess_fit/object")

dir.create("3-data_analysis/monkey_feature_analysis/fuzzy_clustering_loess_data",
           recursive = TRUE)

setwd("3-data_analysis/monkey_feature_analysis/fuzzy_clustering_loess_data")

object

object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::count(mode)

###remove samples without age
object <-
  object %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(!is.na(age))

####clustering
library(Mfuzz)

object[1, , drop = TRUE] %>%
  unlist() %>%
  density() %>%
  plot()

######cluster all variables using all ages
setwd(r4projects::get_project_wd())
dir.create("3-data_analysis/monkey_feature_analysis/fuzzy_clustering_loess_data/all_age")
setwd("3-data_analysis/monkey_feature_analysis/fuzzy_clustering_loess_data/all_age")

object@sample_info$age %>% range

object@sample_info %>%
  ggplot(aes(age)) +
  geom_histogram(binwidth = 1, color = "black")

object@sample_info$age %>% range

age_index <-
  seq(0, 21, by = 1) %>%
  purrr::map(function(x) {
    data.frame(
      from = x,
      to = x + 1,
      number =  sum(object@sample_info$age > x &
                      object@sample_info$age <= x + 1)
    )
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame() %>%
  dplyr::filter(number > 0)

temp_data <-
  log(object@expression_data + 1, 2)

temp_data_mean <-
  seq_len(nrow(age_index)) %>%
  purrr::map(function(i) {
    idx <- as.numeric(age_index[i,])
    temp_data[, which(object@sample_info$age > idx[1] &
                        object@sample_info$age <= idx[2]), drop = FALSE] %>%
      apply(1, mean)
  }) %>%
  dplyr::bind_cols() %>%
  as.data.frame()

colnames(temp_data_mean) <-
  paste(age_index$from, age_index$to, sep = "_")

rownames(temp_data_mean) <- object@variable_info$variable_id

idx = 3

object[idx, , drop = TRUE] %>%
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

cluster <- 15

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
  mfrow = c(4, 4),
  time.labels = time,
  new.window = FALSE
)

library(ComplexHeatmap)

###
cluster_color <-
  ggsci::pal_jama()(n = 7)[1:15]

names(cluster_color) <- as.character(1:15)

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

# ggsave(plot,
#        filename = "cluster_center_plot.pdf",
#        width = 7,
#        height = 7)

centers <- c$centers

cluster_info <-
  data.frame(
    variable_id = names(c$cluster),
    c$membership,
    cluster = c$cluster,
    stringsAsFactors = FALSE
  ) %>%
  arrange(cluster)

# openxlsx::write.xlsx(
#   x = cluster_info,
#   file = "cluster_info.xlsx",
#   asTable = TRUE,
#   overwrite = TRUE
# )

####plot for each cluster
idx <- 1

temp_data <-
  data.s %>% as.data.frame() %>%
  t() %>% as.data.frame()

# for (idx in 1:15) {
#   cat(idx, " ")
#   cluster_data <-
#     cluster_info %>%
#     dplyr::filter(cluster == idx) %>%
#     dplyr::select(1, 1 + idx, cluster)
# 
#   colnames(cluster_data)[2] <- c("membership")
# 
#   cluster_data <-
#     cluster_data %>%
#     dplyr::filter(membership > 0.5)
# 
#   path <- paste("cluster", idx, sep = "_")
#   dir.create(path)
# 
#   openxlsx::write.xlsx(
#     cluster_data,
#     file = file.path(path, paste("cluster", idx, ".xlsx", sep = "")),
#     asTable = TRUE,
#     overwrite = TRUE
#   )
# 
#   temp_center <-
#     centers[idx, , drop = TRUE] %>%
#     data.frame(time = names(.),
#                value = .,
#                stringsAsFactors = FALSE) %>%
#     dplyr::mutate(time = factor(time, levels = time)) %>%
#     dplyr::mutate(time_point = time)
# 
#   temp_center$time <-
#     temp_center$time %>%
#     as.character() %>%
#     stringr::str_split(pattern = "_") %>%
#     lapply(function(x) {
#       mean(as.numeric(x))
#     }) %>%
#     unlist()
# 
#   temp <-
#     temp_data[cluster_data$variable_id,] %>%
#     data.frame(
#       membership = cluster_data$membership,
#       .,
#       stringsAsFactors = FALSE,
#       check.names = FALSE
#     ) %>%
#     tibble::rownames_to_column(var = "variable_id") %>%
#     tidyr::pivot_longer(
#       cols = -c(variable_id, membership),
#       names_to = "time",
#       values_to = "value"
#     ) %>%
#     dplyr::mutate(time = factor(time, levels = unique(time))) %>%
#     dplyr::mutate(time_point = time)
# 
#   temp$time <-
#     temp$time %>%
#     as.character() %>%
#     stringr::str_split(pattern = "_") %>%
#     lapply(function(x) {
#       mean(as.numeric(x))
#     }) %>%
#     unlist()
# 
#   plot <-
#     temp %>%
#     ggplot(aes(time, value, group = variable_id)) +
#     geom_hline(yintercept = 0) +
#     geom_line(aes(color = membership), alpha = 0.7) +
#     theme_bw() +
#     theme(
#       legend.position = c(0, 1),
#       legend.justification = c(0, 1),
#       panel.grid.minor = element_blank(),
#       axis.title = element_text(size = 13),
#       axis.text = element_text(size = 12),
#       axis.text.x = element_text(
#         angle = 45,
#         hjust = 1,
#         vjust = 1,
#         size = 12
#       ),
#       panel.background = element_rect(fill = "transparent", color = NA),
#       plot.background = element_rect(fill = "transparent", color = NA),
#       legend.background = element_rect(fill = "transparent", color = NA)
#     ) +
#     labs(
#       x = "",
#       y = "Z-score",
#       title = paste(
#         "Cluster ",
#         idx,
#         " (",
#         nrow(cluster_data),
#         " metabolic features)",
#         sep = ""
#       )
#     ) +
#     geom_line(
#       mapping = aes(time, value, group = 1),
#       data = temp_center,
#       size = 2
#     ) +
#     # scale_color_gradientn(colours = viridis::cividis(n = 10))
#     scale_x_continuous(breaks = c(temp_center$time), labels = temp_center$time_point)
# 
#   plot
# 
#   ggsave(
#     plot,
#     filename = file.path(path, paste("cluster", idx, ".pdf", sep = "")),
#     width = 8,
#     height = 7
#   )
# 
# }

dim(cluster_data)

table(cluster_info$cluster)

cluster_info <-
  unique(cluster_info$cluster) %>%
  purrr::map(function(x) {
    temp <-
      cluster_info %>%
      dplyr::filter(cluster == x) %>%
      dplyr::select(variable_id, paste0("X", x), cluster)
    colnames(temp)[2] <- "membership"
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
