no_source()

library(tidymass)
library(tidyverse)

rm(list = ls())

setwd(r4projects::get_project_wd())
source("1-code/tools.R")
load("3-data_analysis/monkey_metabolomics_data_preparation/metabolite/object")

setwd("3-data_analysis/monkey_demongraphics_data/")

sample_info <-
  extract_sample_info(object) %>%
  dplyr::filter(group == "Subject")

library(kinship2)

temp_data1 <-
  sample_info %>%
  # dplyr::filter(mother != "unknown" | Father != "unknown") %>%
  dplyr::filter(mother %in% sample_info$sample_id |
                  Father %in% sample_info$sample_id) %>%
  dplyr::select(-c(group, class, dob, injection.order))

miss_id1 <-
  unique(temp_data1$Father)[which(!unique(temp_data1$Father) %in% temp_data1$sample_id)]

miss_id2 <-
  unique(temp_data1$mother)[which(!unique(temp_data1$mother) %in% temp_data1$sample_id)]

temp_data2 <-
  data.frame(
    sample_id = miss_id1,
    sex = "M",
    mother = NA,
    Father = NA,
    age = NA
  )

temp_data3 <-
  data.frame(
    sample_id = miss_id2,
    sex = "F",
    mother = NA,
    Father = NA,
    age = NA
  )

temp_data <-
  rbind(temp_data1,
        temp_data2,
        temp_data3)

df <-
  data.frame(
    id = temp_data$sample_id,
    sex = ifelse(temp_data$sex == "F", "female", "male"),
    dadid = temp_data$Father,
    momid = temp_data$mother,
    famid = 1
  )

foo <-
  pedigree(
    id = df$id,
    dadid = df$dadid,
    momid = df$momid,
    sex = df$sex,
    # relation = relation1,
    famid = df$famid
  )

ped <- foo['1']
plot(ped, cex = 0.3)


#####circos plot
library(circlize)
sex <- sample_info$sex
age <- sample_info$age

n <- nrow(sample_info)

df <-
  data.frame(
    factors = sample_info$sample_id,
    x = 1,
    y = 1,
    sample_info,
    stringsAsFactors = TRUE
  )

circos.par(
  "track.height" = 0.2,
  start.degree = 90,
  clock.wise = TRUE,
  gap.after = c(rep(0, nrow(df) - 1), 90),
  cell.padding = c(0, 0, 0, 0)
)

circos.initialize(factors = df$factors,
                  x = df$x,
                  xlim = c(0.5, 1.5))


## sex
temp_sex <- df$sex
temp_sex[is.na(temp_sex)] <- "grey"
temp_sex[temp_sex == "F"] <- sex_color["F"]
temp_sex[temp_sex == "M"] <- sex_color["M"]

circos.track(
  factors = df$factors,
  # x = df$x,
  y = df$y,
  ylim = c(0, 1),
  bg.border = "black",
  # bg.col = NA,
  track.height = 0.1,
  panel.fun = function(x, y) {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    #text direction (dd) and adjusmtents (aa)
    theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
    dd <-
      ifelse(theta < 90 || theta > 270, "clockwise", "reverse.clockwise")
    aa = c(0.5, 1)
    # if(theta < 90 || theta > 270)  aa = c(0, 0.5)
    #plot country labels
    circos.text(
      x = 1,
      y = 2,
      labels = name,
      facing = "clockwise",
      niceFacing = TRUE,
      cex = 0.8
      # adj = aa
    )
    
    circos.rect(
      xleft = xlim[1],
      ybottom = ylim[1],
      xright = xlim[2],
      ytop = ylim[2],
      col = temp_sex[i],
      bg.border = "black"
    )
  }
)

##mother age
range(df$age, na.rm = TRUE)
temp_value <- df$age

circos.track(
  factors = df$factors,
  # x = df$x,
  y = temp_value,
  ylim = c(0, 1.1 * max(temp_value, na.rm = TRUE)),
  bg.border = "black",
  # bg.col = NA,
  track.height = 0.25,
  panel.fun = function(x, y) {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    circos.yaxis(
      side = "left",
      at = c(0,
             round((
               min(temp_value, na.rm = TRUE) + max(temp_value, na.rm = TRUE)
             ) / 2, 2),
             round(max(
               temp_value, na.rm = TRUE
             ), 2)),
      sector.index = get.all.sector.index()[1],
      labels.cex = 0.4,
      labels.niceFacing = FALSE
    )
    
    circos.lines(
      x = mean(xlim, na.rm = TRUE),
      y =  temp_value[i],
      pch = 16,
      cex = 8,
      type = "h",
      col = ggsci::pal_aaas()(n = 10)[4],
      lwd = 2
    )
    
    circos.points(
      x = mean(xlim),
      y =  temp_value[i],
      pch = 16,
      cex = 0.8,
      col = ggsci::pal_aaas()(n = 10)[4]
    )
  }
)

####links for family tree
temp_data <-
  sample_info %>%
  dplyr::filter(mother %in% sample_info$sample_id)

for (i in 1:nrow(sample_info)) {
  cat(i, " ")
  if (is.na(match(sample_info$mother[i], sample_info$sample_id))) {
    next
  } else{
    circos.link(
      sector.index1 = sample_info$mother[i],
      point1 = 1,
      sector.index2 = sample_info$sample_id[i],
      point2 = 1,
      directional = 1,
      arr.width = 0.1,
      arr.length = 0.2
    )
  }
  
}

temp_data <-
  sample_info %>%
  dplyr::filter(Father %in% sample_info$sample_id)

for (i in 1:nrow(sample_info)) {
  cat(i, " ")
  if (is.na(match(sample_info$Father[i], sample_info$sample_id))) {
    next
  } else{
    circos.link(
      sector.index1 = sample_info$Father[i],
      point1 = 1,
      sector.index2 = sample_info$sample_id[i],
      point2 = 1,
      directional = 1,
      arr.width = 0.1,
      arr.length = 0.2,
      col = ggsci::pal_jama()(n = 5)[2]
    )
  }
  
}

#####age
age <-
  df$age
library(gghalves)
plot_age <-
  age %>%
  data.frame(class = "class", value = .) %>%
  ggplot(aes(x = class, y = value)) +
  gghalves::geom_half_dotplot(
    binaxis = "y",
    stackdir = "down",
    shape = 21,
    fill = ggsci::pal_aaas()(n = 10)[4],
    size = 0.6,
    binwidth = 0.2
  ) +
  gghalves::geom_half_boxplot(side = "l", fill = alpha("skyblue", 0.5)) +
  gghalves::geom_half_violin(side = "r", fill = alpha("skyblue", 0.5)) +
  theme_bw() +
  labs(x = "", y = "") +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
plot_age
ggsave(plot_age,
       filename = "plot_age.pdf",
       width = 3,
       height = 10)

sex <-
  df$sex

plot_sex <-
  sex %>%
  data.frame(class = "class", value = .) %>%
  dplyr::mutate(value = factor(value, levels = c("F", "M"))) %>%
  ggplot(aes(x = class)) +
  geom_bar(
    aes(fill = value),
    color = "black",
    position = "stack",
    show.legend = FALSE,
    width = 2
  ) +
  scale_fill_manual(values = sex_color) +
  theme_bw() +
  labs(x = "", y = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

plot_sex

ggsave(plot_sex,
       filename = "plot_sex.pdf",
       width = 1.2,
       height = 10)

plot <-
  sample_info %>%
  ggplot(aes(age)) +
  ggplot2::geom_histogram(binwidth = 1, color = "black", aes(fill = sex)) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  theme_base +
  labs(x = "Age (years)", y = "Frequency") +
  scale_fill_manual(values = sex_color)
plot
# ggsave(plot,
#        filename = "distributation_plot.pdf",
#        width = 7,
#        height = 3)


####family tree using ggraph
####father
library(igraph)
library(ggraph)
library(tidygraph)

temp_data <-
  sample_info %>%
  dplyr::select(sample_id, sex, mother, Father, age) %>%
  dplyr::rename(father = Father)

miss_mother_id <-
  unique(temp_data$mother)[!unique(temp_data$mother) %in% temp_data$sample_id]

miss_father_id <-
  unique(temp_data$father)[!unique(temp_data$father) %in% temp_data$sample_id]

node_data <-
  rbind(
    data.frame(
      sample_id = miss_mother_id,
      sex = "F",
      mother = NA,
      father = NA,
      age = NA
    ),
    data.frame(
      sample_id = miss_father_id,
      sex = "M",
      mother = NA,
      father = NA,
      age = NA
    ),
    temp_data
  ) %>%
  dplyr::rename(node = sample_id) %>%
  dplyr::distinct(node, .keep_all = TRUE) %>%
  dplyr::filter(node != "unknown")

edge_data <-
  temp_data %>%
  dplyr::select(from = father,
                to = sample_id) %>%
  dplyr::filter(from != "unknown")

node_data <-
  node_data %>%
  dplyr::filter(node %in% unique(c(edge_data$from, edge_data$to)))

edge_data <-
  edge_data %>%
  dplyr::filter(from %in% node_data$node & to %in% node_data$node)

node_data <-
  node_data %>%
  dplyr::filter(node %in% unique(c(edge_data$from, edge_data$to)))

node_data <-
  node_data %>%
  dplyr::mutate(in_dataset =
                  case_when(is.na(age) ~ "No",
                            TRUE ~ "Yes"))
node_data$age[is.na(node_data$age)] <- 5
temp_graph <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = TRUE)

plot <-
  ggraph(temp_graph, layout = "dendrogram") +
  geom_edge_diagonal(arrow = arrow(length = unit(2, 'mm'),
                                   type = "closed")) +
  geom_node_point(aes(fill = sex,
                      size = age,
                      shape = in_dataset)) +
  scale_shape_manual(values = c("Yes" = 21, "No" = 22)) +
  scale_fill_manual(values = sex_color) +
  scale_size_continuous(range = c(1, 5)) +
  ggraph::theme_graph() +
  geom_node_text(
    aes(label = node),
    angle = 45,
    hjust = 1,
    vjust = 1,
    size = 2
  )
extrafont::loadfonts()
plot
# ggsave(plot,
#        filename = "father_tree.pdf",
#        width = 7,
#        height = 3)


####mohter
library(igraph)
library(ggraph)
library(tidygraph)

temp_data <-
  sample_info %>%
  dplyr::select(sample_id, sex, mother, Father, age) %>%
  dplyr::rename(father = Father)

miss_mother_id <-
  unique(temp_data$mother)[!unique(temp_data$mother) %in% temp_data$sample_id]

miss_father_id <-
  unique(temp_data$father)[!unique(temp_data$father) %in% temp_data$sample_id]

node_data <-
  rbind(
    data.frame(
      sample_id = miss_mother_id,
      sex = "F",
      mother = NA,
      father = NA,
      age = NA
    ),
    data.frame(
      sample_id = miss_father_id,
      sex = "M",
      mother = NA,
      father = NA,
      age = NA
    ),
    temp_data
  ) %>%
  dplyr::rename(node = sample_id) %>%
  dplyr::distinct(node, .keep_all = TRUE) %>%
  dplyr::filter(node != "unknown")

edge_data <-
  temp_data %>%
  dplyr::select(from = mother,
                to = sample_id) %>%
  dplyr::filter(from != "unknown")

node_data <-
  node_data %>%
  dplyr::filter(node %in% unique(c(edge_data$from, edge_data$to)))

edge_data <-
  edge_data %>%
  dplyr::filter(from %in% node_data$node & to %in% node_data$node)

node_data <-
  node_data %>%
  dplyr::filter(node %in% unique(c(edge_data$from, edge_data$to)))

node_data <-
  node_data %>%
  dplyr::mutate(in_dataset =
                  case_when(is.na(age) ~ "No",
                            TRUE ~ "Yes"))
node_data$age[is.na(node_data$age)] <- 5

temp_graph <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = TRUE)

plot <-
  ggraph(temp_graph, layout = "tree") +
  geom_edge_diagonal(arrow = arrow(length = unit(2, 'mm'),
                                   type = "closed")) +
  geom_node_point(aes(fill = sex,
                      size = age,
                      shape = in_dataset)) +
  scale_shape_manual(values = c("Yes" = 21, "No" = 22)) +
  scale_fill_manual(values = sex_color) +
  scale_size_continuous(range = c(1, 5)) +
  ggraph::theme_graph() +
  geom_node_text(
    aes(label = node),
    angle = 45,
    hjust = 1,
    vjust = 1,
    size = 2
  )
plot
extrafont::loadfonts()

# ggsave(plot,
#        filename = "mother_tree.pdf",
#        width = 7,
#        height = 6)

