no_source()
####This should be run in workstation

library(tidymass)
library(tidyverse)

rm(list = ls())

setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")

load("3_data_analysis/monkey_metabolomics_data_preparation/metabolite/object")

setwd("3_data_analysis/monkey_aging")

####PCA analysis

pca_object <-
  object %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(class == "Subject") %>%
  `+`(1) %>%
  log(2) %>%
  scale %>%
  massstat::run_pca()

plot <-
  massstat::pca_score_plot(
    object = object %>%
      activate_mass_dataset(what = "sample_info") %>%
      dplyr::filter(class == "Subject"),
    pca_object = pca_object,
    color_by = "age",
    frame = FALSE
  ) +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 9, name = "Reds")) +
  ggrepel::geom_text_repel(aes(label = round(age, 2)))

plot

plot1 <-
  pca_object$x %>%
  as.data.frame() %>%
  dplyr::select(PC1) %>%
  tibble::rownames_to_column(var = "sample_id") %>%
  dplyr::left_join(object@sample_info[, c("sample_id", "age", "sex")],
                   by = "sample_id") %>%
  ggplot(aes(x = PC1, y = age)) +
  geom_bar(
    aes(fill = sex),
    stat = "identity",
    width = 0.1,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = sex_color) +
  theme_base +
  labs(x = "", y = "Age (years)") +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.05))) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

library(patchwork)

plot <-
  plot1 + plot + plot_layout(ncol = 1, heights = c(1, 5))

plot
# ggsave(plot, filename = "pca_age_plot.pdf",
#        width = 9, height = 7)

####heatmap to show the different heatmap
library(ComplexHeatmap)

temp_object <-
  object %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(class == "Subject") %>%
  dplyr::filter(!is.na(age)) %>%
  `+`(1) %>%
  log(2) %>%
  scale()

library(circlize)

cor_fun <-
  circlize::colorRamp2(
    colors = RColorBrewer::brewer.pal(n = 9, name = "Reds"),
    breaks = seq(
      min(extract_sample_info(temp_object)$age),
      max(extract_sample_info(temp_object)$age),
      length.out = 9
    )
  )

h1 <-
  HeatmapAnnotation(
    age = extract_sample_info(temp_object)$age,
    sex = extract_sample_info(temp_object)$sex,
    col = list(sex = sex_color,
               age = cor_fun)
  )

variable_info <-
  extract_variable_info(temp_object)

load(
  here::here(
    "3_data_analysis/monkey_metabolomics_metabolite_annotation/hmdb_ms1.rda"
  )
)

variable_info <-
  variable_info %>%
  dplyr::left_join(hmdb_ms1@spectra.info[, c("HMDB.ID", "Super_class", "Class")],
                   by = "HMDB.ID")

h2 <-
  rowAnnotation(
    Super_class = variable_info$Super_class,
    col = list(Super_class = metabolite_super_class_color)
  )

plot <-
  ComplexHeatmap::Heatmap(
    matrix = temp_object@expression_data,
    name = "z-score",
    show_column_names = FALSE,
    show_row_names = FALSE,
    row_names_gp = gpar(cex = 0.2),
    column_names_gp = gpar(cex = 0.2),
    top_annotation = h1,
    right_annotation = h2,
    clustering_distance_columns = "euclidean",
    clustering_method_columns = "ward.D",
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "ward.D"
  )

plot <- ggplotify::as.ggplot(plot)
plot
# ggsave(plot, filename = "all_heatmap.pdf", width = 10, height = 7)

######linear mixed model to find the metabolite that changing with aging
temp_object <-
  object %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(class == "Subject" & !is.na(age)) %>%
  `+`(1) %>%
  log(2) %>%
  scale()

###cor
# cor_data <-
# seq_len(nrow(temp_object)) %>%
#   purrr::map(function(i){
#     cat(i, " ")
#     value <-
#       as.numeric(unlist(temp_object[i,,drop = TRUE]))
#     sample_info <-
#       temp_object@sample_info
#     temp_data <-
#       data.frame(sample_info, value)
#     temp_data$sex[temp_data$sex == "M"] <- 1
#     temp_data$sex[temp_data$sex == "F"] <- 0
#     temp_data$sex <- as.numeric(temp_data$sex)
#     result <-
#       lm(formula = value ~ sex, data = temp_data)
#     temp_data$value <- result$residuals
#
#     cor_result <-
#     cor.test(temp_data$value, temp_data$age, method = "spearman")
#     data.frame(cor_p = cor_result$p.value,
#                spearman_cor = cor_result$estimate)
#   }) %>%
#   dplyr::bind_rows() %>%
#   as.data.frame()
# save(cor_data, file = "cor_data")

load("cor_data")

temp_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  mutate(cor_p = cor_data$cor_p,
         spearman_cor = cor_data$spearman_cor)

temp_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(cor_p_adjust = p.adjust(cor_p, method = "fdr"))

# ######linear mixed model
# library(lme4)
# linear_mix_model_data <-
# seq_len(nrow(temp_object)) %>%
#   purrr::map(function(i){
#     cat(i, " ")
#     value <-
#       as.numeric(unlist(temp_object[i,,drop = TRUE]))
#     sample_info <-
#       temp_object@sample_info
#     temp_data <-
#       data.frame(sample_info, value)
#     temp_data$sex[temp_data$sex == "M"] <- 1
#     temp_data$sex[temp_data$sex == "F"] <- 0
#     temp_data$sex <- as.numeric(temp_data$sex)
#
#     lm_result <-
#     glm(formula = value ~ age + sex, data = temp_data)
#
#     lm_result <-
#     lm_result %>%
#       broom::tidy()
#
#     data.frame(lm_p = lm_result$p.value[2],
#                coefficient = lm_result$estimate[2])
#   }) %>%
#   dplyr::bind_rows() %>%
#   as.data.frame()
# save(linear_mix_model_data, file = "linear_mix_model_data")
load("linear_mix_model_data")

temp_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  mutate(lm_p = linear_mix_model_data$lm_p,
         coefficient = linear_mix_model_data$coefficient)

temp_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(lm_p_adjust = p.adjust(lm_p, method = "fdr"))


###volcano plot
variable_info <-
  extract_variable_info(temp_object)

top10_up_marker_name <-
  variable_info %>%
  dplyr::filter(cor_p_adjust < 0.05 & spearman_cor > 0) %>%
  dplyr::arrange(desc(spearman_cor)) %>%
  head(11) %>%
  dplyr::pull(Compound.name)

top10_down_marker_name <-
  variable_info %>%
  dplyr::filter(cor_p_adjust < 0.05 & spearman_cor < 0) %>%
  dplyr::arrange(desc(abs(spearman_cor))) %>%
  head(10) %>%
  dplyr::pull(Compound.name)

volcano_plot <-
  variable_info %>%
  mutate(
    marker = case_when(
      cor_p_adjust < 0.05 & spearman_cor > 0 ~ "Up",
      cor_p_adjust < 0.05 &
        spearman_cor < 0 ~ "Down",
      TRUE ~ "No"
    )
  ) %>%
  ggplot(aes(spearman_cor, -log(cor_p_adjust, 10))) +
  geom_point(aes(size = -log(lm_p_adjust, 10),
                 color = marker),
             alpha = 0.7) +
  theme_base +
  scale_color_manual(values = marker_color) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = -log(0.05, 10), linetype = 2) +
  labs(x = "Spearman correlation", y = "-log10(FDR)") +
  ggrepel::geom_text_repel(aes(
    label = ifelse(Compound.name %in% top10_up_marker_name,
                   Compound.name, NA)
  ), size = 3) +
  ggrepel::geom_text_repel(aes(
    label = ifelse(Compound.name %in% top10_down_marker_name,
                   Compound.name, NA)
  ), size = 3)

volcano_plot
# ggsave(volcano_plot, filename = "volcano_plot.pdf", width = 9, height = 7)

####differential expressional metabolites
dim(variable_info)
aging_markers <-
  variable_info %>%
  dplyr::filter(cor_p_adjust < 0.05)

pca_object <-
  object %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(class == "Subject") %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(variable_id %in% aging_markers$variable_id) %>%
  `+`(1) %>%
  log(2) %>%
  scale %>%
  massstat::run_pca()

plot <-
  massstat::pca_score_plot(
    object = object %>%
      activate_mass_dataset(what = "sample_info") %>%
      dplyr::filter(class == "Subject"),
    pca_object = pca_object,
    color_by = "age",
    frame = FALSE
  ) +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 9, name = "Reds")) +
  ggrepel::geom_text_repel(aes(label = round(age, 2)))

plot

plot1 <-
  pca_object$x %>%
  as.data.frame() %>%
  dplyr::select(PC1) %>%
  tibble::rownames_to_column(var = "sample_id") %>%
  dplyr::left_join(object@sample_info[, c("sample_id", "age", "sex")],
                   by = "sample_id") %>%
  ggplot(aes(x = PC1, y = age)) +
  geom_bar(
    aes(fill = sex),
    stat = "identity",
    width = 0.1,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = sex_color) +
  theme_base +
  labs(x = "", y = "Age (years)") +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.05))) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

library(patchwork)

plot <-
  plot1 + plot + plot_layout(ncol = 1, heights = c(1, 5))
plot
# ggsave(plot, filename = "pca_age_plot_with_markers.pdf",
#        width = 9, height = 7)


####heatmap to show the different heatmap
library(ComplexHeatmap)

temp_object2 <-
  temp_object %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(class == "Subject") %>%
  dplyr::filter(!is.na(age)) %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(cor_p_adjust < 0.05)

library(circlize)
cor_fun <-
  circlize::colorRamp2(
    colors = RColorBrewer::brewer.pal(n = 9, name = "Reds"),
    breaks = seq(
      min(extract_sample_info(temp_object2)$age),
      max(extract_sample_info(temp_object2)$age),
      length.out = 9
    )
  )


ha = HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = 2:3),
                                              height = unit(4, "cm")))
v = rnorm(50)

h1 <-
  HeatmapAnnotation(age = extract_sample_info(temp_object2)$age,
                    col = list(age = cor_fun))

variable_info <-
  extract_variable_info(temp_object2)

variable_info <-
  variable_info %>%
  dplyr::left_join(hmdb_ms1@spectra.info[, c("HMDB.ID", "Super_class", "Class")],
                   by = "HMDB.ID")

lipid_class = variable_info$Class
lipid_class[variable_info$Super_class != "Lipids and lipid-like molecules"] <-
  NA
h2 <-
  rowAnnotation(
    Super_class = variable_info$Super_class,
    lipid_class = lipid_class,
    spearman_cor = anno_points(variable_info$spearman_cor),
    linear_coef = anno_points(variable_info$coefficient),
    col = list(Super_class = metabolite_super_class_color,
               lipid_class = lipid_class_color)
  )

plot <-
  ComplexHeatmap::Heatmap(
    matrix = temp_object2@expression_data,
    name = "z-score",
    show_column_names = FALSE,
    show_row_names = FALSE,
    row_names_gp = gpar(cex = 0.2),
    column_names_gp = gpar(cex = 0.2),
    top_annotation = h1,
    right_annotation = h2,
    clustering_distance_columns = "euclidean",
    clustering_method_columns = "ward.D",
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "ward.D",
    border = TRUE,
    km = 2,
    column_split = 2
  )

idx <-
  ComplexHeatmap::column_order(plot)

age = extract_sample_info(temp_object2)$age

# age_sample_info <-
# 1:2 %>%
#   purrr::map(function(i){
#     temp_idx <- idx[[i]]
#     data.frame(extract_sample_info(temp_object2)[temp_idx,],
#                age_class = i)
#   }) %>%
#   dplyr::bind_rows() %>%
#   as.data.frame()
#
# save(age_sample_info, file = "age_sample_info")
load("age_sample_info")

age_sample_info$age
age_sample_info$age_class

library(pROC)
pROC_obj <-
  pROC::roc(
    age_sample_info$age_class,
    age_sample_info$age,
    smoothed = TRUE,
    # arguments for ci
    ci = TRUE,
    ci.alpha = 0.9,
    stratified = FALSE,
    # arguments for plot
    plot = TRUE,
    auc.polygon = TRUE,
    max.auc.polygon = TRUE,
    grid = TRUE,
    print.auc = TRUE,
    show.thres = TRUE
  )

pROC_obj$thresholds[which.max(pROC_obj$sensitivities +
                                pROC_obj$specificities)]

sens.ci <- ci.se(pROC_obj)

library(plotROC)

rocplot <-
  ggplot(age_sample_info, aes(m = age, d = age_class)) +
  geom_roc(n.cuts = 20, labels = FALSE)
rocplot

rocplot <-
  rocplot +
  style_roc(theme = theme_grey) +
  theme_bw()

rocplot
# ggsave(rocplot,
#        filename = "rocplot.pdf",
#        width = 5, height = 2)



age_distributation_plot <-
  age_sample_info %>%
  dplyr::mutate(age_class = as.character(age_class)) %>%
  dplyr::mutate(age_class = case_when(age_class == "1" ~ "young",
                                      age_class == "2" ~ "old")) %>%
  dplyr::mutate(age_class = factor(age_class, levels = c("young", "old"))) %>%
  ggplot(aes(x = age_class, y = age)) +
  geom_boxplot(aes(fill = age_class),
               show.legend = FALSE,
               outlier.shape = NA) +
  geom_jitter(
    aes(color = age_class),
    shape = 21,
    fill = "white",
    show.legend = FALSE,
    size = 3
  ) +
  scale_fill_manual(values = young_color) +
  scale_color_manual(values = young_color) +
  theme_base +
  labs(x = "", y = "Age (years)") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# ggsave(age_distributation_plot,
#        filename = "age_distributation_plot.pdf",
#        width = 5, height = 2)

plot <- ggplotify::as.ggplot(plot)
plot

# ggsave(plot, filename = "marker_heatmap.pdf", width = 10, height = 5.5)




#######class of up and down metabolites
temp_data <-
  temp_object2 %>%
  extract_variable_info %>%
  dplyr::filter(spearman_cor > 0) %>%
  dplyr::left_join(hmdb_ms1@spectra.info[, c("HMDB.ID", "Super_class")],
                   by = "HMDB.ID") %>%
  dplyr::select(Super_class) %>%
  dplyr::filter(!is.na(Super_class)) %>%
  dplyr::count(Super_class) %>%
  dplyr::distinct() %>%
  dplyr::mutate(rate = round(n * 100 / sum(n), 2)) %>%
  dplyr::arrange(rate) %>%
  dplyr::mutate(Super_class = factor(Super_class, levels = Super_class))

library(ggrepel)

temp_data2 <-
  temp_data %>%
  mutate(
    csum = rev(cumsum(rev(rate))),
    pos = rate / 2 + lead(csum, 1),
    pos = if_else(is.na(pos), rate / 2, pos)
  )

library(forcats)

plot <-
  ggplot(temp_data,
         aes(
           x = "" ,
           y = rate,
           fill = fct_inorder(Super_class)
         )) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = metabolite_super_class_color) +
  geom_label_repel(
    data = temp_data2,
    aes(y = pos,
        label = paste0(rate, "%")),
    size = 4.5,
    nudge_x = 1,
    show.legend = FALSE
  ) +
  guides(fill = guide_legend(title = "Super class")) +
  theme_void()

plot

# ggsave(plot,
#        filename = "up_metabolites_super_class.pdf", width = 10, height = 7)




temp_data <-
  temp_object2 %>%
  extract_variable_info %>%
  dplyr::filter(spearman_cor < 0) %>%
  dplyr::left_join(hmdb_ms1@spectra.info[, c("HMDB.ID", "Super_class")],
                   by = "HMDB.ID") %>%
  dplyr::select(Super_class) %>%
  dplyr::filter(!is.na(Super_class)) %>%
  dplyr::count(Super_class) %>%
  dplyr::distinct() %>%
  dplyr::mutate(rate = round(n * 100 / sum(n), 2)) %>%
  dplyr::arrange(rate) %>%
  dplyr::mutate(Super_class = factor(Super_class, levels = Super_class))

library(ggrepel)

temp_data2 <-
  temp_data %>%
  mutate(
    csum = rev(cumsum(rev(rate))),
    pos = rate / 2 + lead(csum, 1),
    pos = if_else(is.na(pos), rate / 2, pos)
  )

library(forcats)

plot <-
  ggplot(temp_data,
         aes(
           x = "" ,
           y = rate,
           fill = fct_inorder(Super_class)
         )) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = metabolite_super_class_color) +
  geom_label_repel(
    data = temp_data2,
    aes(y = pos,
        label = paste0(rate, "%")),
    size = 4.5,
    nudge_x = 1,
    show.legend = FALSE
  ) +
  guides(fill = guide_legend(title = "Super class")) +
  theme_void()

plot

# ggsave(plot,
#        filename = "down_metabolites_super_class.pdf", width = 10, height = 7)




#####pathway enrichment
library(metpath)
###kegg
data("kegg_hsa_pathway", package = "metpath")
kegg_hsa_pathway
#get the class of pathways
pathway_class <-
  metpath::pathway_class(kegg_hsa_pathway)

remain_idx <-
  pathway_class %>%
  purrr::map(function(x) {
    any(stringr::str_detect(x, "Disease"))
  }) %>%
  unlist() %>%
  `!`() %>%
  which()

pathway_database <-
  kegg_hsa_pathway[remain_idx]

kegg_id <-
  temp_object2 %>%
  extract_variable_info %>%
  dplyr::filter(spearman_cor > 0) %>%
  pull(KEGG.ID)

kegg_id <-
  kegg_id[!is.na(kegg_id)] %>%
  unique() %>%
  stringr::str_split(pattern = "\\{\\}") %>%
  unlist() %>%
  unique()

result_kegg_up <-
  enrich_kegg(
    query_id = kegg_id,
    query_type = "compound",
    id_type = "KEGG",
    pathway_database = pathway_database,
    p_cutoff = 0.05,
    p_adjust_method = "BH",
    threads = 3
  )

dir.create("pathway_enrichment")
save(result_kegg_up, file = "pathway_enrichment/result_kegg_up")


kegg_id <-
  temp_object2 %>%
  extract_variable_info %>%
  dplyr::filter(spearman_cor < 0) %>%
  pull(KEGG.ID)

kegg_id <-
  kegg_id[!is.na(kegg_id)] %>%
  unique() %>%
  stringr::str_split(pattern = "\\{\\}") %>%
  unlist() %>%
  unique()

result_kegg_down <-
  enrich_kegg(
    query_id = kegg_id,
    query_type = "compound",
    id_type = "KEGG",
    pathway_database = pathway_database,
    p_cutoff = 0.05,
    p_adjust_method = "BH",
    threads = 3
  )

save(result_kegg_down, file = "pathway_enrichment/result_kegg_down")
load("pathway_enrichment/result_kegg_down")
load("pathway_enrichment/result_kegg_up")

metpath::enrich_bar_plot(result_kegg_up)
metpath::enrich_bar_plot(result_kegg_down, x_axis = "p_value", cutoff = 0.05)

###HMDB
data("hmdb_pathway", package = "metpath")

#get the class of pathways
pathway_class <-
  metpath::pathway_class(hmdb_pathway)

remain_idx <-
  pathway_class %>%
  unlist() %>%
  stringr::str_detect("Metabolic;primary_pathway") %>%
  which()

pathway_database <-
  hmdb_pathway[remain_idx]

hmdb_id <-
  temp_object2 %>%
  extract_variable_info %>%
  dplyr::filter(spearman_cor > 0) %>%
  pull(HMDB.ID)

hmdb_id <-
  hmdb_id[!is.na(hmdb_id)] %>%
  unique() %>%
  stringr::str_split(pattern = "\\{\\}") %>%
  unlist() %>%
  unique()

result_hmdb_up <-
  enrich_hmdb(
    query_id = hmdb_id,
    query_type = "compound",
    id_type = "HMDB",
    pathway_database = pathway_database,
    p_cutoff = 0.05,
    p_adjust_method = "BH",
    threads = 3
  )


save(result_hmdb_up, file = "pathway_enrichment/result_hmdb_up")


hmdb_id <-
  temp_object2 %>%
  extract_variable_info %>%
  dplyr::filter(spearman_cor < 0) %>%
  pull(HMDB.ID)

hmdb_id <-
  hmdb_id[!is.na(hmdb_id)] %>%
  unique() %>%
  stringr::str_split(pattern = "\\{\\}") %>%
  unlist() %>%
  unique()

result_hmdb_down <-
  enrich_hmdb(
    query_id = hmdb_id,
    query_type = "compound",
    id_type = "HMDB",
    pathway_database = pathway_database,
    p_cutoff = 0.05,
    p_adjust_method = "BH",
    threads = 3
  )

save(result_hmdb_down, file = "pathway_enrichment/result_hmdb_down")
load("pathway_enrichment/result_hmdb_down")
load("pathway_enrichment/result_hmdb_up")

metpath::enrich_bar_plot(result_hmdb_up,
                         x_axis = "p_value",
                         cutoff = 0.05,
                         top = 10)

metpath::enrich_bar_plot(result_hmdb_down,
                         x_axis = "p_value",
                         cutoff = 0.05)

metpath::enrich_bar_plot(result_kegg_up,
                         x_axis = "p_value",
                         cutoff = 0.05)

metpath::enrich_bar_plot(result_kegg_down,
                         x_axis = "p_value",
                         cutoff = 0.05)

hmdb_up <-
  result_hmdb_up@result %>%
  dplyr::filter(p_value < 0.05)

kegg_up <-
  result_kegg_up@result %>%
  dplyr::filter(p_value < 0.05) %>%
  dplyr::filter(mapped_number > 1)

hmdb_down <-
  result_hmdb_down@result %>%
  dplyr::filter(p_value < 0.05)

kegg_down <-
  result_kegg_down@result %>%
  dplyr::filter(p_value < 0.05) %>%
  dplyr::filter(mapped_number > 1)


intersect(hmdb_up$pathway_name, hmdb_down$pathway_name)

intersect(kegg_up$pathway_name, kegg_down$pathway_name)


####network to show the enriched pathways (kegg)
edge_data_up <-
  seq_len(nrow(kegg_up)) %>%
  purrr::map(function(i) {
    id <-
      unique(stringr::str_split(kegg_up$mapped_id[i], ";")[[1]])
    data.frame(KEGG.ID = id) %>%
      dplyr::left_join(variable_info[, c("KEGG.ID", "variable_id")],
                       by = c("KEGG.ID")) %>%
      dplyr::mutate(pathway_name = kegg_up$pathway_name[i])
  }) %>%
  dplyr::bind_rows()

edge_data_down <-
  seq_len(nrow(kegg_down)) %>%
  purrr::map(function(i) {
    id <-
      unique(stringr::str_split(kegg_down$mapped_id[i], ";")[[1]])
    data.frame(KEGG.ID = id) %>%
      dplyr::left_join(variable_info[, c("KEGG.ID", "variable_id")],
                       by = c("KEGG.ID")) %>%
      dplyr::mutate(pathway_name = kegg_down$pathway_name[i])
  }) %>%
  dplyr::bind_rows()


edge_data <-
  rbind(edge_data_up,
        edge_data_down)

edge_data <-
  edge_data %>%
  dplyr::select(-c(KEGG.ID)) %>%
  dplyr::rename(from = pathway_name, to = variable_id)  %>%
  dplyr::left_join(variable_info[, c("variable_id", "spearman_cor", "cor_p_adjust")],
                   by = c("to" = "variable_id"))


node_data2 <-
  data.frame(node = unique(edge_data$to)) %>%
  dplyr::left_join(variable_info[, c(
    "variable_id",
    "Compound.name",
    "Super_class",
    "Class",
    "spearman_cor",
    "cor_p_adjust"
  )],
  by = c("node" = "variable_id")) %>%
  dplyr::mutate(node_class = "metabolite")


node_data1 <-
  data.frame(
    node = unique(edge_data$from),
    Compound.name = NA,
    Super_class = NA,
    Class = NA,
    spearman_cor = NA,
    cor_p_adjust = NA,
    node_class = "pathway"
  )

number_info <-
  seq_along(node_data1$node) %>%
  purrr::map(function(i) {
    total_edge_num <-
      sum(edge_data$from == node_data1$node[i])
    positive_edge_num <-
      edge_data %>%
      dplyr::filter(from == node_data1$node[i] &
                      spearman_cor > 0)  %>%
      nrow()
    
    negtive_edge_num <-
      edge_data %>%
      dplyr::filter(from == node_data1$node[i] &
                      spearman_cor < 0)  %>%
      nrow()
    data.frame(total_edge_num, positive_edge_num, negtive_edge_num)
  }) %>%
  dplyr::bind_rows()

node_data1 <-
  data.frame(node_data1, number_info)

node_data2$total_edge_num <-
  node_data2$positive_edge_num <-
  node_data2$negtive_edge_num <- NA

node_data <-
  rbind(node_data1,
        node_data2)

library(ggraph)
library(igraph)
library(tidygraph)

graph <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = FALSE) %>%
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

library(igraph)
library(graphlayouts)

xy =
  ggraph::create_layout(graph, layout = "fr")

V(graph)$x <- xy$x
V(graph)$y <- xy$y
library(scatterpie)

plot <-
  ggraph(
    graph,
    layout = 'manual',
    x = V(graph)$x,
    y = V(graph)$y,
    circular = FALSE
  ) +
  geom_edge_link(aes(color = spearman_cor,
                     width = -log(cor_p_adjust, 10)),
                 alpha = 1,
                 show.legend = TRUE) +
  scale_edge_color_gradient2(
    low = unname(marker_color["Down"]),
    mid = "white",
    high = unname(marker_color["Up"])
  ) +
  scale_edge_width(range = c(0.4, 1)) +
  geom_scatterpie(
    data = as_data_frame(graph, "vertices") %>% filter(node_class == "pathway"),
    cols = c("positive_edge_num", "negtive_edge_num")
  )  +
  scale_fill_manual(values = c(
    "positive_edge_num" = unname(marker_color["Up"]),
    "negtive_edge_num" = unname(marker_color["Down"])
  )) +
  ggnewscale::new_scale_fill() +
  ggraph::geom_node_point(
    aes(size = Degree,
        fill = spearman_cor),
    data = as_data_frame(graph, "vertices") %>%
      filter(node_class == "metabolite"),
    shape = 22
  ) +
  scale_fill_gradient2(low = marker_color["Down"],
                       mid = "white",
                       high = marker_color["Up"]) +
  scale_size_continuous(range = c(2, 8)) +
  coord_equal() +
  shadowtext::geom_shadowtext(
    aes(
      x = x,
      y = y,
      label = ifelse(node_class %in% "pathway",
                     node, NA)
    ),
    size = 3.5,
    check_overlap = TRUE,
    color = "black",
    bg.colour = "white",
    show.legend = FALSE
  ) +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA)
  )

plot

ggsave(plot,
       filename = "enriched_pathway.pdf",
       width = 9,
       height = 7)


######
temp_data <-
  rbind(
    kegg_up %>%
      dplyr::select(pathway_name, p_value, mapped_number) %>%
      dplyr::mutate(class = "Up"),
    kegg_down %>%
      dplyr::select(pathway_name, p_value, mapped_number) %>%
      dplyr::mutate(class = "Down") %>%
      dplyr::mutate(mapped_number = mapped_number * -1)
  )

plot <-
  temp_data %>%
  dplyr::arrange(mapped_number) %>%
  dplyr::mutate(pathway_name = factor(pathway_name, levels = unique(pathway_name))) %>%
  ggplot(aes(mapped_number, pathway_name)) +
  geom_vline(xintercept = 0) +
  geom_bar(width = 0.1, stat = "identity", aes(color = class)) +
  geom_point(aes(size = -log(p_value, 10),
                 fill = class),
             shape = 21)  +
  labs(x = "Mapped metabolite number",
       y = "") +
  theme_base +
  scale_color_manual(values = marker_color) +
  scale_fill_manual(values = marker_color) +
  scale_size_continuous(range = c(1.5, 6)) +
  theme(axis.text.y = element_text(size = 8))


ggsave(plot,
       filename = "enriched_plot_barplot.pdf",
       width = 7,
       height = 7)
