no_source()

library(tidymass)
library(tidyverse)

rm(list = ls())

setwd(r4projects::get_project_wd())
source("1-code/tools.R")

load("3-data_analysis/monkey_metabolomics_data_preparation/object_hilic_pos")
load("3-data_analysis/monkey_metabolomics_data_preparation/object_hilic_neg")
load("3-data_analysis/monkey_metabolomics_data_preparation/object_rplc_pos")
load("3-data_analysis/monkey_metabolomics_data_preparation/object_rplc_neg")

###data quality assessment
load("3-data_analysis/monkey_metabolomics_data_preparation/metabolite/object")

setwd("3-data_analysis/monkey_metabolomics_data_overview")

# object_rplc_pos %>%
#   activate_mass_dataset(what = "sample_info") %>%
#   dplyr::filter(class != "Blank") %>%
#   massqc::massqc_report(path = "RPLC_positive_data_quality")
#
#
# object_rplc_neg %>%
#   activate_mass_dataset(what = "sample_info") %>%
#   dplyr::filter(class != "Blank") %>%
#   massqc::massqc_report(path = "RPLC_negative_data_quality")
#
# object_hilic_pos %>%
#   activate_mass_dataset(what = "sample_info") %>%
#   dplyr::filter(class != "Blank") %>%
#   massqc::massqc_report(path = "HILIC_positive_data_quality")
#
# object_hilic_neg %>%
#   activate_mass_dataset(what = "sample_info") %>%
#   dplyr::filter(class != "Blank") %>%
#   massqc::massqc_report(path = "HILIC_negative_data_quality")
#
# object %>%
#   activate_mass_dataset(what = "sample_info") %>%
#   dplyr::filter(class != "Blank") %>%
#   massqc::massqc_report(path = "data_quality")

object@annotation_table$KEGG.ID
object@annotation_table$HMDB.ID

object_rplc_pos@annotation_table %>%
  dplyr::filter(Level == 1) %>%
  dim()

object_rplc_neg@annotation_table %>%
  dplyr::filter(Level == 1) %>%
  dim()

load(
  here::here(
    "3-data_analysis/monkey_metabolomics_metabolite_annotation/hmdb_ms1.rda"
  )
)

sample_info <-
  object@sample_info %>%
  dplyr::arrange(injection.order)

sample_info$row <-
  sample_info$injection.order %/% 12 + 1

sample_info$row[sample_info$injection.order %% 12 == 0] <-
  sample_info$row[sample_info$injection.order %% 12 == 0] - 1

# sample_info$row <- rev(sample_info$row)

sample_info$column <-
  sample_info$injection.order %% 12

sample_info$column[sample_info$column == 0] <- 12

plot <-
  sample_info %>%
  dplyr::mutate(
    row = factor(as.character(row), levels = as.character(rev(sort(
      unique(row)
    )))),
    column = factor(as.character(column), levels = as.character(sort(unique(
      column
    ))))
  ) %>%
  data.frame(column, row, injection.order) %>%
  ggplot(aes(column, row)) +
  geom_tile(color = "black", aes(fill = class)) +
  geom_text(aes(label = injection.order)) +
  theme_void() +
  labs(x = "", y = "") +
  scale_fill_manual(values = sample_class_color)

plot

ggsave(plot,
       filename = "sample_injection_order.pdf",
       width = 9,
       height = 7)

plot <-
  sample_info %>%
  dplyr::filter(class == "Subject") %>%
  ggplot(aes(injection.order, age)) +
  geom_bar(stat = "identity", aes(fill = sex)) +
  scale_fill_manual(values = sex_color) +
  theme_base +
  labs(x = "Injection order", y = "Age")
plot
ggsave(plot,
       filename = "injection_order_age_sex.pdf",
       width = 9,
       height = 7)

plot <-
  metid::summary_annotation_table(object = object)

plot

ggsave(plot,
       filename = "annotation_summary.pdf",
       width = 7,
       height = 7)

variable_info <-
  extract_variable_info(object)

variable_info <-
  variable_info %>%
  dplyr::filter(!is.na(HMDB.ID)) %>%
  dplyr::mutate(HMDB.ID =
                  stringr::str_split(HMDB.ID, "\\{\\}") %>%
                  purrr::map(function(x) {
                    x[1]
                  }) %>%
                  unlist()) %>%
  dplyr::left_join(hmdb_ms1@spectra.info[, c("HMDB.ID", "Kingdom", "Super_class", "Class", "Sub_class")], by = "HMDB.ID")

####donuts plot to show the class of the metabolites
temp_data <-
  hmdb_ms1@spectra.info %>%
  dplyr::filter(stringr::str_detect(Biospecimen_locations, "Blood")) %>%
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
  ggplot(temp_data, aes(
    x = "" ,
    y = rate,
    fill = fct_inorder(Super_class)
  )) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = metabolite_super_class_color) +
  geom_label_repel(
    data = temp_data2,
    aes(y = pos, label = paste0(rate, "%")),
    size = 4.5,
    nudge_x = 1,
    show.legend = FALSE
  ) +
  guides(fill = guide_legend(title = "Super class")) +
  theme_void()

plot

# ggsave(plot, filename = "all_hmdb_blood_metabolite_super_class.pdf", width = 10, height = 7)




####donuts plot to show the class of the metabolites
temp_data <-
  variable_info %>%
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
  ggplot(temp_data, aes(
    x = "" ,
    y = rate,
    fill = fct_inorder(Super_class)
  )) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = metabolite_super_class_color) +
  geom_label_repel(
    data = temp_data2,
    aes(y = pos, label = paste0(rate, "%")),
    size = 4.5,
    nudge_x = 1,
    show.legend = FALSE
  ) +
  guides(fill = guide_legend(title = "Super class")) +
  theme_void()

plot

# ggsave(plot, filename = "annoated_metabolite_super_class.pdf", width = 10, height = 7)


####lipids
temp_data <-
  variable_info %>%
  dplyr::filter(!is.na(Super_class)) %>%
  dplyr::filter(Super_class == "Lipids and lipid-like molecules")  %>%
  dplyr::select(Class) %>%
  dplyr::count(Class) %>%
  dplyr::distinct() %>%
  dplyr::mutate(rate = round(n * 100 / sum(n), 2)) %>%
  dplyr::arrange(rate) %>%
  dplyr::mutate(Class = factor(Class, levels = Class))

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
  ggplot(temp_data, aes(x = "" , y = rate, fill = fct_inorder(Class))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = lipid_class_color) +
  geom_label_repel(
    data = temp_data2,
    aes(y = pos, label = paste0(rate, "%")),
    size = 4.5,
    nudge_x = 1,
    show.legend = FALSE
  ) +
  guides(fill = guide_legend(title = "Class")) +
  theme_void()

plot

# ggsave(plot, filename = "lipid_class.pdf", width = 10, height = 7)


####match metabolites to pathways
library(metpath)
data("hmdb_pathway")
#get the class of pathways
pathway_class =
  metpath::pathway_class(hmdb_pathway)

remain_idx = which(unlist(pathway_class) == "Metabolic;primary_pathway")

remain_idx
hmdb_pathway =
  hmdb_pathway[remain_idx]

pathway_matched_result <-
  purrr::map(
    1:length(hmdb_pathway@pathway_name),
    .f = function(i) {
      intersect_metabolites <-
        intersect(variable_info$HMDB.ID,
                  hmdb_pathway@compound_list[[i]]$HMDB.ID) %>%
        unique()
      metabolite_number <-
        length(intersect_metabolites)
      metabolite_coverage <-
        metabolite_number / nrow(hmdb_pathway@compound_list[[i]])
      data.frame(
        number = metabolite_number,
        metabolite_coverage,
        matched_metabolites = paste(intersect_metabolites, collapse = "{}")
      )
    }
  ) %>%
  dplyr::bind_rows() %>%
  as.data.frame()

pathway_matched_result$pathway_id <- hmdb_pathway@pathway_id
pathway_matched_result$pathway_name <- hmdb_pathway@pathway_name

openxlsx::write.xlsx(
  pathway_matched_result,
  file = "pathway_matched_result.xlsx",
  asTable = TRUE,
  overwrite = TRUE
)


plot <-
  pathway_matched_result %>%
  dplyr::filter(number > 0) %>%
  dplyr::arrange(desc(number)) %>%
  head(10) %>%
  dplyr::mutate(metabolite_coverage = metabolite_coverage * 100) %>%
  dplyr::mutate(pathway_name =
                  factor(pathway_name, levels = rev(pathway_name))) %>%
  ggplot(aes(x = number, y = pathway_name)) +
  geom_bar(stat = "identity",
           aes(fill = metabolite_coverage),
           color = "black") +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 9, name = "Reds")) +
  theme_base +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.05)))  +
  labs(x = "Mapped metabolite number", y = "")
plot
# ggsave(plot, filename = "mapped_pathways.pdf", width = 10, height = 7)
