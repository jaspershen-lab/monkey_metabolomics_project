group_features4fmsea <-
  function(object_list,
           rt.match.tol = 10,
           threads = 3) {
    ###check object_list
    if (missing(object_list)) {
      stop("object_list is missing.")
    }
    
    message("Checking object_list...")
    check4fmsea(object_list = object_list)
    message("Done.")
    
    message("Combing annotation_table...")
    
    annotation_table <-
      object_list %>%
      purrr::map(function(x) {
        annotation_table <-
          x@annotation_table
        
        variable_info <-
          x@variable_info
        
        setdiff_name <-
          setdiff(colnames(variable_info), colnames(annotation_table))
        
        variable_info <-
          variable_info[, c("variable_id", setdiff_name)]
        
        annotation_table <-
          annotation_table %>%
          dplyr::left_join(variable_info, by = "variable_id") %>%
          dplyr::mutate(variable_id_old = variable_id) %>%
          dplyr::mutate(variable_id = paste(variable_id, column, polarity, sep = '_'))
        annotation_table
      }) %>%
      dplyr::bind_rows()
    
    message("Done.")
    
    rt_table <-
      annotation_table %>%
      dplyr::select(variable_id, mz, rt, mode) %>%
      dplyr::distinct(variable_id, .keep_all = TRUE) %>%
      plyr::dlply(.variables = plyr::.(mode))
    
    message("Get compound classes...")
    
    future::plan(future::multisession, workers = threads)
    
    compound_class <-
      rt_table %>%
      furrr::future_map(function(x) {
        module <-
          stats::dist(x$rt, method = "manhattan") %>%
          stats::hclust() %>%
          stats::cutree(h = rt.match.tol)
        data.frame(
          variable_id = x$variable_id,
          mode = x$mode,
          compound_class = paste(x$mode, module, sep = "_")
        )
      }, .progress = TRUE) %>%
      dplyr::bind_rows() %>%
      as.data.frame()
    message("")
    message("Done.")
    
    annotation_table <-
      annotation_table %>%
      dplyr::left_join(compound_class,
                       by = c("variable_id", "mode")) %>%
      dplyr::mutate(compound_class = paste(HMDB.ID, compound_class, sep = "_")) %>%
      dplyr::arrange(compound_class, rt)
    
    rm(list = "compound_class")
    gc()
    
    message("Score compound class...")
    
    annotation_table <-
      unique(annotation_table$HMDB.ID) %>%
      furrr::future_map(
        .f = function(temp_id) {
          # cat(temp_id, " ")
          # pb$tick()
          x <-
            annotation_table %>%
            dplyr::filter(HMDB.ID == temp_id)
          
          x <-
            unique(x$compound_class) %>%
            purrr::map(function(y) {
              z <-
                x %>%
                dplyr::filter(compound_class == y)
              z$score <- score_compound_class(z)
              z
            }) %>%
            dplyr::bind_rows()
          x
        },
        .progress = TRUE
      ) %>%
      dplyr::bind_rows()
    return(annotation_table)
  }


#' @title group features according to RT
#' @description group features according to RT
#' @author Xiaotao Shen
#' \email{shenxt@@stanford.edu}
#' @param rt rt
#' @param rt.tol rt.tol
#' @export
## group features according to RT
group_features_rt <-
  function(rt,
           rt.tol = 10) {
    rt_class <-
      lapply(rt, function(x) {
        which(rt >= x & rt < x + rt.tol)
      })
    
    rt_class <-
      lapply(seq_along(rt_class)[-1], function(i) {
        setdiff(rt_class[[i]], unlist(rt_class[1:(i - 1)]) %>% unique())
      }) %>%
      `c`(rt_class[1], .)
    
    rt_class <-
      rt_class[which(lapply(rt_class, length) != 0)]
    
    names(rt_class) <- seq_along(rt_class)
    
    rt_class <-
      purrr::map2(
        .x = rt_class,
        .y = names(rt_class),
        .f = function(x, y) {
          data.frame(rt = rt[x],
                     class = y,
                     stringsAsFactors = FALSE)
        }
      ) %>%
      do.call(rbind, .)
    rownames(rt_class) <- NULL
    return(rt_class)
  }


#' @title score_compound_class
#' @description score_compound_class
#' @author Xiaotao Shen
#' \email{shenxt@@stanford.edu}
#' @param compound_class compound_class
#' @export
score_compound_class <-
  function(compound_class) {
    score <- 0
    
    ###only feature is identified by level 1
    if (any(compound_class$Level == 1)) {
      score <- 500
      return(score)
    }
    
    ##This should be optimized in the future
    
    ###if positive adduct are +H, score add 50
    if (any(compound_class$Adduct == "(M+H)+")) {
      score <- score + 50
    }
    
    if (any(compound_class$Adduct == "(M+H)+" &
            compound_class$isotopes != "[M]")) {
      score <- score + 20
    }
    
    if (any(compound_class$polarity == "positive" &
            compound_class$Adduct != "(M+H)+")) {
      score <- score + 20
    }
    
    if (any(
      compound_class$polarity == "positive" &
      compound_class$Adduct != "(M+H)+" &
      compound_class$isotopes != "[M]"
    )) {
      score <- score + 10
    }
    
    ##if any positive peak with MS2 spectra matching Level 2, add 100
    if (any(compound_class$Level == 2 &
            compound_class$polarity == "positive")) {
      score <- score + 100
    }
    
    if (any(compound_class$Adduct == "(M-H)-")) {
      score <- score + 50
    }
    
    if (any(compound_class$Adduct == "(M-H)-" &
            compound_class$isotopes != "[M]")) {
      score <- score + 20
    }
    
    if (any(compound_class$polarity == "negative" &
            compound_class$Adduct != "(M-H)-")) {
      score <- score + 20
    }
    
    
    if (any(
      compound_class$polarity == "negative" &
      compound_class$Adduct != "(M-H)-" &
      compound_class$isotopes != "[M]"
    )) {
      score <- score + 10
    }
    
    if (any(compound_class$Level == 2 &
            compound_class$polarity == "negative")) {
      score <- score + 100
    }
    
    return(score)
  }






# #score distribution
# score_rule <-
#   data.frame(
#     rule =
#       c(
#         "Level 2 positive",
#         "Adduct M+H",
#         "Isotope (Adduct M+H)",
#         "Other positive adduct",
#         "Isotope (Other positive adduct)",
#         "Level 2 negative",
#         "Adduct M-H",
#         "Isotope (Adduct M-H)",
#         "Other negative adduct",
#         "Isotope (Other negative adduct)"
#       ),
#     score = c(100, 50, 20, 20, 10, 100, 50, 20, 20, 10),
#     stringsAsFactors = FALSE
#   )
#
# 3 * 3 * 3 * 3
#
# 2 * 3 * 3 * 2 * 3 * 3
#
# test1 <- c(0, 1)
# test2 <- list(c(1, 1), c(1, 0), c(0, 0))
#
# # comb <- vector(mode = "list", length = 3^4)
# comb <- NULL
#
# for (i in 1:2) {
#   for (j in 1:3) {
#     for (k in 1:3) {
#       for (l in 1:2) {
#         for (m in 1:3) {
#           for (n in 1:3) {
#             comb <-
#               c(comb, list(c(
#                 test1[[i]], test2[[j]], test2[[k]], test1[[l]], test2[[m]], test2[[n]]
#               )))
#           }
#         }
#       }
#     }
#   }
# }
#
# comb <-
#   comb %>%
#   do.call(cbind, .)
#
# remove_idx <-
#   apply(comb, 2, function(x) {
#     all(x == 0)
#   }) %>%
#   which()
#
# comb <- comb[, -remove_idx]
#
# score <- score_rule$score * comb
# score <- score %>% colSums()
#
# idx <- order(score, decreasing = TRUE)
#
# score <- score[idx]
#
# comb <- comb[, idx]
#
# score_rule <-
#   data.frame(score_rule, comb, stringsAsFactors = FALSE)
#
# plot1 <-
#   data.frame(index = factor(paste('X', 1:length(score), sep = ""),
#                             levels = paste('X', 1:length(score), sep = "")),
#              score,
#              stringsAsFactors = TRUE) %>%
#   ggplot(aes(index, score)) +
#   geom_point(
#     stat = "identity",
#     aes(color = as.character(score)),
#     shape = 16,
#     show.legend = FALSE
#   ) +
#   geom_segment(aes(
#     x = index,
#     y = 0,
#     xend = index,
#     yend = score,
#     color = as.character(score)
#   ),
#   show.legend = FALSE) +
#   theme_classic() +
#   scale_y_continuous(expand = expansion(mult = c(0, .05))) +
#   labs(x = "", y = "Score") +
#   theme(
#     axis.title.x = element_blank(),
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank(),
#     axis.title.y = element_text(size = 10),
#     axis.text.y = element_text(size = 10),
#     plot.background = element_rect(fill = "transparent", color = NA),
#     panel.background = element_rect(fill = "transparent", color = NA),
#     legend.background = element_rect(fill = "transparent", color = NA),
#     plot.margin = unit(c(0, 0, 0, 0), "pt")
#   )
#
# plot1
#
# plot2 <-
#   score_rule %>%
#   dplyr::select(-score) %>%
#   tidyr::pivot_longer(cols = -rule,
#                       names_to = "index",
#                       values_to = "value") %>%
#   dplyr::mutate(rule = factor(
#     rule,
#     levels =
#       c(
#         "Level 2 positive",
#         "Adduct M+H",
#         "Isotope (Adduct M+H)",
#         "Other positive adduct",
#         "Isotope (Other positive adduct)",
#         "Level 2 negative",
#         "Adduct M-H",
#         "Isotope (Adduct M-H)",
#         "Other negative adduct",
#         "Isotope (Other negative adduct)"
#       ) %>% rev()
#   )) %>%
#   ggplot(aes(index, rule)) +
#   geom_tile(aes(fill = as.character(value)),
#             color = "#FF6F00FF",
#             show.legend = FALSE) +
#   theme_bw() +
#   labs(x = "", y = "") +
#   scale_y_discrete(expand = expansion(mult = c(0, 0))) +
#   scale_fill_manual(values = c("0" = "white", "1" = "#008EA0FF")) +
#   theme(
#     panel.grid = element_blank(),
#     axis.title.x = element_blank(),
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank(),
#     axis.text.y = element_text(size = 10),
#     plot.margin = unit(c(0, 0, 0, 0), "pt")
#   )
#
# plot2
#
#
# plot3 <-
#   score_rule %>%
#   dplyr::select(rule, score) %>%
#   dplyr::mutate(rule = factor(
#     rule,
#     levels =
#       c(
#         "Level 2 positive",
#         "Adduct M+H",
#         "Isotope (Adduct M+H)",
#         "Other positive adduct",
#         "Isotope (Other positive adduct)",
#         "Level 2 negative",
#         "Adduct M-H",
#         "Isotope (Adduct M-H)",
#         "Other negative adduct",
#         "Isotope (Other negative adduct)"
#       ) %>% rev()
#   )) %>%
#   ggplot(aes(score, rule)) +
#   geom_bar(stat = "identity", color = "black", fill = "grey") +
#   theme_bw() +
#   labs(x = "Score", y = "") +
#   scale_y_discrete(expand = expansion(mult = c(0, 0))) +
#   scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
#   theme(
#     panel.grid = element_blank(),
#     axis.title.y = element_blank(),
#     axis.text.y = element_blank(),
#     axis.ticks.y = element_blank(),
#     axis.text.x = element_text(size = 10),
#     plot.margin = unit(c(0, 0, 0, 0), "pt")
#   )
#
# plot3
#
# library(patchwork)
#
# plot <-
#   (plot1 + plot_spacer() + plot_layout(widths = c(6, 1))) /
#   (plot2 + plot3 + plot_layout(widths = c(6, 1)))
#
# plot
#
# setwd(r4projects::get_project_wd())
#
# setwd("3-data_analysis/monkey_fMSEA_analysis_rplc/4_group_features/")
#
# ggsave(
#   plot,
#   filename = "score_rule.pdf",
#   width = 15,
#   height = 7,
#   bg = "transparent"
# )


mutate_compound_class <-
  function(variable_info, cor_data) {
    cor_data <-
      cor_data %>%
      dplyr::filter(from %in% variable_info$variable_id &
                      to %in% variable_info$variable_id)
    
    rt_error <-
      abs(variable_info$rt[match(cor_data$from, variable_info$variable_id)] -
            variable_info$rt[match(cor_data$to, variable_info$variable_id)])
    
    cor_data$rt_error <-
      rt_error
    
    # calculate standard error
    # standard_error <-
    #   sd(cor_data$rt_error)/sqrt(nrow(cor_data))
    rt_standard_error <- 5
    cor_standard_error <- 0.4
    
    cor_data$rt_score <-
      exp((-(cor_data$rt_error) ^ 2) / (2 * (rt_standard_error ^ 2)))
    
    cor_data$rt_score[cor_data$rt_error > 10] <- 0
    
    di <- 1 - cor_data$correlation / 2
    
    cor_data$cor_score <-
      cor_data$correlation
    
    cor_data$total_score <-
      sqrt(cor_data$rt_score * cor_data$cor_score)
    
    edge_data <-
      cor_data %>%
      dplyr::filter(total_score > 0)
    
    node_data <-
      data.frame(node = variable_info$variable_id)
    
    graph_data <-
      tidygraph::tbl_graph(nodes = node_data,
                           edges = edge_data,
                           directed = FALSE)
    
    membership <-
      igraph::cluster_fast_greedy(
        graph = graph_data,
        weights = igraph::edge_attr(graph = graph_data, name = "total_score")
      ) %>%
      igraph::membership() %>%
      as.character()
    
    # graph_data <-
    # graph_data %>%
    #   tidygraph::activate(what = "nodes") %>%
    #   dplyr::mutate(membership = membership)
    #
    # graph_data %>%
    #   ggraph::ggraph(layout = "kk") +
    #   ggraph::geom_edge_link() +
    #   # ggraph::scale_edge_colour_gradient(low = "blue", high = "red") +
    #   ggraph::geom_node_point(aes(fill = membership), shape = 21,
    #                           size = 5) +
    #   ggraph::geom_node_text(aes(label = membership))
    data.frame(variable_info, compound_class = membership)
  }
