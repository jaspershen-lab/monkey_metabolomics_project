plot_whole <-
  function(annotation_table,
           pathway) {
    all_hmdb_id <-
      seq_along(slot(pathway, "compound_list")) %>%
      purrr::map(function(i) {
        data.frame(
          pathway_name = slot(pathway, "pathway_name")[i],
          HMDB.ID = slot(pathway, "compound_list")[[i]]$HMDB.ID
        )
      }) %>%
      dplyr::bind_rows()
    
    temp_annotation_table <-
      annotation_table %>%
      dplyr::filter(HMDB.ID %in% all_hmdb_id$HMDB.ID)
    
    all_hmdb_id <-
      all_hmdb_id %>%
      dplyr::filter(HMDB.ID %in% temp_annotation_table$HMDB.ID)
    
    temp_data <-
      temp_annotation_table %>%
      dplyr::select(variable_id,
                    HMDB.ID,
                    score,
                    spearman_cor,
                    compound_class,
                    Level) %>%
      dplyr::mutate(Level = as.character(Level))
    
    
    ###network
    edge_data <-
      all_hmdb_id %>%
      dplyr::rename(from = pathway_name,
                    to = HMDB.ID)
    
    node_data <-
      rbind(
        data.frame(node = all_hmdb_id$pathway_name,
                   class = "pathway"),
        data.frame(node = all_hmdb_id$HMDB.ID,
                   class = "metabolite")
      ) %>%
      dplyr::distinct(.keep_all = TRUE)
    
    graph <-
      tidygraph::tbl_graph(nodes = node_data, edges = edge_data) %>%
      dplyr::mutate(Degree = tidygraph::centrality_degree(mode = 'all'))
    
    library(ggraph)
    library(igraph)
    library(tidygraph)
    V(graph)$type <- bipartite_mapping(graph)$type
    
    coords <-
      ggraph::create_layout(graph, layout = "bipartite") %>%
      dplyr::select(x, y, node, class)
    
    plot1 <-
      ggraph(graph,
             layout = 'manual',
             x = coords$x,
             y = coords$y) +
      ggraph::geom_edge_link() +
      ggraph::geom_node_point(aes(color = class),
                              show.legend = FALSE) +
      scale_color_manual(values = c("pathway" = "red",
                                    "metabolite" = "black")) +
      ggraph::theme_graph() +
      scale_x_continuous(expand = expansion(mult = c(0.04, 0.02))) +
      theme(plot.margin = margin(0, 0, 0, 0, "cm"))
    
    compound_order <-
      coords %>%
      dplyr::filter(class == "metabolite") %>%
      dplyr::arrange(x) %>%
      dplyr::pull(node)
    
    plot2 <-
      temp_data %>%
      dplyr::arrange(spearman_cor) %>%
      dplyr::mutate(variable_id = factor(variable_id, levels = unique(variable_id))) %>%
      dplyr::mutate(HMDB.ID = factor(HMDB.ID, levels = compound_order)) %>%
      ggplot(aes(HMDB.ID, variable_id)) +
      # geom_tile(aes(fill = score)) +
      geom_point(aes(size = score, color = Level),
                 alpha = 0.6) +
      theme_bw() +
      labs(x = "Metabolites", y = "") +
      scale_x_discrete(expand = expansion(mult = c(0.02, 0.02))) +
      scale_y_discrete(expand = expansion(mult = c(0.02, 0.02))) +
      theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "cm")
      ) +
      scale_color_manual(values = c(
        "1" = ggsci::pal_tron()(n = 10)[1],
        "2" = ggsci::pal_tron()(n = 10)[2],
        "3" = "black"
      )) +
      scale_size_continuous(range = c(0.1, 3))
    
    plot3 <-
      temp_data %>%
      dplyr::mutate(color = case_when(spearman_cor > 0 ~ "pos",
                                      spearman_cor <= 0 ~ "neg")) %>%
      dplyr::arrange(spearman_cor) %>%
      dplyr::mutate(variable_id = factor(variable_id, levels = unique(variable_id))) %>%
      dplyr::mutate(HMDB.ID = factor(HMDB.ID, levels = compound_order)) %>%
      ggplot(aes(spearman_cor, variable_id)) +
      geom_point(
        alpha = 0.6,
        aes(color = color),
        size = 0.5,
        show.legend = FALSE
      ) +
      scale_color_manual(values = c(
        "pos" = ggsci::pal_d3()(n = 10)[4],
        "neg" = ggsci::pal_d3()(n = 10)[1]
      )) +
      geom_vline(xintercept = 0) +
      theme_bw() +
      scale_x_discrete(expand = expansion(mult = c(0.02, 0.02))) +
      scale_y_discrete(expand = expansion(mult = c(0.02, 0.02))) +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "cm")
      ) +
      labs(x = "Correlation", y = "Metabolic features")
    
    library(patchwork)
    
    plot <-
      (patchwork::plot_spacer() + plot1 + plot_layout(widths = c(1, 5))) /
      (plot3 + plot2 + plot_layout(widths = c(1, 5))) +
      plot_layout(heights = c(1, 5))
    
    plot
    
  }