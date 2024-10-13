library(tidyverse)
library(ggplot2)

metabolite_super_class <-
  c(
    "Acetylides",
    "Alkaloids and derivatives",
    "Allenes",
    "Benzenoids",
    "Homogeneous metal compounds",
    "Homogeneous non-metal compounds",
    "Hydrocarbon derivatives",
    "Hydrocarbons",
    "Lignans, neolignans and related compounds",
    "Lipids and lipid-like molecules",
    "Miscellaneous inorganic compounds",
    "Mixed metal/non-metal compounds",
    "Nucleosides, nucleotides, and analogues",
    "Organic 1,3-dipolar compounds",
    "Organic acids and derivatives",
    "Organic compounds",
    "Organic nitrogen compounds",
    "Organic oxygen compounds",
    "Organic Polymers",
    "Organic salts",
    "Organohalogen compounds",
    "Organoheterocyclic compounds",
    "Organometallic compounds",
    "Organophosphorus compounds",
    "Phenylpropanoids and polyketides"
  )

metabolite_super_class_color <-
  colorRampPalette(colors = ggsci::pal_lancet()(n = 9)[1:9])(n = length(metabolite_super_class))

names(metabolite_super_class_color) <- metabolite_super_class


lipid_class <- c(
  "Fatty Acyls",
  "Glycerolipids",
  "Glycerophospholipids",
  "Prenol lipids",
  "Sphingolipids",
  "Steroids and steroid derivatives"
)

lipid_class_color <- RColorBrewer::brewer.pal(n = 6, name = "Dark2")
names(lipid_class_color) <- lipid_class

sex_color <-
  c("F" = ggsci::pal_aaas()(n = 5)[2],
    "M" = ggsci::pal_aaas()(n = 5)[1])

library(ggplot2)

theme_base <-
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    panel.grid.minor = element_blank()
  )


sample_class_color <-
  c(
    "Blank" = "#6F99ADFF",
    "QC" = "#BC3C29FF",
    "Subject" = "#E18727FF"
  )

marker_color <-
  c(
    "Up" = ggsci::pal_futurama()(n = 9)[6],
    "Down" = ggsci::pal_futurama()(n = 9)[3],
    "No" = ggsci::pal_futurama()(n = 9)[9]
  )


young_color <-
  c(young = ggsci::pal_tron()(n = 9)[2],
    old = ggsci::pal_material()(n = 9)[9])

# from <- "108F"
# to <- "26B"
#
# get_family_score(graph = graph,
#                  from = "72X",
#                  to = "DBL7")


get_node_is_children_from_two_nodes <-
  function(graph,
           child_node_idx,
           parient_node_idx1,
           parient_node_idx2) {
    edge1 <-
      igraph::shortest_paths(
        graph = graph,
        from = parient_node_idx1,
        to = child_node_idx,
        mode = "out"
      )$vpath[[1]]
    
    edge2 <-
      igraph::shortest_paths(
        graph = graph,
        from = parient_node_idx2,
        to = child_node_idx,
        mode = "out"
      )$vpath[[1]]
    
    if (length(edge1) != 0 & length(edge2) != 0) {
      return(TRUE)
    } else{
      return(FALSE)
    }
    
  }




get_parents <-
  function(id) {
    result <-
      id %>%
      purrr::map(function(x) {
        mother <- df[df$id == x, "mother"]
        father <- df[df$id == x, "father"]
        
        id2 <-
          c(mother, father)
        
        id2 <- id2[id2 != "unknown"]
        
        if (length(id2) == 0) {
          return(NULL)
        }
        data.frame(id1 = id,
                   id2 = id2,
                   level = 1)
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    
    if (nrow(result) == 0) {
      return(NULL)
    }
    
    result %>%
      dplyr::arrange(id1)
  }

get_ancestors <-
  function(id, level) {
    temp1 <-
      get_parents(id)
    
    if (!is.null(temp1)) {
      temp2 <-
        get_parents(unique(temp1$id2))
    }else{
      temp2 <- NULL
    }
    
    if (!is.null(temp2)) {
      temp3 <-
        get_parents(unique(temp2$id2))
    }else{
      temp3 <- NULL
    }
    
    if (!is.null(temp3)) {
      temp4 <-
        get_parents(unique(temp3$id2))
    } else{
      temp4 <- NULL
    }
    
    if (!is.null(temp4)) {
      temp5 <-
        get_parents(unique(temp4$id2))
    } else{
      temp5 <- NULL
    }
    
    if (!is.null(temp5)) {
      temp6 <-
        get_parents(unique(temp5$id2))
    } else{
      temp6 <- NULL
    }
    
    if (level == 1) {
      return(temp1)
    }
    
    if (level == 2) {
      if (is.null(temp2)) {
        return(NULL)
      }
      
      result <-
        temp1[, 1:2] %>%
        dplyr::left_join(temp2[, 1:2], by = c("id2" = "id1")) %>%
        dplyr::select(-id2) %>%
        dplyr::rename(id2 = "id2.y") %>%
        dplyr::mutate(level = 2) %>%
        dplyr::filter(!is.na(id1) & !is.na(id2))
      
      return(result)
    }
    
    if (level == 3) {
      if (is.null(temp3)) {
        return(NULL)
      }
      
      result <-
        temp1[, 1:2] %>%
        dplyr::left_join(temp2[, 1:2], by = c("id2" = "id1")) %>%
        dplyr::select(-id2) %>%
        dplyr::rename(id2 = "id2.y") %>%
        dplyr::left_join(temp3[, 1:2], by = c("id2" = "id1")) %>%
        dplyr::select(-id2) %>%
        dplyr::rename(id2 = "id2.y") %>%
        dplyr::mutate(level = 3) %>%
        dplyr::filter(!is.na(id1) & !is.na(id2))
      
      return(result)
    }
    
    if (level == 4) {
      if (is.null(temp4)) {
        return(NULL)
      }
      
      
      result <-
        temp1[, 1:2] %>%
        dplyr::left_join(temp2[, 1:2], by = c("id2" = "id1")) %>%
        dplyr::select(-id2) %>%
        dplyr::rename(id2 = "id2.y") %>%
        dplyr::left_join(temp3[, 1:2], by = c("id2" = "id1")) %>%
        dplyr::select(-id2) %>%
        dplyr::rename(id2 = "id2.y") %>%
        dplyr::left_join(temp4[, 1:2], by = c("id2" = "id1")) %>%
        dplyr::select(-id2) %>%
        dplyr::rename(id2 = "id2.y") %>%
        dplyr::mutate(level = 4) %>%
        dplyr::filter(!is.na(id1) & !is.na(id2))
      
      return(result)
      
      
    }
    
    if (level == 5) {
      if (is.null(temp5)) {
        return(NULL)
      }
      
      result <-
        temp1[, 1:2] %>%
        dplyr::left_join(temp2[, 1:2], by = c("id2" = "id1")) %>%
        dplyr::select(-id2) %>%
        dplyr::rename(id2 = "id2.y") %>%
        dplyr::left_join(temp3[, 1:2], by = c("id2" = "id1")) %>%
        dplyr::select(-id2) %>%
        dplyr::rename(id2 = "id2.y") %>%
        dplyr::left_join(temp4[, 1:2], by = c("id2" = "id1")) %>%
        dplyr::select(-id2) %>%
        dplyr::rename(id2 = "id2.y") %>%
        dplyr::left_join(temp5[, 1:2], by = c("id2" = "id1")) %>%
        dplyr::select(-id2) %>%
        dplyr::rename(id2 = "id2.y") %>%
        dplyr::mutate(level = 5) %>%
        dplyr::filter(!is.na(id1) & !is.na(id2))
      
      return(result)
      
    }
    
  }



find_relationship_internal <-
  function(df, id1, id2) {
    temp1 <-
      get_ancestors(id = id1, level = 1)
    
    temp2 <-
      get_ancestors(id = id1, level = 2)
    
    temp3 <-
      get_ancestors(id = id1, level = 3)
    
    temp4 <-
      get_ancestors(id = id1, level = 4)
    
    temp5 <-
      get_ancestors(id = id1, level = 5)
    
    temp6 <-
      get_ancestors(id = id1, level = 6)
    
    temp <-
      rbind(temp1, temp2, temp3, temp4, temp5, temp6)
    
    if (is.null(temp)) {
      return(data.frame(
        from = id1,
        to = id2,
        family_distance = Inf
      ))
    }
    
    result <-
      temp[temp$id2 == id2, , drop = FALSE] %>%
      dplyr::rename(from = id1,
                    to = id2,
                    family_distance = level)
    
    if (nrow(result) == 0) {
      return(data.frame(
        from = id1,
        to = id2,
        family_distance = Inf
      ))
    }
    
    return(result)
    
  }

find_relationship <-
  function(df, id1, id2) {
    result1 <-
      find_relationship_internal(df = df, id1 = id1, id2 = id2)
    result2 <-
      find_relationship_internal(df = df, id1 = id2, id2 = id1) %>%
      dplyr::rename(from1 = to, to1 = from) %>%
      dplyr::rename(from = from1, to = to1) %>%
      dplyr::select(from, to, family_distance)
    result <-
      rbind(result1, result2) %>%
      dplyr::arrange(family_distance) %>%
      dplyr::slice_head(n = 1)
    result
  }


family_distance_color <-
  c(
    "1" = ggsci::pal_gsea()(n = 10)[10],
    "2" = ggsci::pal_gsea()(n = 10)[7],
    "3" = ggsci::pal_gsea()(n = 10)[6],
    "4" = ggsci::pal_gsea()(n = 10)[1]
  )

convert_mass_dataset2pmd <-
  function(object) {
    data <-
      as.matrix(massdataset::extract_expression_data(object = object))
    group <- massdataset::extract_sample_info(object) %>%
      dplyr::select(sample_id, group) %>%
      dplyr::rename(sample_name = sample_id, sample_group = group) %>%
      as.data.frame()
    list(
      data = data,
      group = group,
      mz = as.numeric(massdataset::extract_variable_info(object)$mz),
      rt = as.numeric(massdataset::extract_variable_info(object)$rt)
    )
  }




polarity_color <-
  c(
    "positive" = "#FC4E07",
    "negative" = "#00AFBB",
    "mix" = "#E7B800"
  )



optimize_loess_span <-
  function(x, y, span_range = seq(0.2, 0.6, 0.1)) {
    span_rmse <-
      purrr::map(span_range, function(span) {
        # cat(span, " ")
        temp_data =
          data.frame(x, y)
        
        prediction <-
          purrr::map(
            2:(nrow(temp_data) - 1),
            .f = function(idx) {
              temp_result =
                loess(formula = y ~ x,
                      data = temp_data[-idx, ],
                      span = span)
              prediction =
                try(predict(object = temp_result, newdata = temp_data[idx, -2, drop = FALSE]))
              
              if (class(prediction) == "try-error") {
                data.frame(real = temp_data$y[idx], prediction = NA)
              } else{
                data.frame(real = temp_data$y[idx],
                           prediction = as.numeric(prediction))
              }
            }
          ) %>%
          dplyr::bind_rows() %>%
          dplyr::filter(!is.na(prediction))
        
        if (all(is.na(prediction$prediction))) {
          temp_rmse = NA
        } else{
          temp_rmse = sqrt(sum((
            prediction$real - prediction$prediction
          ) ^ 2) / nrow(prediction))
        }
        
        data.frame(span = span, rmse = temp_rmse)
      }) %>%
      dplyr::bind_rows()
    
    # span_rmse
    plot <-
      data.frame(x, y) %>%
      ggplot(aes(x, y)) +
      geom_point(size = 5) +
      # geom_line() +
      base_theme
    
    span_rmse =
      span_rmse %>%
      dplyr::filter(!is.na(rmse))
    idx = which.min(span_rmse$rmse)
    # for(i in 1:nrow(span_rmse)){
    plot =
      plot +
      geom_smooth(
        se = FALSE,
        span = span_rmse$span[idx],
        color = ggsci::pal_lancet()(n = 9)[idx]
      )
    # }
    
    plot =
      plot +
      ggplot2::ggtitle(label = paste("Span: ", span_rmse$span[idx])) +
      theme(title = element_text(colour = ggsci::pal_lancet()(n = 9)[idx]))
    
    list(span_rmse, plot)
  }


base_theme =
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"),
    strip.text = element_text(size = 12)
  )


get_go_similarity_matrix <-
  function(go_id,
           ont = c("CC", "BP", "MF"),
           db = 'org.Hs.eg.db',
           measure = c("Wang", "Resnik", "Rel", "Jiang", "Lin", "TCSS"),
           remove_orphan_terms = FALSE,
           sim_cutoff = 0.7) {
    if (length(go_id) == 0) {
      return(data.frame(
        go_id1 = character(),
        go_id2 = character(),
        sim = numeric()
      ))
    }
    measure <- match.arg(measure)
    ont <- match.arg(ont)
    message("Calculating...")
    sim_matrix <-
      simplifyEnrichment::GO_similarity(
        go_id = go_id,
        ont = ont,
        db = db,
        measure = measure,
        remove_orphan_terms
      ) %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "go_id1") %>%
      tidyr::pivot_longer(cols = -go_id1,
                          names_to = "go_id2",
                          values_to = "sim") %>%
      dplyr::filter(go_id1 != go_id2) %>%
      dplyr::filter(sim > sim_cutoff)
    
    name <- apply(sim_matrix, 1, function(x) {
      paste(sort(x[1:2]), collapse = "_")
    })
    
    sim_matrix <-
      sim_matrix %>%
      dplyr::mutate(name = name) %>%
      dplyr::arrange(name) %>%
      dplyr::distinct(name, .keep_all = TRUE) %>%
      dplyr::select(-name) %>%
      dplyr::mutate(ont = ont)
    message("Done.")
    sim_matrix
  }


merge_go_network <-
  function(go_information,
           go_similarity_matrix,
           path = ".",
           layout = c("fr", "kk"),
           mark_hull = FALSE) {
    library(ggraph)
    library(igraph)
    library(openxlsx)
    library(ggforce)
    layout <- match.arg(layout)
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    go_similarity_matrix <-
      as.data.frame(go_similarity_matrix)
    
    go_information <-
      as.data.frame(go_information)
    
    edge_data <-
      rbind(go_similarity_matrix) %>%
      dplyr::rename(from = go_id1, to = go_id2)
    
    node_data <-
      rbind(go_information) %>%
      as.data.frame() %>%
      dplyr::select(ID, everything()) %>%
      dplyr::rename(node = ID)
    
    temp_graph <-
      tidygraph::tbl_graph(nodes = node_data,
                           edges = edge_data,
                           directed = FALSE) %>%
      dplyr::mutate(degree = tidygraph::centrality_degree())
    
    subnetwork <-
      igraph::cluster_edge_betweenness(graph = temp_graph, weights = abs(edge_attr(temp_graph, "sim")))
    cluster <-
      as.character(membership(subnetwork)) %>%
      purrr::map(function(x) {
        if (sum(x == as.character(membership(subnetwork))) == 1) {
          return("Other")
        } else{
          return(x)
        }
      }) %>%
      unlist()
    
    new_cluster <-
      purrr::map(cluster, function(x) {
        paste("Module", match(x, unique(cluster)[unique(cluster) != "Other"]))
      }) %>%
      unlist()
    
    new_cluster[new_cluster == "Module NA"] <- "Other"
    
    temp_graph <-
      temp_graph %>%
      tidygraph::mutate(module = factor(new_cluster, levels = stringr::str_sort(unique(new_cluster), numeric = TRUE)))
    
    ###clustered different GO terms
    result_go <-
      igraph::vertex_attr(temp_graph) %>%
      do.call(cbind, .) %>%
      as.data.frame() %>%
      dplyr::mutate(p.adjust = as.numeric(p.adjust)) %>%
      dplyr::arrange(ONTOLOGY, module, p.adjust)
    
    
    ###output networ data
    wb <- createWorkbook()
    modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roman")
    addWorksheet(wb, sheetName = "Node data")
    addWorksheet(wb, sheetName = "Edge data")
    
    freezePane(
      wb = wb,
      sheet = 1,
      firstRow = TRUE,
      firstCol = TRUE
    )
    
    freezePane(
      wb = wb,
      sheet = 2,
      firstRow = TRUE,
      firstCol = TRUE
    )
    
    writeDataTable(
      wb,
      sheet = 1,
      x = result_go,
      colNames = TRUE,
      rowNames = FALSE
    )
    
    writeDataTable(
      wb,
      sheet = 2,
      x = edge_data,
      colNames = TRUE,
      rowNames = FALSE
    )
    
    saveWorkbook(wb, file.path(path, "network_data.xlsx"), overwrite = TRUE)
    
    ###plot to show the clusters of GO terms
    cluster_label1 <-
      igraph::as_data_frame(temp_graph, what = "vertices") %>%
      # dplyr::filter(module != "Other") %>%
      dplyr::group_by(module) %>%
      dplyr::filter(p.adjust == min(p.adjust)) %>%
      dplyr::pull(Description)
    
    cluster_label2 <-
      igraph::as_data_frame(temp_graph, what = "vertices") %>%
      dplyr::filter(module == "Other") %>%
      pull(Description)
    
    cluster_label <-
      c(cluster_label1, cluster_label2)
    
    center_go_id <-
      result_go$node[result_go$Description %in% cluster_label]
    
    save(center_go_id, file = file.path(path, "center_go_id"))
    
    #####base network 1: network with all the Go terms
    plot1_base_network <-
      temp_graph %>%
      ggraph(layout = layout, circular = FALSE) +
      geom_edge_link(
        aes(width = sim),
        strength = 1,
        color = "black",
        alpha = 1,
        show.legend = TRUE
      ) +
      geom_node_point(
        aes(fill = module, size = -log(p.adjust, 10)),
        shape = 21,
        alpha = 1,
        show.legend = TRUE
      ) +
      guides(fill = guide_legend(ncol = 1)) +
      scale_edge_width_continuous(range = c(0.1, 2)) +
      scale_size_continuous(range = c(1, 7)) +
      ggraph::theme_graph() +
      theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.position = "right",
        legend.background = element_rect(fill = "transparent", color = NA)
      )
    
    
    #####base network 2: network with only the merged network
    plot2_base_network <-
      temp_graph %>%
      tidygraph::filter(module != "Other") %>%
      ggraph(layout = layout, circular = FALSE) +
      geom_edge_link(
        aes(width = sim),
        strength = 1,
        color = "black",
        alpha = 1,
        show.legend = TRUE
      ) +
      geom_node_point(
        aes(fill = module, size = -log(p.adjust, 10)),
        shape = 21,
        alpha = 1,
        show.legend = TRUE
      ) +
      guides(fill = guide_legend(ncol = 1)) +
      scale_edge_width_continuous(range = c(0.1, 2)) +
      scale_size_continuous(range = c(1, 7)) +
      ggraph::theme_graph() +
      theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.position = "right",
        legend.background = element_rect(fill = "transparent", color = NA)
      )
    
    
    #####base network 3: network with only the unmerged network
    plot3_base_network <-
      temp_graph %>%
      tidygraph::filter(module == "Other") %>%
      ggraph(layout = layout, circular = FALSE) +
      geom_edge_link(
        aes(width = sim),
        strength = 1,
        color = "black",
        alpha = 1,
        show.legend = TRUE
      ) +
      geom_node_point(
        aes(fill = module, size = -log(p.adjust, 10)),
        shape = 21,
        alpha = 1,
        show.legend = TRUE
      ) +
      guides(fill = guide_legend(ncol = 1)) +
      scale_edge_width_continuous(range = c(0.1, 2)) +
      scale_size_continuous(range = c(1, 7)) +
      ggraph::theme_graph() +
      theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.position = "right",
        legend.background = element_rect(fill = "transparent", color = NA)
      )
    
    ####network with all the Go terms
    plot1_base_network_label <-
      plot1_base_network +
      geom_node_text(aes(
        x = x,
        y = y,
        label = ifelse(Description %in% cluster_label, Description, NA)
      ),
      size = 3,
      repel = TRUE)
    
    plot2_base_network_label <-
      plot2_base_network +
      geom_node_text(aes(
        x = x,
        y = y,
        label = ifelse(Description %in% cluster_label, Description, NA)
      ),
      size = 3,
      repel = TRUE)
    
    plot3_base_network_label <-
      plot2_base_network +
      geom_node_text(aes(
        x = x,
        y = y,
        label = ifelse(Description %in% cluster_label, Description, NA)
      ),
      size = 3,
      repel = TRUE)
    
    extrafont::loadfonts()
    
    ggplot2::ggsave(
      plot1_base_network,
      filename =
        file.path(path, "plot1_base_network.pdf"),
      width = 9,
      height = 7
    )
    
    ggplot2::ggsave(
      plot2_base_network,
      filename =
        file.path(path, "plot2_base_network.pdf"),
      width = 9,
      height = 7
    )
    
    ggplot2::ggsave(
      plot3_base_network,
      filename =
        file.path(path, "plot3_base_network.pdf"),
      width = 9,
      height = 7
    )
    
    ggplot2::ggsave(
      plot1_base_network_label,
      filename =
        file.path(path, "plot1_base_network_label.pdf"),
      width = 9,
      height = 7
    )
    
    ggplot2::ggsave(
      plot2_base_network_label,
      filename =
        file.path(path, "plot2_base_network_label.pdf"),
      width = 9,
      height = 7
    )
    
    ggplot2::ggsave(
      plot3_base_network_label,
      filename =
        file.path(path, "plot3_base_network_label.pdf"),
      width = 9,
      height = 7
    )
    
  }





