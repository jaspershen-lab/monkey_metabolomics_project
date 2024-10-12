run_msea <-
  function(feature_table,
           metabolite_set,
           exponent = 1,
           permutation_num = 1000,
           min_size = 5,
           max_size = 1000,
           pvalue_cutoff = 0.2,
           p_adjust_method = "fdr",
           seed = FALSE,
           verbose = TRUE,
           ...) {
    ##filter pathways according to pathway size
    feature_table <-
      feature_table %>%
      dplyr::arrange(dplyr::desc(condition))
    
    slot(metabolite_set, "compound_list") <-
      slot(metabolite_set, "compound_list") %>%
      lapply(function(x) {
        if (nrow(x) == 0) {
          return(x)
        }
        x %>%
          dplyr::filter(HMDB.ID %in% feature_table$HMDB.ID)
      })
    
    set_length <-
      lapply(metabolite_set@compound_list, function(x) {
        nrow(x)
      }) %>%
      unlist()
    
    remain_idx <-
      which(set_length >= min_size & set_length <= max_size)
    
    if (length(remain_idx) == 0) {
      return(NULL)
    }
    
    metabolite_set <-
      metabolite_set[remain_idx]
    
    ##find where are the metabolites in the feature table
    match_idx <-
      slot(metabolite_set, name = "compound_list") %>%
      purrr::map(function(x) {
        lapply(x$HMDB.ID, function(y)
          which(y == feature_table$HMDB.ID))
      })
    
    purrr::map2(
      .x = match_idx[[1]],
      .y = seq_len(length(match_idx[[1]])),
      .f = function(x, y) {
        if (length(x) == 0) {
          return(NULL)
        }
        temp <- 
        data.frame(x,
                   y,
                   score = feature_table$score[x],
                   compound_class = feature_table$compound_class[x])
        compound_class_index <-
          data.frame(compound_class = unique(temp$compound_class)) %>% 
          dplyr::mutate(compound_class_index = as.character(seq_along(unique(temp$compound_class))))
        temp <-
          temp %>% 
          dplyr::left_join(compound_class_index, by = "compound_class")
        temp
      }
    ) %>%
      dplyr::bind_rows() %>%
      # dplyr::filter(y == 2) %>%
      ggplot(aes(y, x)) +
      geom_point(aes(size = score, color = compound_class_index),
                 show.legend = FALSE) +
      theme_bw() +
      labs(x = "Metabolites", y = "Metabolic features")
    
    ##for the metabolites in the feature table,
    ##all the peaks should be remained, other peaks only remain one metabolite,
    ##either one is fine
    idx <-
      sort(unique(unlist(match_idx)))
    
    ##for the metabolites that are in feature tables, all the matched
    ###features should be remained
    ###the isotope
    feature_table1 <-
      feature_table[idx,] %>%
      dplyr::filter(isotopes == "[M]")
    
    ##for each compound_class, if the peaks are in different condition_class,
    ##then it should be classed into
    ##different classes and then score again
    ##
    ##and if there are [M+H] or [M-H] in the compound_class, only remain them
    #     feature_table1 <-
    # feature_table1
    #     dplyr::arrange(compound_class) %>%
    #     dplyr::mutate(fc_class =
    #                     dplyr::case_when(fc >= 0 ~ "increase",
    #                                      TRUE ~ "decrease")) %>%
    #     plyr::dlply(.variables = .(compound_class)) %>%
    #     purrr::map(
    #       .f = function(x) {
    #         if (length(unique(x$fc_class)) != 1) {
    #           new_score <-
    #             x %>%
    #             plyr::dlply(.variables = .(fc_class)) %>%
    #             purrr::map(
    #               .f = function(y) {
    #                 score_peak_group(peak_group = y)
    #               }
    #             ) %>%
    #             unlist()
    #           x$score[x$fc_class == "decrease"] <- new_score[1]
    #           x$score[x$fc_class == "increase"] <- new_score[1]
    #           if (length(unique(new_score)) != 1) {
    #             x <-
    #               x %>%
    #               dplyr::filter(score != min(new_score))
    #           }
    #         }
    #
    #         if(any(x$Adduct == "(M+H)+") | any(x$Adduct == "(M-H)-")){
    #           x <-
    #             x %>%
    #             dplyr::filter(Adduct == "(M+H)+" | Adduct == "(M-H)-")
    #         }
    #       }
    #     ) %>%
    #     dplyr::bind_rows() %>%
    #     dplyr::select(-fc_class) %>%
    #     dplyr::arrange(dplyr::desc(fc)) %>%
    #     dplyr::filter(isotope == "[M]")
    #
    #
    # ##for each HMDB.ID, only remain the compound_class with the max score
    #     feature_table1 <-
    #       feature_table1 %>%
    #       plyr::dlply(.variables = .(HMDB.ID)) %>%
    #       purrr::map(function(x) {
    #         x %>%
    #           dplyr::filter(score == max(score))
    #       }) %>%
    #       dplyr::bind_rows() %>%
    #       dplyr::arrange(desc(fc))
    
    
    ##for the metabolites are not in the feature set
    feature_table2 <-
      feature_table[-idx,] %>%
      dplyr::filter(isotopes == "[M]") %>%
      dplyr::filter(!variable_id %in% feature_table1$variable_id) %>%
      dplyr::distinct(variable_id, .keep_all = TRUE)
    
    feature_table <-
      rbind(feature_table1,
            feature_table2) %>%
      dplyr::arrange(dplyr::desc(condition))
    
    if (verbose) {
      message("Calculating observed enrichment scores...")
    }
    
    # browser()
    ######this is a tranditional GSEA algorithm
    observed_info <-
      slot(object = metabolite_set,
           name = "compound_list") %>%
      purrr::map(function(x) {
        get_msea_score(
          metabolite_set_list = x$HMDB.ID,
          feature_table = feature_table,
          exponent = exponent
        )
      })
    
    observed_score <-
      observed_info %>%
      purrr::map(function(x) {
        x$enrichment_score
      }) %>%
      unlist()
    
    if (verbose) {
      message("Calculating permutation scores...")
    }
    
    if (seed) {
      seeds <- sample.int(permutation_num)
    }
    
    permutation_scores <-
      BiocParallel::bplapply(seq_len(permutation_num), function(i) {
        if (seed) {
          set.seed(seeds[i])
        }
        
        permute_msea_enrichment_score(
          feature_table = feature_table,
          metabolite_set = metabolite_set,
          exponent = exponent
        )
      })
    
    
    permutation_scores <-
      do.call(what = "cbind",
              args = permutation_scores)
    
    rownames(permutation_scores) <- names(metabolite_set)
    
    pos.m <- apply(permutation_scores, 1, function(x)
      mean(x[x >= 0]))
    
    neg.m <- apply(permutation_scores, 1, function(x)
      abs(mean(x[x < 0])))
    
    normalize_es <- function(enrichment_score, pos.m, neg.m) {
      s <- sign(enrichment_score)
      m <- numeric(length(enrichment_score))
      m[s == 1] <- pos.m[s == 1]
      m[s == -1] <- neg.m[s == -1]
      enrichment_score / m
    }
    
    NES <- normalize_es(observed_score, pos.m, neg.m)
    
    permutation_scores <-
      apply(permutation_scores,
            2,
            normalize_es,
            pos.m = pos.m,
            neg.m = neg.m)
    
    if (class(permutation_scores) == "numeric") {
      permutation_scores <- matrix(permutation_scores, nrow = 1)
      rownames(permutation_scores) <- names(metabolite_set)
    }
    
    if (verbose)
      message("calculating p values...")
    
    p_values <- sapply(seq_along(observed_score), function(i) {
      if (is.na(NES[i])) {
        NA
      } else if (NES[i] >= 0) {
        (sum(permutation_scores[i, ] >= NES[i]) + 1) / (sum(permutation_scores[i,] >= 0) + 1)
      } else {
        # NES[i] < 0
        (sum(permutation_scores[i, ] <= NES[i]) + 1) / (sum(permutation_scores[i,] < 0) +
                                                          1)
      }
      
    })
    
    p_adjust <- p.adjust(p_values, method = p_adjust_method)
    
    qvalues <- calculate_qvalue(pvals = p_values)
    
    gs.name <- names(metabolite_set)
    
    description <- stringr::str_split(gs.name, pattern = ";") %>%
      lapply(function(x) {
        x[2]
      }) %>%
      unlist()
    
    id <- stringr::str_split(gs.name, pattern = ";") %>%
      lapply(function(x) {
        x[1]
      }) %>%
      unlist()
    
    params <- list(
      pvalue_cutoff = pvalue_cutoff,
      permutation_num = permutation_num,
      p_adjust_method = p_adjust_method,
      exponent = exponent,
      min_size = min_size,
      max_size = max_size
    )
    
    
    res <- data.frame(
      ID = id,
      Description = description,
      setSize = sapply(metabolite_set, length),
      enrichmentScore = observed_score,
      NES = NES,
      pvalue = p_values,
      p_adjust = p_adjust,
      qvalues = qvalues,
      stringsAsFactors = FALSE
    )
    
    
    res <-
      res %>%
      dplyr::filter(!is.na(pvalue)) %>%
      dplyr::filter(pvalue <= pvalue_cutoff) %>%
      dplyr::filter(p_adjust <= pvalue_cutoff) %>%
      dplyr::arrange(pvalue)
    
    if (nrow(res) == 0) {
      message("no term enriched under specific pvalue_cutoff...")
      return(
        new(
          "metPathMSEA",
          result     = res,
          metabolite_set   = metabolite_set,
          feature_list   = feature_list,
          params     = params
        )
      )
    }
    
    row.names(res) <- res$ID
    observed_info <-
      observed_info[sapply(res$ID, function(x)
        grep(x, names(observed_info)))]
    
    if (verbose)
      message("leading edge analysis...")
    
    ledge <- get_leading_edge(observed_info)
    
    res$rank <- ledge$rank
    res$leading_edge <- ledge$leading_edge
    res$core_enrichment <-
      sapply(ledge$core_enrichment, paste0, collapse = '/')
    
    
    if (verbose)
      message("done...")
    
    feature_list <- feature_table$condition
    names(feature_list) <- feature_table$HMDB.ID
    new(
      "metPathMSEA",
      result     = res,
      metabolite_set   = metabolite_set,
      feature_list   = feature_list,
      permutation_scores = permutation_scores,
      params     = params
    )
  }

#########Note, lots of code are from Guangchuang Yu, credits should be his.


##' @name metPathMSEA-class
##' @aliases metPathMSEA-class
##'   show,metPathMSEA-method summary,metPathMSEA-method
##'
##' @docType class
##' @slot result MSEA anaysis result.
##' @slot metaboliteSets metaboliteSets
##' @slot feature_list order rank feature_list
##' @slot permScores permutation scores
##' @slot params parameters
##' @exportClass metPathMSEA
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
##' @keywords classes
setClass(
  "metPathMSEA",
  representation   = representation(
    result          = "data.frame",
    metabolite_set        = "list",
    feature_list        = "numeric",
    permutation_scores      = "matrix",
    params          = "list"
  )
)


###----------------------------------------------------------------------------
# metabolite_set <- kegg_hsa_compound_pathway
# feature_list <- feature_list

####filter_metabolite_set
#' @title filter_metabolite_set
#' @description filter_metabolite_set
#' @author Xiaotao Shen
#' \email{shenxt@@stanford.edu}
#' @param metabolite_set metabolite_set
#' @param feature_table feature_table
#' @param exponent exponent
#' @param min_size min_size
#' @param max_size max_size
#' @export

filter_metabolite_set <- function(metabolite_set,
                                  feature_table,
                                  min_size = 5,
                                  max_size = 1000) {
  len <-
    lapply(metabolite_set@compound_list, function(x) {
      length(intersect(rownames(x), feature_table$HMDB.ID))
    }) %>%
    unlist() %>%
    unname()
  
  remain_idx <-
    which(len > min_size & len < max_size)
  
  metabolite_set <- metabolite_set[remain_idx]
  return(metabolite_set)
}


####get_msea_score
#' @title get_msea_score
#' @description get_msea_score
#' @author Xiaotao Shen
#' \email{shenxt@@stanford.edu}
#' @param metabolite_set metabolite_set
#' @param feature_table feature_table
#' @param exponent exponent
#' @param fortify fortify
#' @export
get_msea_score <-
  function(feature_table,
           metabolite_set_list,
           exponent = 1,
           fortify = FALSE) {
    ###################################################################
    ##    feature_list                                               ##
    ##                                                               ##
    ## 1. Rank order the N metabolites in D to form L = { g_1, ... , g_N}  ##
    ##    according to the correlation, r(g_j)=r_j,                  ##
    ##    of their expression profiles with C.                       ##
    ##                                                               ##
    ###################################################################
    
    ###################################################################
    ##    exponent                                                   ##
    ##                                                               ##
    ## An exponent p to control the weight of the step.              ##
    ##   When p = 0, Enrichment Score ( ES(S) ) reduces to           ##
    ##   the standard Kolmogorov-Smirnov statistic.                  ##
    ##   When p = 1, we are weighting the metabolites in S                 ##
    ##   by their correlation with C normalized                      ##
    ##   by the sum of the correlations over all of the metabolites in S.  ##
    ##                                                               ##
    ###################################################################
    
    ## metabolites defined in metabolite_set should appear in feature_table.
    metabolite_set_list <-
      intersect(metabolite_set_list, feature_table$HMDB.ID)
    
    ##N the length of feature_list
    num_feature_list <- nrow(feature_table)
    # num_feature_list <-
    # length(unique(feature_table$HMDB.ID))
    
    ##num_metabolite_set the length of feature set
    num_metabolite_set <- length(metabolite_set_list)
    
    p_hits <- p_miss <- numeric(num_feature_list)
    hits <-
      feature_table$HMDB.ID %in% metabolite_set_list ## logical
    
    p_hits[hits] <-
      (abs(feature_table$condition[hits]) * (feature_table$score[hits] / 100)) ^ exponent
    
    # plot(abs(feature_table$condition[hits]),
    #      (abs(feature_table$condition[hits]) * (feature_table$score[hits] / 100)))
    #
    #   score <- feature_table$score[hits]
    #   compound_class <- feature_table$compound_class[hits]
    #
    #   temp_data <- data.frame(
    #     index = 1:sum(hits),
    #     p_hits = p_hits[hits],
    #     feature_table[hits, ],
    #     stringsAsFactors = FALSE
    #   )
    #
    # idx <- 17
    #   temp_data %>%
    #     ggplot(aes(index, p_hits)) +
    #     geom_point() +
    #     geom_point(
    #       aes(x = index, y = p_hits, color = compound_class),
    #       show.legend = FALSE,
    #       data = temp_data %>% filter(HMDB.ID == unique(HMDB.ID)[idx])
    #     ) +
    #     ggrepel::geom_label_repel(
    #       aes(
    #         x = index,
    #         y = p_hits,
    #         label = paste(compound_class, score, Adduct, isotope, sep = ":"),
    #         color = compound_class
    #       ),
    #       show.legend = FALSE,
    #       data = temp_data %>% filter(HMDB.ID == unique(HMDB.ID)[idx])
    #     ) +
    #     ggsci::scale_color_aaas() +
    #     theme_bw()
    
    # write.csv(feature_table[hits,], "test.csv", row.names = FALSE)
    
    NR <- sum(p_hits)
    p_hits <- cumsum(p_hits / NR)
    
    p_miss[!hits] <- 1 / (num_feature_list - num_metabolite_set)
    p_miss <- cumsum(p_miss)
    
    running_enrichment_score <- p_hits - p_miss
    
    # data.frame(index = 1:length(p_hits),
    #            p_hits, p_miss, running_enrichment_score,
    #            stringsAsFactors = FALSE) %>%
    #   tidyr::pivot_longer(cols = -index, names_to = "class", values_to = "value") %>%
    #   ggplot(aes(x = index, y = value, color = class)) +
    #   geom_point()
    
    ## enrichment_score is the maximum deviation from zero of p_hits-p_miss
    max_enrichment_score <- max(running_enrichment_score)
    min_enrichment_score <- min(running_enrichment_score)
    
    if (abs(max_enrichment_score) > abs(min_enrichment_score)) {
      enrichment_score <- max_enrichment_score
    } else {
      enrichment_score <- min_enrichment_score
    }
    
    df <- data.frame(
      x = seq_along(running_enrichment_score),
      running_enrichment_score = running_enrichment_score,
      position = as.integer(hits)
    )
    
    # df %>%
    #   ggplot(aes(x, running_enrichment_score)) +
    #   geom_point(aes(color = as.character(position)))
    
    if (fortify) {
      return(df)
    }
    
    df$metabolite = feature_table$HMDB.ID
    result <-
      list(enrichment_score = enrichment_score,
           running_enrichment_score = df)
    return(result)
  }


####permute_feature_table
#' @title permute_feature_table
#' @description permute_feature_table
#' @author Xiaotao Shen
#' \email{shenxt@@stanford.edu}
#' @param feature_table feature_table
#' @export

permute_feature_table <-
  function(feature_table) {
    perm.idx <- sample.int(nrow(feature_table))
    feature_table_permutation <- feature_table
    feature_table_permutation$HMDB.ID <-
      feature_table_permutation$HMDB.ID[perm.idx]
    return(feature_table_permutation)
  }

####permute_msea_enrichment_score
#' @title permute_msea_enrichment_score
#' @description permute_feature_table
#' @author Xiaotao Shen
#' \email{shenxt@@stanford.edu}
#' @param feature_table feature_table
#' @param metabolite_set metabolite_set
#' @param exponent exponent
#' @export

permute_msea_enrichment_score <-
  function(feature_table,
           metabolite_set,
           exponent = 1) {
    feature_table_permutation <-
      permute_feature_table(feature_table)
    
    result <-
      seq_len(length(metabolite_set@compound_list)) %>%
      purrr::map(function(i) {
        get_msea_score(
          feature_table = feature_table_permutation,
          metabolite_set_list = metabolite_set@compound_list[[i]]$HMDB.ID,
          exponent = exponent
        )$enrichment_score
      }) %>%
      unlist()
    return(result)
  }



####calculate_qvalue
#' @title calculate_qvalue
#' @description calculate_qvalue
#' @author Xiaotao Shen
#' \email{shenxt@@stanford.edu}
#' @param pvals pvals
#' @export
calculate_qvalue <- function(pvals) {
  if (length(pvals) == 0)
    return(numeric(0))
  
  qobj <- tryCatch(
    qvalue::qvalue(pvals, lambda = 0.05,
                   pi0.method = "bootstrap"),
    error = function(e)
      NULL
  )
  
  if (class(qobj) == "qvalue") {
    qvalues <- qobj$qvalues
  } else {
    qvalues <- NA
  }
  return(qvalues)
}


###----------------------------------------------------------------------------
####calculate_qvalue
#' @title get_leading_edge
#' @description get_leading_edge
#' @author Xiaotao Shen
#' \email{shenxt@@stanford.edu}
#' @param observed_info observed_info
#' @export
get_leading_edge <- function(observed_info) {
  core_enrichment <- lapply(observed_info, function(x) {
    running_enrichment_score <- x$running_enrichment_score
    running_enrichment_score <-
      running_enrichment_score[running_enrichment_score$position == 1, ]
    enrichment_score <- x$enrichment_score
    if (enrichment_score >= 0) {
      i <- which.max(running_enrichment_score$running_enrichment_score)
      leading_metabolite <- running_enrichment_score$metabolite[1:i]
    } else {
      i <- which.min(running_enrichment_score$running_enrichment_score)
      leading_metabolite <-
        running_enrichment_score$metabolite[-c(1:(i - 1))]
    }
    return(leading_metabolite)
  })
  
  rank <- sapply(observed_info, function(x) {
    running_enrichment_score <- x$running_enrichment_score
    enrichment_score <- x$enrichment_score
    if (enrichment_score >= 0) {
      rr <- which.max(running_enrichment_score$running_enrichment_score)
    } else {
      i <- which.min(running_enrichment_score$running_enrichment_score)
      rr <- nrow(running_enrichment_score) - i + 1
    }
    return(rr)
  })
  
  tags <- sapply(observed_info, function(x) {
    running_enrichment_score <- x$running_enrichment_score
    running_enrichment_score <-
      running_enrichment_score[running_enrichment_score$position == 1, ]
    enrichment_score <- x$enrichment_score
    if (enrichment_score >= 0) {
      i <- which.max(running_enrichment_score$running_enrichment_score)
      res <- i / nrow(running_enrichment_score)
    } else {
      i <- which.min(running_enrichment_score$running_enrichment_score)
      res <-
        (nrow(running_enrichment_score) - i + 1) / nrow(running_enrichment_score)
    }
    return(res)
  })
  
  ll <- sapply(observed_info, function(x) {
    running_enrichment_score <- x$running_enrichment_score
    enrichment_score <- x$enrichment_score
    if (enrichment_score >= 0) {
      i <- which.max(running_enrichment_score$running_enrichment_score)
      res <- i / nrow(running_enrichment_score)
    } else {
      i <- which.min(running_enrichment_score$running_enrichment_score)
      res <-
        (nrow(running_enrichment_score) - i + 1) / nrow(running_enrichment_score)
    }
    return(res)
  })
  
  N <- nrow(observed_info[[1]]$running_enrichment_score)
  setSize <-
    sapply(observed_info, function(x)
      sum(x$running_enrichment_score$position))
  signal <- tags * (1 - ll) * (N / (N - setSize))
  
  tags <- paste0(round(tags * 100), "%")
  ll <- paste0(round(ll * 100), "%")
  signal <- paste0(round(signal * 100), "%")
  leading_edge <-
    paste0('tags=', tags, ", list=", ll, ", signal=", signal)
  
  res <- list(
    rank = rank,
    tags = tags,
    list = ll,
    signal = signal,
    leading_edge = leading_edge,
    core_enrichment = core_enrichment
  )
  return(res)
}



##' show method for \code{metPathMSEA} instance
##'
##' @name show
##' @docType methods
##' @rdname show-methods
##' @title show method
##' @return message
##' @importFrom methods show
##' @exportMethod show
##' @usage show(object)
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
setMethod("show", signature(object = "metPathMSEA"),
          function (object) {
            params <- object@params
            cat("#\n# Gene Set Enrichment Analysis\n#\n")
            cat("#...@feature_list", "\t")
            str(object@feature_list)
            cat("#...nPerm", "\t", params$permutation_num, "\n")
            cat(
              "#...pvalues adjusted by",
              paste0("'", params$p_adjust_method, "'"),
              paste0("with cutoff <", params$pvalueCutoff),
              "\n"
            )
            cat(paste0("#...", nrow(object@result)), "enriched terms found\n")
            str(object@result)
          })



##' visualize analyzing result of MSEA
##'
##' plotting function for mseaResult
##' @title msea_plot
##' @rdname mseaplot
##' @param x object of msea result
##' @param metaboliteSetID metaboliteSet ID
##' @param by one of "running_enrichment_score" or "position"
##' @param title plot title
##' @param ... additional parameters
##' @return ggplot2 object
##' @export
##' @examples
##' library(DOSE)
##' data(metaboliteList)
##' x <- gseDO(metaboliteList)
##' mseaplot(x, metaboliteSetID=1)
setGeneric(name = "msea_plot",
           function(x,
                    metabolite_set_idx = 1,
                    by = "all",
                    title = "",
                    color = 'black',
                    color.line = "#8DD3C7",
                    color.vline = "#FB8072",
                    ...) {
             standardGeneric("msea_plot")
           })


##' @rdname msea_plot
##' @exportMethod msea_plot
setMethod(f = "msea_plot", signature(x = "metPathMSEA"),
          function (x,
                    metabolite_set_idx = 1,
                    by = "all",
                    title = "",
                    color = 'black',
                    color.line = "#8DD3C7",
                    color.vline = "#FB8072",
                    ...) {
            msea_plot.metPathMSEA(
              x,
              metabolite_set_idx = metabolite_set_idx,
              by = by,
              title = title,
              color = color,
              color.line = color.line,
              color.vline = color.vline,
              ...
            )
          })

##' @rdname msea_plot
##' @param color color of line segments
##' @param color.line color of running enrichment score line
##' @param color.vline color of vertical line which indicating the maximum/minimal running enrichment score
##' @return ggplot2 object
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 geom_linerange
##' @importFrom ggplot2 geom_line
##' @importFrom ggplot2 geom_vline
##' @importFrom ggplot2 geom_hline
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 xlim
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 ggplotGrob
##' @importFrom ggplot2 geom_segment
##' @importFrom ggplot2 ggplot_gtable
##' @importFrom ggplot2 ggplot_build
##' @importFrom ggplot2 ggtitle
##' @importFrom ggplot2 element_text
##' @importFrom ggplot2 rel
##' @importFrom cowplot plot_grid
##' @author Guangchuang Yu

msea_plot.metPathMSEA <-
  function (x,
            metabolite_set_idx = 1,
            by = "all",
            title = "",
            color = 'black',
            color.line = "green",
            color.vline = "#FA5860",
            ...) {
    if (is.null(x)) {
      return(NULL)
    }
    by <-
      match.arg(by, c("running_enrichment_score", "preranked", "all"))
    gs_data <- get_gs_info(x, metabolite_set_idx)
    
    p <- ggplot(gs_data, aes_(x = ~ x)) +
      theme_dose() +
      xlab("Position in the ranked list of features")
    
    if (by == "running_enrichment_score" || by == "all") {
      p.res <-
        p + geom_linerange(aes_(ymin =  ~ ymin, ymax =  ~ ymax), color = color)
      p.res <-
        p.res + geom_line(aes_(y = ~ running_enrichment_score),
                          color = color.line,
                          size =
                            1)
      enrichmentScore <-
        x@result[metabolite_set_idx, "enrichmentScore"]
      
      es.df <-
        data.frame(es = which.min(abs(
          p$data$running_enrichment_score - enrichmentScore
        )))
      
      p.res <-
        p.res + geom_vline(
          data = es.df,
          aes_(xintercept = ~ es),
          colour = color.vline,
          linetype = "dashed"
        )
      p.res <- p.res + ylab("Running enrichment score")
      p.res <- p.res + geom_hline(yintercept = 0)
    }
    
    if (by == "preranked" || by == "all") {
      df2 <- data.frame(x = which(p$data$position == 1))
      df2$y <- p$data$feature_list[df2$x]
      p.pos <-
        p + geom_segment(data = df2,
                         aes_(
                           x =  ~ x,
                           xend =  ~ x,
                           y =  ~ y,
                           yend = 0
                         ),
                         color = color)
      p.pos <-
        p.pos + ylab("Ranked list metric") + xlim(0, length(p$data$feature_list))
    }
    if (by == "running_enrichment_score")
      return(p.res + ggtitle(title))
    if (by == "preranked")
      return(p.pos + ggtitle(title))
    
    p.pos <-
      p.pos + xlab(NULL) + theme(axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank())
    p.pos <- p.pos + ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5, size = rel(2)))
    cowplot::plot_grid(p.pos, p.res, ncol = 1, align = "v")
  }



get_gs_info <- function(object,
                        metabolite_set_idx = 1) {
  feature_list <- object@feature_list
  # if (is.numeric(metabolite_set_idx))
  #   metabolite_set_idx <- object@result[metabolite_set_idx, "ID"]
  metabolite_set <- object@metabolite_set[[metabolite_set_idx]]
  exponent <- object@params[["exponent"]]
  df <-
    get_msea_score(feature_list, metabolite_set, exponent, fortify = TRUE)
  df$ymin = 0
  df$ymax = 0
  pos <- df$position == 1
  h <- diff(range(df$running_enrichment_score)) / 20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$feature_list <- feature_list
  return(df)
}


##' ggplot theme of DOSE
##'
##' @title theme_dose
##' @param font.size font size
##' @return ggplot theme
##' @importFrom ggplot2 theme_bw
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 element_text
##' @importFrom ggplot2 margin
##' @examples
##' library(ggplot2)
##' qplot(1:10) + theme_dose()
##' @export
theme_dose <- function(font.size = 13) {
  theme_bw() +
    theme(
      axis.text.x = element_text(
        colour = "black",
        size = 12,
        vjust = 1
      ),
      axis.text.y = element_text(
        colour = "black",
        size = 12,
        hjust = 1
      ),
      axis.title = element_text(
        margin = margin(10, 5, 0, 0),
        color = "black",
        size = 13
      ),
      axis.title.y = element_text(angle = 90)
    )
}
