remove_redundancy4fmsea <-
  function(annotation_table,
           score_cutoff = 20,
           threads = 3) {
    ###remove the compound classes with score <= 20
    annotation_table <-
      annotation_table %>%
      dplyr::filter(score > score_cutoff)
    
    ###remove some compound classes based on adducts
    message("Removing some compound classess according to adducts...")
    
    compound_class <-
      annotation_table %>%
      plyr::dlply(.variables = plyr::.(compound_class))
    
    adduct_table <-
      annotation_table %>%
      dplyr::filter(Level != 3) %>%
      dplyr::filter(isotopes == "[M]") %>%
      select(Adduct, polarity) %>%
      dplyr::count(Adduct, polarity)
    
    #####remove some compound classes
    remain_index <-
      purrr::map(compound_class, function(x) {
        if (any(
          x$Adduct %in% c(
            adduct_table$Adduct,
            "(M+H)+",
            "(M+Na)+",
            "(M+K)+",
            "(M-H)-",
            "(M+Cl)-"
          )
        ) | nrow(x) > 1) {
          return(TRUE)
        } else{
          return(FALSE)
        }
      }) %>%
      unlist() %>%
      which()
    
    compound_class <-
      compound_class[remain_index]
    
    annotation_table <-
      compound_class %>%
      data.table::rbindlist()
    
    annotation_table <-
      remove_redundancy(annotation_table = annotation_table,
                        threads = threads)
    
    message("Done.")
    
    annotation_table
  }




##for one compound, if one compound class with score > 100, then remove the
##compound class with score <= 20
remove_redundancy <-
  function(annotation_table,
           threads = 3) {
    redundancy_diff = c(-1, -1)
    future::plan(future::multisession, workers = threads)
    message("Removing redundancy...")
    while (any(redundancy_diff < 0)) {
      # cat("i", " ")
      before_redundancy <-
        calculate_redundancy(annotation_table = annotation_table)
      message("Removing some compound classes...")
      annotation_table <-
        annotation_table %>%
        plyr::dlply(.variables = plyr::.(HMDB.ID)) %>%
        furrr::future_map(function(x) {
          ##if any compound class with Level 1 annotation,
          ##only remain this compound class
          if (any(x$score >= 500)) {
            x <-
              x %>%
              dplyr::filter(score >= 500)
            return(x)
          }
          # if (any(x$Level == 2) & any(x$score >= 200)) {
          #   x <-
          #     x %>%
          #     dplyr::filter(score >= 200)
          # }
          if (any(x$score >= 200)) {
            x <-
              x %>%
              dplyr::filter(score >= 100)
          }
          
          if (any(x$score >= 100)) {
            x <-
              x %>%
              dplyr::filter(score > 30)
          }
          
          return(x)
        }, .progress = TRUE) %>%
        dplyr::bind_rows() %>%
        as.data.frame()
      
      ##for one peak, if it has a annotation with score > 100, then remove other
      ##annotations with score <= 20
      message("")
      message("Removing some annotations...")
      annotation_table <-
        annotation_table %>%
        plyr::dlply(.variables = plyr::.(variable_id)) %>%
        furrr::future_map(function(x) {
          ##if any compound class with Level 1 annotation,
          ##only remain this compound class
          if (any(x$Level == 1)) {
            x <-
              x %>%
              dplyr::filter(Level == 1)
            return(x)
          }
          
          if (any(x$Level == 2)) {
            x <-
              x %>%
              dplyr::filter(Level == 2)
            return(x)
          }
          
          if (any(x$score >= 200)) {
            x <-
              x %>%
              dplyr::filter(score >= 100)
          }
          
          if (any(x$score >= 100)) {
            x <-
              x %>%
              dplyr::filter(score > 30)
          }
          return(x)
        }, .progress = TRUE) %>%
        dplyr::bind_rows() %>%
        as.data.frame()
      
      ###re-calculated confidence score for each compound class.
      message("")
      message("Rescoring...")
      annotation_table <-
        annotation_table %>%
        plyr::dlply(.variables = plyr::.(compound_class)) %>%
        furrr::future_map(function(x) {
          score <- score_compound_class(x)
          x$score <- score
          x
        }, .env_globals = TRUE) %>%
        dplyr::bind_rows()
      
      after_redundancy <-
        calculate_redundancy(annotation_table = annotation_table)
      
      redundancy_diff = after_redundancy - before_redundancy
    }
    
    ###remove some annotation based on variable_id and HMDB.ID
    temp <-
      data.table::data.table(annotation_table)[, .N, by = .(variable_id, HMDB.ID)]
    
    temp <-
      temp[which(temp$N > 1),]
    # browser()
    if (nrow(temp) > 0) {
      remove_table <-
        seq_len(nrow(temp)) %>%
        purrr::map(function(i) {
          temp_table <-
            annotation_table %>%
            dplyr::filter(variable_id == temp$variable_id[i],
                          HMDB.ID == temp$HMDB.ID[i]) %>%
            dplyr::arrange(dplyr::desc(RT.error))
          temp_table %>%
            head(nrow(temp_table) - 1)
        }) %>%
        dplyr::bind_rows()
      
      annotation_table <-
        annotation_table %>%
        dplyr::anti_join(remove_table, by = colnames(annotation_table))
    }
    
    return(annotation_table)
    
  }





#' @title score_peak_group
#' @description score_peak_group
#' @author Xiaotao Shen
#' \email{shenxt@@stanford.edu}
#' @param peak_group peak_group
#' @export
score_peak_group <- function(peak_group) {
  score <- 0
  
  ##This should be optimized in the future
  
  ###if positive adduct are +H, score add 50
  if (any(peak_group$Adduct == "(M+H)+")) {
    score <- score + 50
  }
  
  if (any(peak_group$Adduct == "(M+H)+" &
          peak_group$isotope != "[M]")) {
    score <- score + 20
  }
  
  if (any(peak_group$polarity == "positive" &
          peak_group$Adduct != "(M+H)+")) {
    score <- score + 20
  }
  
  if (any(
    peak_group$polarity == "positive" & peak_group$Adduct != "(M+H)+" &
    peak_group$isotope != "[M]"
  )) {
    score <- score + 10
  }
  
  if (any(peak_group$Adduct == "(M-H)-")) {
    score <- score + 50
  }
  
  if (any(peak_group$Adduct == "(M-H)-" &
          peak_group$isotope != "[M]")) {
    score <- score + 20
  }
  
  if (any(peak_group$polarity == "negative" &
          peak_group$Adduct != "(M-H)-")) {
    score <- score + 20
  }
  
  
  if (any(
    peak_group$polarity == "negative" & peak_group$Adduct != "(M-H)-" &
    peak_group$isotope != "[M]"
  )) {
    score <- score + 10
  }
  
  return(score)
}
