# object_list <- list(object_rplc_pos,
#                     object_rplc_neg)

annotate_metabolite4fmsea <-
  function(object_list,
           databases_list,
           use_default_database = TRUE,
           same_lc_with_inhouse_library = TRUE,
           rt.match.tol = 30,
           threads = 3,
           candidate.num = 10000,
           remove_fragment_intensity_cutoff = 0.01,
           path = ".") {
    ###check object_list
    if (missing(object_list)) {
      stop("object_list is missing.")
    }
    message("Checking object_list...")
    check4fmsea(object_list = object_list)
    message("Done.")
    
    if (missing(databases_list)) {
      use_default_database <- TRUE
    } else{
      use_default_database <- FALSE
    }
    
    if (use_default_database) {
      #####load database
      message("Downloading databases...")
      load(
        here::here(
          "3-data_analysis/monkey_fMSEA_analysis_rplc/2_metabolite_annotation/database/hmdb_ms2.rda"
        )
      )
      load(
        here::here(
          "3-data_analysis/monkey_fMSEA_analysis_rplc/2_metabolite_annotation/database/massbank_ms2.rda"
        )
      )
      load(
        here::here(
          "3-data_analysis/monkey_fMSEA_analysis_rplc/2_metabolite_annotation/database/metlin_ms2.rda"
        )
      )
      load(
        here::here(
          "3-data_analysis/monkey_fMSEA_analysis_rplc/2_metabolite_annotation/database/mona_ms2.rda"
        )
      )
      load(
        here::here(
          "3-data_analysis/monkey_fMSEA_analysis_rplc/2_metabolite_annotation/database/mpsnyder_hilic_ms2.rda"
        )
      )
      load(
        here::here(
          "3-data_analysis/monkey_fMSEA_analysis_rplc/2_metabolite_annotation/database/mpsnyder_rplc_ms2.rda"
        )
      )
      load(
        here::here(
          "3-data_analysis/monkey_fMSEA_analysis_rplc/2_metabolite_annotation/database/nist_ms2.rda"
        )
      )
      load(
        here::here(
          "3-data_analysis/monkey_fMSEA_analysis_rplc/2_metabolite_annotation/database/hmdb_ms1.rda"
        )
      )
      
      ###only remain the metabolites have HMDB.ID, CAS.ID or KEGG.ID
      hmdb_ms2 <-
        hmdb_ms2 %>%
        dplyr::filter(!is.na(HMDB.ID))
      
      massbank_ms2 <-
        massbank_ms2 %>%
        dplyr::filter(!is.na(HMDB.ID))
      
      metlin_ms2 <-
        metlin_ms2 %>%
        dplyr::filter(!is.na(HMDB.ID))
      
      mpsnyder_hilic_ms2 <-
        mpsnyder_hilic_ms2 %>%
        dplyr::filter(!is.na(HMDB.ID))
      
      mpsnyder_rplc_ms2 <-
        mpsnyder_rplc_ms2 %>%
        dplyr::filter(!is.na(HMDB.ID))
      
      mona_ms2 <-
        mona_ms2 %>%
        dplyr::filter(!is.na(HMDB.ID))
      
      nist_ms2 <-
        nist_ms2 %>%
        dplyr::filter(!is.na(HMDB.ID))
      
      hmdb_ms1 <-
        hmdb_ms1 %>%
        dplyr::filter(!is.na(HMDB.ID)) %>%
        dplyr::filter(!is.na(Formula)) %>%
        dplyr::filter(nchar(Formula) > 4) %>%
        dplyr::filter(!stringr::str_detect(Formula, "Cu|Mg|Al|Si|Ca"))
      message("Done.")
      
      #####metabolite annotation
      ###mpsnyder_rplc
      object_list <-
        seq_along(object_list) %>%
        purrr::map(function(i) {
          x <- object_list[[i]]
          column <-
            slot(x, "variable_info")$column %>%
            unique()
          polarity <-
            slot(x, "variable_info")$polarity %>%
            unique()
          
          message(column, " ", polarity)
          
          temp_output_path <-
            file.path(path,
                      "metabolite_annotation",
                      paste(column, polarity, sep = "_"))
          
          dir.create(temp_output_path,
                     showWarnings = FALSE,
                     recursive = TRUE)
          
          if ("x" %in% dir(temp_output_path)) {
            load(file.path(temp_output_path, "x"))
            message("Use previouse data.")
            return(x)
          }
          
          ###In house library
          message("In house library...")
          
          if (column == "RPLC") {
            x <-
              metid::annotate_metabolites_mass_dataset(
                object = x,
                polarity = polarity,
                rt.match.tol = ifelse(same_lc_with_inhouse_library, rt.match.tol, 100000),
                column = ifelse(column == "RPLC", "rp", "hilic"),
                candidate.num = candidate.num,
                database = mpsnyder_rplc_ms2,
                threads = threads,
                remove_fragment_intensity_cutoff = remove_fragment_intensity_cutoff
              )
          } else{
            x <-
              metid::annotate_metabolites_mass_dataset(
                object = x,
                polarity = polarity,
                rt.match.tol = ifelse(same_lc_with_inhouse_library, rt.match.tol, 100000),
                column = ifelse(column == "RPLC", "rp", "hilic"),
                candidate.num = candidate.num,
                database = mpsnyder_hilic_ms2,
                threads = threads,
                remove_fragment_intensity_cutoff = remove_fragment_intensity_cutoff
              )
          }
          
          ##MassBank
          message("massbank_ms2...")
          x <-
            metid::annotate_metabolites_mass_dataset(
              object = x,
              polarity = polarity,
              rt.match.tol = rt.match.tol,
              column = ifelse(column == "RPLC", "rp", "hilic"),
              candidate.num = candidate.num,
              database = massbank_ms2,
              threads = threads,
              remove_fragment_intensity_cutoff = remove_fragment_intensity_cutoff
            )
          
          ###METLIN
          message("metlin_ms2...")
          x <-
            metid::annotate_metabolites_mass_dataset(
              object = x,
              polarity = polarity,
              rt.match.tol = rt.match.tol,
              column = ifelse(column == "RPLC", "rp", "hilic"),
              candidate.num = candidate.num,
              database = metlin_ms2,
              threads = threads,
              remove_fragment_intensity_cutoff = remove_fragment_intensity_cutoff
            )
          
          ###MONA
          message("mona_ms2...")
          x <-
            metid::annotate_metabolites_mass_dataset(
              object = x,
              polarity = polarity,
              rt.match.tol = rt.match.tol,
              column = ifelse(column == "RPLC", "rp", "hilic"),
              candidate.num = candidate.num,
              database = mona_ms2,
              threads = threads,
              remove_fragment_intensity_cutoff = remove_fragment_intensity_cutoff
            )
          
          ###NIST
          message("nist_ms2...")
          x <-
            metid::annotate_metabolites_mass_dataset(
              object = x,
              polarity = polarity,
              rt.match.tol = rt.match.tol,
              column = ifelse(column == "RPLC", "rp", "hilic"),
              candidate.num = candidate.num,
              database = nist_ms2,
              threads = threads,
              remove_fragment_intensity_cutoff = remove_fragment_intensity_cutoff
            )
          
          ###HMDB MS2
          message("hmdb_ms2...")
          x <-
            metid::annotate_metabolites_mass_dataset(
              object = x,
              polarity = polarity,
              rt.match.tol = rt.match.tol,
              column = ifelse(column == "RPLC", "rp", "hilic"),
              candidate.num = candidate.num,
              database = hmdb_ms2,
              threads = threads,
              remove_fragment_intensity_cutoff = remove_fragment_intensity_cutoff
            )
          
          ###HMDB MS1
          message("hmdb_ms1...")
          x <-
            metid::annotate_metabolites_mass_dataset(
              object = x,
              polarity = polarity,
              rt.match.tol = 100000,
              column = ifelse(column == "RPLC", "rp", "hilic"),
              candidate.num = candidate.num,
              database = hmdb_ms1,
              threads = threads,
              remove_fragment_intensity_cutoff = remove_fragment_intensity_cutoff
            )
          save(x, file = file.path(temp_output_path, "x"))
          x
        })
    } else{
      databases_list <-
        databases_list %>%
        purrr::map(function(x) {
          x %>%
            dplyr::filter(!is.na(HMDB.ID))
        })
      
      object_list <-
        seq_along(object_list) %>%
        purrr::map(function(i) {
          x <- object_list[[i]]
          column <-
            slot(x, "variable_info")$column %>%
            unique()
          polarity <-
            slot(x, "variable_info")$polarity %>%
            unique()
          
          message(column, " ", polarity)
          
          temp_output_path <-
            file.path(path,
                      "metabolite_annotation",
                      paste(column, polarity, sep = "_"))
          
          dir.create(temp_output_path,
                     showWarnings = FALSE,
                     recursive = TRUE)
          
          if ("x" %in% dir(temp_output_path)) {
            load(file.path(temp_output_path, "x"))
            message("Use previouse data.")
            return(x)
          }
          
          ###In house library
          for (i in seq_along(databases_list)) {
            temp_database <-
              databases_list[[i]]
            message(
              paste(
                temp_database@database.info$Source,
                temp_database@database.info$Version,
                sep = "_"
              )
            )
            
            x <-
              metid::annotate_metabolites_mass_dataset(
                object = x,
                polarity = polarity,
                rt.match.tol = rt.match.tol,
                column = ifelse(column == "RPLC", "rp", "hilic"),
                candidate.num = candidate.num,
                database = temp_database,
                threads = threads,
                remove_fragment_intensity_cutoff = remove_fragment_intensity_cutoff
              )
          }
          
          save(x, file = file.path(temp_output_path, "x"))
          x
        })
      
      
    }
    return(object_list)
  }
