annotate_isotopes4fmsea <-
  function(object_list,
           threads = 3,
           rt_tol = 5,
           mz_tol = 15,
           intensity_tol = 0.5,
           max_isotope = 2) {
    ###check object_list
    if (missing(object_list)) {
      stop("object_list is missing.")
    }
    
    message("Checking object_list...")
    check4fmsea(object_list = object_list)
    message("Done.")
    
    message("Loading database...")
    load(
      here::here(
        "3-data_analysis/monkey_fMSEA_analysis_rplc/2_metabolite_annotation/database/hmdb_ms1.rda"
      )
    )
    message("Done.")
    
    message("Isotope annotation...")
    
    object_list <-
      seq_along(object_list) %>%
      purrr::map(function(i) {
        temp_data <-
          object_list[[i]]
        column <-
          slot(temp_data, "variable_info")$column %>%
          unique()
        polarity <-
          slot(temp_data, "variable_info")$polarity %>%
          unique()
        message("-----")
        message(column, " ", polarity)
        message("-----")
        
        ###Remove some redundancy
        message("Remove some redundancy...")
        pb <-
          progress::progress_bar$new(total = length(unique(
            temp_data@variable_info$variable_id
          )))
        
        slot(object = temp_data, name = "annotation_table") <-
          slot(object = temp_data, name = "annotation_table") %>%
          plyr::dlply(.variables = plyr::.(variable_id)) %>%
          purrr::map(function(x) {
            pb$tick()
            if (any(x$Level == 1)) {
              x <-
                x %>%
                dplyr::filter(Level == 1) %>%
                dplyr::group_by(HMDB.ID) %>%
                dplyr::filter(SS == max(SS)) %>%
                head(1) %>%
                dplyr::ungroup()
            }
            
            # if (any(x$Level == 2)) {
            #   x <-
            #     x %>%
            #     dplyr::filter(Level == 2) %>%
            #     dplyr::group_by(HMDB.ID) %>%
            #     dplyr::filter(SS == max(SS)) %>%
            #     head(1) %>%
            #     dplyr::ungroup()
            # }
            return(x)
          }) %>%
          dplyr::bind_rows() %>%
          as.data.frame()
        
        annotation_table <-
          temp_data@annotation_table %>%
          dplyr::left_join(hmdb_ms1@spectra.info[, c("HMDB.ID", "Formula")],
                           by = "HMDB.ID")
        
        variable_info <-
          temp_data@variable_info
        
        annotation_table <-
          annotation_table %>%
          dplyr::left_join(variable_info[, c("variable_id", "mz", "rt", "mean_intensity")],
                           by = "variable_id")
        
        annotation_table <-
          annotation_table %>%
          dplyr::filter(!is.na(Formula)) %>%
          dplyr::filter(nchar(Formula) > 4) %>%
          dplyr::filter(!stringr::str_detect(Formula, "Cu|Mg|Al|Si|Ca"))
        
        new_annotation_table <-
          annotation_table %>%
          dplyr::distinct(variable_id, Formula, Adduct, .keep_all = TRUE)
        
        # pb <-
        #   progress::progress_bar$new(total = nrow(new_annotation_table))
        
        future::plan(future::multisession, workers = threads)
        
        message("")
        message("Isotope annotation for each compound...")
        
        isotope_table <-
          seq_len(nrow(new_annotation_table)) %>%
          furrr::future_map(function(j) {
            # pb$tick()
            temp <-
              annotate_isotope(
                formula = new_annotation_table$Formula[j],
                adduct = new_annotation_table$Adduct[j],
                feature_mz = new_annotation_table$mz[j],
                feature_rt = new_annotation_table$rt[j],
                feature_intensity = new_annotation_table$mean_intensity[j],
                variable_info = variable_info,
                rt_tol = rt_tol,
                mz_tol = mz_tol,
                intensity_tol = intensity_tol,
                max_isotope = max_isotope
              )
            if (!is.null(temp)) {
              temp <-
                data.frame(
                  original_variable_id = new_annotation_table$variable_id[j],
                  Adduct = new_annotation_table$Adduct[j],
                  Formula = new_annotation_table$Formula[j],
                  temp
                )
            }
            return(temp)
            
          }, .progress = TRUE)
        
        isotope_table <-
          isotope_table %>%
          data.table::rbindlist()
        
        if (nrow(isotope_table) == 0) {
          isotope_table <- NULL
        } else{
          isotope_table <-
            annotation_table[, c(
              "variable_id",
              "Formula",
              "Adduct",
              "Compound.name",
              "CAS.ID",
              "HMDB.ID",
              "KEGG.ID",
              "Lab.ID",
              "Database"
            )] %>%
            dplyr::left_join(
              isotope_table,
              by = c("variable_id" = "original_variable_id",
                     "Formula", "Adduct")
            ) %>%
            dplyr::select(-c(feature_mz, feature_rt, variable_id)) %>%
            dplyr::rename(RT.error = rt.error,
                          variable_id = variable_id.y) %>%
            dplyr::filter(!is.na(isotopes)) %>%
            dplyr::mutate(
              ms2_files_id = NA,
              ms2_spectrum_id = NA,
              mz.match.score = NA,
              RT.match.score = NA,
              CE = NA,
              SS = NA,
              Total.score = NA,
              Level = 3
            ) %>%
            dplyr::select(
              variable_id,
              ms2_files_id,
              ms2_spectrum_id,
              Compound.name,
              isotopes,
              CAS.ID,
              HMDB.ID,
              KEGG.ID,
              Lab.ID,
              Adduct,
              mz.error,
              mz.match.score,
              RT.error,
              RT.match.score,
              CE,
              SS,
              Total.score,
              Database,
              Level,
              Formula,
              mz,
              rt
            )
        }
        
        annotation_table <-
          annotation_table %>%
          dplyr::mutate(isotopes = "[M]") %>%
          dplyr::select(
            variable_id,
            ms2_files_id,
            ms2_spectrum_id,
            Compound.name,
            isotopes,
            CAS.ID,
            HMDB.ID,
            KEGG.ID,
            Lab.ID,
            Adduct,
            mz.error,
            mz.match.score,
            RT.error,
            RT.match.score,
            CE,
            SS,
            Total.score,
            Database,
            Level,
            Formula,
            mz,
            rt
          )
        
        # dim(annotation_table)
        # dim(isotope_table)
        
        annotation_table <-
          rbind(annotation_table,
                isotope_table)
        
        slot(object = temp_data, name = "annotation_table") <-
          annotation_table
        rm(
          list = c(
            "annotation_table",
            "isotope_table",
            "new_annotation_table",
            "variable_info"
          )
        )
        temp_data
      })
    
    return(object_list)
    
  }

















####KEGG database
#' @title annotate_isotope
#' @description Annotate isotopes for known formula compounds.
#' @author Xiaotao Shen
#' \email{shenxt@@stanford.edu}
#' @param formula formula
#' @param adduct adduct
#' @param charge charge
#' @param feature_mz feature_mz
#' @param feature_rt feature_rt
#' @param feature_intensity feature_intensity
#' @param variable_info variable_info
#' @param rt_tol rt_tol
#' @param max_isotope max_isotope
#' @export

annotate_isotope <-
  function(formula = "C6H6N2O",
           adduct = "(M+H)+",
           charge = 1,
           feature_mz = 123.0553,
           feature_rt = 76.367,
           feature_intensity = 159517274,
           ## peak information
           variable_info,
           ## other parameters
           rt_tol = 5,
           mz_tol = 15,
           intensity_tol = 0.5,
           ##relative
           max_isotope = 2) {
    adduct <-
      adduct %>%
      stringr::str_extract("\\(.+\\)") %>%
      stringr::str_replace("\\(", "") %>%
      stringr::str_replace("\\)", "")
    formula1 <-
      masstools::sum_formula(formula = formula,
                             adduct = adduct)
    ###should be fix latter
    if (is.na(formula1)) {
      formula1 <- formula
    }
    
    molecule <-
      tryCatch(
        Rdisop::getMolecule(formula = formula1,
                            maxisotopes = max_isotope + 4),
        error = function(e)
          NULL
      )
    
    if (is.null(molecule)) {
      return(NULL)
    }
    
    isotopes <-
      t(Rdisop::getIsotope(molecule = molecule)) %>%
      as.data.frame()
    
    colnames(isotopes) <- c("mz", "intensity")
    
    exact_mass <- Rdisop::getMass(molecule)
    
    isotope <-
      seq_len(nrow(isotopes)) - which(isotopes$mz == exact_mass)
    
    isotope[isotope > 0] <-
      paste0("+", isotope[isotope > 0])
    
    isotopes$isotope <-
      paste0("[M", isotope, "]")
    
    
    # rownames(isotopes) <-
    #   c("[M]", paste("[M", "+", c(1:(nrow(
    #     isotopes
    #   ) - 1)), "]", sep = ""))
    #
    # isotopes <-
    #   data.frame(isotopes,
    #              rownames(isotopes),
    #              stringsAsFactors = FALSE)
    #
    # colnames(isotopes) <- c("mz", "intensity", "isotope")
    
    isotopes$intensity <-
      isotopes$intensity * feature_intensity / max(isotopes$intensity)
    
    accurate_mz <- feature_mz
    
    isotopes <-
      isotopes %>%
      dplyr::filter(isotope != "[M0]") %>%
      head(max_isotope)
    
    ###rt filtering
    if (any(colnames(variable_info) == "rt.error")) {
      variable_info <-
        variable_info %>%
        dplyr::collect(-rt.error)
    }
    variable_info <-
      variable_info %>%
      dplyr::mutate(rt.error = abs(feature_rt - rt)) %>%
      dplyr::filter(rt.error <= rt_tol)
    
    if (nrow(variable_info) == 0) {
      return(NULL)
    }
    
    iso_info <-
      seq_len(nrow(isotopes)) %>%
      purrr::map(function(i) {
        peak_mz_error <-
          abs(isotopes$mz[i] - variable_info$mz) * 10 ^ 6 / ifelse(isotopes$mz[i] >= 400, isotopes$mz[i], 400)
        intensity_error <-
          abs(isotopes$intensity[i] - variable_info$mean_intensity) / isotopes$intensity[i]
        idx <-
          which(peak_mz_error <= mz_tol &
                  intensity_error <= intensity_tol)
        if (length(idx) == 0) {
          return(NULL)
        }
        
        if (length(idx) > 0) {
          idx <- idx[which.min(intensity_error[idx])]
          return(
            data.frame(
              variable_id = variable_info$variable_id[idx],
              mz.error = peak_mz_error[idx],
              isotope = isotopes$isotope[i]
            )
          )
        }
        
      }) %>%
      dplyr::bind_rows() %>%
      as.data.frame()
    
    if (nrow(iso_info) == 0) {
      return(NULL)
    }
    
    if (!"[M+1]" %in% iso_info$isotope) {
      return(NULL)
    }
    
    colnames(iso_info) <- c("variable_id",
                            "mz.error",
                            "isotopes")
    
    if (!"[M+2]" %in% iso_info$isotopes) {
      iso_info <- iso_info[1, , drop = FALSE]
    }
    
    iso_info <-
      iso_info %>%
      dplyr::left_join(variable_info[, c("variable_id", "mz", "rt", "rt.error")],
                       by = "variable_id") %>%
      dplyr::mutate(feature_mz = feature_mz, feature_rt = feature_rt) %>%
      dplyr::select(feature_mz,
                    feature_rt,
                    variable_id,
                    isotopes,
                    mz,
                    rt,
                    mz.error,
                    rt.error)
    # rm(list = c("peak_mz", "peak_rt", "peak.int", "cor"))
    # gc()
    return(iso_info)
  }
