do_fmsea <-
  function(object_list,
           databases_list,
           use_default_database = TRUE,
           rt.match.tol = 30,
           threads = 3,
           remove_fragment_intensity_cutoff = 0.01,
           isotope_match_rt_tol = 10,
           isotope_match_mz_tol = 15,
           isotope_match_intensity_tol = 0.5,
           max_isotope = 2,
           grouping_rt_tol = 10,
           score_cutoff = 20,
           candidate.num = 10000,
           path = ".") {
    ###check data
    if (missing(object_list)) {
      stop("object_list is missing.")
    }
    
    message("Checking object_list...")
    check4fmsea(object_list = object_list)
    message("Done.")
    
    ####Metabolite annotation
    message()
    message("Annotating metabolites...")
    dir.create(path, showWarnings = FALSE, recursive = TRUE)
    metabolite_annotation_path <-
      file.path(path, "metabolite_annotation")
    dir.create(metabolite_annotation_path)

    if (any(dir(metabolite_annotation_path) == "object_list")) {
      message("Downloading previous dataset.")
      load(file.path(metabolite_annotation_path, "object_list"))
    } else{
      object_list <-
        annotate_metabolite4fmsea(
          object_list = object_list,
          use_default_database = use_default_database,
          databases_list = databases_list,
          rt.match.tol = rt.match.tol,
          threads = threads,
          remove_fragment_intensity_cutoff = remove_fragment_intensity_cutoff,
          path = path,
          candidate.num = candidate.num
        )
      
      save(object_list,
           file = file.path(metabolite_annotation_path, "object_list"))
    }
    
    ####Isotope annotation
    message()
    message("Isotope annotation...")
    dir.create(path, showWarnings = FALSE, recursive = TRUE)
    isotope_annotation_path <-
      file.path(path, "isotope_annotation")
    dir.create(isotope_annotation_path,
               showWarnings = FALSE,
               recursive = TRUE)
    
    if (any(dir(isotope_annotation_path) == "object_list")) {
      message("Downloading previous dataset.")
      load(file.path(isotope_annotation_path, "object_list"))
    } else{
      object_list <-
        annotate_isotopes4fmsea(
          object_list = object_list,
          threads = threads,
          rt_tol = isotope_match_rt_tol,
          mz_tol = isotope_match_mz_tol,
          intensity_tol = isotope_match_intensity_tol,
          max_isotope = max_isotope
        )
      
      save(object_list,
           file = file.path(isotope_annotation_path, "object_list"))
    }
    
    ####Grouping features
    message()
    message("Grouping features...")
    feature_grouping_path <-
      file.path(path, "feature_grouping")
    dir.create(feature_grouping_path,
               showWarnings = FALSE,
               recursive = TRUE)
    
    if (any(dir(feature_grouping_path) == "annotation_table")) {
      message("Downloading previous dataset.")
      load(file.path(feature_grouping_path, "annotation_table"))
    } else{
      annotation_table <-
        group_features4fmsea(
          object_list = object_list,
          rt.match.tol = grouping_rt_tol,
          threads = threads
        )
      
      save(annotation_table,
           file = file.path(feature_grouping_path, "annotation_table"))
    }
    
    ####Removing redundancy
    message()
    message("Removing redundancy...")
    redundancy_removing_path <-
      file.path(path, "redundancy_removing")
    dir.create(redundancy_removing_path,
               showWarnings = FALSE,
               recursive = TRUE)
    
    if (any(dir(redundancy_removing_path) == "annotation_table")) {
      message("Downloading previous dataset.")
      load(file.path(redundancy_removing_path, "annotation_table"))
    } else{
      annotation_table <-
        remove_redundancy4fmsea(
          annotation_table = annotation_table,
          score_cutoff = score_cutoff,
          threads = threads
        )
      
      save(annotation_table,
           file = file.path(redundancy_removing_path, "annotation_table"))
    }
    
    message()
    message("All done.")
    
  }