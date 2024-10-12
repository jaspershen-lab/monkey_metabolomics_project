no_source()

setwd(r4projects::get_project_wd())

library(tidymass)

setwd("D:/shenxt/project/metabolomics_aging/hmp_data/metabolomics/HILIC/MS1/mzxml/POS")


massprocesser::process_data(path = ".", 
                            polarity = "positive", 
                            ppm = 15, 
                            peakwidth = c(10, 60),
                            threads = 10, 
                            output_tic = FALSE, 
                            output_bpc = FALSE, 
                            min_fraction = 0.5,
                            group_for_figure = "QC")


setwd("D:/shenxt/project/metabolomics_aging/hmp_data/metabolomics/HILIC/MS1/mzxml/NEG")


massprocesser::process_data(path = ".", 
                            polarity = "negative", 
                            ppm = 15, 
                            peakwidth = c(10, 60),
                            threads = 10, 
                            output_tic = FALSE, 
                            output_bpc = FALSE, 
                            min_fraction = 0.5,
                            group_for_figure = "QC")
