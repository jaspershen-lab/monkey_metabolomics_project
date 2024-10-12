####This should be run in workstation

library(tidymass)

###RPLC pos
setwd(
  "D:/shenxt/project/metabolomics_aging/monkey_data/metabolomics/R_Macaque/mzxml/RPLC/pos"
)

massprocesser::process_data(
  path = ".",
  polarity = 'positive',
  ppm = 15,
  peakwidth = c(5, 30),
  snthresh = 10,
  prefilter = c(3, 500),
  fitgauss = FALSE,
  integrate = 2,
  mzdiff = 0.01,
  noise = 500,
  threads = 20,
  binSize = 0.025,
  bw = 5,
  output_tic = TRUE,
  output_bpc = TRUE,
  output_rt_correction_plot = TRUE,
  min_fraction = 0.5,
  fill_peaks = FALSE,
  group_for_figure = "QC"
)


###RPLC neg
setwd(
  "D:/shenxt/project/metabolomics_aging/monkey_data/metabolomics/R_Macaque/mzxml/RPLC/neg"
)

massprocesser::process_data(
  path = ".",
  polarity = 'negative',
  ppm = 15,
  peakwidth = c(5, 30),
  snthresh = 10,
  prefilter = c(3, 500),
  fitgauss = FALSE,
  integrate = 2,
  mzdiff = 0.01,
  noise = 500,
  threads = 20,
  binSize = 0.025,
  bw = 5,
  output_tic = TRUE,
  output_bpc = TRUE,
  output_rt_correction_plot = TRUE,
  min_fraction = 0.5,
  fill_peaks = FALSE,
  group_for_figure = "QC"
)






###HILIC pos
setwd(
  "D:/shenxt/project/metabolomics_aging/monkey_data/metabolomics/R_Macaque/mzxml/HILIC/pos"
)

massprocesser::process_data(
  path = ".",
  polarity = 'positive',
  ppm = 15,
  peakwidth = c(10, 60),
  snthresh = 10,
  prefilter = c(3, 500),
  fitgauss = FALSE,
  integrate = 2,
  mzdiff = 0.01,
  noise = 500,
  threads = 20,
  binSize = 0.025,
  bw = 5,
  output_tic = TRUE,
  output_bpc = TRUE,
  output_rt_correction_plot = TRUE,
  min_fraction = 0.5,
  fill_peaks = FALSE,
  group_for_figure = "QC"
)


###RPLC neg
setwd(
  "D:/shenxt/project/metabolomics_aging/monkey_data/metabolomics/R_Macaque/mzxml/HILIC/neg"
)

massprocesser::process_data(
  path = ".",
  polarity = 'negative',
  ppm = 15,
  peakwidth = c(10, 60),
  snthresh = 10,
  prefilter = c(3, 500),
  fitgauss = FALSE,
  integrate = 2,
  mzdiff = 0.01,
  noise = 500,
  threads = 20,
  binSize = 0.025,
  bw = 5,
  output_tic = FALSE,
  output_bpc = FALSE,
  output_rt_correction_plot = TRUE,
  min_fraction = 0.5,
  fill_peaks = FALSE,
  group_for_figure = "QC"
)
