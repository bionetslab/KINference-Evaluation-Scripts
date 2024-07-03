pars <- list(
  intensity.path = './data/Wilkes2015/intensity_measurements.tsv',
  output.path = './results/Wilkes2015/',
  disease_condition = 'MCF7-G2'
)

source('./src/r_src_scripts.R')
run_pipeline_withGivenLog2FC(
  intensity.path = pars$intensity.path,
  output.path = pars$output.path,
  disease_condition = pars$disease_condition,
  gamma = 1.5
)