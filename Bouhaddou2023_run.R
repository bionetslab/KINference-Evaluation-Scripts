source('./src/r_src_scripts.R')

pars <- list(
  intensity.path = './data/Bouhaddou2023/intensity_measurements.tsv',
  output.path = './results/Bouhaddou2023/',
  disease_condition = 'VIC_10h',
  control_condition = 'Mock_10h'
)

run_pipeline(
  intensity.path = pars$intensity.path,
  output.path = pars$output.path,
  disease_condition = pars$disease_condition,
  control_condition = pars$control_condition,
  compute_CORR_filter = F,
  m = 9,
  gamma = 3.5
)
