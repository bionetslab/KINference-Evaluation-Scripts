source('./src/r_src_scripts.R')

pars <- list(
  x1.path = './data/Bouhaddou2023/intensities_VIC10h.tsv',
  x0.path = './data/Bouhaddou2023/intensities_Mock10h.tsv',
  output.path = './results/Bouhaddou2023/VIC10h_vs_Mock10h',
  output.id = 'VIC10h_vs_Mock10h',
  species = 'Human',
  paired_samples = F,
  log_intensities = T, 
  # all default values except for m = 9
  alpha = 0.9, 
  n = 15,
  beta = 0.4,
  gamma = 1.0,
  delta = 0.8,
  epsilon = 0.05,
  m = 9
)

run_KINference(
  x1.path = pars$x1.path,
  x0.path = pars$x0.path,
  output.path = pars$output.path,
  output.id = pars$output.id,
  species = pars$species,
  paired_samples = pars$paired_samples,
  log_intensities = pars$log_intensities, 
  alpha = pars$alpha,
  n = pars$n,
  beta = pars$beta,
  gamma = pars$gamma,
  delta = pars$delta,
  epsilon = pars$epsilon,
  m = pars$m
)
