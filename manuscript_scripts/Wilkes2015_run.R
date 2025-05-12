source('./src/r_src_scripts.R')

pars <- list(
  f.path = './data/Wilkes2015/intensities_log2fcTransformedVector.tsv',
  output.path = './results/Wilkes2015/',
  output.id = 'MCF7G2_vs_Mock',
  species = 'Human',
  # all default values except for m = 9
  log_intensities = F,
  alpha = 0.9, 
  n = 15,
  beta = 0.4,
  gamma = 1.0
)

run_KINference(
  f.path = pars$f.path,
  output.path = pars$output.path,
  output.id = pars$output.id,
  species = pars$species,
  alpha = pars$alpha,
  n = pars$n,
  beta = pars$beta,
  gamma = pars$gamma
)
