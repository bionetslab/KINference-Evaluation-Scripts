library(Kinference)

tested_alpha <- c(0.85, 0.9, 0.95)
tested_n <- c(20, 15, 10)
tested_beta <- c(0.0, 0.2, 0.4, 0.6)
tested_gamma <- c(0.0, 0.5, 1.0, 1.5, 2.0, 2.5)
tested_delta <- c(0.7, 0.8, 0.9)

# Wilkes et al. 2015
f.path = '../data/Wilkes2015/intensities_log2fcTransformedVector.tsv'
run_Hyperparameter_evaluation(
   f.path=f.path, 
   alpha=tested_alpha, 
   n=tested_n, 
   beta=tested_beta, 
   gamma=tested_gamma, 
   output.id='Wilkes2015', 
   output.path='../results/Hyperparameters/Wilkes2015/'
)

# Bouhaddou et al. 2023
x0.path = '../data/Bouhaddou2023/intensities_Mock10h.tsv'
x1.path = '../data/Bouhaddou2023/intensities_VIC10h.tsv'
run_Hyperparameter_evaluation(
    x0.path=x0.path,
    x1.path=x1.path,
    paired.samples=T,
    alpha=tested_alpha, 
    n=tested_n, 
    apply.CORR = TRUE,
    beta=tested_beta, 
    gamma=tested_gamma, 
    delta=tested_delta,
    m=9,
    output.id='Bouhaddou2023', 
    output.path='../results/Hyperparameters/Bouhaddou2023/'
)
