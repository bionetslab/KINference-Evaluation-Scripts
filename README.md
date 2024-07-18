# KINference-Evaluation-Scripts

## Installation instruction

1) Create a conda environment (mamba is recommended):
```
mamba create -n 'KINference' -c conda-forge r-base python
mamba activate KINference
```
2) Install `pcst_fast`
```
pip install pcst_fast
```
3) Install required R packages (in an R terminal):
```
install.packages(c('tidyverse', 'data.table', 'BiocManager', 'devtools', 'argparse'))
devtools::install_github("evocellnet/funscoR")
BiocManager::install(c('OmnipathR', 'UniProt.ws'))
```

## Run the code for result generation:
Note: You have to run the code from the terminal and not inside an R terminal as the PCST calculation using the python implementatiion of pcst_fast will not work otherwise.
- Wilkes 2015 et al. (https://doi.org/10.1073/pnas.142334411)
```
Rscript Wilkes2015_run.R
```
- Bouhaddou 2023 et al. (https://doi.org/10.1016/j.cell.2023.08.026)
```
Rscript Bouhaddou2023_run.R
```