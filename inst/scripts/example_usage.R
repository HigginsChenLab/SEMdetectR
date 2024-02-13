devtools::install_github("yaromar/SEMdetectR")
library(SEMdetectR)



DNA_methylation <- read.csv("path/to/your/data.csv") # Example data loading must be beta values; rows = subjects; columns = probes (with Illumina 450k names if using RF models)

SEM_results <- detectSEM(DNA_methylation, num_cores = 4, rf = FALSE, cluster = FALSE)
