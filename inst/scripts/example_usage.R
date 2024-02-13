# Ensure necessary libraries are installed and loaded
if (!requireNamespace("SEMdetectR", quietly = TRUE)) {
  devtools::install_github("HigginsChenLab/SEMdetectR") # Replace with devtools::install_github("username/YourPackageName") if hosted on GitHub
}

library(SEMdetectR)


# Example call to your main function
DNA_methylation <- read.csv("path/to/your/data.csv") # Example data loading must be beta values; rows = subjects; columns = probes (with Illumina 450k names if using RF models)
#mean imputation

DNA_methylation <- tDNAmFHS[,290000:300000]
SEM_results <- detectSEM(DNA_methylation, num_cores = 4, rf = TRUE, cluster = TRUE)
betas <- tDNAmFHS[,290000:300000]

SEM_results1 <- CountSEM_rf(DNA_methylation, model_features = model_features, cell_counts = cell_counts)
SEM_results1 <- as.data.frame(SEM_results1)


summary(lm(SEM_results1$hypoSEM ~ SEM_results$hypoSEM))
summary(lm(SEM_results1$hyperSEM ~ SEM_results$hyperSEM))
summary(lm(SEM_results1$hypoSEM ~ SEM_results$all_probes.hypoSEM))
summary(lm(SEM_results1$hyperSEM ~ SEM_results$all_probes.hyperSEM))
