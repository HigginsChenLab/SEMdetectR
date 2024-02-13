#' Random Forest Model used to filter out unreliable hypoSEMs
#'
#' This model is used internally by the package's function CountSEM_rf. It was trained to predict whether a hypoSEM detected by the IQR-based method would have been detected in a technical replicate.
#' This model was trained using the randomForest package with
#' parameters: mtry = 8, ntree = 245
#'
#' @name treeModelData_hypo.rda
#' @format An object of class `randomForest` as defined by the randomForest package.
#' @source The model was trained on the GSE55763 dataset of whole blood technical replicates. SEMs were detected with the IQR-based method in both replicates and annotated with Illumina manifest and inferred cell counts (Houseman).
#' @references Markov, Y., Levine, M., & Higgins-Chen, A. T. (2023). Stochastic Epigenetic Mutations: Reliable Detection and Associations with Cardiovascular Aging. bioRxiv. 10.1101/2023.12.12.571149
"treeModelData_hypo"
