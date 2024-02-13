utils::globalVariables(c("probe_annotations", "treeModelData_hyper", "treeModelData_hypo"))
#' Calculate Statistics for Each probe in a DNA methylation dataframe
#'
#' Computes mean, quartiles, skewness, standard deviation, minimum, maximum, and kurtosis.
#' Can operate in parallel to speed up computations.
#'
#' @param betas A data frame or matrix of numerical data.
#' @param num_cores The number of cores to use for parallel computation.
#' @return A data frame with computed statistics for each column.
#'
#' @keywords internal
#'
#' @importFrom parallel makeCluster clusterEvalQ parLapply
#' @importFrom moments skewness kurtosis
#' @importFrom stats quantile sd
#'
#' @note This function is intended for internal use within the package.
GetStats<- function(betas, num_cores=1) {

  if(num_cores > 1){
    stats <- data.frame(
      mean = colMeans(betas)
    )

    # Set up parallel computing
    cl <- parallel::makeCluster(num_cores)
    # Ensure worker nodes have access to 'moments' package functions
    parallel::clusterEvalQ(cl, library(moments))

    # Define a function to compute statistics for a single column
    computeStats <- function(col){
      c(q1 = stats::quantile(col, probs = .25),
        q3 = stats::quantile(col, probs = .75),
        skewness = moments::skewness(col),
        stdev = stats::sd(col),
        min = min(col),
        max = max(col),
        kurtosis = moments::kurtosis(col))
    }

    # Apply the computeStats function in parallel across columns
    otherStats <- parallel::parLapply(cl, betas, computeStats)

    # Clean up parallel computing resources
    parallel::stopCluster(cl)
    rm(cl)
    gc()
    # Combine computed statistics into a single data frame
    otherStatsDf <- do.call(rbind, otherStats)
    stats <- cbind(stats, otherStatsDf)
  }
  else{
    # Compute statistics without parallel processing
    stats <- data.frame(mean=colMeans(betas),
                        q1=apply(betas,2,stats::quantile,probs=.25),
                        q3=apply(betas,2,stats::quantile,probs=.75),
                        skewness=apply(betas, 2, moments::skewness),
                        stdev=apply(betas,2,stats::sd),
                        min=apply(betas,2,min),
                        max=apply(betas, 2, max),
                        kurtosis=apply(betas, 2,moments::kurtosis)
    )
  }

  # Add interquartile range and range calculations
  stats$iqr <- stats$q3 - stats$q1
  stats$range <- stats$max - stats$min

  return(stats)
}


#' Cluster Columns Based on Their Statistics
#'
#' Performs k-means clustering on the scaled statistics of probes in DNA methylation dataframe,
#' aiming to group columns into 3 clusters based on their statistical properties.
#'
#' @param probe_stats A data frame of statistics for each probe, as returned by `GetStats`.
#' @return An object of class `kmeans` representing the clustering result.
#' @keywords internal
#' @importFrom stats kmeans
#' @note This function is intended for internal use within the package.
GetClusters <- function(probe_stats){
  clusters <- stats::kmeans(scale(probe_stats), centers = 3)
  return(clusters)
}


#' Annotate Probe Statistics with Clustering Information
#'
#' This function merges probe statistics with probe annotations and clustering information,
#' then annotates the merged data with additional probe characteristics.
#'
#' @param probe_stats A dataframe containing statistics calculated by GetStats() for each probe.
#' @param probe_clusters A clustering object with clusters for each probe based on methylation status calculated by GetClusters().
#' @keywords internal
#' @importFrom dplyr mutate_at case_when if_else vars mutate %>%
#' @importFrom stringr str_count
#' @return A dataframe of annotated probe statistics including clustering information.
AnnotateStats <- function(probe_stats, probe_clusters){
  # Assuming probe_annotations is available in the package's data directory

  # Merging probe statistics with global probe annotations
  # This adds additional details like Color, SourceSeq, etc.
  probe_stats <- merge(probe_stats, probe_annotations, by="row.names", all.x=TRUE)

  # Counting occurrences of "CG" in the SourceSeq to identify CpG sites
  probe_stats$nCpG <- stringr::str_count(probe_stats$SourceSeq, "CG")

  # Renaming columns for RF models
  colnames(probe_stats)[3] <- 'q1'
  colnames(probe_stats)[4] <- 'q3'

  # This block replaces empty values with "none" and transforms categorical columns into factors
  probe_stats <- probe_stats %>%
    dplyr::mutate_at(dplyr::vars(Color, DMR, Enhancer, regFeature, DHS),
              ~ dplyr::if_else(. == "", "none", .)) %>%
    dplyr::mutate(SBE_rs = dplyr::case_when(is.na(SBE_rs) ~ "none",
                              !is.na(SBE_rs) ~ "snp")) %>%
    dplyr::mutate_at(dplyr::vars(Color, SBE_rs, DMR, Enhancer, regFeature, location, DHS),
              as.factor)

  # Merging with clustering information
  probe_stats <- merge(probe_stats, probe_clusters$cluster, by.x="Row.names", by.y='row.names', all.x=TRUE)
  colnames(probe_stats)[22] <- 'clusters_repl.cluster' #cluster labels should be named clusters_repl.cluster for the model
  probe_stats$clusters_repl.cluster <- as.factor(probe_stats$clusters_repl.cluster)

  # Setting row names and removing unnecessary columns
  rownames(probe_stats) <- probe_stats$Name
  probe_stats <- probe_stats[, !(colnames(probe_stats) %in% c('Row.names', 'SourceSeq', 'Name'))]

  #ensure that the categorical variables have all levels requred by rf models
  probe_stats <- dplyr::mutate(probe_stats,
                        Color = factor(Color, levels = c("Grn", "none", "Red")),
                        DMR = factor(DMR, levels = c("CDMR", "DMR", "none", "RDMR")),
                        Enhancer = factor(Enhancer, levels = c("none", "TRUE")),
                        DHS = factor(DHS, levels = c("none", "TRUE")),
                        SBE_rs = factor(SBE_rs, levels = c("none", "snp")),
                        regFeature = factor(regFeature, levels = c("Gene_Associated", "Gene_Associated_Cell_type_specific", "none",
                                                                  "NonGene_Associated", "NonGene_Associated_Cell_type_specific",
                                                                  "Promoter_Associated", "Promoter_Associated_Cell_type_specific",
                                                                  "Unclassified", "Unclassified_Cell_type_specific")),
                        location = factor(location, levels = c("Island", "N_Shelf", "N_Shore", "OpenSea", "S_Shelf", "S_Shore")),
                        clusters_repl.cluster = factor(clusters_repl.cluster, levels = c("1", "2", "3")))

  return(probe_stats)
}


#' Detect Stochastic Epigenetic Mutations (SEMs) Using IQR Method
#'
#' This function detects hyper and hypo stochastic epigenetic mutations (SEMs) in DNA methylation data
#' by applying an interquartile range (IQR) based method.
#'
#' @param data A data frame or matrix with methylation data, samples as rows and CpG sites as columns.
#' @param reference Optional reference data frame or matrix to calculate quantiles for SEM detection.
#'        If not provided, `data` itself is used as reference.
#' @return A list with two elements: 'hyperSEM' and 'hypoSEM', vectors containing counts of hyper and hypo SEMs per sample.
#' @keywords internal
#############################################################ORIGINAL SEM METHOD
#' @importFrom stats quantile
#' @note Intended for internal package use.
CountSEM <- function(data, reference=data) {
  # Initialize vectors to hold counts of SEMs
  countsHyper = integer(nrow(data))
  countsHypo = integer(nrow(data))

  # Loop through each CpG site to detect SEMs
  for(i in 1:ncol(data)){
    probe.name <- colnames(data)[i]
    # Calculate the lower and upper quartiles from the reference data
    lowerq = stats::quantile(reference[, probe.name])[2]
    upperq = stats::quantile(reference[, probe.name])[4]

    # Calculate the interquartile range
    iqr = upperq - lowerq
    if(iqr == 0){next} # Skip if IQR is zero

    # Define SEM detection thresholds
    extreme.threshold.upper = (iqr * 3) + upperq
    extreme.threshold.lower = lowerq - (iqr * 3)

    # Detect hyper and hypo SEMs based on thresholds
    resultHyper <- which(data[,i] > extreme.threshold.upper)
    resultHypo <- which(data[,i] < extreme.threshold.lower)

    # Increment SEM counts for affected samples
    for(j in resultHyper){
      countsHyper[j] = countsHyper[j] + 1
    }
    for(j in resultHypo){
      countsHypo[j] = countsHypo[j] + 1
    }
  }
  rm(list=c("lowerq", "upperq", "iqr", "extreme.threshold.upper", "extreme.threshold.lower", "j", "i"))
  return(list("hyperSEM"=countsHyper, "hypoSEM"=countsHypo))
}


#' Detect Stochastic Epigenetic Mutations (SEMs) Using Random Forest Models
#'
#' Applies a Random Forest model to refine the detection of stochastic epigenetic mutations (SEMs)
#' identified by the IQR method, considering additional features and cell counts.
#'
#' @param data A data frame or matrix with methylation data, samples as rows and CpG sites as columns.
#' @param model_features A data frame containing annotations for each feature.
#' @param cell_counts A data frame containing estimated cell counts for each sample by estimate.cell.proportions().
#' @return A list with two elements: 'hyperSEM' and 'hypoSEM', vectors containing counts of hyper and hypo SEMs per sample.
#' @keywords internal
#' @import randomForest
#' @importFrom dplyr mutate
#' @note This function is intended for internal use and requires pre-calculated statistics and cell count estimates.
CountSEM_rf <- function(data, model_features, cell_counts) {
  # Initialize vectors to hold counts of SEMs
  countsHyper = integer(nrow(data))
  countsHypo = integer(nrow(data))

  # Define features used by the Random Forest models
  features_hyper <- c("skewness", "kurtosis", "q3", "range", "stdev", "iqr", "clusters_repl.cluster", "Enhancer", "DHS", "regFeature", "DMR")
  features_hypo <- c("min", "max", "mean", "nCpG", "iqr", "range", "stdev", "clusters_repl.cluster", "SBE_rs", "Color", "location")

  # Loop over each CpG site to detect SEMs using the Random Forest model
  for(i in 1:ncol(data)){
    probe.name <- colnames(data)[i]

    stats_row <- model_features[probe.name, ]

    iqr <- stats_row[["iqr"]]

    if(iqr == 0){next} # Skip if IQR is zero
    lowerq <- stats_row[["q1"]]
    upperq <- stats_row[["q3"]]

    # Define SEM detection thresholds
    extreme.threshold.upper = (iqr * 3) + upperq    #3IQR
    extreme.threshold.lower = lowerq - (iqr * 3)   #3IQR

    # Detect SEMs and extract relevant cell counts
    resultHyper <- which(data[,i] > extreme.threshold.upper)
    resultHypo <- which(data[,i] < extreme.threshold.lower)

    cellsHyper <- cell_counts[resultHyper, c("CD8T", "Bcell", "NK", "CD4T")]
    cellsHypo <-  cell_counts[resultHypo, c("CD8T", "Bcell", "CD4T")]

    # Calculate distances from threshold and use Random Forest models for prediction
    IQRhyper <- (data[resultHyper, i] - upperq) / iqr
    IQRhypo <- (data[resultHypo, i] - lowerq) / iqr

    if(length(IQRhyper) > 0){

      predictionHyper <- predict(treeModelData_hyper,
                                 newdata = dplyr::mutate((stats_row[features_hyper])[rep(1, each = length(IQRhyper)),],
                                                  beta=IQRhyper,
                                                  CD8T=cellsHyper[["CD8T"]],
                                                  Bcell=cellsHyper[["Bcell"]],
                                                  NK=cellsHyper[["NK"]],
                                                  CD4T=cellsHyper[["CD4T"]]))

      countsHyper[resultHyper] <- countsHyper[resultHyper] + (as.integer(predictionHyper) - 1)
    }
    if(length(IQRhypo) > 0){

      predictionHypo <- predict(treeModelData_hypo,
                                newdata = dplyr::mutate((stats_row[features_hypo])[rep(1, each = length(IQRhypo)),],
                                                 beta=IQRhypo,CD8T=cellsHypo[["CD8T"]],
                                                 Bcell=cellsHypo[["Bcell"]],
                                                 CD4T=cellsHypo[["CD4T"]]))

      countsHypo[resultHypo] <- countsHypo[resultHypo] + (as.integer(predictionHypo) - 1)
    }

  }

  rm(list=c("lowerq", "upperq", "iqr", "probe.name", "extreme.threshold.upper", "extreme.threshold.lower", "resultHyper", "resultHypo", "IQRhyper", "IQRhypo", "i"))
  return(list("hyperSEM"=countsHyper, "hypoSEM"=countsHypo))
}
