utils::globalVariables(c("probe_annotations"))
#' @title Detect SEMs in DNA Methylation Beta Values
#' @details This function serves as the main entry point for detecting stochastic epigenetic mutations (SEMs)
#' in DNA methylation data using either an IQR-based method or a more sophisticated Random Forest (RF)-based method.
#' It also offers the capability to cluster probes based on their methylation status
#' and analyze SEMs accordingly. The function is designed to be flexible, allowing for parallel processing to speed up
#' computations on large datasets. User interaction is required for confirming proceeding after probe clustering when
#' running interactively. For automated scripts, consider wrapping calls to this function in a way that handles user prompts programmatically.
#' @importFrom plyr mapvalues
#' @importFrom parallel detectCores makeCluster clusterEvalQ clusterExport parLapply stopCluster
#' @importFrom graphics hist
#'
#' @param betas DNA methylation dataframe with beta values. Rows are subjects, columns are probes. Columns must have probe names from the Illumina 450k array if using the RF-based detection method.
#' @param num_cores Number of cores to be used. Default is 1 (no parallel processing).
#' @param rf Logical (TRUE/FALSE) indicating whether to use the RF-based method. The default is FALSE (using the IQR-based method).
#' @param probes An optional list of character elements containing the names of the probes in 'betas' that should be analyzed. Default is to analyze all probes present in 'betas'.
#' @param cluster Logical (TRUE/FALSE) indicating whether to detect cluster probes into unmethylated, intermediate, and methylated probes based on their methylation status in 'betas' (later will incorporate a reference for blood) and group SEMs based on the type of probe they originate from. Clustering is required for RF-based models.
#'
#' @return A dataframe containing different types of SEM counts. Subjects are in organized in the same order as in \code{betas}.
#' @export
################################################################################################
detectSEM <- function(betas, num_cores=1, rf=FALSE, probes=NULL, cluster=FALSE) {
  ################################################################################################
  # Check if 'betas' is a dataframe
  if (!is.data.frame(betas)) {
    stop("'betas' must be a data frame with subjects as rows and CpGs as columns.")
  }

  if(all(betas >= 0 & betas <= 1)) {
    cat("All DNA methylation values are within the range [0, 1].\n")
  } else {
    stop("Some elements in 'betas' are outside the range [0, 1]. Possibly NA values?\n")
  }

  if(rf & ncol(betas) <= 2e+05) {
    warning("Houseman et al. suggest 'you use at least 2x10^5 probes' for cell count inference (required by the RF models)")
  }

  # Check if 'betas' has enough observations to calculate stats
  if (nrow(betas) < 2) {
    stop(paste("'betas' must have at least 2 observations (preferably more)."))
  } else if(nrow(betas) < 50) {
    warning("We recommend to have at least 50 subjects for reliable SEM calls.")
  }

  # Ensure rf is either TRUE or FALSE
  if (!is.logical(rf) || length(rf) != 1) {
    stop("'rf' must be either TRUE or FALSE.")
  }

  # Ensure cluster is either TRUE or FALSE
  if (!is.logical(cluster) || length(cluster) != 1) {
    stop("'cluster' must be either TRUE or FALSE.")
  }

  # Check if 'probes' is a character vector
  if (!is.null(probes) && !is.character(probes)) {
    stop("If specifying a list of probes, all elements must be of character type.")
  }

  if ((rf || cluster) && ncol(betas) < 3) {
    stop("There must be at least 3 probes for clustering (clustering is required if running the RF model).")
  }

  # Ensure num_cores is a positive integer >= 1
  if (!is.numeric(num_cores) || num_cores < 1 || num_cores != round(num_cores)) {
    stop("'num_cores' must be a positive integer >= 1.")
  }

  set.seed(1)

  ########################################################################



  if(num_cores > 1){
    #check whether the system has that many cores
    avaliable_cores <- parallel::detectCores()
    if (num_cores > avaliable_cores) {
      warning(paste("Specified num_cores exceeds available cores. Using", avaliable_cores-2, "instead."))
      num_cores <- avaliable_cores-2
    }
    rm(avaliable_cores)
  }


  #cluster probes by methylation status
  if (rf || cluster) {

    cat("Calculating probe statistics...\n")
    if (num_cores > 1) {
      #get statistics of probes
      probe_stats <- GetStats(betas=betas, num_cores=num_cores)
    } else {
      probe_stats <- GetStats(betas=betas, num_cores=1)
    }

    cat("Clustering probes by methylation status...\n")
    probe_clusters <- GetClusters(probe_stats)

    cluster_means <- probe_clusters$centers[,'mean']

    # Order the cluster means of means and use the order to rank them
    # This will give ranks such that the lowest mean gets rank 1, and so on
    ranked_clusters <- order(cluster_means)

    #map cluster labels so that 1=unmethylated, 2=intermediate, 3=methylated
    probe_clusters$cluster <- plyr::mapvalues(probe_clusters$cluster, c(ranked_clusters), c(1,2,3))


    # First, plot the histograms
    graphics::hist(probe_stats[probe_clusters$cluster == 1, "mean"], main="Unmethylated probes", xlab="beta means of probes")
    # You might want to pause between plots, depending on how you're executing this (e.g., in a script vs. interactively in RStudio)
    Sys.sleep(2) # Pause for 2 seconds (adjust timing as needed)

    graphics::hist(probe_stats[probe_clusters$cluster == 2, "mean"], main="Intermediate probes", xlab="beta means of probes")
    Sys.sleep(2) # Pause

    graphics::hist(probe_stats[probe_clusters$cluster == 3, "mean"], main="Methylated probes", xlab="beta means of probes")
    Sys.sleep(2) # Pause

    # Initialize a variable to store the user's response
    valid_response <- FALSE
    while (!valid_response) {
      response <- readline(prompt = "Here're the results of clustering. Please inspect the 3 histograms in the Plots pannel. Do you want to proceed with the analysis? (yes/no): ")
      response <- tolower(response)  # Normalize the response to lowercase

      if (response == "yes") {
        cat("Proceeding with the analysis...\n")
        valid_response <- TRUE  # Set flag to TRUE to exit the loop
      } else if (response == "no") {
        rm(probe_clusters)
        cat("Analysis terminated by the user.\n")
        return(NULL)  # Terminate the function early
      } else {
        cat("Invalid response. Please answer 'yes' or 'no'.\n")
        # The loop will continue, giving the user another chance to respond
      }
    }
    remove(valid_response, response, cluster_means, ranked_clusters)

    if (!is.null(probes)) {
      probe_stats <- probe_stats[probes, ]
    }
  }

  if(rf) {
    #now we add probe annotations to the model_features dataframe
    cat("Creating probe annotations for the RF models...\n")

    tryCatch({
      model_features <- AnnotateStats(probe_stats, probe_clusters, array=array)
    }, error = function(e) {
      cat("Error in AnnotateStats: ", e$message, "\n")
      return()  # Exiting the main function due to the error
  })
    
    
    betas <- betas[, colnames(betas) %in% probe_annotations$Name]
    if(!is.null(probes)){
      probes <- intersect(probes, probe_annotations$Name)
    }
    cat("Estimating cell counts...\n")
    betas_transpose <- t(betas)
    cell_counts <- estimate.cell.proportions(betas_transpose)
    cell_counts <- as.data.frame(cell_counts)
    rm(betas_transpose)
  }

  if (!is.null(probes)) {
    betas <- betas[, colnames(betas) %in% probes]
  } else {
    probes <- colnames(betas)
  }

  cat("Cleaning memory...\n")
  gc()



  if(num_cores > 1){
    cat("Preparing multiple cores...\n")
    cl <- parallel::makeCluster(num_cores)

    if(rf){
      parallel::clusterEvalQ(cl, library(randomForest))
      parallel::clusterEvalQ(cl, library(dplyr))
      parallel::clusterExport(cl=cl, varlist = c("CountSEM_rf", "betas", "model_features", "cell_counts", 'treeModelData_hypo', 'treeModelData_hyper'), envir = environment())

      if(!cluster) {
        cat("Detecting SEMs with RF models...\n")
        indices <- 1:ncol(betas)

        results <- parallel::parLapply(cl, indices, function(index) {
          data_chunk <- betas[, index, drop = FALSE]  # Keep it as a dataframe
          model_feature_rows <- model_features[index, , drop = FALSE]  # Corresponding model features
          CountSEM_rf(data_chunk, model_feature_rows, cell_counts)
        })

        SEM_results <- Reduce(function(x, y) {
          list(hyperSEM = x$hyperSEM + y$hyperSEM, hypoSEM = x$hypoSEM + y$hypoSEM)
        }, results)
        SEM_results <- as.data.frame(SEM_results)

        parallel::stopCluster(cl)
        rm(indices, cl, results, cell_counts)
        gc()

        cat("DONE\n")

        return(SEM_results)

      } else {
        # if user requested clustering with RF
        indices_unmethylated <- intersect(names(probe_clusters$cluster)[names(probe_clusters$cluster) %in% probes],
                                          names(probe_clusters$cluster[probe_clusters$cluster == 1]))
        indices_intermediate <- intersect(names(probe_clusters$cluster)[names(probe_clusters$cluster) %in% probes],
                                          names(probe_clusters$cluster[probe_clusters$cluster == 2]))
        indices_methylated <- intersect(names(probe_clusters$cluster)[names(probe_clusters$cluster) %in% probes],
                                        names(probe_clusters$cluster[probe_clusters$cluster == 3]))

        cat("Detecting SEMs in unmethylated probes with RF models...\n")
        results_unmethylated <- parallel::parLapply(cl, indices_unmethylated, function(index) {
          data_chunk <- betas[, index, drop = FALSE]  # Keep it as a dataframe
          model_feature_rows <- model_features[index, , drop = FALSE]  # Corresponding model features
          CountSEM_rf(data_chunk, model_feature_rows, cell_counts)
        })
        cat("Cleaning memory...\n")
        gc()

        cat("Detecting SEMs in intermediate probes with RF models...\n")
        results_intermediate <- parallel::parLapply(cl, indices_intermediate, function(index) {
          data_chunk <- betas[, index, drop = FALSE]  # Keep it as a dataframe
          model_feature_rows <- model_features[index, , drop = FALSE]  # Corresponding model features
          CountSEM_rf(data_chunk, model_feature_rows, cell_counts)
        })
        cat("Cleaning memory...\n")
        gc()

        cat("Detecting SEMs in methylated probes with RF models...\n")
        results_methylated <- parallel::parLapply(cl, indices_methylated, function(index) {
          data_chunk <- betas[, index, drop = FALSE]  # Keep it as a dataframe
          model_feature_rows <- model_features[index, , drop = FALSE]  # Corresponding model features
          CountSEM_rf(data_chunk, model_feature_rows, cell_counts)
        })
        cat("Cleaning memory...\n")
        gc()

        SEM_results_unmethylated <- Reduce(function(x, y) {
          list(hyperSEM = x$hyperSEM + y$hyperSEM, hypoSEM = x$hypoSEM + y$hypoSEM)
        }, results_unmethylated)

        SEM_results_intermediate <- Reduce(function(x, y) {
          list(hyperSEM = x$hyperSEM + y$hyperSEM, hypoSEM = x$hypoSEM + y$hypoSEM)
        }, results_intermediate)

        SEM_results_methylated <- Reduce(function(x, y) {
          list(hyperSEM = x$hyperSEM + y$hyperSEM, hypoSEM = x$hypoSEM + y$hypoSEM)
        }, results_methylated)

        results_DF <- data.frame(unmethylated_probes=SEM_results_unmethylated,
                                 intermediate_probes=SEM_results_intermediate,
                                 methylated_probes=SEM_results_methylated)
        results_DF$all_probes.hypoSEM <- results_DF$unmethylated_probes.hypoSEM + results_DF$intermediate_probes.hypoSEM + results_DF$methylated_probes.hypoSEM
        results_DF$all_probes.hyperSEM <- results_DF$unmethylated_probes.hyperSEM + results_DF$intermediate_probes.hyperSEM + results_DF$methylated_probes.hyperSEM

        cat("DONE\n")
        parallel::stopCluster(cl)
        rm(indices_unmethylated, indices_intermediate, indices_methylated, cl, results_unmethylated, results_intermediate, results_methylated, SEM_results_unmethylated, SEM_results_intermediate, SEM_results_methylated, cell_counts)
        gc()
        return(results_DF)
      }

    } else {
      parallel::clusterExport(cl=cl, varlist = c("betas", "CountSEM"), envir = environment())

      if(!cluster){
        cat("Detecting SEMs with the original method...\n")

        indices <- 1:ncol(betas)

        results <- parallel::parLapply(cl, indices, function(index) {
          data_chunk <- betas[, index, drop = FALSE]  # Keep it as a dataframe
          CountSEM(data_chunk)
        })

        SEM_results <- Reduce(function(x, y) {
          list(hyperSEM = x$hyperSEM + y$hyperSEM, hypoSEM = x$hypoSEM + y$hypoSEM)
        }, results)
        SEM_results <- as.data.frame(SEM_results)

        parallel::stopCluster(cl)
        rm(indices, cl, results)
        gc()
        cat("DONE\n")

        return(SEM_results)

      } else {
        #clusters
        indices_unmethylated <- intersect(names(probe_clusters$cluster)[names(probe_clusters$cluster) %in% probes],
                                          names(probe_clusters$cluster[probe_clusters$cluster == 1]))
        indices_intermediate <- intersect(names(probe_clusters$cluster)[names(probe_clusters$cluster) %in% probes],
                                          names(probe_clusters$cluster[probe_clusters$cluster == 2]))
        indices_methylated <- intersect(names(probe_clusters$cluster)[names(probe_clusters$cluster) %in% probes],
                                        names(probe_clusters$cluster[probe_clusters$cluster == 3]))

        cat("Detecting SEMs in unmethylated probes with the original method...\n")
        results_unmethylated <- parallel::parLapply(cl, indices_unmethylated, function(index) {
          data_chunk <- betas[, index, drop = FALSE]  # Keep it as a dataframe
          CountSEM(data_chunk)
        })
        cat("Cleaning memory...\n")
        gc()

        cat("Detecting SEMs in intermediate probes with the original method...\n")
        results_intermediate <- parallel::parLapply(cl, indices_intermediate, function(index) {
          data_chunk <- betas[, index, drop = FALSE]  # Keep it as a dataframe
          CountSEM(data_chunk)
        })
        cat("Cleaning memory...\n")
        gc()

        cat("Detecting SEMs in methylated probes with the original method...\n")
        results_methylated <- parallel::parLapply(cl, indices_methylated, function(index) {
          data_chunk <- betas[, index, drop = FALSE]  # Keep it as a dataframe
          CountSEM(data_chunk)
        })
        cat("Cleaning memory...\n")
        gc()


        SEM_results_unmethylated <- Reduce(function(x, y) {
          list(hyperSEM = x$hyperSEM + y$hyperSEM, hypoSEM = x$hypoSEM + y$hypoSEM)
        }, results_unmethylated)

        SEM_results_intermediate <- Reduce(function(x, y) {
          list(hyperSEM = x$hyperSEM + y$hyperSEM, hypoSEM = x$hypoSEM + y$hypoSEM)
        }, results_intermediate)

        SEM_results_methylated <- Reduce(function(x, y) {
          list(hyperSEM = x$hyperSEM + y$hyperSEM, hypoSEM = x$hypoSEM + y$hypoSEM)
        }, results_methylated)


        results_DF <- data.frame(unmethylated_probes=SEM_results_unmethylated,
                                 intermediate_probes=SEM_results_intermediate,
                                 methylated_probes=SEM_results_methylated)
        results_DF$all_probes.hypoSEM <- results_DF$unmethylated_probes.hypoSEM + results_DF$intermediate_probes.hypoSEM + results_DF$methylated_probes.hypoSEM
        results_DF$all_probes.hyperSEM <- results_DF$unmethylated_probes.hyperSEM + results_DF$intermediate_probes.hyperSEM + results_DF$methylated_probes.hyperSEM

        parallel::stopCluster(cl)
        rm(indices_unmethylated, indices_intermediate, indices_methylated, cl, results_unmethylated, results_intermediate, results_methylated, SEM_results_unmethylated, SEM_results_intermediate, SEM_results_methylated)
        gc()
        cat("DONE\n")

        return(results_DF)

      }
    }
  } else {
    #single core
    if(rf) {
      if(!cluster) {
        cat("Detecting SEMs with the RF models...\n")
        SEM_results <- CountSEM_rf(betas, model_features, cell_counts)
        SEM_results <- as.data.frame(SEM_results)

        rm(cell_counts)
        cat("DONE\n")
        return(SEM_results)

      } else {
        indices_unmethylated <- intersect(names(probe_clusters$cluster)[names(probe_clusters$cluster) %in% probes],
                                          names(probe_clusters$cluster[probe_clusters$cluster == 1]))
        indices_intermediate <- intersect(names(probe_clusters$cluster)[names(probe_clusters$cluster) %in% probes],
                                          names(probe_clusters$cluster[probe_clusters$cluster == 2]))
        indices_methylated <- intersect(names(probe_clusters$cluster)[names(probe_clusters$cluster) %in% probes],
                                        names(probe_clusters$cluster[probe_clusters$cluster == 3]))

        cat("Detecting SEMs in unmethylated probes with the RF models...\n")
        SEM_results_unmethylated <- CountSEM_rf(betas[, indices_unmethylated, drop = FALSE],
                                                model_features[indices_unmethylated, , drop = FALSE],
                                                cell_counts)
        cat("Cleaning memory...\n")
        gc()

        cat("Detecting SEMs in intermediate probes with the RF models...\n")
        SEM_results_intermediate <- CountSEM_rf(betas[, indices_intermediate, drop = FALSE],
                                                model_features[indices_intermediate, , drop = FALSE],
                                                cell_counts)
        cat("Cleaning memory...\n")
        gc()

        cat("Detecting SEMs in methylated probes with the RF models...\n")
        SEM_results_methylated <- CountSEM_rf(betas[, indices_methylated, drop = FALSE],
                                              model_features[indices_methylated, , drop = FALSE],
                                              cell_counts)
        cat("Cleaning memory...\n")
        gc()


        results_DF <- data.frame(unmethylated_probes=SEM_results_unmethylated,
                                 intermediate_probes=SEM_results_intermediate,
                                 methylated_probes=SEM_results_methylated)
        results_DF$all_probes.hypoSEM <- results_DF$unmethylated_probes.hypoSEM + results_DF$intermediate_probes.hypoSEM + results_DF$methylated_probes.hypoSEM
        results_DF$all_probes.hyperSEM <- results_DF$unmethylated_probes.hyperSEM + results_DF$intermediate_probes.hyperSEM + results_DF$methylated_probes.hyperSEM

        rm(indices_unmethylated, indices_intermediate, indices_methylated, SEM_results_unmethylated, SEM_results_intermediate, SEM_results_methylated, cell_counts)

        cat("DONE\n")
        return(results_DF)
      }
    } else {
      if(!cluster) {
        cat("Detecting SEMs with the original method...\n")
        SEM_results <- CountSEM(betas)
        SEM_results <- as.data.frame(SEM_results)

        cat("DONE\n")
        return(SEM_results)

      } else {
        indices_unmethylated <- intersect(names(probe_clusters$cluster)[names(probe_clusters$cluster) %in% probes],
                                          names(probe_clusters$cluster[probe_clusters$cluster == 1]))
        indices_intermediate <- intersect(names(probe_clusters$cluster)[names(probe_clusters$cluster) %in% probes],
                                          names(probe_clusters$cluster[probe_clusters$cluster == 2]))
        indices_methylated <- intersect(names(probe_clusters$cluster)[names(probe_clusters$cluster) %in% probes],
                                        names(probe_clusters$cluster[probe_clusters$cluster == 3]))

        cat("Detecting SEMs in unmethylated probes with the original method...\n")
        SEM_results_unmethylated <- CountSEM(betas[, indices_unmethylated, drop = FALSE])
        cat("Cleaning memory...\n")
        gc()

        cat("Detecting SEMs in intermediate probes with the original method...\n")
        SEM_results_intermediate <- CountSEM(betas[, indices_intermediate, drop = FALSE])
        cat("Cleaning memory...\n")
        gc()

        cat("Detecting SEMs in methylated probes with the original method...\n")
        SEM_results_methylated <- CountSEM(betas[, indices_methylated, drop = FALSE])
        cat("Cleaning memory...\n")
        gc()

        results_DF <- data.frame(unmethylated_probes=SEM_results_unmethylated,
                                 intermediate_probes=SEM_results_intermediate,
                                 methylated_probes=SEM_results_methylated)
        results_DF$all_probes.hypoSEM <- results_DF$unmethylated_probes.hypoSEM + results_DF$intermediate_probes.hypoSEM + results_DF$methylated_probes.hypoSEM
        results_DF$all_probes.hyperSEM <- results_DF$unmethylated_probes.hyperSEM + results_DF$intermediate_probes.hyperSEM + results_DF$methylated_probes.hyperSEM

        rm(indices_unmethylated, indices_intermediate, indices_methylated, SEM_results_unmethylated, SEM_results_intermediate, SEM_results_methylated)

        cat("DONE\n")
        return(results_DF)
      }
    }
  }
}
