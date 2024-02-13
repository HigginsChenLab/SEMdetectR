#' Illumina Probe Annotations Dataframe
#'
#' A dataframe containing probe annotations used as features in the Random Forest models.
#' Includes information on Name, Color, Source Sequence, presence of SNP at SBE, overlap with a differentially methylated region, Enhancer, regulatory features, DNase I hypersensitive sites, and type of genomic location of the probe.
#' @name probe_annotations.rda
#' @format A data frame with 485512 rows and 9 columns:
#' \describe{
#'   \item{\code{Name}}{Name of Illumina probe.}
#'   \item{\code{Color}}{Color of Illumina probe.}
#'   \item{\code{SourceSeq}}{Nucleotide sequence.}
#'   \item{\code{SBE_rs}}{SNPs at single-base extension.}
#'   \item{\code{DMR}}{Differentially methylated region.}
#'   \item{\code{Enhancer}}{Enhancer.}
#'   \item{\code{regFeature}}{Regulatory features.}
#'   \item{\code{DHS}}{DNase I hypersensitive sites.}
#'   \item{\code{location}}{Genomic location.}
#' }
#' @source Taken from the Illumina manifest for the 450k array
"probe_annotations"
