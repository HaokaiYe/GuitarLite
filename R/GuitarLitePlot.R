#' Generate Metagene and Transcript Region Plots
#'
#' @description \code{GuitarLitePlot} creates visualizations for metagene distributions and transcript region topology,
#' specifically for depicting the distribution of m6A modification sites across the transcript regions.
#'
#' @details \code{GuitarLitePlot} supports input in the form of GRanges or GRangesList and integrates transcript annotation (TxDb) to generate metagene plots for 5' UTR, CDS, and 3' UTR regions.
#' The function calculates the relative positions of sites within each region by dividing the distance from the 5' end of the region by the total width of the region. A site at the start of the region has a relative position of 0, and a site at the end has a relative position of 1.
#' The metagene plot is generated using a density plot that visualizes the distribution of these relative positions.
#'
#' @param x A \code{\link{GRanges}} or \code{\link{GRangesList}} object containing the genomic locations of features.
#' @param txdb A \code{\link{TxDb}} object containing transcript annotation information. This can be obtained from Bioconductor or using \code{\link{makeTxDbFromGFF}}.
#' @param labels Optional. A \code{character} vector specifying labels for each element in \code{x}. The order of \code{labels} must match the order of elements in \code{x}. If not provided, labels are derived from the names of \code{x}.
#' @param colors Optional. A \code{character} vector specifying colors for each label. The order of \code{colors} must match the order of elements in \code{x}. If not provided, colors are generated automatically.
#' @param region_weights Optional. A \code{numeric} vector specifying the relative weights for the 5' UTR, CDS, and 3' UTR regions. If not provided, they are calculated based on the transcript regions.
#' @param show.legend Logical. Whether to display the legend. Defaults to \code{TRUE}.
#' @param type A \code{character} specifying the type of plot: \code{"density"} or \code{"line"}. Defaults to \code{"density"}.
#'
#' @return A \code{\link[patchwork]{patchwork}} object containing the metagene and transcript region plots.
#'
#' @import ggplot2
#' @import patchwork
#' @importFrom GenomicRanges mapToTranscripts
#' @importFrom GenomicFeatures fiveUTRsByTranscript cdsBy threeUTRsByTranscript
#'
#'
#' @examples
#'
#' # Example 1: Visualizing a single sample using a GRanges object
#'
#' library(GuitarLite)
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#'
#' # Load the transcript annotation database
#' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#'
#' # Load a single sample as a GRanges object
#' x <- readRDS(system.file("extdata", "grg.rds", package = "GuitarLite"))
#'
#' # Generate the topology plot for the single sample
#' topology_plot <- GuitarLitePlot(x, txdb)
#'
#' # Display the plot
#' print(topology_plot)
#'
#'
#' # Example 2: Visualizing multiple samples using a GRangesList object
#'
#' # Load multiple samples as a GRangesList object
#' x <- readRDS(system.file("extdata", "grl.rds", package = "GuitarLite"))
#'
#' # Define region weights (optional, default is automatically calculated)
#' region_weights <- c(1/3, 1/3, 1/3)
#'
#' # Define labels and colors for each sample
#' labels <- c("True m6A", "False positive")
#' colors <- c("#B1182D", "#2464AB")
#'
#' # Generate the topology plot for multiple samples with custom settings
#' topology_plot <- GuitarLitePlot(x, txdb, region_weights = region_weights, labels = labels, colors = colors)
#'
#' # Display the plot
#' print(topology_plot)
#'
#'
#' @export

GuitarLitePlot <- function(x, txdb, labels = NULL, colors = NULL, region_weights = NULL, show.legend = TRUE, type = c("density", "line")) {
  if (!is(x, "GRanges") && !is(x, "GRangesList")) {
    stop("Input x must be either a GRanges or GRangesList object.")
  }
  type <- match.arg(type)

  message("## Processing input and generating metagene data ...")

  if (is(x, "GRangesList")) {
    labels <- if (is.null(labels)) names(x) else labels
    if (length(labels) != length(x)) stop("The length of 'labels' must match the length of the GRangesList.")

    colors <- if (is.null(colors)) scales::hue_pal()(length(labels)) else colors
    if (length(colors) != length(x)) stop("The length of 'colors' must match the length of the GRangesList.")

    metagene_list <- lapply(seq_along(x), function(i) {
      metagene(x[[i]], labels[i], txdb, region_weights)
    })
    metagene_df <- do.call(rbind, lapply(metagene_list, function(item) item$metagene_df))

  } else if (is(x, "GRanges")) {
    labels <- if (is.null(labels)) "Sample" else labels
    colors <- if (is.null(colors)) scales::hue_pal()(1) else colors

    metagene_list <- metagene(x, labels, txdb, region_weights)
    metagene_df <- metagene_list$metagene_df
  }

  message("## Plotting metagene distributions and transcript regions ...")

  region_weights <- if (is.null(region_weights)) {
    if (is(x, "GRangesList")) metagene_list[[1]]$region_weights else metagene_list$region_weights
  } else {
    region_weights
  }

  metagene_df$class <- factor(metagene_df$class, levels = labels)

  metagene_plot <- plot_metagene(metagene_df, region_weights, labels, colors, show.legend = show.legend, type = type)
  tx_plot <- plot_tx(region_weights)

  topology_plot <- metagene_plot / tx_plot + plot_layout(heights = c(10, 1))
  return(topology_plot)
}
