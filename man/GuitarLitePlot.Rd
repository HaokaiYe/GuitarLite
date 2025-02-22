% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GuitarLitePlot.R
\name{GuitarLitePlot}
\alias{GuitarLitePlot}
\title{Generate Metagene and Transcript Region Plots}
\usage{
GuitarLitePlot(
  x,
  txdb,
  labels = NULL,
  colors = NULL,
  region_weights = NULL,
  show.legend = TRUE,
  type = c("density", "line")
)
}
\arguments{
\item{x}{A \code{\link{GRanges}} or \code{\link{GRangesList}} object containing the genomic locations of features.}

\item{txdb}{A \code{\link{TxDb}} object containing transcript annotation information. This can be obtained from Bioconductor or using \code{\link{makeTxDbFromGFF}}.}

\item{labels}{Optional. A \code{character} vector specifying labels for each element in \code{x}. The order of \code{labels} must match the order of elements in \code{x}. If not provided, labels are derived from the names of \code{x}.}

\item{colors}{Optional. A \code{character} vector specifying colors for each label. The order of \code{colors} must match the order of elements in \code{x}. If not provided, colors are generated automatically.}

\item{region_weights}{Optional. A \code{numeric} vector specifying the relative weights for the 5' UTR, CDS, and 3' UTR regions. If not provided, they are calculated based on the transcript regions.}

\item{show.legend}{Logical. Whether to display the legend. Defaults to \code{TRUE}.}

\item{type}{A \code{character} specifying the type of plot: \code{"density"} or \code{"line"}. Defaults to \code{"density"}.}
}
\value{
A \code{\link[patchwork]{patchwork}} object containing the metagene and transcript region plots.
}
\description{
\code{GuitarLitePlot} creates visualizations for metagene distributions and transcript region topology,
specifically for depicting the distribution of m6A modification sites across the transcript regions.
}
\details{
\code{GuitarLitePlot} supports input in the form of GRanges or GRangesList and integrates transcript annotation (TxDb) to generate metagene plots for 5' UTR, CDS, and 3' UTR regions.
The function calculates the relative positions of sites within each region by dividing the distance from the 5' end of the region by the total width of the region. A site at the start of the region has a relative position of 0, and a site at the end has a relative position of 1.
The metagene plot is generated using a density plot that visualizes the distribution of these relative positions.
}
\examples{

# Example 1: Visualizing a single sample using a GRanges object

library(GuitarLite)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Load the transcript annotation database
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Load a single sample as a GRanges object
x <- readRDS(system.file("extdata", "grg.rds", package = "GuitarLite"))

# Generate the topology plot for the single sample
topology_plot <- GuitarLitePlot(x, txdb)

# Display the plot
print(topology_plot)


# Example 2: Visualizing multiple samples using a GRangesList object

# Load multiple samples as a GRangesList object
x <- readRDS(system.file("extdata", "grl.rds", package = "GuitarLite"))

# Define region weights (optional, default is automatically calculated)
region_weights <- c(1/3, 1/3, 1/3)

# Define labels and colors for each sample
labels <- c("True m6A", "False positive")
colors <- c("#B1182D", "#2464AB")

# Generate the topology plot for multiple samples with custom settings
topology_plot <- GuitarLitePlot(x, txdb, region_weights = region_weights, labels = labels, colors = colors)

# Display the plot
print(topology_plot)


}
