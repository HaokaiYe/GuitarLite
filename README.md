# GuitarLite

### Introduction
**GuitarLite** is a lightweight R package designed for visualizing the distribution of genomic features, such as m6A modification sites, across transcript regions including the 5' UTR, CDS, and 3' UTR. By aligning these features into transcript regions and calculating their relative positions, **GuitarLite** generates metagene plots that depict the distribution of relative positions of features within each region. This package supports input in the form of `GRanges` or `GRangesList` objects and provides customizable options for labels, colors, and plot types to enhance data visualization and analysis.

### Installation
To install **GuitarLite** from GitHub, use the following command in the R console:

```r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("HaokaiYe/GuitarLite")
```

### Usage
First, load the package into R.
``` r
library(GuitarLite)
```

Example 1: Visualizing a single sample using a GRanges object
``` r
# Load the transcript annotation database
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Load a single sample as a GRanges object
x <- readRDS(system.file("extdata", "grg.rds", package = "GuitarLite"))

# Generate the topology plot for the single sample
topology_plot <- GuitarLitePlot(x, txdb)

# Display the plot
print(topology_plot)
```


Example 2: Visualizing multiple samples using a GRangesList object
``` r
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
```
