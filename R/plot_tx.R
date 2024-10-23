plot_tx <- function(region_weights) {
  tx_plot <- ggplot() +
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 0), linetype = "solid", size = 1.5) +
    geom_rect(aes(xmin = region_weights[1], xmax = sum(region_weights[1:2]), ymin = -0.05, ymax = 0.05), fill = "#7e7e7e", color = "black", size = 0.5) +
    annotate("text", x = region_weights[1] / 2, y = -0.15, label = "5'UTR", size = 3.6) +
    annotate("text", x = region_weights[2] / 2 + region_weights[1], y = -0.15, label = "CDS", size = 3.6) +
    annotate("text", x = region_weights[3] / 2 + sum(region_weights[1:2]), y = -0.15, label = "3'UTR", size = 3.6) +
    ylim(-0.25, 0.05) +
    theme_void()

  return(tx_plot)
}
