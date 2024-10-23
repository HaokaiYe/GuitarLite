plot_metagene <- function(metagene_df, region_weights, labels, colors, show.legend = TRUE, type = c("line", "density")) {
  type <- match.arg(type)

  metagene_plot <- ggplot()

  if (type == "line") {
    metagene_plot <- metagene_plot +
      geom_line(data = metagene_df, aes(x = relative_pos, y = ..density.., color = class), stat = "density")
  } else if (type == "density") {
    metagene_plot <- metagene_plot +
      geom_density(data = metagene_df, aes(x = relative_pos, color = class, fill = class), alpha = 0.3)
  }

  metagene_plot <- metagene_plot +
    geom_vline(xintercept = c(region_weights[1], sum(region_weights[1:2])), linetype = 5) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    ylab("Density") +
    scale_color_manual(values = colors, breaks = labels, labels = labels) +
    scale_fill_manual(values = colors, breaks = labels, labels = labels) +
    theme_classic() +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(color = "black", size = 10),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(color = "black", size = 12),
          legend.background = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(color = "black", size = 10, hjust = 0))

  if (!show.legend) {
    metagene_plot <- metagene_plot + theme(legend.position = "none")
  }

  return(metagene_plot)
}
