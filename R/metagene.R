relative_pos_on_region <- function(x, region) {
  region_map <- mapToTranscripts(x, region)
  region_width <- sum(width(region))[region_map$transcriptsHits]
  start_on_region <- start(region_map)
  return(start_on_region / region_width)
}

metagene <- function(x, x_label, txdb, region_weights = NULL) {
  u5bytx <- fiveUTRsByTranscript(txdb)
  cdsbytx <- cdsBy(txdb, by = "tx")
  u3bytx <- threeUTRsByTranscript(txdb)

  if (is.null(region_weights)) {
    UTR5_width <- mean(sum(width(u5bytx)))
    CDS_width <- mean(sum(width(cdsbytx)))
    UTR3_width <- mean(sum(width(u3bytx)))
    sum_width <- UTR5_width + CDS_width + UTR3_width
    region_weights <- c(UTR5_width / sum_width, CDS_width / sum_width, UTR3_width / sum_width)
  }

  utr5_pos <- relative_pos_on_region(x, u5bytx) * region_weights[1]
  cds_pos <- (relative_pos_on_region(x, cdsbytx) * region_weights[2]) + region_weights[1]
  utr3_pos <- (relative_pos_on_region(x, u3bytx) * region_weights[3]) + sum(region_weights[1:2])

  metagene_df <- data.frame(
    relative_pos = c(utr5_pos, cds_pos, utr3_pos),
    tx_region = factor(rep(c("5'UTR", "CDS", "3'UTR"),
                           c(length(utr5_pos), length(cds_pos), length(utr3_pos))),
                       levels = c("5'UTR", "CDS", "3'UTR")),
    class = x_label
  )

  return(list(metagene_df = metagene_df, region_weights = region_weights))
}
