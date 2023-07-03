msg <- function(...,
                prefix = "readthis>",
                collapse = "",
                col = "steelblue3",
                verbose = TRUE) {
  msg <- str_c(list(prefix, ...), collapse = collapse)
  if (verbose) {
    cli::cat_line(msg, col = col)
  }
}


use_chrom_naming_convention <- function(tbl, chrom_convention) {
  if (chrom_convention == "UCSC") {
    use_UCSC_chrom_names(tbl)
  } else if (chrom_convention == "NCBI") {
    use_NCBI_chrom_names(tbl)
  } else if (chrom_convention == "keep") {
    tbl
  } else {
    stop("chrom_convention should be on of: UCSC/NCBI/keep")
  }
}


use_UCSC_chrom_names <- function(tbl) {
  tbl |>
    mutate(
      chrom = if_else(
        str_detect(.data$chrom, "^chr"),
        as.character(.data$chrom),
        str_c("chr", .data$chrom)
      )
    )
}


use_NCBI_chrom_names <- function(tbl) {
  tbl |>
    mutate(
      chrom = if_else(
        str_detect(.data$chrom, "^chr"),
        str_replace(.data$chrom, "^chr", ""),
        as.character(.data$chrom)
      )
    )
}
