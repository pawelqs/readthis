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



get_files <- function(path, pattern, sample_id_pattern) {
  all_files <- list.files(path, full.names = TRUE)
  files <- grep(pattern, all_files, value = TRUE)
  names(files) <- str_extract(files, pattern = sample_id_pattern)
  files
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


# analyze_path <- function(path) {
#   if () {
#     "MANY_FILES"
#   }
# }


is_single_dir <- function(path) {
  if (length(path) > 1) {
    FALSE
  } else {
    file.exists(path) & dir.exists(path)
  }
}


is_single_file <- function(path) {
  if (length(path) > 1) {
    FALSE
  } else {
    file.exists(path) & (!dir.exists(path))
  }
}


is_list_of_files <- function(path) {
  if (length(path) == 1) {
    FALSE
  } else {
    all(map_lgl(path, is_single_file))
  }
  # length(path) > 1 && all(is_file(path))
}
