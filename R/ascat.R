## ------------------------------- Export --------------------------------------

#' Read ASCAT CNV calls
#'
#' Reads the CNV variant calls and sample statistics from
#' [ASCAT](https://github.com/VanLoo-lab/ascat) CNV caller
#'
#' @param path Can be either:
#'   1) path to a single csv file with ASCAT CNV calls
#'   2) tibble with sample_id, cnvs, and (optionally) sample_statistics columns
#'      containing sample_ids and paths.
#'   3) directory containing multiple ASCAT files with "*.csv" and
#'      "*.samplestatistics.*" names.
#'
#'   If path is a single csv file, sampleID can be passed in sample_id argument.
#'   If not provided, path will be used in the output sample_id column.
#'   If path is a directory, sample names will be inferred using str_expr() and
#'   sample_id_pattern.
#' @param sample_statistics Path to the ".samplestatistics.* file, used when path
#'   leads to a single csv file
#' @param sample_id Sample ID, used with path is a sigle file
#' @param sample_id_pattern Pattern used to extract sample IDs from the filenames
#'   when path is a directory
#' @param chrom_convention UCSC/NCBI/keep
#'
#' @name ascat
NULL



#' @rdname ascat
#' @details
#' read_ascat_csv() ...
#' @export
read_ascat_files <- function(path,
                             sample_statistics = NULL,
                             sample_id = path,
                             sample_id_pattern = "(?<=\\/)[:alnum:]*(?=\\.)",
                             chrom_convention = "UCSC") {
  if (is_single_file(path)) {
    cnvs <- read_ascat_cnvs(path) |>
      mutate(sample_id = sample_id, .before = "chrom")
    if (!is.null(sample_statistics) && is_single_file(sample_statistics)) {
      stats <- read_ascat_samplestatistics(sample_statistics) |>
        mutate(sample_id = sample_id, .before = "normal_contamination")
    } else {
      stats <- empty_ascat_samplestatistics()
    }
  } else if (is.data.frame(path)) {
    cnvs <- path |>
      select("sample_id", "csv") |>
      deframe() |>
      map(read_ascat_cnvs) |>
      bind_rows(.id = "sample_id")
    if (is.null(path[["sample_statistics"]])) {
      stats <- empty_ascat_samplestatistics()
    } else {
      stats <- path |>
        select("sample_id", "sample_statistics") |>
        deframe() |>
        keep(\(x) !is.na(x)) |>
        map(read_ascat_samplestatistics) |>
        bind_rows(.id = "sample_id")
    }
  } else if (is_single_dir(path)) {
    csv_files <- get_files(path, ".csv", sample_id_pattern)
    cnvs <- csv_files |>
      map(read_ascat_cnvs) |>
      bind_rows(.id = "sample_id")

    stat_files <- get_files(path, "samplestatistics", sample_id_pattern)
    stats <- stat_files |>
      map(read_ascat_samplestatistics) |>
      bind_rows(.id = "sample_id")
  }

  ascat <- list(
    cnvs = use_chrom_naming_convention(cnvs, chrom_convention),
    sample_statistics = stats
  )
  structure(ascat, class = c("cevo_ASCAT"))
}



## ----------------------------- Functions ------------------------------------


read_ascat_cnvs <- function(path) {
  csv_cols <- c("i", "chrom", "start", "end", "normal_cn_total", "normal_cn_minor", "total_cn", "minor_cn")
  cnvs <- read_csv(path, col_names = csv_cols, col_types = "dcdddddd") |>
    mutate(major_cn = .data$total_cn - .data$minor_cn) |>
    select(
      "chrom", "start", "end",
      "total_cn", "major_cn", "minor_cn",
      normal_cn = "normal_cn_total"
    )

  cnvs
}



read_ascat_samplestatistics <- function(path) {
  sample_statistics <- path %>%
    read_delim(delim = " ", col_names = c("var", "val"), col_types = "cc") |>
    deframe() |>
    as_tibble_row() |>
    rename(
      normal_contamination = "NormalContamination",
      ploidy = "Ploidy",
      goodness_of_fit = "goodnessOfFit",
      gender_chr = "GenderChr",
      gender_chr_found = "GenderChrFound"
    ) |>
    mutate(across(
      c("normal_contamination", "ploidy", "rho", "psi", "goodness_of_fit"),
      as.double
    ))
  sample_statistics
}



empty_ascat_samplestatistics <- function() {
  tibble(
    normal_contamination = double(),
    ploidy = double(),
    rho = double(),
    psi = double(),
    goodness_of_fit = double(),
    gender_chr = character(),
    gender_chr_found = character()
  )
}
