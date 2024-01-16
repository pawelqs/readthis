## ------------------------------- Export --------------------------------------

#' Read FACETS CNA calls
#'
#' Reads the variant calls from [FACETS](https://github.com/mskcc/facets/tree/master)
#' CNA caller
#'
#' @param path Can be either:
#'   1) path to a single file, sample ID can be passed using sample_id argument
#'   2) vector of file paths, element names will be used as sample IDs
#'   3) directory containing multiple FACETS .csv files, sample IDs will be
#'      guessed from the file names
#' @param sample_id Sample ID, used with path is a single file
#' @param chrom_convention UCSC/NCBI/keep
#'
#' @examples
#' library(readthis)
#'
#' file1 <- system.file("extdata", "FACETS", "S1.csv", package = "readthis")
#' read_facets_cnas(file1)
#'
#' file2 <- system.file("extdata", "FACETS", "S2.csv", package = "readthis")
#' files <- c(S1 = file1, S2 = file2)
#' read_facets_cnas(files)
#'
#' dir <- system.file("extdata", "FACETS", package = "readthis")
#' read_facets_cnas(dir)
#'
#' @name facets
NULL


#' @rdname facets
#' @details
#' read_facets_csv() reads a single csv file. If sample_id is not provided,
#'   file path is used
#' @export
read_facets_cnas <- function(path,
                             sample_id = path,
                             chrom_convention = "UCSC") {
  if (is_single_file(path)) {
    cnas <- read_facets_csv(path) |>
      mutate(sample_id = sample_id, .before = "chrom")
  } else if (is_list_of_files(path)) {
    cnas <- map(path, read_facets_csv) |>
      bind_rows(.id = "sample_id")
  } else if (is_single_dir(path)) {
    files <- list.files(path, full.names = TRUE)
    sample_ids <- files |>
      str_split(pattern = "/") |>
      map_chr(last) |>
      str_replace(".csv", "")
    names(files) <- sample_ids
    cnas <- files |>
      set_names(sample_ids) |>
      map(read_facets_csv) |>
      bind_rows(.id = "sample_id")
  }

  cnas |>
    use_chrom_naming_convention(chrom_convention)
}


## ----------------------------- Functions ------------------------------------


read_facets_csv <- function(file, chrom_convention = "UCSC") {
  cnas <- read_csv(file, col_types = "cddddddddiidiidd") |>
    prepare_FACETS_columns()
  class(cnas) <- c("cevo_FACETS", class(cnas))
  cnas
}


prepare_FACETS_columns <- function(tbl) {
  tbl |>
    rename(
      total_cn = "tcn.em",
      minor_cn = "lcn.em",
      log_ratio = "cnlr.median"
    ) |>
    mutate(
      major_cn = .data$total_cn - .data$minor_cn
    ) |>
    select(
      "chrom", "start", "end",
      "total_cn", "major_cn", "minor_cn", "log_ratio",
      everything()
    )
}

