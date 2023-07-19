## ------------------------------- Export --------------------------------------

#' Read Strelka variants
#'
#' Reads the variant calls from Strelka2 Small Variant Caller
#'
#'
#' @inheritParams mutect
#'
#' @examples
#' library(readthis)
#'
#' file1 <- system.file("extdata", "Strelka", "S1.somatic.snvs.vcf.gz", package = "readthis")
#' read_strelka_somatic_snvs(file1)
#'
#' file2 <- system.file("extdata", "Strelka", "S2.somatic.snvs.vcf.gz", package = "readthis")
#' files <- c(S1 = file1, S2 = file2)
#' read_strelka_somatic_snvs(files, verbose = FALSE)
#'
#' dir <- system.file("extdata", "Strelka", package = "readthis")
#' read_strelka_somatic_snvs(dir)
#'
#' @name strelka
NULL



#' @rdname strelka
#' @export
read_strelka_somatic_snvs <- function(path,
                                      sample_ids = "drop_first",
                                      PASS_only = TRUE,
                                      patient_id_pattern = "(?<=\\/)[:alnum:]*(?=\\.)",
                                      chrom_convention = "UCSC",
                                      verbose = TRUE) {
  if (is_single_file(path)) {
    snvs <- read_strelka_vcf(path, sample_ids, PASS_only, verbose)
  } else if (is_list_of_files(path)) {
    snvs <- map(path, read_strelka_vcf, sample_ids, PASS_only, verbose) |>
      bind_rows(.id = "patient_id")
  } else if (is_single_dir(path)) {
    files <- list.files(path, full.names = TRUE)
    patient_ids <- files |>
      str_extract(pattern = patient_id_pattern)
    names(files) <- patient_ids
    snvs <- files |>
      map(read_strelka_vcf, sample_ids, PASS_only, verbose) |>
      bind_rows(.id = "patient_id")
  }

  class(snvs) <- c("cevo_Strelka", class(snvs))
  snvs |>
    use_chrom_naming_convention(chrom_convention)
}



## ------------------------------ Functions ------------------------------------


read_strelka_vcf <- function(file,
                             sample_ids = NULL,
                             # overwrite_sample_id = NULL,
                             PASS_only = TRUE,
                             verbose = TRUE) {
  tidyvcf <- read_vcf(file, PASS_only, verbose = verbose)

  if (sample_ids == "drop_first") {
    sample_ids <- unique(tidyvcf$dat$Indiv)
    sample_ids <- sample_ids[-1]
    msg("Guessing sample_ids: ", str_c(sample_ids, collapse = ", "), verbose = verbose)
  } else if (sample_ids == "all") {
    sample_ids <- unique(tidyvcf$dat$Indiv)
  }

  snvs <- tidyvcf$dat |>
    filter(.data$Indiv %in% sample_ids) |>
    calc_strelka_somatic_VAF() |>
    select(
      sample_id = "Indiv",
      chrom = "CHROM", pos = "POS",
      ref = "REF", alt = "ALT",
      ref_reads = "ref_tier1",
      alt_reads = "alt_tier1",
      VAF = "VAF", DP = "DP"
    )

  # if (!is.null(overwrite_sample_id)) {
  #   snvs$sample_id <- overwrite_sample_id
  #   msg("sample_id set to: ", overwrite_sample_id, verbose = verbose)
  # }

  snvs
}


# Based on Strelka README:
# https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/README.md#somatic
calc_strelka_somatic_VAF <- function(tbl) {
  tbl |>
    mutate(
      gt_ref = case_match(
        .data$REF,
        "T" ~ .data$gt_TU,
        "C" ~ .data$gt_CU,
        "G" ~ .data$gt_GU,
        "A" ~ .data$gt_AU
      ),
      gt_alt = case_match(
        .data$ALT,
        "T" ~ .data$gt_TU,
        "C" ~ .data$gt_CU,
        "G" ~ .data$gt_GU,
        "A" ~ .data$gt_AU
      )
    ) |>
    separate("gt_ref", into = c("ref_tier1", "ref_tier2")) |>
    separate("gt_alt", into = c("alt_tier1", "alt_tier2")) |>
    mutate(across(c("ref_tier1", "alt_tier1"), parse_integer)) |>
    mutate(VAF = .data$alt_tier1 / ( .data$alt_tier1 + .data$ref_tier1))
}

