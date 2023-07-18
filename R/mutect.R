## ------------------------------- Export --------------------------------------

#' Read Mutect variants
#'
#' Reads the variant calls from Mutect2 somatic variants caller
#'
#' @param path Can be either:
#'   1) path to a single file
#'   2) vector of file paths, element names will be used as patient IDs
#'   3) directory containing multiple Mutect .vcf files, patient IDs will be
#'      guessed from the file names (should follow convention: <patient_ID>.XXX.XXX.vcf)
#' @param sample_ids Either:
#'   1) "drop_first"
#'   2) "all"
#'   3) ID(s) of the selected tumor
#'      samples. Default: "drop_first"
#' @param PASS_only Keep FILTER == PASS variants only?
#' @param patient_id_pattern If path is a dir only: pattern for str_extract()
#'   that should be used to extract the patient_id from the filenames
#' @param chrom_convention UCSC/NCBI/keep
#' @param extract_VEP_fields If VCF file contains VEP annotations, following
#'   fields will be extracted: Variant_Classification, impact,
#'   gene_symbol and entrez_id (epxerimental, not tested)
#' @param verbose Verbose?
#'
#' @name mutect
NULL



#' @rdname mutect
#' @export
read_mutect_snvs <- function(path,
                             sample_ids = "drop_first",
                             PASS_only = TRUE,
                             patient_id_pattern = "(?<=\\/)[:alnum:]*(?=\\.)",
                             chrom_convention = "UCSC",
                             extract_VEP_fields = FALSE,
                             verbose = TRUE) {
  if (is_single_file(path)) {
    snvs <- read_mutect_vcf(path, sample_ids, PASS_only, extract_VEP_fields, verbose)
  } else if (is_list_of_files(path)) {
    snvs <- map(path, read_mutect_vcf, sample_ids, PASS_only, extract_VEP_fields, verbose) |>
      bind_rows(.id = "patient_id")
  } else if (is_single_dir(path)) {
    files <- list.files(path, full.names = TRUE)
    patient_ids <- files |>
      str_extract(pattern = patient_id_pattern)
    names(files) <- patient_ids
    snvs <- files |>
      map(read_mutect_vcf, sample_ids, PASS_only, extract_VEP_fields, verbose) |>
      bind_rows(.id = "patient_id")
  }

  snvs |>
    use_chrom_naming_convention(chrom_convention)
}


## ------------------------------ Functions ------------------------------------


read_mutect_vcf <- function(file,
                            sample_ids = NULL,
                            PASS_only = TRUE,
                            extract_VEP_fields = FALSE,
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
    separate("gt_AD", into = c("ref_reads", "alt_reads"), convert = TRUE) |>
    mutate(VAF = .data$alt_reads / (.data$alt_reads + .data$ref_reads)) |>
    select(
      sample_id = "Indiv",
      chrom = "CHROM", pos = "POS",
      ref = "REF", alt = "ALT",
      "FILTER",
      "ref_reads", "alt_reads",
      "VAF", AF = "gt_AF", DP = "gt_DP",
      all_of("CSQ")
    )

  if (extract_VEP_fields && "CSQ" %in% names(snvs)) {
    snvs <- extract_VEP_fields(snvs)
  }

  class(snvs) <- c("Mutect_tbl", class(snvs))

  snvs
}



extract_VEP_fields <- function(tbl) {
  tbl |>
    mutate(
       CSQ_split  = str_split(.data$CSQ, pattern = "\\|"),
       Variant_Classification = map_chr(.data$CSQ_split, 2),
       impact = map_chr(.data$CSQ_split, 3),
       gene_symbol = map_chr(.data$CSQ_split, 4),
       entrez_id = map_chr(.data$CSQ_split, 5)
    ) |>
    select(-"CSQ_split")
}
