## ------------------------------- Export --------------------------------------

#' Read Strelka variants
#'
#' Reads the variant calls from Strelka2 Small Variant Caller
#'
#' @param file VCF file to read
#' @param sample_id ID(s) of the tumor sample to use. The last ID will be used
#'   if none provided
#' @param PASS_only Keep FILTER == PASS variants only?
#' @param chrom_convention UCSC/NCBI/keep
#' @param verbose Verbose?
#'
#' @name strelka
NULL



#' @rdname strelka
#' @export
read_strelka_somatic_snvs <- function(file,
                                      sample_id = NULL,
                                      PASS_only = TRUE,
                                      chrom_convention = "UCSC",
                                      verbose = TRUE) {
  tidyvcf <- read_vcf(file, PASS_only, verbose = verbose)

  if (is.null(sample_id)) {
    sample_ids <- unique(tidyvcf$dat$Indiv)
    sample_id <- last(sample_ids)
    msg("Guessing sample_id: ", sample_id, verbose = verbose)
  }

  snvs <- tidyvcf$dat |>
    filter(.data$Indiv %in% sample_id) |>
    calc_strelka_somatic_VAF() |>
    select(
      sample_id = "Indiv",
      chrom = "CHROM", pos = "POS",
      ref = "REF", alt = "ALT",
      ref_reads = "ref_tier1",
      alt_reads = "alt_tier1",
      VAF = "VAF", DP = "DP"
    ) |>
    use_chrom_naming_convention(chrom_convention)

  class(snvs) <- c("cevo_Strelka", class(snvs))
  snvs
}



## ------------------------------ Functions ------------------------------------


# Read VCF into a tidy tibble format
read_vcf <- function(file, PASS_only = TRUE, ..., verbose = FALSE) {
  vcf <- vcfR::read.vcfR(file, verbose = verbose)
  if (PASS_only) {
    vcf <- vcf[vcfR::getFILTER(vcf) %in% c("PASS", NA)]
  }
  vcfR::vcfR2tidy(vcf, single_frame = TRUE, ..., verbose = verbose)
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

