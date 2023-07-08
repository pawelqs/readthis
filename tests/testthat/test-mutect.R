
test_that("read_mutect_snvs() works with single file", {
  path <- test_path("Mutect", "S1.Mutect2.filter.pass.phased.annot.vcf")

  res <- read_mutect_snvs(path, verbose = FALSE)
  exp_chrs <- str_c("chr", c(1, 1, 2, 2, "X", "X", "Y", "Y"))
  exp_colnames <- c("sample_id", "chrom", "pos", "ref", "alt", "FILTER", "ref_reads", "alt_reads", "VAF", "AF", "DP", "CSQ")

  expect_s3_class(res, "Mutect_tbl")
  expect_named(res, exp_colnames)
  expect_equal(res$chrom, exp_chrs)
  expect_equal(nrow(res), 8)
})



test_that("read_mutect_snvs() works with directory", {
  path <- test_path("Mutect")

  res <- read_mutect_snvs(path, verbose = FALSE)
  exp_chrs <- str_c("chr", rep(c(1, 1, 2, 2, "X", "X", "Y", "Y"), 2))
  exp_colnames <- c(
    "patient_id", "sample_id", "chrom", "pos", "ref", "alt", "FILTER",
    "ref_reads", "alt_reads", "VAF", "AF", "DP", "CSQ")

  expect_s3_class(res, "Mutect_tbl")
  expect_named(res, exp_colnames)
  expect_equal(res$chrom, exp_chrs)
  expect_equal(nrow(res), 16)
})


test_that("read_mutect_snvs() works with list of files", {
  path <- c(
    test_path("Mutect", "S1.Mutect2.filter.pass.phased.annot.vcf"),
    test_path("Mutect", "S2.Mutect2.filter.pass.phased.annot.vcf")
  )

  res <- read_mutect_snvs(path, verbose = FALSE)
  exp_chrs <- str_c("chr", rep(c(1, 1, 2, 2, "X", "X", "Y", "Y"), 2))
  exp_colnames <- c(
    "patient_id", "sample_id", "chrom", "pos", "ref", "alt", "FILTER",
    "ref_reads", "alt_reads", "VAF", "AF", "DP", "CSQ")

  expect_s3_class(res, "Mutect_tbl")
  expect_named(res, exp_colnames)
  expect_equal(res$chrom, exp_chrs)
  expect_equal(nrow(res), 16)
})


test_that("read_mutect_snvs() extracts VEP fields", {
  path <- test_path("Mutect", "S1.Mutect2.filter.pass.phased.annot.vcf")

  res <- read_mutect_snvs(path, extract_VEP_fields = TRUE, verbose = FALSE)
  exp_colnames <- c(
    "sample_id", "chrom", "pos", "ref", "alt", "FILTER",
    "ref_reads", "alt_reads", "VAF", "AF", "DP", "CSQ",
    "Variant_Classification", "impact", "gene_symbol", "entrez_id"
    )

  expect_named(res, exp_colnames)
  expect_equal(nrow(res), 8)
})

