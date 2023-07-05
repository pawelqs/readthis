test_that("read_vcf() with PASS_only = FALSE works correctly", {
  file <- test_path("Strelka", "S1.somatic.snvs.vcf")
  res <- read_vcf(file, PASS_only = FALSE)
  expect_s3_class(res$dat, "tbl")
  expect_s3_class(res$meta, "tbl")
  expect_equal(res$dat$POS, rep(1:12, each = 2))
})


test_that("read_vcf() with PASS_only = TRUE works correctly", {
  file <- test_path("Strelka", "S1.somatic.snvs.vcf")
  res <- read_vcf(file, PASS_only = TRUE)
  expect_s3_class(res$dat, "tbl")
  expect_s3_class(res$meta, "tbl")
  expect_equal(res$dat$POS, rep(c(1, 3, 5, 7, 8, 9, 10, 11, 12), each = 2))
})



test_that("read_vcf() works with gzipped VCFs", {
  file <- test_path("Strelka", "S1.somatic.snvs.vcf.gz")
  res <- read_vcf(file, PASS_only = TRUE)
  expect_s3_class(res$dat, "tbl")
  expect_s3_class(res$meta, "tbl")
  expect_equal(res$dat$POS, rep(c(1, 3, 5, 7, 8, 9, 10, 11, 12), each = 2))
})


test_that("calc_strelka_somatic_VAF() works correctly", {
  tbl <- tibble(
    REF = c("T", "A"),
    ALT = c("G", "C"),
    gt_TU = c("90,99", "0,0"),
    gt_CU = c("0,0", "35,40"),
    gt_GU = c("10,25", "0,0"),
    gt_AU = c("0,0", "65,100")
  )
  expected_VAFs <- c(0.1, 0.35)
  res <- calc_strelka_somatic_VAF(tbl)
  expect_equal(res$VAF, expected_VAFs)
})


test_that("read_strelka_somatic_snvs() works correctly", {
  file <- test_path("Strelka", "S1.somatic.snvs.vcf")

  res <- read_strelka_somatic_snvs(file, verbose = FALSE)
  exp_chrs <- str_c("chr", c(1, 2, 3, 4, 5, 6, 7, "X", "Y"))
  exp_colnames <- c("sample_id", "chrom", "pos", "ref", "alt", "ref_reads", "alt_reads", "VAF", "DP")

  expect_s3_class(res, "tbl")
  expect_named(res, exp_colnames)
  expect_equal(res$chrom, exp_chrs)
})

