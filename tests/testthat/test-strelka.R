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


test_that("read_strelka_somatic_snvs() works with single file", {
  path <- test_path("Strelka", "S1.somatic.snvs.vcf.gz")

  res <- read_strelka_somatic_snvs(path, verbose = FALSE)
  exp_chrs <- str_c("chr", c(1, 2, 3, 4, 5, 6, 7, "X", "Y"))
  exp_colnames <- c("sample_id", "chrom", "pos", "ref", "alt", "ref_reads", "alt_reads", "VAF", "DP")

  expect_s3_class(res, "cevo_Strelka")
  expect_s3_class(res, "tbl")
  expect_named(res, exp_colnames)
  expect_equal(res$chrom, exp_chrs)
})


test_that("read_strelka_somatic_snvs() works with list of files", {
  path <-c(
    S1 = test_path("Strelka", "S1.somatic.snvs.vcf.gz"),
    S2 = test_path("Strelka", "S2.somatic.snvs.vcf.gz")
  )

  res <- read_strelka_somatic_snvs(path, verbose = FALSE)
  exp_chrs <- str_c("chr", c(1, 2, 3, 4, 5, 6, 7, "X", "Y"))
  exp_colnames <- c("patient_id", "sample_id", "chrom", "pos", "ref", "alt", "ref_reads", "alt_reads", "VAF", "DP")

  expect_s3_class(res, "cevo_Strelka")
  expect_s3_class(res, "tbl")
  expect_named(res, exp_colnames)
  expect_equal(res$chrom, rep(exp_chrs, times = 2))
})


test_that("read_strelka_somatic_snvs() works with dir", {
  path <- test_path("Strelka")

  res <- read_strelka_somatic_snvs(path, verbose = FALSE)
  exp_chrs <- str_c("chr", c(1, 2, 3, 4, 5, 6, 7, "X", "Y"))
  exp_colnames <- c("patient_id", "sample_id", "chrom", "pos", "ref", "alt", "ref_reads", "alt_reads", "VAF", "DP")

  expect_s3_class(res, "cevo_Strelka")
  expect_s3_class(res, "tbl")
  expect_named(res, exp_colnames)
  expect_equal(res$chrom, rep(exp_chrs, times = 2))
})

