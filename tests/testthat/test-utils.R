test_that("use_UCSC_chrom_names() adds chr prefix if absent", {
  tbl <- tibble(chrom = 1:3)
  res <- use_UCSC_chrom_names(tbl)
  expect_equal(res$chrom, c("chr1", "chr2", "chr3"))
})


test_that("use_UCSC_chrom_names() does not add chr prefix if present", {
  tbl <- tibble(chrom = str_c("chr", 1:3))
  res <- use_UCSC_chrom_names(tbl)
  expect_equal(res$chrom, c("chr1", "chr2", "chr3"))
})


test_that("use_NCBI_chrom_names() removes chr prefix if present", {
  tbl <- tibble(chrom = str_c("chr", 1:3))
  res <- use_NCBI_chrom_names(tbl)
  expect_equal(res$chrom, as.character(1:3))
})
