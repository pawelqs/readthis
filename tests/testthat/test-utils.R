test_that("get_files() works", {
  tbl <- test_path("ASCAT")

  csv_files <- get_files(tbl, ".csv", "(?<=\\/)[:alnum:]*(?=\\.)")
  expected_csv_files <- c(S1 = "ASCAT/S1.csv", S2 = "ASCAT/S2.csv")
  expect_equal(csv_files, expected_csv_files)

  stat_files <- get_files(tbl, "samplestatistics", "(?<=\\/)[:alnum:]*(?=\\.)")
  expected_stat_files <- c(S1 = "ASCAT/S1.samplestatistics.txt", S2 = "ASCAT/S2.samplestatistics.txt")
  expect_equal(stat_files, expected_stat_files)
})


## --------------------- chrom names conventions -------------------------------

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


## --------------------- analyze path -------------------------------


test_that("directory paths are recognized correctly", {
  path <- test_path("FACETS")
  expect_equal(is_single_dir(path), TRUE)
  expect_equal(is_single_file(path), FALSE)
  expect_equal(is_list_of_files(path), FALSE)
})


test_that("file paths are recognized correctly", {
  path <- test_path("FACETS", "S1.csv")
  expect_equal(is_single_dir(path), FALSE)
  expect_equal(is_single_file(path), TRUE)
  expect_equal(is_list_of_files(path), FALSE)
})


test_that("many paths are recognized correctly", {
  dir <- test_path("FACETS")
  paths <- list.files(dir, full.names = TRUE)
  expect_equal(is_single_dir(paths), FALSE)
  expect_equal(is_single_file(paths), FALSE)
  expect_equal(is_list_of_files(paths), TRUE)
})
