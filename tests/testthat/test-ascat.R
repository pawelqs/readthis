test_that("read_ascat_csv() works", {
  path <- test_path("ASCAT", "S1.csv")
  res <- read_ascat_cnvs(path)
  expect_s3_class(res, "tbl")
  expect_equal(dim(res), c(10, 7))
})



test_that("read_ascat_samplestatistics() works", {
  path <- test_path("ASCAT", "S1.samplestatistics.txt")
  res <- read_ascat_samplestatistics(path)
  expected_names <- c("normal_contamination", "ploidy", "rho", "psi", "goodness_of_fit", "gender_chr", "gender_chr_found")
  expect_s3_class(res, "tbl")
  expect_named(res, expected_names)
})



test_that("read_ascat_files() works with cnvs only", {
  path <- test_path("ASCAT", "S1.csv")
  res <- read_ascat_files(path)
  expect_s3_class(res, "cevo_ASCAT")
  expect_equal(dim(res$cnvs), c(10, 8))
  expect_equal(unique(res$cnvs$sample_id), c("ASCAT/S1.csv"))
})



test_that("read_ascat_files() works with cnvs and sample statistics", {
  path <- test_path("ASCAT", "S1.csv")
  sample_statistics <- test_path("ASCAT", "S1.samplestatistics.txt")
  res <- read_ascat_files(path, sample_statistics)
  expect_s3_class(res, "cevo_ASCAT")
  expect_equal(dim(res$cnvs), c(10, 8))
  expect_equal(dim(res$sample_statistics), c(1, 8))
  expect_equal(unique(res$cnvs$sample_id), c("ASCAT/S1.csv"))
  expect_equal(unique(res$sample_statistics$sample_id), c("ASCAT/S1.csv"))
})



test_that("read_ascat_files() works with tibble of files", {
  dir <- test_path("ASCAT")
  files <- list.files(dir, full.names = TRUE)
  path <- tibble(
    sample_id = c("S1", "S2"),
    csv = grep("csv", files, value = TRUE),
    sample_statistics = grep("samplestatistics", files, value = TRUE)
  )
  res <- read_ascat_files(path)
  expect_s3_class(res, "cevo_ASCAT")
  expect_equal(dim(res$cnvs), c(20, 8))
  expect_equal(dim(res$sample_statistics), c(2, 8))
  expect_equal(unique(res$sample_statistics$sample_id), c("S1", "S2"))
})



test_that("read_ascat_files() works with tibble of files with NAs", {
  dir <- test_path("ASCAT")
  files <- list.files(dir, full.names = TRUE)
  path <- tibble(
    sample_id = c("S1", "S2"),
    csv = grep("csv", files, value = TRUE),
    sample_statistics = c(grep("samplestatistics", files, value = TRUE)[1], NA_character_)
  )
  res <- read_ascat_files(path)
  expect_s3_class(res, "cevo_ASCAT")
  expect_equal(dim(res$cnvs), c(20, 8))
  expect_equal(dim(res$sample_statistics), c(1, 8))
  expect_equal(unique(res$cnvs$sample_id), c("S1", "S2"))
  expect_equal(unique(res$sample_statistics$sample_id), "S1")
})



test_that("read_ascat_files() works with directory", {
  path <- test_path("ASCAT")
  res <- read_ascat_files(path)
  expect_s3_class(res, "cevo_ASCAT")
  expect_equal(dim(res$cnvs), c(20, 8))
  expect_equal(dim(res$sample_statistics), c(2, 8))
  expect_equal(unique(res$sample_statistics$sample_id), c("S1", "S2"))
})
