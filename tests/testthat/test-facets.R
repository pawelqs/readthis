test_that("read_facets_cnvs() works with single file", {
  path <- test_path("FACETS", "S1.csv")
  res <- read_facets_cnvs(path)
  expect_s3_class(res, "cevo_FACETS")
  expect_equal(dim(res), c(64, 18))
  expect_equal(unique(res$sample_id), path)
})


test_that("read_facets_cnvs() works with dir path file", {
  path <- test_path("FACETS")
  res <- read_facets_cnvs(path)
  expect_s3_class(res, "cevo_FACETS")
  expect_equal(dim(res), c(128, 18))
  expect_equal(unique(res$sample_id), c("S1", "S2"))
})


test_that("read_facets_cnvs() works with list of files", {
  path <- test_path("FACETS")
  path <- list.files(path, full.names = TRUE) |>
    set_names(c("S1", "S2"))
  res <- read_facets_cnvs(path)
  expect_s3_class(res, "cevo_FACETS")
  expect_equal(dim(res), c(128, 18))
  expect_equal(unique(res$sample_id), c("S1", "S2"))
})
