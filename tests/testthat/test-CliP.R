

test_that("recognize_clip_sample_dirs() recognizes samples in CliP results", {
  dir <- test_path("CliP")
  res <- recognize_clip_sample_dirs(dir)
  expected <- tibble(
    sample_id = str_c("Sample", LETTERS[1:2]),
    sample_dir = file.path(test_path("CliP"), sample_id, "results")
  )
  expect_identical(res, expected)
})


test_that("list_clip_files() finds all lambda results", {
  dir <- test_path("CliP")
  sample_dirs <- recognize_clip_sample_dirs(dir)
  res <- list_clip_files(sample_dirs, best_only = FALSE)
  # write_tsv(res, "clip_files.tsv")
  expected <- read_tsv(test_path("clip_files.tsv"), col_types = "ccccl")
  expect_identical(res, expected)
})


test_that("list_clip_files() finds best lambda results", {
  dir <- test_path("CliP")
  sample_dirs <- recognize_clip_sample_dirs(dir)
  res <- list_clip_files(sample_dirs)
  expected <- read_tsv(test_path("clip_files.tsv"), col_types = "ccccl") |>
    filter(best_lambda)
  expect_identical(res, expected)
})


