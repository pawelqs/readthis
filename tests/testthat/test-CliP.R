
# # Prepare clip files
# clip_files <- list.files(test_path("CliP"), recursive = TRUE, full.names = TRUE)
# mut_assignments <- grep("mutation_assignments", clip_files, value = TRUE)
# map(
#   mut_assignments,
#   function(file) {
#     file |>
#       read_tsv() |>
#       mutate(position = row_number()) |>
#       write_tsv(file)
#   }
# )



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


test_that("read_clip_best_lambda() reads best lambda results correctly", {
  dir <- test_path("CliP")
  res <- read_clip_all(dir)

  expect_s3_class(res$mutation_assignments, "tbl")
  expect_named(
    res$mutation_assignments,
    c("sample_id", "chrom", "pos", "cluster_index", "lambda", "best_lambda")
  )
  expect_equal(nrow(res$mutation_assignments), 8497 * 4)

  expect_s3_class(res$subclonal_structure, "tbl")
  expect_named(
    res$subclonal_structure,
    c("sample_id", "cluster_index", "num_SNV", "cellular_prevalence", "lambda", "best_lambda")
  )
  expect_equal(nrow(res$subclonal_structure), 35)
})


test_that("read_clip_all_wide() pivots mutation assignments correctly", {
  dir <- test_path("CliP")
  res <- read_clip_all_wide(dir)

  expect_s3_class(res$mutation_assignments, "tbl")
  expect_named(
    res$mutation_assignments,
    c("sample_id", "chrom", "pos", "0.01", "0.03", "0.1", "0.2")
  )
  expect_equal(nrow(res$mutation_assignments), 8497)
})
