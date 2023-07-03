## ------------------------------- Export --------------------------------------

#' Read CliP results
#'
#' [CliP](https://github.com/wwylab/CliP) is an algorithm for clonal structure
#' identification through penalizing pairwise differences by Wenyi Wang Lab
#' at MD Anderson Cancer Center in Houston.
#'
#' @param dir directory in which to seek for the CliP results
#'
#' @examples
#' dir <- system.file("CliP", package = "readthis")
#' read_clip_best_lambda(dir)
#'
#' @name clip
NULL



#' @describeIn clip Read CliP results for all lambdas
#' @export
read_clip_all <- function(dir) {
  sample_dirs <- recognize_clip_sample_dirs(dir)
  files <- list_clip_files(sample_dirs, best_only = FALSE)
  read_clip_files(files)
}



#' @describeIn clip Read CliP results for all lambdas, in the wider format
#' @export
read_clip_all_wide <- function(dir) {
  sample_dirs <- recognize_clip_sample_dirs(dir)
  files <- list_clip_files(sample_dirs, best_only = FALSE)

  res <- read_clip_files(files)
  res$mutation_assignments <- res$mutation_assignments |>
    select(-"best_lambda") |>
    pivot_wider(names_from = "lambda", values_from = "cluster_index")

  res
}



#' @describeIn clip Read CliP results for best lambda only
#' @export
read_clip_best_lambda <- function(dir) {
  sample_dirs <- recognize_clip_sample_dirs(dir)
  files <- list_clip_files(sample_dirs)
  read_clip_files(files)
}


## ------------------------------- Functions --------------------------------------

recognize_clip_sample_dirs <- function(dir) {
  dirs <- list.dirs(dir) |>
    str_replace_all("//", "/")
  best_lambda_dirs <- grep("Best_lambda", dirs, value = TRUE)

  path_divided <- best_lambda_dirs |>
    str_split("/") |>
    map(\(elements) set_names(elements, str_c("element", seq_along(elements)))) |>
    map(as_tibble_row) |>
    bind_rows()

  sample_names_mask <- path_divided |>
    map_lgl(~ length(unique(.x)) > 1)
  sample_ids <- unlist(path_divided[sample_names_mask])
  names(sample_ids) <- NULL

  tibble(
    sample_id = sample_ids, # if there is one sample only, this column does not exist
    sample_dir = str_replace(best_lambda_dirs, "/Best_lambda", "")
  )
}



list_clip_files <- function(sample_dirs, best_only = TRUE) {
  single_sample <- !"sample_id" %in% names(sample_dirs)
  if (single_sample) {
    sample_dirs$sample_id <- "fake_id"
  }

  files <- sample_dirs |>
    rowwise("sample_id") |>
    reframe(
      path = list.files(.data$sample_dir, recursive = TRUE, full.names = TRUE)
    ) |>
    mutate(
      best_lambda = str_detect(.data$path, "Best_lambda"),
      lambda = str_extract(.data$path, "(?<=lam)[0-9.]*(?=.txt)"),
      file_type = case_when(
        str_detect(.data$path, "mutation_assignments") ~ "mutation_assignments",
        str_detect(.data$path, "subclonal_structure") ~ "subclonal_structure",
        TRUE ~ NA_character_
      )
    ) |>
    pivot_wider(names_from = "file_type", values_from = "path")

  best_lambdas <- files |>
    filter(.data$best_lambda) |>
    select("sample_id", "lambda", "best_lambda")

  files <- files |>
    filter(!.data$best_lambda) |>
    select(-"best_lambda") |>
    left_join(best_lambdas, by = c("sample_id", "lambda")) |>
    replace_na(list(best_lambda = FALSE)) |>
    arrange("sample_id", "lambda")

  if (single_sample) {
    files$sample_id <- NULL
  }

  if (best_only) {
    filter(files, .data$best_lambda)
  } else {
    files
  }
}



read_clip_files <- function(files) {
  mutation_assignments <- files |>
    rowwise(any_of("sample_id"), "lambda", "best_lambda") |>
    reframe(read_tsv(mutation_assignments, show_col_types = FALSE)) |>
    rename(chrom = "chromosome_index", pos = "position") |>
    mutate(chrom = str_c("chr", .data$chrom)) |>
    select(any_of("sample_id"), "chrom", "pos", "cluster_index", "lambda", "best_lambda")

  subclonal_structure <- files |>
    rowwise(any_of("sample_id"), "lambda", "best_lambda") |>
    reframe(read_tsv(subclonal_structure, show_col_types = FALSE)) |>
    select(any_of("sample_id"), "cluster_index", "num_SNV", "cellular_prevalence", "lambda", "best_lambda")

  lst(mutation_assignments, subclonal_structure)
}

