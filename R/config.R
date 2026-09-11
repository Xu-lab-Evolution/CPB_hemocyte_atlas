# Shared paths and helper functions for the CPB hemocyte scRNA-seq workflow.

get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (!length(file_arg)) return(NULL)
  normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
}

project_from_env <- Sys.getenv("CPB_PROJECT_DIR", unset = "")
if (nzchar(project_from_env)) {
  PROJECT_DIR <- normalizePath(project_from_env, mustWork = TRUE)
} else {
  script_path <- get_script_path()
  PROJECT_DIR <- if (!is.null(script_path) && basename(script_path) == "run_analysis.R") {
    dirname(script_path)
  } else {
    normalizePath(getwd(), mustWork = TRUE)
  }
}

INPUT_DIR <- normalizePath(
  Sys.getenv("CPB_INPUT_DIR", unset = PROJECT_DIR),
  mustWork = TRUE
)
OUTPUT_DIR <- Sys.getenv(
  "CPB_OUTPUT_DIR",
  unset = file.path(PROJECT_DIR, "results")
)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUTPUT_DIR <- normalizePath(OUTPUT_DIR, mustWork = TRUE)

WORKDIR <- PROJECT_DIR
OUT <- OUTPUT_DIR
SEED <- 1234L

input_file <- function(name, directory = FALSE) {
  path <- file.path(INPUT_DIR, name)
  exists <- if (directory) dir.exists(path) else file.exists(path)
  if (!exists) {
    type <- if (directory) "directory" else "file"
    stop("Missing input ", type, ": ", path, call. = FALSE)
  }
  normalizePath(path, mustWork = TRUE)
}

require_files <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop("Missing required files: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  invisible(paths)
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(path)
}
