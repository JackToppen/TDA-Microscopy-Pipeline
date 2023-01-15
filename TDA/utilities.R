

get_csvs <- function(path) {
    # get list of discretized CSVs recursively
    files <- list.files(path, pattern="\\.csv$", include.dirs=TRUE, recursive=TRUE)
}


concat_path <- function(path, filename) {
  # make sure path has file separator at the end
  n <- nchar(path)
  if (substr(path, n, n) != .Platform$file.sep) {
    path <- paste0(path, .Platform$file.sep)
  }

  # return correct joined path to file
  return(paste0(path, filename))
}

