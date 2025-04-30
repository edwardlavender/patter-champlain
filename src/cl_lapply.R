# Iteratively read iteration[[.file]] and apply a function(.sim, file_contents)
cl_lapply_iteration_file <- function(.iteration, .file, .fun, .read = qs::qread, ...) {
  
  # Checks
  iteration <- copy(.iteration)
  stopifnot(rlang::has_name(iteration, .file))
  iteration[, dot_file := iteration[[.file]]]
  
  #  Focus on iteration rows for which dot_file exists
  n0        <- nrow(iteration)
  iteration <- iteration[file.exists(dot_file), ]
  n1        <- nrow(iteration)
  if (n0 != n1) {
    message(glue::glue("Focusing on {n1}/{n0} iteration rows."))
  }
  stopifnot(nrow(n1) > 0L)
  
  # Iterate over iteration rows, read the file and apply the function
  cl_lapply(split(iteration, seq_len(nrow(iteration))), function(.sim) {
    output <- .read(.sim$dot_file)
    .fun(.sim, output)
  }, ...)
  
}
