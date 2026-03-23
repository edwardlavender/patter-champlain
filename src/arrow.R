write_feather_compressed <- function(x, sink, ...) {
  arrow::write_feather(x, sink, compression = "zstd", compression_level = 9, ...)
}