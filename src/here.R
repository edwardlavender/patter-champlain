# here::here() wrappers

here_data_raw_mf <- function(...) {
  here_data_raw("mfutia", "model_comparison", ...)
}

here_fig_sim <- function(...) {
  here_fig("sim", ...)
}

here_fig_real <- function(...) {
  here_fig("real", ...)
}

here_input <- function(...) {
  here::here("data", "input", ...)
}

here_input_sim <- function(...) {
  here_input("sim", ...)
}

here_input_real <- function(...) {
  here_input("real", ...)
}

here_output <- function(...) {
  here::here("data", "output", ...)
}

here_output_sim <- function(...) {
  here_output("sim", ...)
}

here_output_sim_main <- function(...) {
  here_output("sim", "main", ...)
}

here_output_real <- function(...) {
  here_output("real", ...)
}

here_output_real_main <- function(...) {
  here_output("real", "main", ...)
}
