switch_function <- function(analysis = c("sim", "real"), fun_sim, fun_real) {
  analysis <- match.arg(analysis)
  switch(
    analysis,
    sim  = fun_sim,
    real = fun_real,
    stop("`analysis` must be 'sim' or 'real'")
  )
}

switch_constructor_ac_analysis <- function(analysis) {
  switch_function(analysis, constructor_ac_sim, constructor_ac_real)
}

switch_here_input_analysis <- function(analysis) {
  switch_function(analysis, here_input_sim, here_input_real)
}

switch_here_output_analysis_main <- function(analysis) {
  switch_function(analysis, here_output_sim_main, here_output_real_main)
}