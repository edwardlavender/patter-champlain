switch_here_input_analysis_subanalysis <- function(analysis = c("sim", "real"), subanalysis = "main") {
  analysis    <- match.arg(analysis)
  subanalysis <- match.arg(subanalysis)
  switch(analysis, 
         sim = function(...) here_input(analysis, subanalysis, ...), 
         real = function(...) here_input(analysis, subanalysis, ...))
}

switch_here_fig_analysis_subanalysis <- function(analysis = c("sim", "real"), subanalysis = "main") {
  analysis    <- match.arg(analysis)
  subanalysis <- match.arg(subanalysis)
  switch(analysis, 
         sim = function(...) here_fig(analysis, subanalysis, ...), 
         real = function(...) here_fig(analysis, subanalysis, ...))
}
