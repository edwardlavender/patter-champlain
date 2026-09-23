switch_here_input_analysis_subanalysis <-
  function(analysis = c("sim", "real", "validation"), subanalysis = "main") {
  analysis    <- match.arg(analysis)
  subanalysis <- match.arg(subanalysis)
  switch(analysis, 
         sim = function(...) here_input(analysis, subanalysis, ...), 
         real = function(...) here_input(analysis, subanalysis, ...), 
         validation = function(...) here_input(analysis, subanalysis, ...))
}

switch_here_output_analysis_subanalysis <- 
  function(analysis = c("sim", "real", "validation"), subanalysis = "main") {
  analysis    <- match.arg(analysis)
  subanalysis <- match.arg(subanalysis)
  switch(analysis, 
         sim = function(...) here_output(analysis, subanalysis, ...), 
         real = function(...) here_output(analysis, subanalysis, ...),
         validation = function(...) here_output(analysis, subanalysis, ...))
}

switch_here_fig_analysis_subanalysis <- 
  function(analysis = c("sim", "real", "validation"), subanalysis = "main") {
  analysis    <- match.arg(analysis)
  subanalysis <- match.arg(subanalysis)
  switch(analysis, 
         sim = function(...) here_fig(analysis, subanalysis, ...), 
         real = function(...) here_fig(analysis, subanalysis, ...), 
         validation = function(...) here_fig(analysis, subanalysis, ...))
}
