if (!patter:::os_linux() | (patter:::os_linux() & !patter:::julia_session())) {
  
  # Compute area spanned by 95 % of particles via 2D histogram
  # (Computation via ks::kde() and ks::contourSizes() is too slow)
  particle_hr <- function(.map, .coord, .percentage = 0.75, .summarise = TRUE) {
    
    # Compute area of grid cell (assuming UTM grid)
    A <- prod(terra::res(.map))
    
    # Compute area spanned by 95 % of particles for each t
    areas <- 
      .coord |>
      lazy_dt() |>
      select(timestep, x, y) |>
      # filter(timestep %in% seq(1, max(timestep), by = 100)) |> 
      # Discretise coordinates
      mutate(id = terra::cellFromXY(.map, cbind(.data$x, .data$y)),
             x = terra::xFromCell(.map, .data$id),
             y = terra::yFromCell(.map, .data$id)) |>
      # Assign equal weights (marks)
      group_by(.data$timestep) |>
      mutate(mark = 1 / n()) |>
      ungroup() |>
      # Calculate the total weight of each location within time steps (2D histogram)
      group_by(.data$timestep, .data$id) |>
      summarise(mark = sum(.data$mark)) |>
      ungroup() |>
      # Calculate for each timestep the area containing 95 % of probability mass
      # * Sum weights
      # * Identify the number of cells required to put us above 95 % probability mass * A
      # * This gives the area in units of A (m^2)
      group_by(.data$timestep) |> 
      arrange(desc(mark), .by_group = TRUE) |> 
      mutate(cmark = cumsum(mark)) |>
      summarise(area = which(cmark >= .percentage)[1] * A) |>
      ungroup() 
    
    # Optionally summarise areas
    if (!.summarise) {
      return(as.data.table(areas))
    } else {
      # Compute mean area over all time steps
      areas |> 
        summarise(area_mean = mean(area)) |>
        ungroup() |> 
        pull(area_mean)
    }
    
  }
  
}