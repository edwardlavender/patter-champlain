# Discretise coordinates for kde estimation
# * This code is based on .map_coord() and .map_marks()
# * We return a time series with the distribution for each time step
kde_coord <- function(.map, .coord) {
  .coord |>
    lazy_dt() |>
    # Discretise coordinates for speed
    mutate(id = terra::cellFromXY(.map, cbind(.data$x, .data$y)),
           x = terra::xFromCell(.map, .data$id),
           y = terra::yFromCell(.map, .data$id)) |>
    # Assign equal weights (marks)
    group_by(.data$timestep) |>
    mutate(mark = 1 / n()) |>
    ungroup() |>
    # Calculate the total weight of each location within time steps
    group_by(.data$timestep, .data$id) |>
    mutate(mark = sum(.data$mark)) |>
    slice(1L) |>
    ungroup() |>
    as.data.table()
}

# Compute area spanned by 95 % of particles via kde
# * Use binned = TRUE for speed
hr_area <- function(x, y) {
  fhat <- ks::kde(cbind(x, y), 
                  binned = TRUE, 
                  compute.cont = FALSE, 
                  approx.cont = FALSE)
  # plot(fhat, cont = 95)
  ks::contourSizes(fhat, cont = 95, approx = TRUE)
}
