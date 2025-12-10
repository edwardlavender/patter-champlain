# Check if spatial dependencies are loaded (linux helper)
expect_no_geospatial <- function() {
  if (any(c("terra", "sf") %in% loadedNamespaces())) {
    warning("Geospatial dependencies loaded!")
  }
  invisible(NULL)
}

# Read a map with terra, off Linux
rast_map <- function(file) {
  if (!patter:::os_linux()) {
    map <- terra::rast(file)
  } else {
    map <- NULL
  }
  map
}

if (!patter:::os_linux() | (patter:::os_linux() & !patter:::julia_session())) {
  
  # Convert shapefile to SpatRaster
  as_SpatRaster <- function(.x, .simplify = NULL, .utm = NULL,
                            .res, .field, .touches = TRUE, 
                            .plot = FALSE, ...) {
    
    # Translate sf objects to SpatVectors
    if (inherits(.x, "sf")) {
      .x <- terra::vect(.x)
    }
    
    # (optional) Simplify SpatVector 
    # > This improves speed
    if (!is.null(.simplify)) {
      .x <- terra::simplifyGeom(.x, tolerance = .simplify)
    }
    
    # Translate to UTM
    if (!is.null(.utm)) {
      .x <- terra::project(.x, .utm)
    }
    
    # Define map_value
    .x$map_value <- .x[[.field]]
    
    # Define blank raster for rasterization 
    r <- terra::rast(terra::ext(.x), 
                     crs = terra::crs(.x),
                     res = .res)
    
    # Rasterise the SpatVector 
    map <- terra::rasterize(.x, r, field = "map_value", touches = .touches, ...)
    
    # (optional) Plot 
    if (.plot) {
      terra::plot(map)
      terra::lines(.x)
    }
    
    # Return list
    list(SpatVector = .x, SpatRaster = map)
    
  }
  
  # Read all
  readAll <- function(x) {
    x |> 
      terra::wrap() |> 
      terra::unwrap()
  }
  
}

# Compute the smallest absolute rotation between two angles
# Source: wahoo-flapper
# abs_angle_difference <- function(a1, a2) {
#   angle_delta <- abs((a1 - a2) %% (2 * pi))
#   pmin(angle_delta, 2 * pi - angle_delta)
# }

signed_angle_difference <- function(a1, a2) {
  ((a1 - a2 + pi) %% (2*pi)) - pi
}

# Compute turning angles for the VPS data
# Modified from wahoo-flapper
if (!patter:::os_linux() | (patter:::os_linux() & !patter:::julia_session())) {
  
  calc_angle_vps <- function(.data) {
    
    # Checks
    check_names(.data, c("Longitude", "Latitude"))
    .data <- as.data.table(.data)
    
    # Define UTM coordinates from 'Longitude' and 'Latitude'
    xy <- 
      cbind(.data$Longitude, .data$Latitude) |>
      terra::vect(crs = "WGS84") |> 
      terra::project("EPSG:3175") |> 
      terra::geom() |> 
      as.data.table()
    .data[, x := xy$x]
    .data[, y := xy$y]
    
    # Compute turning angles
    .data |> 
      lazy_dt(immutable = FALSE) |>
      mutate(dx = x - lag(x), 
             dy = y - lag(y), 
             heading = atan2(dy, dx), 
             heading = ((pi / 2 - heading + pi) %% (2 * pi)) - pi,
             turning_angle = signed_angle_difference(heading, lag(heading))) |> 
      pull(turning_angle)
    
  }

}
