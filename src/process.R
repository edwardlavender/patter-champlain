# Select detection data that meet minimum quality criteria

# cf. patter-flapper criteria:
# - individuals must be detected in at least two different weeks 
# - on a total of seven days in a given month
# - (NB: this study included depth observations)

# patter-champlain criteria development:
# - Trout mobility is reasonable relative to the size of the study area
# - It would take a trout one day at top speed to cross the study area
# - Even at slower speeds, locations are uncertain after relatively short times
# - Without regular detections, locations are uncertain

# map  <- terra::rast(here_input("map.tif"))
# pars <- qs::qread(here_input("pars-patter.qs"))
# ydim <- terra::ext(map)[4] - terra::ext(map)[3] # 174928.5 m
# ydim / (24 * 60/2 * pars$mobility[1])           # 1 day to cross area

filter_detections <- function(detections, plot = TRUE) {
  
  # Check inputs
  proj.build::check_inherits(detections, "data.table")
  proj.build::check_names(detections, "unit_id")
  detections <- copy(detections)
  
  # Compute the proportion of days per time block with detections
  unit_id   <- NULL
  durations <- 
    detections |> 
    group_by(unit_id) |> 
    summarise(
      duration = as.numeric(difftime(max(timestamp), min(timestamp), units = "days")), 
      ndays = length(unique(lubridate::floor_date(timestamp, "days"))),
      pdays = ndays / duration
    ) |> 
    as.data.table()
  
  # (Interactive) Examine the proportion of days with detections
  if (plot) {
    # quantile(durations$pdays)
    plot(ecdf(durations$pdays), xlim = c(0, 1))
  }

  # Select individual/month combinations with detections X % of days
  # * Note that unit_id is _not_ redefined
  unit_ids   <- durations$unit_id[durations$pdays >= 0.75]
  n0         <- length(unique(detections$unit_id))
  n1         <- length(unique(unit_ids))
  detections <- detections[unit_id %in% unit_ids, ]
  cat(glue::glue("{n1} / {n0} individual/time block(s) ('unit_id(s)') retained.
                 
                 "))
  
  # Return filtered detection data.table
  detections
  
}