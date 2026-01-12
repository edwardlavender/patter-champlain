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

filter_detections <- function(.detections) {
  
  # Check inputs
  proj.build::check_inherits(.detections, "data.table")
  detections <- copy(.detections)
  
  # Define {individual_id}/{time_id} blocks
  unit_id_tmp <- NULL
  detections[, unit_id_tmp := paste(individual_id, time_id)]
  
  # Approach (1): Filter based on maximum gap duration
  # * Compute maximum gap along start/end of time series between sequential detections
  maxgaps <-
    detections |>
    group_by(unit_id_tmp) |>
    # Define timeline start + end and times of detections 
    mutate(timeline = list(sort(unique(
      c(timestamp, 
        min(time_id),
        min(time_id) + lubridate::as.period("1 month") - 60 * 2))))) |>
    # Define duration of max gap in detections
    summarise(max_gap = as.numeric(max(diff(timeline[[1]], units = "days")))) |>
    ungroup() |> 
    as.data.table()
  
  # Approach (2): Filter based on proportion of days per time block with detections
  # * Compute proportion of days per time block with detections 
  pdays <-
    detections |>
    group_by(unit_id_tmp) |>
    summarise(
      duration = as.numeric(difftime(lubridate::ceiling_date(max(timestamp), "months"),
                                     lubridate::floor_date(min(timestamp), "months"),
                                     units = "days")),
      ndays = length(unique(lubridate::floor_date(timestamp, "days"))),
      pdays = ndays / duration
    ) |>
    as.data.table()

  # Select individual/month combinations 
  unit_ids_maxgap <- maxgaps$unit_id_tmp[maxgaps$max_gap <= 7]
  unit_ids_pdays  <- pdays$unit_id_tmp[pdays$pdays >= 0.5]
  unit_ids        <- intersect(unit_ids_maxgap, unit_ids_pdays)
  unit_ids        <- unit_ids_maxgap
  n0              <- length(unique(detections$unit_id_tmp))
  n1              <- length(unique(unit_ids))
  detections      <- detections[unit_id_tmp %in% unit_ids, ]
  cat(glue::glue("{n1} / {n0} individual/time block(s) ('unit_id(s)') retained.
                 
                 "))
  
  # Return filtered detection data.table
  detections[, unit_id_tmp := NULL]
  detections
  
}