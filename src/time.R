# Season definitions (Futia et al. 2024)
# * Winter: December 1–March 31
# * Spring (transitional season): April 1–May 31
# * Summer (stratified season): June 1–September 30
# * Fall (spawning season): October 1–November 30)

# Assign seasonal labels for a vector of time stamps
season_factor <- function(x) {
  stopifnot(inherits(x, c("Date", "POSIXct")))
  m <- lubridate::month(x)
  m <- dplyr::case_when(
    m %in% c(12, 1, 2, 3) ~ "winter",
    m %in% c(4, 5)        ~ "spring",
    m %in% 6:9            ~ "summer",
    m %in% c(10, 11)      ~ "fall",
    TRUE                  ~ NA_character_
  )
  m <- factor(m, levels = c("fall", "winter", "spring", "summer"), 
              labels = c("Fall", "Winter", "Spring", "Summer"))
  m
}

# Extract seasonal labels from yyyy-season e.g., 2016-spring
season_factor.ys <- function(x) {
  stopifnot(inherits(x, "character"))
  m <- stringr::str_split_fixed(x, pattern = "-", n = 2L)[, 2]
  m <- factor(m, levels = c("fall", "winter", "spring", "summer"), 
              labels = c("Fall", "Winter", "Spring", "Summer"))
  stopifnot(!is.na(m))
  m
}
