# Season definitions (Futia et al. 2024)
# * Winter: December 1–March 31
# * Spring (transitional season): April 1–May 31
# * Summer (stratified season): June 1–September 30
# * Fall (spawning season): October 1–November 30)

# Assign seasonal labels for a vector of time stamps
season_factor <- function(x) {
  m <- lubridate::month(x)
  m <- dplyr::case_when(
    m %in% c(12, 1, 2, 3) ~ "winter",
    m %in% c(4, 5)        ~ "spring",
    m %in% 6:9            ~ "summer",
    m %in% c(10, 11)      ~ "fall",
    TRUE                  ~ NA_character_
  )
  m <- factor(m, levels = c("winter", "spring", "summer", "fall"), 
              labels = c("Winter", "Spring", "Summer", "Fall"))
  m
}