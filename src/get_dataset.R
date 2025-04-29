# Define the timeline using one of two options
# A) Model movements over the same time period for each month
# B) Model movements over the time period of the observations for each month
# -> We select (B)

get_dataset_timeline <- function(.sim) {
  proj.build::check_names(.sim, "file_detections")
  detections <- qs::qread(.sim$file_detections)
  assemble_timeline(.datasets = list(detections), .step = "2 mins")
}