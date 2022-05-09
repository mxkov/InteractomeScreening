## Simply writes a sessionInfo report.

writeSessionInfo <- function(si.file) {
  writeLines(capture.output(sessionInfo()), si.file)
  message(paste('Session info written to file', si.file))
}
