args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Provide the output file path.")
lines <- c(
  "status=ok",
  paste0("r_version=", R.version.string),
  paste0("platform=", R.version$platform)
)
writeLines(lines, con = args[[1]], useBytes = TRUE)

