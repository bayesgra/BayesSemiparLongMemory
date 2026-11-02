# Migration entry point. The implementation now lives in R/longmemory_ksb.R.
# Source this file from any directory using source('/path/to/longmemory_DP.R').
local({
  frame_files <- lapply(sys.frames(), function(frame) frame$ofile)
  paths <- Filter(function(x) is.character(x) && length(x) == 1L, frame_files)
  if (!length(paths)) stop("Source this file with source(); see README.md")
  here <- dirname(normalizePath(tail(paths, 1L)[[1L]], mustWork = TRUE))
  sys.source(file.path(here, "R", "longmemory_ksb.R"), envir = parent.env(environment()))
})

longmemory_dp_sampler <- function(...) {
  stop(paste("The old DP model has been replaced.",
             "Use longmemory_ksb_sampler(Y, omega, ...) for the paper's model.",
             "Arguments and output have changed; see docs/migration.md."), call. = FALSE)
}
