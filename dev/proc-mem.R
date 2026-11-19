# Peak resident set size of this process, in kB, from /proc/self/status
# (VmHWM).  `bench`'s mem_alloc is NA on R built without memory profiling
# (this host: capabilities("profmem") is FALSE), so this is the only peak the
# dev scripts can report here.
proc_peak_rss_kb <- function() {
  if (!file.exists("/proc/self/status")) {
    return(NA_real_)
  }
  line <- grep("^VmHWM:", readLines("/proc/self/status", warn = FALSE), value = TRUE)
  if (!length(line)) {
    return(NA_real_)
  }
  as.numeric(sub("^VmHWM:\\s+([0-9]+) kB.*$", "\\1", line[[1L]]))
}

fmt_peak_rss <- function(kb) {
  if (!is.finite(kb)) {
    return("n/a")
  }
  sprintf("%.1f MB", kb / 1024)
}
