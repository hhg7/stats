#!/usr/bin/env Rscript
# Time reading an .xlsx file over 10 runs and report memory usage.
#
# Equivalent of benchmark.pl's timing loop, but for read_table on .xlsx:
#   read xlsx -> readxl::read_excel
#
# Reads titanic.more.complete.xlsx by default (first arg overrides).
# Emits xlsx.read.R.tsv, one row per run (run\tseconds).
# Prints memory usage to stdout:
#   - peak R heap use during the run (gc "max used")
#   - in-memory object size of the result
#   - process peak resident set size (VmHWM from /proc/self/status)

# User library holds packages rebuilt for this R version (site-library copies
# were compiled against an older R ABI and fail to load).
user_lib <- path.expand("~/R/libs")
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))
suppressMessages(library(readxl))

save_times <- function(fname, times) {
    out <- data.frame(run = seq_along(times) - 1L, seconds = times)
    write.table(out, fname, sep = "\t", quote = FALSE, row.names = FALSE)
}

human <- function(nbytes) {
    units <- c("B", "KiB", "MiB", "GiB")
    i <- 1L
    while (nbytes >= 1024 && i < length(units)) {
        nbytes <- nbytes / 1024
        i <- i + 1L
    }
    sprintf("%.2f %s", nbytes, units[i])
}

vmhwm_bytes <- function() {
    st <- readLines("/proc/self/status")
    line <- grep("^VmHWM:", st, value = TRUE)
    if (length(line) == 0) return(NA_real_)
    as.numeric(sub("\\D+(\\d+).*", "\\1", line)) * 1024  # kB -> bytes
}

main <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    src <- if (length(args) > 0) args[1] else "titanic.more.complete.xlsx"

    gc(reset = TRUE)
    read_t <- numeric(10)
    for (i in seq_len(10)) {
        t0 <- Sys.time()
        d <- read_excel(src)
        t1 <- Sys.time()
        read_t[i] <- as.numeric(t1) - as.numeric(t0)
    }
    save_times("xlsx.read.R.tsv", read_t)

    # gc() matrix column 6 is the "max used" figure in Mb (col 5 is raw cell
    # counts); sum the Ncells and Vcells rows for peak heap over the run.
    g <- gc()
    peak_mb <- sum(g[, 6])
    obj_bytes <- as.numeric(object.size(d))
    rss <- vmhwm_bytes()

    cat(sprintf("file                : %s\n", src))
    cat(sprintf("shape               : %d rows x %d cols\n", nrow(d), ncol(d)))
    cat(sprintf("median read time    : %.6f s\n", median(read_t)))
    cat(sprintf("object in-memory    : %s\n", human(obj_bytes)))
    cat(sprintf("peak R heap (gc)    : %s\n", human(peak_mb * 1024 * 1024)))
    cat(sprintf("process peak RSS    : %s\n", if (is.na(rss)) "n/a" else human(rss)))
}

main()
