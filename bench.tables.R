#!/usr/bin/env Rscript
# Time read/write of t/HepatitisCdata.csv over 10 runs each.
#   read_table  -> read.csv                     (his read_table -> _parse_csv_file)
#   write_table -> write.csv / write.table(sep = ",")
# Reads t/HepatitisCdata.csv; writes the same data to a temp file in /tmp.
# Emits R.read_table.tsv and R.write_table.tsv, one row per run.
# Optional first arg overrides the input path.

save_times <- function(fname, times) {
    out <- data.frame(run = seq_along(times) - 1L, seconds = times)
    write.table(out, fname, sep = "\t", quote = FALSE, row.names = FALSE)
}

main <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    src <- if (length(args) > 0) args[1] else "t/HepatitisCdata.csv"
    d <- read.csv(src, stringsAsFactors = FALSE)
    path <- tempfile(tmpdir = "/tmp", fileext = ".csv")

    write_t <- numeric(10)
    for (i in seq_len(10)) {
        t0 <- Sys.time()
        write.csv(d, path, row.names = FALSE)
        t1 <- Sys.time()
        write_t[i] <- as.numeric(t1) - as.numeric(t0)
    }
    save_times("R.write_table.tsv", write_t)

    read_t <- numeric(10)
    for (i in seq_len(10)) {
        t0 <- Sys.time()
        invisible(read.csv(src, stringsAsFactors = FALSE))
        t1 <- Sys.time()
        read_t[i] <- as.numeric(t1) - as.numeric(t0)
    }
    save_times("R.read_table.tsv", read_t)

    unlink(path)
}

main()
